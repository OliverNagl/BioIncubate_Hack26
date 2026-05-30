#!/usr/bin/env -S bash -c '"${PIPELINES_REPO:-../..}/.venv/bin/python" "$0" "$@"'
"""Co-fold Cetuximab scFv variants with EGFR Domain III using AF3 + ESMFold2.

Orchestration script modelled on ``../sweep_pipeline.py`` (absl flags +
``init_cli()`` + per-job timing/logging). Both engines consume the *same*
``inputs/<variant>/complex.json`` produced by ``prepare_inputs.py``:

  * AF3    -- driven through the pipelines AF3 *block* via ``engine.run(af3)``,
              exactly as sweep_pipeline.py invokes AF3 (the repo has no raw AF3
              command / SIF; the framework resolves the container + weights).
              The block builds its own AF3 JSON from a structure, so we hand it
              a 2-chain CA-only structure built from the SAME complex.json
              sequences -- the JSON remains the single source of truth.
  * ESMFold2 -- NEW: a self-contained Apptainer container that reads the
              AF3-style complex.json directly (apptainer subprocess).

Usage:
    ./fold_pipeline.py                       # all variants, both engines
    ./fold_pipeline.py --variants Cetuximab --engines esm2
    ./fold_pipeline.py --dry-run
    ./fold_pipeline.py --seeds 42,123
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

# Normalise the documented hyphenated flag to absl's underscore form so both
# `--dry-run` (as advertised) and `--dry_run` work.
sys.argv = ["--dry_run" if a == "--dry-run" else a for a in sys.argv]

from absl import flags

from blocks import Blocks
from pipelines import ExecutionEngine, PipelineDataHandler, init_cli

SCRIPT_DIR = Path(__file__).resolve().parent
INPUTS_DIR = SCRIPT_DIR / "inputs"
RESULTS_DIR = SCRIPT_DIR / "results"

# Container paths -- hardcoded per spec (the only permitted absolute paths).
ESMFOLD2_SIF = "/net/software/containers/users/magnusb/esmfold2.sif"

# ESM perplexity scoring (not folding) -- separate tool for future use:
# apptainer run --nv -B /software/scripts/esm /software/containers/esm.sif \
#   /software/scripts/esm/score_sequences.py --help
# Perplexity 1-20; natural sequences ~6; filter threshold ~6-10 for de novo designs.

AF3_SAMPLES = 5  # diffusion samples per AF3 fold

# ── CLI flags (mirrors sweep_pipeline.py: define, then init_cli) ──────────

flags.DEFINE_list("variants", ["all"],
                  "Comma-separated variant names to fold, or 'all'.")
flags.DEFINE_list("engines", ["all"],
                  "Folding engines: af3 | esm2 | all.")
flags.DEFINE_bool("dry_run", False,
                  "Print the commands that would run; execute nothing.")
flags.DEFINE_list("seeds", ["42"],
                  "Comma-separated model seeds (written to each complex.json).")

init_cli()
FLAGS = flags.FLAGS

(AF3,) = Blocks.getall("af3")

# 3-letter codes so a CA-only AtomArray parses back as a polypeptide chain.
_AA1_TO_3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}


# ── Helpers ──────────────────────────────────────────────────────────────


def discover_variants() -> list[str]:
    """Variant dirs under inputs/ that carry a complex.json."""
    if not INPUTS_DIR.is_dir():
        return []
    return sorted(
        p.name for p in INPUTS_DIR.iterdir()
        if (p / "complex.json").is_file()
    )


def resolve_variants() -> list[str]:
    available = discover_variants()
    requested = [v.strip() for v in FLAGS.variants if v.strip()]
    if not requested or requested == ["all"]:
        return available
    chosen, missing = [], []
    for v in requested:
        (chosen if v in available else missing).append(v)
    for v in missing:
        print(f"[fold] WARN: requested variant {v!r} has no inputs/{v}/complex.json; skipping")
    return chosen


def resolve_engines() -> list[str]:
    req = [e.strip().lower() for e in FLAGS.engines if e.strip()]
    if not req or "all" in req:
        return ["af3", "esm2"]
    return [e for e in req if e in ("af3", "esm2")]


def apply_seeds(complex_json: Path) -> None:
    """Write the requested --seeds into the shared complex.json modelSeeds.

    Keeps ONE complex.json per variant (read by both engines); default
    [42] matches prepare_inputs.py, so this is a no-op unless --seeds changes.
    """
    seeds = [int(s) for s in FLAGS.seeds]
    data = json.loads(complex_json.read_text())
    if data.get("modelSeeds") != seeds:
        data["modelSeeds"] = seeds
        complex_json.write_text(json.dumps(data, indent=2))
        print(f"[fold]   set modelSeeds={seeds} in {complex_json}")


def build_two_chain_structure(complex_json: Path, out_cif: Path) -> Path:
    """CA-only 2-chain structure (NaN coords) from complex.json sequences.

    The pipelines AF3 block builds its own AF3 JSON from a structure's atom
    array (sequence + chain composition); coordinates are irrelevant since AF3
    folds from sequence. Mirrors handler._sequence_to_cif, extended to the two
    proteinChains of the complex (chains A=scFv, B=EGFR D3).
    """
    import numpy as np
    from atomworks.io.utils.io_utils import to_cif_file
    from biotite.structure import AtomArray

    # complex.json is the official alphafold3 dialect: {"protein": {"id","sequence"}}.
    seqs = [
        s["protein"]["sequence"]
        for s in json.loads(complex_json.read_text())["sequences"]
        if "protein" in s
    ]
    residues = []  # (chain_id, res_id, res_name)
    for chain_id, seq in zip("ABCDEFGH", seqs):
        for i, aa in enumerate(seq.upper()):
            if aa in _AA1_TO_3:
                residues.append((chain_id, i + 1, _AA1_TO_3[aa]))

    arr = AtomArray(len(residues))
    for i, (chain_id, res_id, res_name) in enumerate(residues):
        arr.coord[i] = [np.nan, np.nan, np.nan]
        arr.chain_id[i] = chain_id
        arr.res_id[i] = res_id
        arr.res_name[i] = res_name
        arr.atom_name[i] = "CA"
        arr.element[i] = "C"
    out_cif.parent.mkdir(parents=True, exist_ok=True)
    to_cif_file(arr, str(out_cif))
    return out_cif


# ── Engines ──────────────────────────────────────────────────────────────


def run_af3(variant: str, complex_json: Path) -> bool:
    """Fold via the pipelines AF3 block (engine.run), like sweep_pipeline.py."""
    out_dir = RESULTS_DIR / "af3" / variant
    out_dir.mkdir(parents=True, exist_ok=True)
    cif = build_two_chain_structure(complex_json, out_dir / f"{variant}_input.cif")
    handler = PipelineDataHandler.from_paths([str(cif)], name=f"{variant}_complex")

    job_args = {
        # A100-class GPU on gpu-train (sm_70..sm_90 JAX kernels + enough VRAM),
        # mirroring sweep_pipeline.py's AF3 stage.
        "partition": "gpu-train",
        "gres": "gpu:large:1",
        "constraint": "Ampere",
        "mem_gb": 64,
        "rows_per_job": 1,
    }
    print(f"[fold] AF3 {variant}: engine.run(af3, data={cif.name}, "
          f"diffusion_batch_size={AF3_SAMPLES}) -> {out_dir}")
    if FLAGS.dry_run:
        return True

    engine = ExecutionEngine(rundir=out_dir / "_run")
    result = engine.run(AF3(diffusion_batch_size=AF3_SAMPLES), data=handler,
                        job_args=job_args)

    # Surface AF3's native outputs (model CIF + confidences) into out_dir.
    n_linked = 0
    for structures, _ in result:
        for s in structures:
            p = s.get("path")
            if p and Path(p).is_file():
                dest = out_dir / Path(p).name
                if dest != Path(p):
                    dest.write_bytes(Path(p).read_bytes())
                    n_linked += 1
    print(f"[fold] AF3 {variant}: surfaced {n_linked} output file(s) under {out_dir}")
    return n_linked > 0


def run_esm2(variant: str, complex_json: Path) -> bool:
    """Fold via the ESMFold2 Apptainer container (reads complex.json directly)."""
    out_dir = RESULTS_DIR / "esm2" / variant
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "apptainer", "exec", "--nv",
        ESMFOLD2_SIF,
        "esmfold2_predict",
        "--input", str(complex_json),
        "--output-dir", str(out_dir),
        "--full-metrics",   # writes *_confidences.pkl next to each CIF
    ]
    print(f"[fold] ESM2 {variant}: {' '.join(cmd)}")
    if FLAGS.dry_run:
        return True

    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        print(f"[fold] ESM2 {variant}: non-zero exit {proc.returncode}; continuing")
        return False

    # The container writes the pkl (via --full-metrics); we only verify it landed.
    ok = True
    for pattern, label in (
        ("*_model.cif", "model CIF"),
        ("*_summary_confidences.json", "summary confidences"),
        ("*_confidences.pkl", "tensor confidences (pkl)"),
    ):
        hits = list(out_dir.rglob(pattern))
        if not hits:
            print(f"[fold] ESM2 {variant}: WARN missing {label} ({pattern}) under {out_dir}")
            ok = False
        elif label.startswith("tensor"):
            for h in hits:
                print(f"[fold] ESM2 {variant}: pkl -> {h}")
    return ok


ENGINES = {"af3": run_af3, "esm2": run_esm2}


# ── Driver ───────────────────────────────────────────────────────────────


def main() -> int:
    variants = resolve_variants()
    engines = resolve_engines()
    if not variants:
        print("[fold] ERROR: no variants to fold (run prepare_inputs.py first?)",
              file=sys.stderr)
        return 1
    if not engines:
        print(f"[fold] ERROR: no valid engines in {FLAGS.engines}", file=sys.stderr)
        return 1

    print(f"[fold] variants={variants} engines={engines} "
          f"seeds={FLAGS.seeds} dry_run={FLAGS.dry_run}")

    n_total = n_success = 0
    for variant in variants:
        complex_json = INPUTS_DIR / variant / "complex.json"
        apply_seeds(complex_json)
        for engine in engines:
            n_total += 1
            t0 = time.time()
            print(f"\n[fold] === {variant} / {engine} ===  start={time.strftime('%H:%M:%S')}")
            try:
                ok = ENGINES[engine](variant, complex_json)
            except Exception as exc:  # continue on failure
                ok = False
                print(f"[fold] {variant}/{engine}: EXCEPTION {exc}")
            dt = time.time() - t0
            status = "ok" if ok else "FAILED"
            print(f"[fold] === {variant} / {engine} === {status} ({dt:.1f}s)")
            n_success += int(ok)

    print(f"\n[fold] {n_success}/{n_total} jobs completed successfully")
    return 0 if n_success == n_total else 1


if __name__ == "__main__":
    sys.exit(main())
