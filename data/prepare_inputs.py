#!/usr/bin/env python
"""Prepare AF3-style inputs for Cetuximab scFv x EGFR Domain III co-folding.

Standalone input-prep step. Depends only on ``biopython`` + ``requests``.

For every variant in ``cetuximab_variants.fasta`` (headers ``>{name}_VH`` /
``>{name}_VL``) this builds an scFv (``VH-(GGGGS)3-VL``) and writes, per variant:

    inputs/<variant>/scfv.fasta    single-entry FASTA (>{variant}_scfv)
    inputs/<variant>/complex.json  AF3-style JSON shared by AF3 *and* ESMFold2

EGFR Domain III is extracted once from PDB 1YY9 (chain A, residues 1-165 as
deposited -- no renumbering, insertion codes preserved) to:

    data/egfr_domain3.pdb
    data/egfr_domain3.fasta

Run standalone:  python prepare_inputs.py
Exit code 0 on success, non-zero on any failure.
"""
from __future__ import annotations

import sys
from pathlib import Path

# NOTE: sweep_pipeline.py contains no reusable AF3 JSON-builder -- its
# generate_sweep_json() builds RFD3 specs and build_symmetry_strings() builds
# MPNN tie strings. So the AF3-style complex.json builder lives here.

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
INPUTS_DIR = SCRIPT_DIR / "inputs"
FASTA = SCRIPT_DIR / "cetuximab_variants.fasta"

LINKER = "GGGGSGGGGSGGGGS"          # (GGGGS)3
PDB_ID = "1YY9"
EGFR_CHAIN = "A"
EGFR_RESRANGE = (1, 165)            # residues 1-165 as deposited
MODEL_SEEDS = [42]


# ── FASTA parsing ────────────────────────────────────────────────────────


def parse_variants(fasta_path: Path) -> dict[str, dict[str, str]]:
    """Parse ``>{name}_VH`` / ``>{name}_VL`` records into {name: {VH, VL}}.

    Names missing either chain are warned about and dropped downstream.
    """
    from Bio import SeqIO

    variants: dict[str, dict[str, str]] = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        header = record.id
        if "_" not in header:
            print(f"[prepare] WARN: header {header!r} lacks a _VH/_VL suffix; skipping")
            continue
        name, _, chain = header.rpartition("_")
        chain = chain.upper()
        if chain not in ("VH", "VL"):
            print(f"[prepare] WARN: header {header!r} suffix {chain!r} not VH/VL; skipping")
            continue
        variants.setdefault(name, {})[chain] = str(record.seq)
    return variants


# ── EGFR Domain III extraction ───────────────────────────────────────────


def _fetch_pdb(pdb_id: str, dest: Path) -> None:
    """Download a PDB from RCSB to ``dest`` (only if absent upstream)."""
    import requests

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"[prepare] {dest} absent -> fetching {url}")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(resp.text)


def extract_egfr_domain3() -> str:
    """Extract EGFR D3 (1YY9 chain A, res 1-165 as deposited) -> sequence.

    Writes ``data/egfr_domain3.pdb`` (subset, original numbering preserved) and
    ``data/egfr_domain3.fasta``. Returns the one-letter sequence.
    """
    from Bio.PDB import PDBIO, PDBParser, Select
    from Bio.PDB.Polypeptide import index_to_one, is_aa, three_to_index

    pdb_path = DATA_DIR / f"{PDB_ID}.pdb"
    if not pdb_path.is_file():
        _fetch_pdb(PDB_ID, pdb_path)

    structure = PDBParser(QUIET=True).get_structure(PDB_ID, str(pdb_path))
    model = next(structure.get_models())
    if EGFR_CHAIN not in model:
        raise ValueError(f"chain {EGFR_CHAIN} not found in {pdb_path}")

    lo, hi = EGFR_RESRANGE

    class _D3Select(Select):
        """Keep only standard residues lo..hi of the EGFR chain."""

        def accept_chain(self, chain):  # noqa: N802 (biopython API)
            return chain.id == EGFR_CHAIN

        def accept_residue(self, residue):  # noqa: N802
            hetflag, resseq, _icode = residue.id
            return hetflag == " " and lo <= resseq <= hi

    # Build the sequence from the same accepted residues, in deposited order.
    seq_chars: list[str] = []
    for residue in model[EGFR_CHAIN]:
        hetflag, resseq, _icode = residue.id
        if hetflag != " " or not (lo <= resseq <= hi):
            continue
        if not is_aa(residue, standard=True):
            continue
        seq_chars.append(index_to_one(three_to_index(residue.resname)))
    sequence = "".join(seq_chars)
    if not sequence:
        raise ValueError(
            f"no standard residues in {PDB_ID} chain {EGFR_CHAIN} range {lo}-{hi}"
        )

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(DATA_DIR / "egfr_domain3.pdb"), _D3Select())
    (DATA_DIR / "egfr_domain3.fasta").write_text(f">egfr_domain3\n{sequence}\n")
    print(f"[prepare] EGFR D3: {len(sequence)} aa -> data/egfr_domain3.{{pdb,fasta}}")
    return sequence


# ── AF3-style complex.json ───────────────────────────────────────────────


def build_complex_json(name: str, scfv_seq: str, egfr_seq: str) -> dict:
    """AF3-style two-chain co-fold JSON (consumed by both AF3 and ESMFold2).

    SPEC DEVIATION (justified): the prompt specified the AlphaFold-*Server*
    dialect (``{"proteinChain": {"sequence", "count"}}``). The ESMFold2
    container's AF3 reader (``_esmfold2_inputs._parse_af3_sequences``) only
    accepts the official ``alphafold3`` dialect -- each entry must be keyed
    ``protein``/``dna``/``rna``/``ligand`` with an explicit chain ``id`` -- and
    raises "must contain exactly one of protein, dna, rna, or ligand" on the
    Server dialect (verified at runtime, job 2314545). We therefore emit the
    official dialect, which is the canonical AF3 input format and still the
    single shared file for both engines (constraint 3 preserved).
    """
    return {
        "name": f"{name}_complex",
        "modelSeeds": list(MODEL_SEEDS),
        "sequences": [
            {"protein": {"id": "A", "sequence": scfv_seq}},
            {"protein": {"id": "B", "sequence": egfr_seq}},
        ],
        "dialect": "alphafold3",
        "version": 1,
    }


# ── Main ─────────────────────────────────────────────────────────────────


def main() -> int:
    if not FASTA.is_file():
        print(f"[prepare] ERROR: missing {FASTA}", file=sys.stderr)
        return 1

    import json

    variants = parse_variants(FASTA)
    if not variants:
        print(f"[prepare] ERROR: no variants parsed from {FASTA}", file=sys.stderr)
        return 1

    egfr_seq = extract_egfr_domain3()

    rows: list[tuple[str, int, int, int, int]] = []
    for name in sorted(variants):
        chains = variants[name]
        vh, vl = chains.get("VH"), chains.get("VL")
        if not vh or not vl:
            missing = "VH" if not vh else "VL"
            print(f"[prepare] WARN: variant {name!r} missing {missing}; skipping")
            continue

        scfv = vh + LINKER + vl
        out_dir = INPUTS_DIR / name
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "scfv.fasta").write_text(f">{name}_scfv\n{scfv}\n")
        with open(out_dir / "complex.json", "w") as fh:
            json.dump(build_complex_json(name, scfv, egfr_seq), fh, indent=2)
        rows.append((name, len(vh), len(vl), len(scfv), len(egfr_seq)))
        print(f"[prepare] {name}: inputs/{name}/{{scfv.fasta,complex.json}}")

    if not rows:
        print("[prepare] ERROR: no complete variants written", file=sys.stderr)
        return 1

    # Summary table
    header = ("variant", "VH_len", "VL_len", "scFv_len", "EGFR_D3_len")
    width = max(len(header[0]), *(len(r[0]) for r in rows))
    print("\n" + f"{header[0]:<{width}}  " + "  ".join(f"{h:>11}" for h in header[1:]))
    for r in rows:
        print(f"{r[0]:<{width}}  " + "  ".join(f"{v:>11}" for v in r[1:]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
