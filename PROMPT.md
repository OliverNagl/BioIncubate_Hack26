# OPTIMIZED PROMPT — Bioincubate_HACK26
# Target: Claude Opus 4.8 CLI (claude --model claude-opus-4-8)
# Usage:  cat PROMPT.md | claude --model claude-opus-4-8
#         or paste directly into claude.ai

---

## Task

You are working inside the repository `Bioincubate_HACK26/`. Read `CLAUDE.md` first for
full project context, methodology, and structural data definitions.

Your goal is to set up and execute a protein structure analysis + visualization pipeline
for EGFR binder optimization, starting from Cetuximab.

---

## Environment setup

1. Create a Python venv at `.venv/` and install from `requirements.txt`.
   If `requirements.txt` is missing, generate it with these exact packages:
   `biopython>=1.83`, `py3Dmol>=2.1.0`, `nglview>=3.1.0`, `ipywidgets>=8.1.0`,
   `biotite>=0.40.0`, `numpy>=1.26.0`, `pandas>=2.2.0`, `matplotlib>=3.8.0`,
   `requests>=2.31.0`, `jupyterlab>=4.1.0`, `ipykernel>=6.29.0`.

2. Register the venv as a Jupyter kernel named `hack26`:
   `.venv/bin/python -m ipykernel install --user --name hack26 --display-name "Hack26"`

---

## Notebook: `egfr_binder_analysis.ipynb`

The notebook must run top-to-bottom without errors. If it already exists, validate and
fix any cells that would fail. If it does not exist, create it from scratch.

### Required cells in order:

**Cell 0 — Auto-install** (run `pip install -r requirements.txt -q` via subprocess;
print ✅ on success, print stderr on failure; skip if already installed).

**Cell 1 — Fetch structures** into `data/` directory:
- `8HGO.pdb` from `https://files.rcsb.org/download/8HGO.pdb`
- `1YY9.pdb` from `https://files.rcsb.org/download/1YY9.pdb`
- `AF_P00533.pdb` from `https://alphafold.ebi.ac.uk/files/AF-P00533-F1-model_v4.pdb`
- `UNIPROT_P00533.fasta` from `https://rest.uniprot.org/uniprotkb/P00533.fasta`
- Use `requests` with a 60 s timeout; skip if file already exists; print file sizes.

**Cell 2 — Parse & summarise** with `biotite`:
- Load each PDB, filter to canonical amino acids, print chain IDs and residue ranges.
- Identify which chain is EGFR and which are Fab H/L in each structure.
- Assign: `EGFR_CHAINS_1YY9`, `FAB_CHAINS_1YY9`, `EGFR_CHAINS_8HGO`, `FAB_CHAINS_8HGO`.

**Cell 3 — Interface mapping** (1YY9 focus, 4 Å cutoff):
- Compute Fab–EGFR contact residues using pairwise heavy-atom distances.
- Output two DataFrames: `egfr_contacts_1yy9`, `fab_contacts_1yy9`.
- Print residue lists; cross-reference against known hotspots R29, Q411, K465, S440, G441.

**Cell 4 — Visualize 1YY9** with `py3Dmol` (900×500 px, dark background `#1a1a2e`):
- EGFR cartoon in light blue; Fab heavy chain green, light chain orange.
- Interface residues as colored sticks (EGFR red, Fab yellow).
- Add residue number labels for top-5 highest-contact EGFR residues.
- Save screenshot path reference in a markdown cell.

**Cell 5 — Visualize 8HGO** (same style as Cell 4).

**Cell 6 — Visualize AF-P00533**:
- Full EGFR cartoon coloured by pLDDT (`colorscheme: "pLDDT"`).
- Highlight Domain III (residues 340–510) in blue.
- Label key hotspot residues (R29, Q411, K465, S440, G441) as red sticks with text labels.
- Zoom to Domain III region.

**Cell 7 — Sequence analysis**:
- Load EGFR FASTA; report length, signal peptide boundary (aa 1–24), Domain III range.
- Hardcode Cetuximab VH and VL reference sequences (from PDB chains in 8HGO/1YY9).
- Annotate CDR regions with IMGT positions; print CDR sequences.

**Cell 8 — Known mutations table + KD plot**:
- Print a comparison table: Cetuximab vs Cradle (1.21 nM) vs ConvergeAB (315 pM).
- Mutations per variant (VH count, VL count, specific positions).
- Bar chart (log scale y-axis) with dark theme; save as `data/kd_comparison.png`.

**Cell 9 — In-silico mutagenesis**:
- Apply ConvergeAB mutations (`VH T61A, S87A, N88D` | `VL V9A, N32D, N93A`) to
  Cetuximab sequences using a `apply_mutations(seq, [(pos, orig, new), ...])` function.
- Print final mutant sequences; compute and print % sequence identity to parent.
- Export both variants as FASTA to `data/cetuximab_variants.fasta`.

**Cell 10 — Next steps markdown cell** listing:
- ColabFold AF2-M scoring pipeline
- ESM2/ESM3 log-likelihood scoring
- BindCraft hallucination campaign setup
- (See `CLAUDE.md` Section 5 for full workflow)

---

## Constraints

- Do NOT use `nglview` as primary visualizer (widget rendering is unreliable in some
  Jupyter environments). Use `py3Dmol` for all 3D views.
- Do NOT hard-code absolute paths. All paths relative to notebook location.
- Each visualization cell must call `.show()` so the view renders inline.
- Handle fetch failures gracefully: print warning and continue; do not crash notebook.
- Keep all cells idempotent: re-running the full notebook must not duplicate data.

---

## Deliverables

After completing the notebook, confirm:
1. `data/` contains all 4 fetched files.
2. All 10 cells execute without exceptions.
3. `data/kd_comparison.png` exists.
4. `data/cetuximab_variants.fasta` exists with 4 sequences (VH+VL × 2 variants).
5. Print a final summary: structure chain assignments, interface residue counts,
   total mutations applied.
