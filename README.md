# Bioincubate_HACK26 — EGFR Binder Optimization

> Hackathon project: improve Cetuximab binding to EGFR using computational protein design.

## Quick start

```bash
# Clone and enter repo
git clone <repo-url> && cd Bioincubate_HACK26

# Create venv and install
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# Register Jupyter kernel
python -m ipykernel install --user --name hack26 --display-name "Hack26"

# Launch notebook
jupyter lab egfr_binder_analysis.ipynb
```

The notebook **auto-installs dependencies** in Cell 0, so pulling the repo and running
the notebook is sufficient on a fresh machine (Python 3.10+ required).

## Files

| File | Description |
|------|-------------|
| `CLAUDE.md` | Full methodology context (read this first) |
| `PROMPT.md` | Optimized prompt for Claude Opus 4.8 CLI |
| `egfr_binder_analysis.ipynb` | Main analysis notebook |
| `requirements.txt` | Python dependencies |
| `data/` | Fetched structures + outputs (auto-created) |

## Structures analysed

- **8HGO** — Cetuximab Fab + EGFR extracellular domain
- **1YY9** — IMC-11F8 Fab + EGFR Domain III (2.80 Å) — primary binding interface
- **AF-P00533** — AlphaFold2 full EGFR model
- **UniProt P00533** — EGFR canonical sequence

## Key literature

See `CLAUDE.md` Section 8 for full references. Key findings:
- Cradle (framework mutations only) → **1.21 nM** (8.2× improvement)  
- ConvergeAB (6 mutations, zero-shot pLM) → **315 pM** (2.1× better than Cetuximab)
- DSM (masked diffusion design) → ppKd 9.14 in-silico

## Claude CLI usage

```bash
# Run the optimized prompt against the full project
claude --model claude-opus-4-8 < PROMPT.md
```
