# CLAUDE.md — Bioincubate_HACK26: EGFR Binder Optimization

> **Purpose**: Context document for AI agents and human contributors.
> Covers the project goal, structural data, computational methodology,
> and literature-grounded design strategies for improving Cetuximab binding to EGFR.

---

## 1. Project Goal

Design improved variants of **Cetuximab** (anti-EGFR IgG1) with enhanced binding affinity
to the **Epidermal Growth Factor Receptor (EGFR)**, while preserving or improving
developability (solubility, thermal stability, low polyspecificity).

**Primary metric**: equilibrium dissociation constant KD (lower = tighter binding).  
**Baseline**: Cetuximab KD ≈ 560–786 pM (SPR, full IgG1); BLI control in Adaptyv competition ≈ 9.94 nM (cell-free scFv format).

---

## 2. Key Structural Data (Inputs to Notebook)

| ID | Description | Resolution | Relevance |
|----|-------------|-----------|-----------|
| `8HGO` | Cetuximab Fab + EGFR extracellular domain | — | Full paratope–epitope complex |
| `1YY9` | Fab (IMC-11F8) + EGFR Domain III | 2.80 Å | **Primary binding interface** — Domain III contacts |
| `UniProt: P00533` | EGFR canonical sequence (1210 aa) | — | Full target sequence for design |
| `AF: P00533` | AlphaFold2 model of EGFR | — | Predicted structure for hotspot analysis |

**EGFR binding domain focus**: Domains I and III (leucine-rich L-domains).  
Validated hotspot residues for Domain III binders: **R29, Q411, K465, S440, G441**.  
Cetuximab CDR-H3: `ARALTYYDYEFAY` — preserved in all optimization strategies.

---

## 3. Optimization Strategies — Literature Summary

### 3.1 Framework Region Mutation (Cradle approach) ✅ Experimentally validated
**Source**: Adaptyv EGFR Binder Competition paper (Cotet et al., bioRxiv 2025.04.17.648362)

- Start from Cetuximab scFv; preserve all CDRs (primary binding determinants).
- Mutate **framework regions only** (VH/VL non-CDR positions).
- Hypothesis: framework mutations induce small conformational changes in CDRs → improved interface without disrupting binding mode.
- Result: **1.21 nM KD** (8.2× improvement over Cetuximab control in BLI assay).
- Post-competition best: **339 pM** with ≤19 mutations.
- Tools: MMseqs2 (evolutionary data), uniref30 + colabfold_envdb databases.

### 3.2 Zero-Shot Protein Language Model Maturation (ConvergeAB approach) ✅ Experimentally validated
**Source**: Weiner, bioRxiv 2026.05.05.722890

- Platform: ConvergeAB™ — proprietary pLM with ~100,000 target-aware candidates per campaign.
- Input: Cetuximab VH + VL sequences paired with EGFR target sequence. **No task-specific training.**
- Filtering pipeline: binding affinity prediction + docking + thermal stability + solubility predictors.
- Result: **315 pM KD** (mean; 2.1× better than Cetuximab, 4.4× better than Cradle).
- 6 mutations total: `VH T61A`, `VH S87A`, `VH N88D`, `VL V9A`, `VL N32D`, `VL N93A`.
- Global Cα RMSD vs Cetuximab: **0.15 Å** — near-identical fold.
- Mechanism: paratope preorganization (Ala substitutions reduce conformational heterogeneity) + electrostatic optimization (Asn→Asp).
- Developability: DLS radius 5.85 nm (< Cetuximab 6.80 nm), PSR scores indistinguishable from negative control.

### 3.3 Diffusion Sequence Modeling (DSM approach) ✅ In-silico validated
**Source**: Hallee et al., arXiv 2506.08293

- Framework: **Diffusion Sequence Model (DSM)** — extends ESM2 with masked diffusion (LLaDA framework).
- Key advantage: unified representation learning AND generative design in one model.
- Fine-tuned variant **DSMppi**: conditioned on target sequence to generate binders.
- Template-based design: randomly mask 0–100% of known binder → regenerate → score with Synteract2 (ppKd).
- EGFR results: conditional ppKd 8.05 (success rate 0.59%), unconditional ppKd 8.25 (success rate 3.55%).
- Best EGFR design: ppKd 9.14–9.18 (template ppKd 8.92) with 45% of template masked.
- Benchmark: **BenchBB** (Bench-tested Binder Benchmark) — standardized targets including EGFR, IL-7Rα, PD-L1, BHRF1, SpCas9, BBF-14, MBP.

### 3.4 De Novo Design Tools (competition context)
**Source**: Adaptyv competition paper

| Tool | Approach | EGFR hit rate | Notes |
|------|----------|--------------|-------|
| **BindCraft** | AF2-M backpropagation hallucination | 82 nM best de novo | Beta-sheet designs outperformed helical |
| **RFdiffusion** | Diffusion on RoseTTAFold backbone | ~20% general hit rate | Most popular; combine with ProteinMPNN |
| **ProteinMPNN** | Inverse folding (sequence for backbone) | — | Use soluble weights for higher expression |
| **ESM2/ESM3** | Language model scoring & generation | — | ESM3 log-prob correlates with binding strength |

---

## 4. Computational Metrics for Ranking & Filtering

Based on Adaptyv competition analysis (ROC AUC scores):

### Expression prediction (binder alone):
| Metric | ROC AUC | Interpretation |
|--------|---------|----------------|
| `composition_GLU` (% glutamate) | **0.77** | Higher → better expression |
| `composition_LYS` (% lysine) | **0.73** | Higher → better expression |
| Rosetta `fa_intra_sol_xover4_per_aa` | 0.68 | Solvation energy |
| `ss_prop_alpha_helix` | 0.67 | Helix content |
| ipTM | 0.58 | Weak predictor of expression |
| ESM2 PLL | 0.44 | Poor predictor of expression |

### Binding prediction (complex with EGFR):
| Metric | ROC AUC |
|--------|---------|
| ESM2 PLL (normalized by length) | **0.72** |
| Foldseek `fident_pdb` | 0.68 |
| pLDDT | 0.66 |
| ipTM | 0.64 |

**Key insight**: ipTM and iPAE are useful screens but **not strong predictors** of experimental KD.
ESM3/ESM-C log-probabilities show stronger correlation to binding categories than ESM2.

### Interface quality metrics (successful binders vs non-binders):
- `pDockQ / mpDockQ` — relevant for all molecule types
- `int_area` — larger interface → more likely binder
- `binding_energy` — useful but molecule-type-specific
- Successful binders: more hydrogen bonds, better ΔG, more interface contacts, better shape complementarity

---

## 5. Design Workflow for This Project

```
1. STRUCTURAL ANALYSIS
   ├── Fetch 8HGO, 1YY9, P00533 (UniProt + AF)
   ├── Align structures → identify binding interface residues
   ├── Map EGFR Domain III contacts (1YY9 focus: residues within 4Å)
   └── Identify Cetuximab CDR positions (IMGT numbering)

2. SEQUENCE PREPARATION
   ├── Extract Cetuximab VH + VL from PDB/literature
   ├── Format as single-chain scFv (for competition-style scoring)
   └── Identify framework vs CDR regions

3. CANDIDATE GENERATION (choose one or combine)
   ├── Strategy A: Framework mutation scan (Cradle-style)
   │   ├── Tools: ESM2 log-likelihood, ColabFold AF2-M
   │   └── Constrain: CDRs fixed, vary framework positions
   ├── Strategy B: Masked sequence regeneration (DSM-style)
   │   ├── Tools: ESM2/ESM3, random masking 0-50%
   │   └── Score with AlphaFold2-Multimer ipTM/iPAE
   └── Strategy C: Direct pLM-guided evolution
       ├── Tools: ESM2 pseudo-log-likelihood maximization
       └── Iterative MLDE loop

4. FILTERING PIPELINE
   ├── ColabFold AF2-M → ipTM > 0.7, iPAE < 15
   ├── ESM2 PLL (length-normalized) → percentile filter
   ├── Rosetta/DE-STRESS → aggregation, solvation energy
   ├── NetSolP / SoluProt → solubility prediction
   └── Aggrescan3D → aggregation propensity

5. EXPERIMENTAL VALIDATION (future)
   └── BLI / SPR affinity measurement (Adaptyv pipeline)
```

---

## 6. Cetuximab Sequence Reference

**VH CDRs (IMGT)**:
- CDR-H1: `GFSLTNYG` (positions ~27-38)
- CDR-H2: `IWSGGNT` (positions ~56-65)
- CDR-H3: `ARALTYYDYEFAY` (positions ~105-117) ← **primary binding determinant**

**VL CDRs (IMGT)**:
- CDR-L1: `RASQSIGTNIT` (positions ~27-38)
- CDR-L2: `YASESIS` (positions ~56-65)
- CDR-L3: `QQNNNWPTT` (positions ~105-117)

**Key mutations from literature**:
- Cradle winner: VH {5,70,71,73,87,88} + VL {45,49,79,80} — framework only
- ConvergeAB: `VH T61A`, `VH S87A`, `VH N88D`, `VL V9A`, `VL N32D`, `VL N93A`

---

## 7. Environment & Reproducibility

- Python venv with `requirements.txt` — auto-installs on repo pull.
- Notebook: `egfr_binder_analysis.ipynb`
- Key packages: `biotite`, `py3Dmol`, `nglview`, `biopython`, `requests`, `mdanalysis`.
- PDB fetching: RCSB REST API (`https://files.rcsb.org/download/{ID}.pdb`).
- AlphaFold fetching: `https://alphafold.ebi.ac.uk/files/AF-{UNIPROT}-F1-model_v4.pdb`.
- UniProt fetching: `https://rest.uniprot.org/uniprotkb/{ID}.fasta`.

---

## 8. Key References

1. Cotet et al. (2025). *Crowdsourced Protein Design: Lessons From the Adaptyv EGFR Binder Competition.* bioRxiv 2025.04.17.648362.
2. Weiner, I. (2026). *Zero-Shot Design of a Biobetter Cetuximab.* bioRxiv 2026.05.05.722890.
3. Hallee et al. (2025). *Diffusion Sequence Models for Enhanced Protein Representation and Generation.* arXiv 2506.08293.
4. Pacesa et al. (2024). *BindCraft: one-shot design of functional protein binders.* bioRxiv 2024.09.30.615802.
5. Watson et al. (2023). *De novo design of protein structure and function with RFdiffusion.* Nature 620.
6. Lin et al. (2023). *Evolutionary-scale prediction of atomic-level protein structure with a language model.* Science 379.

---

## 9. BenchBB Benchmark Context

The **Bench-tested Binder Benchmark** (proposed in the Adaptyv paper) standardizes evaluation:
- Targets: EGFR, IL-7Rα, PD-L1, BHRF1, SpCas9, BBF-14, MBP
- Methods failing to yield measurable-affinity binders across a practical screening set are unlikely to outperform current SOTA.
- EGFR best known binder for DSM benchmark: Cetuximab-derived, pKd 8.92 (KD ~1.2 nM scFv format).

---

*Last updated: May 2026 | Synthesized from Adaptyv 2025, Converge Bio 2026, Hallee/Synthyra 2025*
