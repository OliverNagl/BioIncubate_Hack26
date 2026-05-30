#!/usr/bin/env python
"""Extract 'Table S3' (Hu225 DMS enrichment ratios) from the supplementary PDF to CSV.

The table lists, per CDR position, the enrichment ratio of every NNK codon (31 codons;
TAG stop excluded). Some cells are blank in the PDF (wild-type codon / unobserved
codons). We assign each number to a codon column by its x-coordinate relative to the
header, so blanks are preserved in the correct column rather than shifting values.

Usage:  python extract_table_s3.py [input.pdf] [output.csv]
"""
import csv
import re
import sys
import pathlib

import pdfplumber

PDF = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "data/Table S3.pdf")
OUT = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else "data/Table_S3.csv")

META_HEADERS = {"position", "WT", "codon", "distance"}
POS_RE = re.compile(r"^V[HL]:[A-Z]\d+[a-z]?$")  # e.g. VH:N31, VL:T97, VH:Y100a
CODON_RE = re.compile(r"^[ACGT]{3}$")


def center(w):
    return (w["x0"] + w["x1"]) / 2


def cluster_rows(words):
    """Group words into text rows keyed by rounded vertical position."""
    rows = {}
    for w in words:
        rows.setdefault(round(w["top"]), []).append(w)
    return [sorted(rows[k], key=lambda w: w["x0"]) for k in sorted(rows)]


def main():
    with pdfplumber.open(str(PDF)) as pdf:
        page = pdf.pages[0]
        rows = cluster_rows(page.extract_words(use_text_flow=False))

    # ── Header: find the row with 'position' and read codon column centers ──
    header = next(r for r in rows if any(w["text"] == "position" for w in r))
    codon_cols = [(w["text"], center(w)) for w in header
                  if w["text"] not in META_HEADERS]
    codons = [c for c, _ in codon_cols]
    centers = [x for _, x in codon_cols]
    assert all(CODON_RE.match(c) for c in codons), "unexpected codon headers"

    # half the smallest inter-column gap = max allowed assignment distance
    gaps = [centers[i + 1] - centers[i] for i in range(len(centers) - 1)]
    tol = min(gaps) / 2.0

    # ── Data rows ──
    records = []
    for r in rows:
        if not r or not POS_RE.match(r[0]["text"]):
            continue
        position = r[0]["text"]
        wt_codon = r[1]["text"]
        distance = r[2]["text"]
        assert CODON_RE.match(wt_codon), f"{position}: bad WT codon {wt_codon!r}"

        values = {c: "" for c in codons}
        for w in r[3:]:
            cx = center(w)
            j = min(range(len(centers)), key=lambda k: abs(centers[k] - cx))
            assert abs(centers[j] - cx) <= tol, (
                f"{position}: value {w['text']!r} @ {cx:.1f} not under any codon "
                f"(nearest {codons[j]} Δ{cx-centers[j]:+.1f})")
            assert values[codons[j]] == "", (
                f"{position}: two values map to {codons[j]}")
            values[codons[j]] = w["text"]
        records.append([position, wt_codon, distance] + [values[c] for c in codons])

    # ── Write CSV ──
    OUT.parent.mkdir(exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["position", "WT_codon", "distance"] + codons)
        writer.writerows(records)

    n_blank = sum(1 for rec in records for v in rec[3:] if v == "")
    print(f"Wrote {OUT}: {len(records)} positions x {len(codons)} codons "
          f"(+3 meta cols), {n_blank} blank cells preserved.")
    print(f"Positions: {records[0][0]} ... {records[-1][0]}")


if __name__ == "__main__":
    main()
