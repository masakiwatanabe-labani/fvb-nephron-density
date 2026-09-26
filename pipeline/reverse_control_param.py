#!/usr/bin/env python3
"""Genome-driven count shifts (D11): B6 reads counted on the personalised genome vs on GRCm39.
The hybrid matrix takes C57BL/6J from GRCm39 and FVB/N from the personalised genome, so any gene whose
B6 count changes when ONLY the genome changes can appear differentially expressed for technical reasons.
Parameterised version of rerun2/G_hybrid_artefact/reverse_control_screen.py (paths as arguments).
Usage: reverse_control_param.py --uncorrected M --arm M --gene-chr F --out-prefix P [--min-reads 60] [--threshold 1.0]
"""
import argparse
import numpy as np, pandas as pd
ap = argparse.ArgumentParser()
ap.add_argument("--uncorrected", required=True); ap.add_argument("--arm", required=True)
ap.add_argument("--gene-chr", required=True); ap.add_argument("--out-prefix", required=True)
ap.add_argument("--min-reads", type=int, default=60); ap.add_argument("--threshold", type=float, default=1.0)
a = ap.parse_args()
U = pd.read_csv(a.uncorrected, sep="\t", index_col=0)
P = pd.read_csv(a.arm, sep="\t", index_col=0)
chrom = pd.read_csv(a.gene_chr, sep="\t").set_index("gene_id")["chr"]
flagged = set()
for tp, tag in [("E13.5", "135"), ("P1", "P1")]:
    b6 = [c for c in U.columns if c.startswith("B6") and tag in c]
    bu, bp = U[b6].sum(axis=1), P[b6].sum(axis=1)
    ok = (bu + bp) >= a.min_reads
    lr = np.log2((bp + 1) / (bu + 1))
    f = lr.index[ok & (lr.abs() > a.threshold)]
    pd.DataFrame({"Geneid": f, "chr": chrom.reindex(f).values,
                  "B6_reads_GRCm39": bu[f].values, "B6_reads_arm": bp[f].values,
                  "log2_shift": lr[f].round(3).values}).to_csv(f"{a.out_prefix}_{tp}.tsv", sep="\t", index=False)
    flagged |= set(f)
    print(f"{tp}: {int(ok.sum())} genes with >= {a.min_reads} reads; {len(f)} flagged (|log2 shift| > {a.threshold})")
print(f"union across stages: {len(flagged)}")
