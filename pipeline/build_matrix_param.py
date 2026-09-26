#!/usr/bin/env python3
"""Assemble a gene x sample count matrix from featureCounts per-sample outputs.
Usage: build_matrix_param.py --counts-dir DIR --pattern 'SAMPLE.<tag>.counts.txt' --out matrix.tsv [--summary out.tsv]"""
import argparse, os
import pandas as pd
ap = argparse.ArgumentParser()
ap.add_argument("--counts-dir", required=True)
ap.add_argument("--pattern", required=True, help="file name with SAMPLE where the sample name goes")
ap.add_argument("--samples", default="B6135-1,B6135-2,B6135-3,B6P1-1,B6P1-2,B6P1-3,FVB135-1,FVB135-2,FVB135-3,FVBP1-1,FVBP1-2,FVBP1-3")
ap.add_argument("--out", required=True); ap.add_argument("--summary", default=None)
a = ap.parse_args()
cols, summ = {}, {}
for s in a.samples.split(","):
    f = os.path.join(a.counts_dir, a.pattern.replace("SAMPLE", s))
    d = pd.read_csv(f, sep="\t", comment="#")
    cols[s] = pd.Series(d.iloc[:, -1].values, index=d["Geneid"].values)
    sf = f + ".summary"
    if os.path.exists(sf):
        t = pd.read_csv(sf, sep="\t", index_col=0)
        summ[s] = t.iloc[:, 0]
m = pd.DataFrame(cols)
m.index.name = "gene_id"
m.to_csv(a.out, sep="\t")
print(f"matrix: {m.shape[0]} genes x {m.shape[1]} samples -> {a.out}")
print(f"total assigned counts per sample (millions):\n{(m.sum()/1e6).round(2).to_string()}")
if a.summary and summ:
    pd.DataFrame(summ).to_csv(a.summary, sep="\t")
    print(f"summary -> {a.summary}")
