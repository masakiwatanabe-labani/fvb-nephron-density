#!/usr/bin/env python3
"""Forward and reverse genome-swap shifts, computed as the same quantity on the same genes.

Results 2.3 / Figure 3F compares a reverse-control quantity (B6 reads counted on the personalised
genome vs on GRCm39, mean |log2| 0.026-0.033) against 1.99/1.88, which is the mean |log2 fold change|
between FVB/N and C57BL/6J in genes selected for having a large strain difference. Those are not the
same quantity, so the "61- to 73-fold" ratio compares a genome-swap shift with a strain fold change.

Here both directions are the same quantity:
    forward  = mean over FVB samples of log2(count on the personalised genome / count on GRCm39)
    reverse  = mean over B6  samples of log2(count on the personalised genome / count on GRCm39)
Normalisation and the expression filter follow bias_qc_param.py:31-35 exactly (CPM on the full
matrix, keep a gene when either arm has mean CPM >= 1, means taken of log2(CPM + 1)).

Usage: forward_vs_reverse.py --matdir DIR --variant-counts BED --outdir DIR
"""
import argparse
import numpy as np, pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--matdir", required=True); ap.add_argument("--variant-counts", required=True)
ap.add_argument("--outdir", required=True); ap.add_argument("--min-cpm", type=float, default=1.0)
a = ap.parse_args()
import os; os.makedirs(a.outdir, exist_ok=True)

load = lambda n: pd.read_csv(f"{a.matdir}/{n}_counts.tsv", sep="\t", index_col=0)
U, A, B = load("uncorrected"), load("condA"), load("condB")
GROUPS = {"E13.5": {"FVB": ["FVB135-1", "FVB135-2", "FVB135-3"], "B6": ["B6135-1", "B6135-2", "B6135-3"]},
          "P1":    {"FVB": ["FVBP1-1", "FVBP1-2", "FVBP1-3"],    "B6": ["B6P1-1", "B6P1-2", "B6P1-3"]}}
cpm = lambda df: df.div(df.sum(axis=0), axis=1) * 1e6

def swap(mat_pers, mat_ref, samples):
    """log2(personalised / reference) for one group of samples, same filter as bias_qc_param.py."""
    cp, cr = cpm(mat_pers)[samples], cpm(mat_ref)[samples]
    keep = (cp.mean(axis=1) >= a.min_cpm) | (cr.mean(axis=1) >= a.min_cpm)
    return np.log2(cp[keep] + 1).mean(axis=1) - np.log2(cr[keep] + 1).mean(axis=1)

def strain_fc(mat_f, mat_b, fs, bs):
    cf, cb = cpm(mat_f)[fs], cpm(mat_b)[bs]
    keep = (cf.mean(axis=1) >= a.min_cpm) | (cb.mean(axis=1) >= a.min_cpm)
    return np.log2(cf[keep] + 1).mean(axis=1) - np.log2(cb[keep] + 1).mean(axis=1)

gvc = pd.read_csv(a.variant_counts, sep="\t", header=None,
                  names=["chr", "start", "end", "gene_id", "n_variants"]).set_index("gene_id")
gvc["variant_density"] = gvc["n_variants"] / ((gvc["end"] - gvc["start"]) / 1000.0)

metrics, quart, scat = [], [], []
for tp, g in GROUPS.items():
    unc = strain_fc(U, U, g["FVB"], g["B6"])
    for arm, M in [("condA", A), ("condB", B)]:
        fwd = swap(M, U, g["FVB"])
        rev = swap(M, U, g["B6"])
        d = pd.DataFrame({"forward": fwd, "reverse": rev}).join(
            pd.DataFrame({"uncorrected_strain_log2FC": unc}), how="inner").join(
            gvc[["n_variants", "variant_density"]], how="left")
        d.to_csv(f"{a.outdir}/forward_reverse_{tp}_{arm}.tsv", sep="\t")
        hi = d[d.uncorrected_strain_log2FC.abs() > 1]
        nov = d[d.n_variants.fillna(0) == 0]
        rho, prho = stats.spearmanr(d.forward, d.reverse)
        for label, sub in [("genome_wide", d), ("high_bias_absFC_gt1", hi), ("no_variant_genes", nov)]:
            metrics.append(dict(timepoint=tp, arm=arm, gene_set=label, n_genes=len(sub),
                                forward_mean_abs=sub.forward.abs().mean(), forward_mean_signed=sub.forward.mean(),
                                reverse_mean_abs=sub.reverse.abs().mean(), reverse_mean_signed=sub.reverse.mean(),
                                forward_over_reverse_abs=sub.forward.abs().mean() / sub.reverse.abs().mean()))
        metrics[-3].update(spearman_fwd_rev=rho, spearman_p=prho)
        v = d.replace([np.inf, -np.inf], np.nan).dropna(subset=["variant_density"])
        v = v[v.variant_density > 0].copy()
        v["quartile"] = pd.qcut(v.variant_density, 4, labels=["Q1_low", "Q2", "Q3", "Q4_high"], duplicates="drop")
        q = v.groupby("quartile", observed=True).apply(
            lambda x: pd.Series({"n": len(x), "density_median": x.variant_density.median(),
                                 "forward_mean_abs": x.forward.abs().mean(), "forward_mean_signed": x.forward.mean(),
                                 "reverse_mean_abs": x.reverse.abs().mean(), "reverse_mean_signed": x.reverse.mean()}),
            include_groups=False).reset_index()
        q.insert(0, "arm", arm); q.insert(0, "timepoint", tp)
        quart.append(q)
        scat.append(d.assign(timepoint=tp, arm=arm))

M = pd.DataFrame(metrics); M.to_csv(f"{a.outdir}/metrics.tsv", sep="\t", index=False)
Q = pd.concat(quart, ignore_index=True); Q.to_csv(f"{a.outdir}/quartiles.tsv", sep="\t", index=False)
pd.concat(scat).to_csv(f"{a.outdir}/per_gene_forward_reverse.tsv.gz", sep="\t", compression="gzip")
pd.set_option("display.width", 200)
print(M.round(4).to_string(index=False)); print(); print(Q.round(4).to_string(index=False))
