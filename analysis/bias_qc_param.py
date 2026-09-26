#!/usr/bin/env python3
"""Mapping-bias QC across arms (uncorrected / condition A / condition B).
Re-implementation of bias_qc.py + bias_qc_final.py with paths as arguments and an optional
contig exclusion applied to ALL arms before library sizes (CPM denominators) are computed.
Logic (CPM, min_cpm=1 filter, log2(CPM+1) group means, |uncorrected|>0.1 for density analysis,
quartiles by pd.qcut) is copied unchanged from the original scripts."""
import argparse, os, json
import numpy as np, pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--matdir", required=True, help="dir with {uncorrected,condA,condB}_counts.tsv")
ap.add_argument("--gene-chr", required=True, help="TSV with gene_id and chr columns")
ap.add_argument("--variant-counts", required=True, help="BED-like: chr start end gene_id n_variants")
ap.add_argument("--exclude-chroms", default="", help="comma-separated, e.g. chrM,chrY")
ap.add_argument("--outdir", required=True)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

gchr = pd.read_csv(a.gene_chr, sep="\t", usecols=["gene_id", "chr"]).set_index("gene_id")["chr"]
excl = [c for c in a.exclude_chroms.split(",") if c]
def load(arm):
    m = pd.read_csv(f"{a.matdir}/{arm}_counts.tsv", sep="\t", index_col=0)
    if excl:
        m = m[~gchr.reindex(m.index).isin(excl).values]
    return m
U, A, Bm = load("uncorrected"), load("condA"), load("condB")
FVB135 = ["FVB135-1", "FVB135-2", "FVB135-3"]; B6135 = ["B6135-1", "B6135-2", "B6135-3"]
FVBP1 = ["FVBP1-1", "FVBP1-2", "FVBP1-3"]; B6P1 = ["B6P1-1", "B6P1-2", "B6P1-3"]
cpm = lambda df: df.div(df.sum(axis=0), axis=1) * 1e6
def log2fc(mf, mb, fs, bs, min_cpm=1.0):
    cf, cb = cpm(mf)[fs], cpm(mb)[bs]
    keep = (cf.mean(axis=1) >= min_cpm) | (cb.mean(axis=1) >= min_cpm)
    return np.log2(cf[keep] + 1).mean(axis=1) - np.log2(cb[keep] + 1).mean(axis=1)

gvc = pd.read_csv(a.variant_counts, sep="\t", header=None, names=["chr", "start", "end", "gene_id", "n_variants"])
gvc["variant_density"] = gvc["n_variants"] / ((gvc["end"] - gvc["start"]) / 1000.0)
gvc = gvc.set_index("gene_id")

summary, metrics = [], {}
for tp, fs, bs in [("E13.5", FVB135, B6135), ("P1", FVBP1, B6P1)]:
    r = {"uncorrected": log2fc(U, U, fs, bs), "condA": log2fc(A, U, fs, bs), "condB": log2fc(Bm, U, fs, bs),
         "reverse_condA": log2fc(A, U, bs, bs), "reverse_condB": log2fc(Bm, U, bs, bs)}
    for arm, s in r.items():
        summary.append(dict(timepoint=tp, arm=arm, n_genes=len(s), median_log2FC=s.median(),
                            mean_abs_log2FC=s.abs().mean(), n_abs_gt1=int((s.abs() > 1).sum())))
    df = pd.DataFrame({f"log2FC_{k}": v for k, v in r.items()})
    df["delta_uncorrected_to_condA"] = df["log2FC_condA"] - df["log2FC_uncorrected"]
    df["delta_uncorrected_to_condB"] = df["log2FC_condB"] - df["log2FC_uncorrected"]
    df.to_csv(f"{a.outdir}/log2FC_{tp}.tsv", sep="\t")
    df = df.join(gvc[["n_variants", "variant_density"]], how="left")
    m = {"n_genes_table": len(df)}
    high = df[df["log2FC_uncorrected"].abs() > 1]
    m["high_bias_n"] = len(high)
    for arm in ["uncorrected", "condA", "condB"]:
        m[f"high_bias_mean_abs_{arm}"] = high[f"log2FC_{arm}"].abs().mean()
        m[f"genomewide_mean_abs_{arm}"] = df[f"log2FC_{arm}"].abs().mean()
        m[f"genomewide_median_{arm}"] = df[f"log2FC_{arm}"].median()
    ab = df[["log2FC_condA", "log2FC_condB"]].dropna()
    m["r_condA_condB"] = np.corrcoef(ab.iloc[:, 0], ab.iloc[:, 1])[0, 1]; m["n_condA_condB"] = len(ab)
    m["reverse_mean_abs_condA"] = df["log2FC_reverse_condA"].abs().mean()
    m["reverse_mean_abs_condB"] = df["log2FC_reverse_condB"].abs().mean()
    v = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["variant_density", "log2FC_uncorrected", "delta_uncorrected_to_condA", "delta_uncorrected_to_condB"])
    ve = v[v["log2FC_uncorrected"].abs() > 0.1].copy()
    m["density_n"] = len(ve)
    for arm in ["condA", "condB"]:
        rho, p = stats.spearmanr(ve["variant_density"], ve[f"delta_uncorrected_to_{arm}"].abs())
        m[f"spearman_rho_{arm}"] = rho; m[f"spearman_p_{arm}"] = p
    ve["density_quartile"] = pd.qcut(ve["variant_density"], 4, labels=["Q1_low", "Q2", "Q3", "Q4_high"], duplicates="drop")
    q = ve.groupby("density_quartile", observed=True)[["delta_uncorrected_to_condA", "delta_uncorrected_to_condB"]].apply(lambda x: x.abs().mean())
    for qq in q.index:
        m[f"quartile_{qq}_condA"] = q.loc[qq, "delta_uncorrected_to_condA"]; m[f"quartile_{qq}_condB"] = q.loc[qq, "delta_uncorrected_to_condB"]
    metrics[tp] = m
    df.to_csv(f"{a.outdir}/log2FC_{tp}_with_density.tsv", sep="\t")
pd.DataFrame(summary).to_csv(f"{a.outdir}/summary_by_arm.tsv", sep="\t", index=False)
pd.DataFrame(metrics).to_csv(f"{a.outdir}/bias_metrics.tsv", sep="\t")
print(pd.DataFrame(metrics).to_string())
print("excluded contigs:", excl, "| genes per arm after exclusion:", len(U), len(A), len(Bm))
