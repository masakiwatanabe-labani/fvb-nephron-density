#!/usr/bin/env python3
"""Quartile means of the correction amount and of the reverse control, on de-duplicated variant density.

The quantity plotted in Figure 3D and quoted in Results 2.3 is
    delta = |log2FC(condA) - log2FC(uncorrected)|
and because log2FC(arm) = FVB(arm) - B6(uncorrected) in bias_qc_param.py:43, the B6 term cancels:
    delta = |FVB(personalised) - FVB(GRCm39)|
which is exactly the "forward" genome-swap shift. The two are the same quantity; they differed only
in the gene filter (Results 2.3 restricts to |log2FC_uncorrected| > 0.1; rerun4 did not).

Both filters are reported here. Density is recomputed from unique variant positions per kb.

Usage: quartile_recompute.py --matdir DIR --variant-counts BED --outdir DIR
"""
import argparse, os
import numpy as np, pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--matdir", required=True); ap.add_argument("--variant-counts", required=True)
ap.add_argument("--outdir", required=True); ap.add_argument("--min-cpm", type=float, default=1.0)
a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)

load = lambda n: pd.read_csv(f"{a.matdir}/{n}_counts.tsv", sep="\t", index_col=0)
U, A, B = load("uncorrected"), load("condA"), load("condB")
GR = {"E13.5": (["FVB135-1","FVB135-2","FVB135-3"], ["B6135-1","B6135-2","B6135-3"]),
      "P1":    (["FVBP1-1","FVBP1-2","FVBP1-3"],    ["B6P1-1","B6P1-2","B6P1-3"])}
cpm = lambda df: df.div(df.sum(axis=0), axis=1) * 1e6

def lfc(mf, mb, fs, bs):
    cf, cb = cpm(mf)[fs], cpm(mb)[bs]
    keep = (cf.mean(axis=1) >= a.min_cpm) | (cb.mean(axis=1) >= a.min_cpm)
    return np.log2(cf[keep] + 1).mean(axis=1) - np.log2(cb[keep] + 1).mean(axis=1)

gvc = pd.read_csv(a.variant_counts, sep="\t", header=None,
                  names=["chr", "start", "end", "gene_id", "n_variants"]).set_index("gene_id")
gvc["variant_density"] = gvc["n_variants"] / ((gvc["end"] - gvc["start"]) / 1000.0)

out, novar, corr = [], [], []
for tp, (fs, bs) in GR.items():
    unc = lfc(U, U, fs, bs)
    for arm, M in [("condA", A), ("condB", B)]:
        d = pd.DataFrame({
            "correction":  lfc(M, U, fs, bs) - unc,          # = FVB(pers) - FVB(ref)
            "reverse":     lfc(M, U, bs, bs),                # B6(pers) - B6(ref)
            "uncorrected": unc,
        }).join(gvc[["n_variants", "variant_density"]], how="left")
        d = d.replace([np.inf, -np.inf], np.nan)
        nv = d[d.n_variants.fillna(0) == 0]
        novar.append(dict(timepoint=tp, arm=arm, n_genes=len(nv.dropna(subset=["correction"])),
                          correction_mean_abs=nv.correction.abs().mean(), correction_mean_signed=nv.correction.mean(),
                          reverse_mean_abs=nv.reverse.abs().mean(), reverse_mean_signed=nv.reverse.mean()))
        for fname, sub in [("manuscript_filter_absFC_gt0.1", d[d.uncorrected.abs() > 0.1]),
                           ("no_extra_filter", d)]:
            v = sub.dropna(subset=["variant_density", "correction", "reverse"])
            v = v[v.variant_density > 0].copy()
            rho_c, p_c = stats.spearmanr(v.variant_density, v.correction.abs())
            rho_r, p_r = stats.spearmanr(v.variant_density, v.reverse.abs())
            rho_fr, _ = stats.spearmanr(v.correction, v.reverse)
            corr.append(dict(timepoint=tp, arm=arm, gene_filter=fname, n_genes=len(v),
                             spearman_density_vs_correction=rho_c, p_correction=p_c,
                             spearman_density_vs_reverse=rho_r, p_reverse=p_r,
                             spearman_correction_vs_reverse=rho_fr))
            v["quartile"] = pd.qcut(v.variant_density, 4, labels=["Q1_low","Q2","Q3","Q4_high"], duplicates="drop")
            for q, g in v.groupby("quartile", observed=True):
                out.append(dict(timepoint=tp, arm=arm, gene_filter=fname, quartile=q, n=len(g),
                                density_median=g.variant_density.median(),
                                correction_mean_abs=g.correction.abs().mean(),
                                correction_mean_signed=g.correction.mean(),
                                reverse_mean_abs=g.reverse.abs().mean(),
                                reverse_mean_signed=g.reverse.mean(),
                                excess_over_reverse=g.correction.abs().mean() - g.reverse.abs().mean()))
Q = pd.DataFrame(out)
NV = pd.DataFrame(novar)
# excess over the no-variant floor, using the matching filter's floor
flo = NV.set_index(["timepoint","arm"])[["correction_mean_abs","reverse_mean_abs"]]
Q["excess_over_novariant_floor"] = [r.correction_mean_abs - flo.loc[(r.timepoint, r.arm), "correction_mean_abs"]
                                    for r in Q.itertuples()]
Q.to_csv(f"{a.outdir}/quartile_means.tsv", sep="\t", index=False)
NV.to_csv(f"{a.outdir}/no_variant_genes.tsv", sep="\t", index=False)
pd.DataFrame(corr).to_csv(f"{a.outdir}/spearman.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(Q[Q.gene_filter == "manuscript_filter_absFC_gt0.1"].round(4).to_string(index=False))
print(); print(NV.round(4).to_string(index=False))
print(); print(pd.DataFrame(corr).round(4).to_string(index=False))
