#!/usr/bin/env python3
"""Figure 3F: forward and reverse genome-swap shifts as the same quantity on the same genes.

Results 2.3 previously compared the reverse-control shift (B6 reads counted on the personalised
genome vs on GRCm39) against 1.99/1.88, the mean |log2 fold change| between FVB/N and C57BL/6J in
genes selected for a large strain difference. Those are not the same quantity, so that comparison
is withdrawn. Here both directions are the same quantity:

    forward = mean over FVB samples of log2(count on the personalised genome / count on GRCm39)
    reverse = mean over B6  samples of log2(count on the personalised genome / count on GRCm39)

Normalisation, the expression filter and the |log2FC_uncorrected| > 0.1 gene filter follow
rerun5/density/quartile_recompute.py:24-47 exactly (CPM on the full matrix, keep a gene when
either arm has mean CPM >= 1, means of log2(CPM + 1)). The gene set is intersected between
condition A and condition B, so the two arms are summarised over one common set: 7,028 genes at
E13.5 and 10,720 at P1, which is what the figure states. quartile_recompute.py reports the same
quantity per variant-density quartile and rounds its summary to four decimals; this script keeps
full precision so that the forward/reverse ratio can be quoted without rounding.

The signed-shift Spearman rho and the no-variant floor are taken from that script's own outputs
(spearman.tsv and no_variant_genes.tsv), not recomputed.

Usage: figure3F_matched_shifts.py --matdir DIR --spearman F --no-variant F --out F
"""
import argparse
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--matdir", required=True, help="rerun3/S2/out/mat")
ap.add_argument("--spearman", required=True, help="rerun5/density/spearman.tsv")
ap.add_argument("--no-variant", required=True, help="rerun5/density/no_variant_genes.tsv")
ap.add_argument("--min-cpm", type=float, default=1.0)
ap.add_argument("--out", required=True)
a = ap.parse_args()

load = lambda n: pd.read_csv(f"{a.matdir}/{n}_counts.tsv", sep="\t", index_col=0)
U, A, B = load("uncorrected"), load("condA"), load("condB")
GR = {"E13.5": (["FVB135-1", "FVB135-2", "FVB135-3"], ["B6135-1", "B6135-2", "B6135-3"]),
      "P1":    (["FVBP1-1", "FVBP1-2", "FVBP1-3"],    ["B6P1-1", "B6P1-2", "B6P1-3"])}
cpm = lambda df: df.div(df.sum(axis=0), axis=1) * 1e6

def lfc(mf, mb, fs, bs):
    cf, cb = cpm(mf)[fs], cpm(mb)[bs]
    keep = (cf.mean(axis=1) >= a.min_cpm) | (cb.mean(axis=1) >= a.min_cpm)
    return np.log2(cf[keep] + 1).mean(axis=1) - np.log2(cb[keep] + 1).mean(axis=1)

rows = []
for tp, (fs, bs) in GR.items():
    unc = lfc(U, U, fs, bs)
    per = {}
    for arm, M in [("condA", A), ("condB", B)]:
        per[arm] = pd.DataFrame({"correction": lfc(M, U, fs, bs) - unc,
                                 "reverse": lfc(M, U, bs, bs),
                                 "uncorrected": unc}).replace([np.inf, -np.inf], np.nan)
    keep = lambda d: d[d.uncorrected.abs() > 0.1].dropna(subset=["correction", "reverse"]).index
    common = keep(per["condA"]).intersection(keep(per["condB"]))
    for arm in ("condA", "condB"):
        s = per[arm].loc[common]
        f, r = s.correction.abs().mean(), s.reverse.abs().mean()
        rows.append(dict(timepoint=tp, arm=arm, n_genes=len(common),
                         forward=f, reverse=r, forward_over_reverse=f / r))

t = pd.DataFrame(rows)
sp = (pd.read_csv(a.spearman, sep="\t").query("density == 'dup'")
        [["timepoint", "arm", "spearman_correction_vs_reverse"]]
        .rename(columns={"spearman_correction_vs_reverse": "spearman_signed"}))
nv = (pd.read_csv(a.no_variant, sep="\t").query("density == 'dup'")
        [["timepoint", "arm", "n", "correction_mean_abs", "reverse_mean_abs"]]
        .rename(columns={"n": "novariant_n", "correction_mean_abs": "novariant_forward",
                         "reverse_mean_abs": "novariant_reverse"}))
t = t.merge(sp, on=["timepoint", "arm"]).merge(nv, on=["timepoint", "arm"])
t.to_csv(a.out, sep="\t", index=False)
pd.set_option("display.float_format", lambda x: f"{x:.10f}")
print(t.to_string(index=False))
