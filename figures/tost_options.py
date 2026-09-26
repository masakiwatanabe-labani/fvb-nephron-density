#!/usr/bin/env python3
"""TOST for FVB/N vs BALB/c males (kidney weight; kidney/body weight ratio), margin +-20% of the BALB/c male mean.
Three variants of the standard error / degrees of freedom:
  current : unpooled SE with df = n1+n2-2 (as implemented in plot_phenotype.py:97-106; internally inconsistent)
  welch   : unpooled SE with Welch-Satterthwaite df            (option 1)
  pooled  : pooled-variance SE with df = n1+n2-2 (Student)      (option 2)
Outputs the difference, margins in native units, both one-sided P, TOST P, 90% and 95% CI."""
import argparse
import numpy as np, pandas as pd
from scipy import stats

def tost(a, b, mode, margin_frac=0.20):
    a, b = np.asarray(a, float), np.asarray(b, float); na, nb = len(a), len(b)
    d = a.mean() - b.mean(); va, vb = a.var(ddof=1), b.var(ddof=1)
    if mode == "pooled":
        sp2 = ((na - 1) * va + (nb - 1) * vb) / (na + nb - 2); se = np.sqrt(sp2 * (1 / na + 1 / nb)); df = na + nb - 2
    else:
        se = np.sqrt(va / na + vb / nb)
        df = na + nb - 2 if mode == "current" else (va / na + vb / nb) ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    bound = margin_frac * b.mean()
    p_lo = 1 - stats.t.cdf((d + bound) / se, df); p_hi = stats.t.cdf((d - bound) / se, df)
    t90, t95 = stats.t.ppf(0.95, df), stats.t.ppf(0.975, df)
    return dict(mode=mode, n_FVB=na, n_BALB=nb, mean_FVB=a.mean(), mean_BALB=b.mean(), diff=d, SE=se, df=df, margin=bound,
                p_lower=p_lo, p_upper=p_hi, TOST_P=max(p_lo, p_hi), CI90_lo=d - t90 * se, CI90_hi=d + t90 * se,
                CI95_lo=d - t95 * se, CI95_hi=d + t95 * se, equivalent_at_alpha05=bool(max(p_lo, p_hi) < 0.05),
                CI90_within_margin=bool(-bound < d - t90 * se and d + t90 * se < bound))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--data", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    df = pd.read_csv(a.data); df["kw_bw_pct"] = df.kidney_weight_g / df.body_weight_g * 100
    m = df[df.sex == "M"]; rows = []
    for label, col, unit in [("kidney_weight", "kidney_weight_g", "g"), ("kw_bw_ratio", "kw_bw_pct", "percentage points")]:
        for mode in ["current", "welch", "pooled"]:
            r = tost(m[m.strain == "FVB"][col], m[m.strain == "BALB"][col], mode); r.update(measure=label, unit=unit); rows.append(r)
    out = pd.DataFrame(rows)[["measure", "unit", "mode", "n_FVB", "n_BALB", "mean_FVB", "mean_BALB", "diff", "SE", "df", "margin", "p_lower", "p_upper",
                              "TOST_P", "CI90_lo", "CI90_hi", "CI95_lo", "CI95_hi", "equivalent_at_alpha05", "CI90_within_margin"]]
    out.to_csv(a.out, sep="\t", index=False); pd.set_option("display.width", 250); print(out.to_string(index=False))
