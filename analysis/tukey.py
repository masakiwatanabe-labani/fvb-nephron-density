"""Tukey HSD for a one-way layout, implemented on scipy's studentized range.
statsmodels is not installed in the current environment; analyze_phenotype.py:81 used
statsmodels.stats.multicomp.pairwise_tukeyhsd. This reproduces that computation."""
import itertools
import numpy as np
from scipy.stats import studentized_range

def tukey_hsd(values, groups):
    values = np.asarray(values, float); groups = np.asarray(groups)
    lv = sorted(set(groups)); k = len(lv)
    n = {g: int((groups == g).sum()) for g in lv}
    m = {g: values[groups == g].mean() for g in lv}
    N = len(values)
    sse = sum(((values[groups == g] - m[g]) ** 2).sum() for g in lv)
    dfe = N - k; mse = sse / dfe
    out = {}
    for a, b in itertools.combinations(lv, 2):
        se = np.sqrt(mse / 2 * (1 / n[a] + 1 / n[b]))
        q = abs(m[a] - m[b]) / se
        out[(a, b)] = float(studentized_range.sf(q, k, dfe))
    return out
