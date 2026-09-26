#!/usr/bin/env python3
"""Figure 5 (cell-composition deconvolution). Extracted from plot_omics.py (Figure 4 block) with all
paths as arguments and the plotted values written out, so the figure no longer depends on a fixed
project layout. Panel content is unchanged.
Usage: plot_figure5.py --deconv-dir DIR --out-prefix P
"""
import argparse
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--deconv-dir", required=True); ap.add_argument("--out-prefix", required=True)
a = ap.parse_args()
C_B6, C_FVB = "#8172B2", "#C44E52"
def style(ax): ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=11.2)

CTS = ["Cap_Mesenchyme", "Nephron_Progenitor", "Ureteric_Bud", "Stromal", "Podocytes",
       "Proximal_Tubule", "Distal_Tubule", "Loop_of_Henle", "Endothelial"]
fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.6))
rows = []

ax = axes[0]
v = pd.read_csv(f"{a.deconv_dir}/pseudobulk_validation.csv", index_col=0)
T = v.iloc[:, :9].values.ravel(); E = v.iloc[:, 9:].values.ravel()
r = float(np.corrcoef(T, E)[0, 1]); mae = float(np.abs(T - E).mean())
ax.scatter(T, E, s=34, color="#4C72B0", alpha=.7, edgecolor="white", linewidth=.5)
ax.plot([0, .5], [0, .5], "k--", lw=.9, alpha=.6)
ax.set_xlabel("True proportion", fontsize=12.6); ax.set_ylabel("Estimated proportion", fontsize=12.6)
ax.set_title("A  Pseudo-bulk validation", fontsize=14, weight="bold"); style(ax)
ax.text(.04, .95, f"r = {r:.3f}\nMAE = {mae:.3f}", transform=ax.transAxes, va="top", fontsize=11.2,
        bbox=dict(boxstyle="round,pad=.3", fc="#f5f5f5", ec="#bbb", lw=.5))
rows += [{"panel": "A", "item": "pseudo-bulk r", "value": round(r, 4)},
         {"panel": "A", "item": "pseudo-bulk MAE", "value": round(mae, 4)}]

for j, tp in enumerate(["E13.5", "P1"]):
    ax = axes[j + 1]
    p = pd.read_csv(f"{a.deconv_dir}/proportions_v2_{tp}.csv")
    present = [c for c in CTS if p[c].sum() > 0.001]
    x = np.arange(len(present)); w = 0.36
    for i, (st, c) in enumerate([("B6", C_B6), ("FVB", C_FVB)]):
        d = p[p.strain == st][present]
        ax.bar(x + (i - 0.5) * w, d.mean(), w, yerr=d.std(), color=c, alpha=.85,
               error_kw=dict(lw=.9, capsize=2.5), label=st)
        for ct in present:
            rows.append({"panel": "BC"[j], "item": f"{tp} {st} {ct}", "value": round(float(d[ct].mean()), 4)})
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("_", " ") for c in present], fontsize=10.1, rotation=38, ha="right")
    ax.set_ylabel("Estimated proportion", fontsize=12.6)
    ax.set_title(f"{'B' if j == 0 else 'C'}  Cell composition: {tp}", fontsize=14, weight="bold")
    style(ax)
    npc = p["Cap_Mesenchyme"] + p["Nephron_Progenitor"]
    b = npc[p.strain == "B6"]; f_ = npc[p.strain == "FVB"]
    pv = float(stats.ttest_ind(f_, b).pvalue)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.42)
    ax.legend(fontsize=11.2, frameon=False, loc="upper left")
    ax.text(.97, .97, f"NPC lineage (Cap Mes. + Neph. Prog.)\n"
                      f"B6 {b.mean():.3f} vs FVB {f_.mean():.3f},  p = {pv:.2f} (n.s.)",
            transform=ax.transAxes, ha="right", va="top", fontsize=9.8,
            bbox=dict(boxstyle="round,pad=.3", fc="#f5f5f5", ec="#bbb", lw=.5))
    rows += [{"panel": "BC"[j], "item": f"{tp} NPC lineage B6", "value": round(float(b.mean()), 4)},
             {"panel": "BC"[j], "item": f"{tp} NPC lineage FVB", "value": round(float(f_.mean()), 4)},
             {"panel": "BC"[j], "item": f"{tp} NPC lineage p", "value": round(pv, 4)}]

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{a.out_prefix}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
pd.DataFrame(rows).to_csv(f"{a.out_prefix}_plotted_values.tsv", sep="\t", index=False)
print(f"saved {a.out_prefix}")
