#!/usr/bin/env python3
"""Figure 7 (method caveats). Extracted from plot_omics2.py (Figure 6 block) with:
  * v2 (2026-09-23): panel A titles are arguments, because on the allele-wise variant set OSR1 is
    enriched (OR 4.10, P = 0.0037) while SIX2 and WT1 are not, so the blanket "not enriched" title
    no longer holds. Defaults reproduce the submitted wording.
  * all /tmp inputs replaced by arguments, scratch files written under --outdir/work;
  * panel C reading background/candidate rates and the binomial P from concordance_summary.tsv
    (B6J-fixed definition) instead of recomputing from the superseded rank_concordant column and hard-coding 10/24 and P.
Panels A and B are otherwise unchanged."""
import argparse, os, subprocess
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
ap = argparse.ArgumentParser()
ap.add_argument("--motifbreakr", required=True); ap.add_argument("--allreg-full", required=True); ap.add_argument("--allreg-bed", required=True)
ap.add_argument("--in-six2", required=True); ap.add_argument("--in-osr1", required=True); ap.add_argument("--in-wt1", required=True)
ap.add_argument("--six2-summits", required=True); ap.add_argument("--concordance-summary", required=True); ap.add_argument("--outdir", required=True)
ap.add_argument("--title-a", default="A  Motif-disruption calls are not enriched\nfor measured TF binding")
ap.add_argument("--title-b", default="B  Insensitive to window size")
ap.add_argument("--title-c", default="C  Four-strain ranking has no\ndiscriminating power")
args = ap.parse_args(); OUT = args.outdir; W = os.path.join(OUT, "work"); os.makedirs(W, exist_ok=True)
C_B6, C_FVB = "#8172B2", "#C44E52"
def style(ax): ax.spines[["top","right"]].set_visible(False); ax.tick_params(labelsize=11.2)
# ==================== Figure 6: 方法論的注意喚起 ====================
mb = pd.read_csv(args.motifbreakr, sep="\t")
allv = pd.read_csv(args.allreg_full, sep="\t")

def enrich(motif, chip_bed):
    pos = set(zip(mb[mb.geneSymbol.str.upper()==motif].seqnames,
                  mb[mb.geneSymbol.str.upper()==motif].start))
    pred = pd.Series([(c,p_) in pos for c,p_ in zip(allv.c, allv.e)])
    ip = pd.read_csv(chip_bed, sep="\t", header=None, names=["c","s","e"])
    inset = set(zip(ip.c, ip.e))
    inchip = pd.Series([(c,e) in inset for c,e in zip(allv.c, allv.e)])
    a, b = inchip[pred], inchip[~pred]
    tbl = [[a.sum(), len(a)-a.sum()], [b.sum(), len(b)-b.sum()]]
    orr, pv = fisher_exact(tbl)
    # log OR の標準誤差 (Woolf)
    with np.errstate(divide="ignore"):
        se = np.sqrt(sum(1/max(x,0.5) for row in tbl for x in row))
    lo, hi = np.exp(np.log(max(orr,1e-9))-1.96*se), np.exp(np.log(max(orr,1e-9))+1.96*se)
    return orr, lo, hi, pv, a.mean(), b.mean(), len(a)

fig, axes = plt.subplots(1, 3, figsize=(15.6, 5.8), gridspec_kw={"wspace": 0.42})

# A: forest plot
ax = axes[0]
res = []
for motif, chip in [("SIX2",args.in_six2),
                    ("OSR1",args.in_osr1),
                    ("WT1",args.in_wt1)]:
    res.append((motif,)+enrich(motif, chip))
for i,(tf,orr,lo,hi,pv,fa,fb,n) in enumerate(res):
    ax.plot([lo,hi],[i,i], color="#444", lw=1.5)
    ax.plot([orr],[i], marker="s", ms=9, color="#4C72B0", zorder=3)
    ax.text(hi*1.08, i, f"OR={orr:.2f}\np={pv:.2f}", fontsize=9.8, va="center")
ax.axvline(1, color=C_FVB, ls="--", lw=1.2)
ax.set_yticks(range(len(res))); ax.set_yticklabels([r[0] for r in res], fontsize=12.6)
ax.set_xscale("log"); ax.set_xlim(0.28, 9)
ax.set_xticks([0.5,1,2,4,8]); ax.set_xticklabels(["0.5","1","2","4","8"], fontsize=11.2)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_ylim(-0.75, len(res)-0.4)
ax.set_xlabel("Odds ratio (95% CI)", fontsize=12.6)
ax.set_title(args.title_a,
             fontsize=13.3, weight="bold")
style(ax)
ax.text(.5,-0.38,"variants in H3K27ac regions (n=5,931), split by whether\nmotifbreakR called a disruption for that TF",
        transform=ax.transAxes, fontsize=9.5, color="#555", ha="center")

# B: 窓幅感度
ax = axes[1]
wins = [100,250,500,1000,2000]
ors, pvs = [], []
for w in wins:
    subprocess.run(["bash","-c",
        f"awk 'BEGIN{{OFS=\"\\t\"}} {{s=$2-{w}; if(s<0)s=0; print $1,s,$3+{w}}}' "
        f"{args.six2_summits} | sort -k1,1 -k2,2n | bedtools merge -i - > {W}/f6w.bed && "
        f"bedtools intersect -a <(sort -k1,1 -k2,2n {args.allreg_bed}) -b {W}/f6w.bed -u > {W}/f6in.bed"],
        capture_output=True)
    o = enrich("SIX2", f"{W}/f6in.bed"); ors.append(o[0]); pvs.append(o[3])
ax.plot(range(len(wins)), ors, "o-", color="#4C72B0", lw=1.8, ms=7)
ax.axhline(1, color=C_FVB, ls="--", lw=1.2)
ax.set_xticks(range(len(wins)))
ax.set_xticklabels([f"±{w}" for w in wins], fontsize=11.2)
ax.set_xlabel("Window around ChIP-seq summit (bp)", fontsize=12.6)
ax.set_ylabel("Odds ratio (SIX2)", fontsize=12.6)
ax.set_ylim(0, 2)
ax.set_title(args.title_b, fontsize=13.3, weight="bold")
style(ax)
for i,(o,pv) in enumerate(zip(ors,pvs)):
    ha = "left" if i == 0 else ("right" if i == len(wins)-1 else "center")
    ax.text(i, o+0.14, f"p={pv:.2f}", ha=ha, fontsize=9.5, color="#555")
ax.set_ylim(0.35, 1.62)

# C: rank concordance の背景率
ax = axes[2]
cs = pd.read_csv(args.concordance_summary, sep="\t").set_index("definition").loc["B6J_fixed_reference (primary)"]
bg = cs.background_pct; cand_rate = cs.candidate_pct
ax.bar([0,1], [bg, cand_rate], color=["#999999", "#4C72B0"], alpha=.85, width=.55)
ax.axhline(bg, color="#999", ls="--", lw=1.1)
ax.set_xticks([0,1])
ax.set_xticklabels([f"All FVB private\nvariants\n(n={int(cs.background_n):,})", f"Prioritised\ncandidates\n(n={int(cs.candidate_n)})"], fontsize=11.2)
ax.set_ylabel("% concordant with phenotype rank", fontsize=12.6)
ax.set_ylim(0, 72)
ax.set_title(args.title_c, fontsize=13.3, weight="bold")
style(ax)
for x,v in zip([0,1],[bg,cand_rate]):
    ax.text(x, v+1.6, f"{v:.1f}%", ha="center", fontsize=12.6, weight="bold")
ax.text(.5,-0.40, f"binomial test vs background: p = {cs.binom_one_sided_P:.2f}\n(inbred genotypes are effectively binary; >50% pass by construction)",
        transform=ax.transAxes, ha="center", fontsize=9.5, color="#444")

fig.tight_layout(rect=[0,0.16,1,0.94])
fig.savefig(f"{OUT}/Figure7_method_caveats_B6J.png", dpi=600, bbox_inches="tight")
fig.savefig(f"{OUT}/Figure7_method_caveats_B6J.pdf", bbox_inches="tight"); plt.close(fig)
pd.DataFrame(res, columns=["TF","OR","CI_lo","CI_hi","P","frac_pred","frac_notpred","n_pred"]).to_csv(f"{OUT}/panelA_values.tsv", sep="\t", index=False)
pd.DataFrame({"window_bp": wins, "OR": ors, "P": pvs}).to_csv(f"{OUT}/panelB_values.tsv", sep="\t", index=False)
print("panel C:", round(bg, 3), round(cand_rate, 3), cs.binom_one_sided_P)
print("Figure7 saved")
