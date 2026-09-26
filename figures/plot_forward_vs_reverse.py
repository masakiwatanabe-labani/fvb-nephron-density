#!/usr/bin/env python3
"""Forward vs reverse genome-swap shift, per gene, coloured by variant density quartile.
Usage: plot_forward_vs_reverse.py --dir DIR --out-prefix P"""
import argparse, glob, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
ap=argparse.ArgumentParser(); ap.add_argument("--dir",required=True); ap.add_argument("--out-prefix",required=True)
a=ap.parse_args()
combos=[("E13.5","condA"),("E13.5","condB"),("P1","condA"),("P1","condB")]
fig,axes=plt.subplots(1,4,figsize=(19.2,5.0),gridspec_kw={"wspace":.3})
COL=["#C8C8C8","#9DB7D4","#4C72B0","#C44E52"]
for ax,(tp,arm) in zip(axes,combos):
    d=pd.read_csv(f"{a.dir}/forward_reverse_{tp}_{arm}.tsv",sep="\t",index_col=0)
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=["forward","reverse","variant_density"])
    d=d[d.variant_density>0].copy()
    d["q"]=pd.qcut(d.variant_density,4,labels=["Q1 low","Q2","Q3","Q4 high"],duplicates="drop")
    for (lab,sub),c in zip(d.groupby("q",observed=True),COL):
        ax.scatter(sub.forward,sub.reverse,s=4,c=c,alpha=.5,linewidths=0,label=f"{lab} (n={len(sub)})")
    lim=np.percentile(np.abs(np.r_[d.forward,d.reverse]),99.5)
    ax.plot([-lim,lim],[lim,-lim],"k--",lw=.8,alpha=.5)
    ax.axhline(0,color="k",lw=.6); ax.axvline(0,color="k",lw=.6)
    ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
    rho,_=stats.spearmanr(d.forward,d.reverse)
    ax.set_xlabel("forward: FVB samples,\n$\\log_2$(personalised / GRCm39)",fontsize=10.5)
    ax.set_ylabel("reverse: B6 samples,\n$\\log_2$(personalised / GRCm39)",fontsize=10.5)
    ax.set_title(f"{tp}, condition {arm[-1]}",fontsize=12.6,weight="bold")
    ax.text(.03,.03,f"Spearman $\\rho$ = {rho:+.2f}\nmean |fwd| {d.forward.abs().mean():.3f}\n"
                    f"mean |rev| {d.reverse.abs().mean():.3f}",transform=ax.transAxes,fontsize=9,va="bottom")
    ax.spines[["top","right"]].set_visible(False)
    if ax is axes[0]: ax.legend(fontsize=8,frameon=False,loc="upper right",markerscale=2.5)
fig.tight_layout()
for e in ("png","pdf"): fig.savefig(f"{a.out_prefix}.{e}",dpi=300,bbox_inches="tight")
print("saved",a.out_prefix)
