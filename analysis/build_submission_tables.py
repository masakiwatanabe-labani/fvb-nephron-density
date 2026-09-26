"""投稿用 Table 1 / Table 3 の生成"""
import pandas as pd, numpy as np, os
from scipy import stats

BASE = "/usr/local/jupyter/FVB_B6_glo"
OUT = f"{BASE}/submission/tables"
os.makedirs(OUT, exist_ok=True)

# ---------- Table 1: 表現型 ----------
df = pd.read_csv(f"{BASE}/data/glomeruli_4strain.csv")
df["glom_per_g_kw"] = df.glomeruli_per_kidney / df.kidney_weight_g
df["glom_per_g_bw"] = df.glomeruli_per_kidney / df.body_weight_g
df["kw_bw_pct"] = df.kidney_weight_g / df.body_weight_g * 100
df["ug_per_glom"] = df.kidney_weight_g / df.glomeruli_per_kidney * 1e6

ORDER = ["BALB", "DBA", "B6", "FVB"]
FULL = {"BALB": "BALB/c", "DBA": "DBA/2", "B6": "C57BL/6J", "FVB": "FVB/N"}

rows = []
for s in ORDER:
    d = df[df.strain == s]
    def ms(col, dec=0):
        return f"{d[col].mean():,.{dec}f} ± {d[col].std():,.{dec}f}"
    rows.append({
        "Strain": FULL[s], "n": len(d),
        "Body weight (g)": ms("body_weight_g", 1),
        "Kidney weight (g)": ms("kidney_weight_g", 3),
        "Kidney/body weight (%)": ms("kw_bw_pct", 3),
        "Glomeruli per kidney": ms("glomeruli_per_kidney"),
        "Glomeruli per g kidney weight": ms("glom_per_g_kw"),
        "Glomeruli per g body weight": ms("glom_per_g_bw"),
        "Kidney weight per glomerulus (ug)": ms("ug_per_glom", 1),
    })
t1 = pd.DataFrame(rows)
t1.to_csv(f"{OUT}/Table1_phenotype_summary.csv", index=False)
print("=== Table 1 ===")
print(t1.to_string(index=False))

# 脚注用の統計
print("\n--- Table 1 脚注用 ---")
from statsmodels.formula.api import ols
import statsmodels.api as sm
aov = sm.stats.anova_lm(ols("glomeruli_per_kidney ~ C(strain)*C(sex)", data=df).fit(), typ=2)
print(f"Two-way ANOVA (glomeruli per kidney): strain F={aov.loc['C(strain)','F']:.1f}, "
      f"P={aov.loc['C(strain)','PR(>F)']:.1e}; sex F={aov.loc['C(sex)','F']:.1f}, "
      f"P={aov.loc['C(sex)','PR(>F)']:.3f}; interaction F={aov.loc['C(strain):C(sex)','F']:.2f}, "
      f"P={aov.loc['C(strain):C(sex)','PR(>F)']:.3f}")
aov2 = sm.stats.anova_lm(ols("glomeruli_per_kidney ~ C(strain)+C(sex)+kidney_weight_g", data=df).fit(), typ=2)
print(f"ANCOVA with kidney weight as covariate: sex P={aov2.loc['C(sex)','PR(>F)']:.3f}")

# ---------- Table 3: 候補遺伝子 ----------
p = pd.read_csv(f"{BASE}/results/phase4/cis_discovery_passed.tsv", sep="\t")
c = pd.read_csv(f"{BASE}/results/phase4/cis_discovery_candidates.tsv", sep="\t").set_index("gene_id")
cols = ["chr","start","end","var_id","g_chr","g_start","g_end","gene_id","dist"]
reg = pd.concat([
    pd.read_csv("/tmp/reg_snps_nearest_gene_v2.bed", sep="\t", header=None, names=cols),
    pd.read_csv("/tmp/reg_indels_nearest_gene_v2.bed", sep="\t", header=None, names=cols),
], ignore_index=True).drop_duplicates(subset=["chr","start","var_id"])

chip = {}
for tf in ["Six2","Six2_BF","Osr1","Wt1","Hoxd11"]:
    h = pd.read_csv(f"/tmp/cis_in_{tf}.bed", sep="\t", header=None, names=["chr","start","end"])
    chip[tf] = set(zip(h.chr, h.start))
reg["tfs"] = reg.apply(lambda r: ";".join(tf for tf,s in chip.items() if (r.chr,r.start) in s), axis=1)
reg_chip = reg[reg.tfs != ""]

# summit までの最短距離
summits = {}
for tf, f in [("Six2","Six2"),("Osr1","Osr1_BF"),("Wt1","Wt1"),("Hoxd11","Hoxd11_BF")]:
    d = pd.read_csv(f"{BASE}/refs/six2_chipseq/{f}.grcm39.summit.bed", sep="\t",
                    header=None, names=["c","s","e","id","score"])
    summits[tf] = d

rows = []
for _, r in p.iterrows():
    v = reg_chip[reg_chip.gene_id == r.gene_id]
    if len(v) == 0: continue
    v0 = v.iloc[0]
    dmin, tfmin = None, None
    for tf, d in summits.items():
        sub = d[(d.c == v0.chr) & (d.e.between(v0.end - 400, v0.end + 400))]
        for _, x in sub.iterrows():
            dd = abs(x.e - v0.end)
            if dmin is None or dd < dmin: dmin, tfmin = dd, tf
    cc = c.loc[r.gene_id]
    rows.append({
        "Gene": r.gene_name,
        "Ensembl ID": r.gene_id,
        "Variant": f"{v0.chr}:{v0.end}",
        "TF peaks containing variant": v0.tfs.replace("Six2_BF","Six2(BioTagFLAG)"),
        "Distance to nearest summit (bp)": dmin,
        "Nearest summit factor": tfmin,
        "log2FC E13.5": round(cc.log2FC_E135, 3) if pd.notna(cc.log2FC_E135) else None,
        "FDR E13.5": f"{cc.padj_E135:.2e}" if pd.notna(cc.padj_E135) else None,
        "log2FC P1": round(cc.log2FC_P1, 3) if pd.notna(cc.log2FC_P1) else None,
        "FDR P1": f"{cc.padj_P1:.2e}" if pd.notna(cc.padj_P1) else None,
        "Mean normalised count": round(r.baseMean, 0),
    })
t3 = pd.DataFrame(rows)
# 最小FDR順
t3["_s"] = t3[["FDR E13.5","FDR P1"]].apply(
    lambda x: min([float(v) for v in x if v is not None] or [1]), axis=1)
t3 = t3.sort_values("_s").drop(columns="_s")
t3.to_csv(f"{OUT}/Table3_prioritised_candidates.csv", index=False)
print(f"\n=== Table 3 ({len(t3)} genes) ===")
print(t3.head(6).to_string(index=False))
