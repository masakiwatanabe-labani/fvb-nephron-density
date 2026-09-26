import pandas as pd

BASE = "/usr/local/jupyter/FVB_B6_glo"
master = pd.read_csv(f"{BASE}/results/phase4/candidate_genes_raw_v3.tsv", sep="\t", index_col=0)

# Original developmental kidney-axis gene set (from the task's own priority list)
KIDNEY_AXIS_GENES = {
    "Ret","Gdnf","Gfra1","Six2","Eya1","Sall1","Wt1","Wnt9b","Wnt11","Bmp7","Grem1",
    "Itga8","Fgf20","Osr1","Hnf1b","Pax2","Pax8","Foxd1","Tcf21","Robo2","Slit2",
    "Etv4","Etv5","Spry1","Cited1","Meox1","Gata3","Lhx1",
}

# CAKUT (Congenital Anomalies of Kidney and Urinary Tract) gene panel, PanelApp Genomics England
# v1.182, "Green" (highest-confidence) genes only, human symbols converted to mouse case convention.
# Source: https://panelapp.genomicsengland.co.uk/panels/234/
# NOTE: this panel is curated for human clinical genetics and is no longer actively maintained by
# PanelApp; used here as a proxy "kidney disease susceptibility gene" set because no mouse strain-level
# QTL for nephron/glomerular number could be found with usable genomic coordinates (see chat log:
# BXD/GeneNetwork, MGI, Mouse Phenome Database, and the FVB congenic Ckdp1/Ckdp2 locus were all checked;
# Ckdp1/Ckdp2 in particular has no published Mb-level coordinates, only a 17-cM interval, so it could
# not be used for a region-restricted analysis).
CAKUT_GENES = {
    "Ace","Actg2","Agt","Agtr1","Anos1","Bnc2","Cep55","Chd7","Chrm3","Chrna3","Ctu2",
    "Dhcr7","Dstyk","Eya1","Fam58a","Fras1","Frem1","Frem2","Gata3","Gli3","Gpc3",
    "Greb1l","Grip1","Haao","Hnf1b","Hoxa13","Hpse2","Itga8","Jag1","Kdm6a","Kmt2d",
    "Kynu","Lifr","Lrig2","Lrp4","Myocd","Nadsyn1","Nipbl","Notch2","Nphp3","Pax2",
    "Plvap","Ren","Ret","Robo1","Sall1","Stra6","Tbc1d1","Tbx18","Tfap2a","Tmem260",
    "Trap1","Wbp11","Zic3","Zmym2",
}

SUSCEPTIBILITY_GENES = KIDNEY_AXIS_GENES | CAKUT_GENES
print(f"Combined susceptibility gene set: {len(SUSCEPTIBILITY_GENES)} genes "
      f"({len(KIDNEY_AXIS_GENES)} developmental-axis + {len(CAKUT_GENES)} CAKUT, "
      f"{len(KIDNEY_AXIS_GENES & CAKUT_GENES)} overlap)")

# relaxed expression-difference threshold (per user approval): FDR<0.1, |log2FC|>0.5 (was >1)
LOG2FC_THRESHOLD = 0.5

def has_expr_diff(row):
    e135 = (pd.notna(row["padj_E135"]) and row["padj_E135"] < 0.1 and
            pd.notna(row["log2FC_E135"]) and abs(row["log2FC_E135"]) > LOG2FC_THRESHOLD)
    p1 = (pd.notna(row["padj_P1"]) and row["padj_P1"] < 0.1 and
          pd.notna(row["log2FC_P1"]) and abs(row["log2FC_P1"]) > LOG2FC_THRESHOLD)
    return e135 or p1

master["is_susceptibility_gene"] = master["gene_name"].isin(SUSCEPTIBILITY_GENES)
master["is_kidney_axis_gene"] = master["gene_name"].isin(KIDNEY_AXIS_GENES)
master["is_cakut_gene"] = master["gene_name"].isin(CAKUT_GENES)
master["has_expr_diff"] = master.apply(has_expr_diff, axis=1)
master["has_regulatory_and_motif"] = master["regulatory_overlap"] & master["motif_disrupted"].notna()

def assign_tier(row):
    if row["is_susceptibility_gene"] and row["has_expr_diff"] and row["has_regulatory_and_motif"]:
        return 1
    if row["is_susceptibility_gene"] and row["has_expr_diff"]:
        return 2
    return 3

master["tier"] = master.apply(assign_tier, axis=1)

n_tier1 = (master["tier"] == 1).sum()
n_tier2 = (master["tier"] == 2).sum()
n_tier3 = (master["tier"] == 3).sum()
print(f"\nTier 1: {n_tier1} genes")
print(f"Tier 2: {n_tier2} genes")
print(f"Tier 3: {n_tier3} genes")

out_cols = [
    "gene_name", "n_private_snps", "n_private_indels", "max_vep_impact",
    "regulatory_overlap", "motif_disrupted",
    "log2FC_E135", "padj_E135", "log2FC_P1", "padj_P1",
    "rank_concordant", "is_kidney_axis_gene", "is_cakut_gene", "tier",
]
out = master[out_cols].copy()
out.index.name = "gene_id"
out["_has_motif"] = master["motif_disrupted"].notna()
out["_n_var"] = master["n_private_snps"] + master["n_private_indels"]
out = out.sort_values(["tier", "_has_motif", "_n_var"], ascending=[True, False, False])
out = out.drop(columns=["_has_motif", "_n_var"])
out.to_csv(f"{BASE}/results/phase4/candidate_genes_v3_ckd_genes.tsv", sep="\t")
print(f"\nWrote candidate_genes_v3_ckd_genes.tsv: {len(out)} total candidate genes "
      f"(threshold: FDR<0.1, |log2FC|>{LOG2FC_THRESHOLD})")

print("\n=== Tier 1 genes ===")
print(out[out["tier"] == 1].to_string())

print("\n=== Tier 2 genes ===")
print(out[out["tier"] == 2][["gene_name","log2FC_E135","padj_E135","log2FC_P1","padj_P1",
                              "n_private_snps","n_private_indels","is_kidney_axis_gene","is_cakut_gene"]].to_string())
