#!/bin/bash
# Regenerate the intermediate files that the analysis scripts read from /tmp, under <outdir> instead.
# Every command is reconstructed from the file contents and verified against the 2026-09-15 snapshot of /tmp
# (see verify_against_snapshot.tsv). Usage: make_intermediates.sh <project_root> <outdir>
set -euo pipefail
ROOT=$1; O=$2; mkdir -p "$O/build_check"
GTF="$ROOT/refs/gencode.vM39.annotation.gtf"
# 1. gene annotation (GTF order), gene_id -> gene_name, and sorted gene BED (0-based start)
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" {match($9,/gene_id "[^"]+"/); id=substr($9,RSTART+9,RLENGTH-10);
  match($9,/gene_name "[^"]+"/); nm=substr($9,RSTART+11,RLENGTH-12); print $1,$4-1,$5,id > "'"$O"'/genes.bed"; print id,nm > "'"$O"'/gene_id2name.tsv"}' "$GTF"
sort -k1,1 -k2,2n "$O/genes.bed" > "$O/genes.sorted.bed"
# 2. all FVB_NJ variants (1-bp BED at POS) and per-gene counts -> variant density for the bias QC
for v in snps indels; do bcftools query -f '%CHROM\t%POS\n' "$ROOT/results/phase1/variants_chr/FVB.$v.chr.vcf.gz"; done \
  | awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2}' > "$O/fvb_variants.bed"
sort -k1,1 -k2,2n "$O/fvb_variants.bed" | bedtools intersect -sorted -a "$O/genes.sorted.bed" -b - -c > "$O/gene_variant_counts.bed"
# 3. private variants as BED (start = POS-1, end = POS-1+len(REF)), id = CHROM:POS:REF:ALT
for v in snps indels; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$ROOT/results/phase1/variants_chr/private_$v.chr.vcf.gz" \
    | awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2-1+length($3),$1":"$2":"$3":"$4}' > "$O/build_check/priv_${v}_chr.bed"
  sort -k1,1 -k2,2n "$O/build_check/priv_${v}_chr.bed" > "$O/build_check/priv_${v}_chr.sorted.bed"
done
# 4. merged H3K27ac regions (E13.5 + P0) on GRCm39, and private variants inside them
cat "$ROOT/refs/atac_gse124804_grcm39/GSM5993850_PO_H3K27ac_peaks.grcm39.bed" "$ROOT/refs/atac_gse124804_grcm39/GSM5993854_E13.5_H3K27ac_peaks.grcm39.bed" \
  | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > "$O/build_check/kidney_h3k27ac_grcm39.sorted.bed"
for v in snps indels; do
  bedtools intersect -u -a "$O/build_check/priv_${v}_chr.sorted.bed" -b "$O/build_check/kidney_h3k27ac_grcm39.sorted.bed" > "$O/build_check/reg_${v}_v2.sorted.bed"
  bedtools closest -a "$O/build_check/reg_${v}_v2.sorted.bed" -b "$O/genes.sorted.bed" -d > "$O/reg_${v}_nearest_gene_v2.bed"
done
cut -f1-3 "$O/build_check/reg_snps_v2.sorted.bed" > "$O/build_check/allreg.bed"
# 5. regulatory SNVs inside TF ChIP summit +-250 bp (used by Figure 7A)
for tf in Six2 Osr1_BF Wt1; do
  awk 'BEGIN{OFS="\t"} {s=$2-250; if(s<0)s=0; print $1,s,$3+250}' "$ROOT/refs/six2_chipseq/$tf.grcm39.summit.bed" | sort -k1,1 -k2,2n | bedtools merge -i - > "$O/build_check/peaks_$tf.bed"
  bedtools intersect -u -a "$O/build_check/allreg.bed" -b "$O/build_check/peaks_$tf.bed" > "$O/build_check/in_$tf.bed"
done
# 6. private variant positions (SNVs + indels)
cat "$O/build_check/priv_snps_chr.bed" "$O/build_check/priv_indels_chr.bed" | cut -f1-3 | sort -k1,1 -k2,2n > "$O/build_check/priv_all_sorted.bed"
# 7. allreg_full.tsv (coordinates + id; the six2_pred column is added by make_allreg_full.py)
python3 "$(dirname "$0")/make_allreg_full.py" --reg "$O/build_check/reg_snps_v2.sorted.bed" \
  --motifbreakr "$ROOT/results/phase1/regulatory_v3/motifbreakr_results.tsv" --out "$O/build_check/allreg_full.tsv"
echo INTERMEDIATES_DONE
