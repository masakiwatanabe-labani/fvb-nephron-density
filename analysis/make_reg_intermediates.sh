#!/bin/bash
# Rebuild the regulatory-variant intermediates from a given pair of private-variant VCFs.
# Steps 3-7 of rerun2/E3/scripts/make_intermediates.sh, with the private VCFs as arguments so the
# allele-wise set can be run through the identical code path as the original set.
# Usage: make_reg_intermediates.sh <root> <priv_snps.chr.vcf.gz> <priv_indels.chr.vcf.gz> <outdir> <genes.sorted.bed>
set -euo pipefail
ROOT=$1; PSNP=$2; PIND=$3; O=$4; GENES=$5
mkdir -p "$O/build_check"

# 3. private variants as BED (start = POS-1, end = POS-1+len(REF)), id = CHROM:POS:REF:ALT
for v in snps indels; do
  src=$([ "$v" = snps ] && echo "$PSNP" || echo "$PIND")
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$src" \
    | awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2-1+length($3),$1":"$2":"$3":"$4}' > "$O/build_check/priv_${v}_chr.bed"
  sort -k1,1 -k2,2n "$O/build_check/priv_${v}_chr.bed" > "$O/build_check/priv_${v}_chr.sorted.bed"
done

# 4. merged H3K27ac regions (E13.5 + P0) on GRCm39, private variants inside them, nearest gene
cat "$ROOT/refs/atac_gse124804_grcm39/GSM5993850_PO_H3K27ac_peaks.grcm39.bed" \
    "$ROOT/refs/atac_gse124804_grcm39/GSM5993854_E13.5_H3K27ac_peaks.grcm39.bed" \
  | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > "$O/build_check/kidney_h3k27ac_grcm39.sorted.bed"
for v in snps indels; do
  bedtools intersect -u -a "$O/build_check/priv_${v}_chr.sorted.bed" -b "$O/build_check/kidney_h3k27ac_grcm39.sorted.bed" \
    > "$O/build_check/reg_${v}_v2.sorted.bed"
  bedtools closest -a "$O/build_check/reg_${v}_v2.sorted.bed" -b "$GENES" -d > "$O/reg_${v}_nearest_gene_v2.bed"
done
cut -f1-3 "$O/build_check/reg_snps_v2.sorted.bed" > "$O/build_check/allreg.bed"

# 5. regulatory SNVs inside TF ChIP summit +-250 bp (Figure 7A)
for tf in Six2 Osr1_BF Wt1; do
  awk 'BEGIN{OFS="\t"} {s=$2-250; if(s<0)s=0; print $1,s,$3+250}' "$ROOT/refs/six2_chipseq/$tf.grcm39.summit.bed" \
    | sort -k1,1 -k2,2n | bedtools merge -i - > "$O/build_check/peaks_$tf.bed"
  bedtools intersect -u -a "$O/build_check/allreg.bed" -b "$O/build_check/peaks_$tf.bed" > "$O/build_check/in_$tf.bed"
done

# 6. private variant positions (SNVs + indels)
cat "$O/build_check/priv_snps_chr.bed" "$O/build_check/priv_indels_chr.bed" | cut -f1-3 | sort -k1,1 -k2,2n \
  > "$O/build_check/priv_all_sorted.bed"

echo "private SNVs            : $(wc -l < "$O/build_check/priv_snps_chr.bed")"
echo "private indels          : $(wc -l < "$O/build_check/priv_indels_chr.bed")"
echo "regulatory SNVs (H3K27ac): $(wc -l < "$O/build_check/reg_snps_v2.sorted.bed")"
echo "regulatory indels        : $(wc -l < "$O/build_check/reg_indels_v2.sorted.bed")"
echo REG_INTERMEDIATES_DONE
