#!/bin/bash
set -uo pipefail
cd /usr/local/jupyter/FVB_B6_glo
samples="B6135-1 B6135-2 B6135-3 B6P1-1 B6P1-2 B6P1-3 FVB135-1 FVB135-2 FVB135-3 FVBP1-1 FVBP1-2 FVBP1-3"
for sample in $samples; do
  echo "[$(date +%T)] STAR mapping: $sample"
  outdir="results/phase2/star_uncorrected/${sample}/"
  mkdir -p "$outdir"
  STAR --genomeDir refs/star_grcm39 \
    --readFilesIn "results/phase0/fastp/${sample}_trimmed_R1.fastq.gz" "results/phase0/fastp/${sample}_trimmed_R2.fastq.gz" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 8 \
    --outFileNamePrefix "$outdir" \
    > "${outdir}star.log" 2>&1
  samtools index "${outdir}Aligned.sortedByCoord.out.bam"
  echo "[$(date +%T)] done: $sample"
done
echo "ALL_STAR_UNCORRECTED_DONE"
