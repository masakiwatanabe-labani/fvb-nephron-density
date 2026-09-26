#!/bin/bash
set -uo pipefail
source /root/miniforge3/etc/profile.d/conda.sh
conda activate fvb_pipeline
mkdir -p results/phase0/fastp
tail -n +2 samplesheet.tsv | while IFS=$'\t' read -r sample strain timepoint sex r1 r2; do
  echo "[$(date +%T)] fastp: $sample"
  fastp \
    -i "$r1" -I "$r2" \
    -o "results/phase0/fastp/${sample}_trimmed_R1.fastq.gz" \
    -O "results/phase0/fastp/${sample}_trimmed_R2.fastq.gz" \
    -j "results/phase0/fastp/${sample}.fastp.json" \
    -h "results/phase0/fastp/${sample}.fastp.html" \
    -w 4 \
    2> "results/phase0/fastp/${sample}.fastp.log"
done
echo "ALL_FASTP_DONE"
