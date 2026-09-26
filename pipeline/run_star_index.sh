#!/bin/bash
set -euo pipefail
cd /usr/local/jupyter/FVB_B6_glo
STAR --runMode genomeGenerate \
  --runThreadN 8 \
  --genomeDir refs/star_grcm39 \
  --genomeFastaFiles refs/GRCm39.primary_assembly.genome.fa \
  --sjdbGTFfile refs/gencode.vM39.annotation.gtf \
  --sjdbOverhang 149 \
  --genomeSAindexNbases 14
echo "STAR_INDEX_DONE"
