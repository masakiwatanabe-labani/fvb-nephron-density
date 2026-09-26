#!/bin/bash
set -euo pipefail
cd /usr/local/jupyter/FVB_B6_glo
STAR --runMode genomeGenerate \
  --runThreadN 6 \
  --genomeDir refs/star_condA \
  --genomeFastaFiles results/phase2/g2g/condA/FVB_snponly.fa \
  --sjdbGTFfile results/phase2/g2g/condA/FVB_snponly.gtf \
  --sjdbOverhang 149 \
  --genomeSAindexNbases 14
echo "STAR_CONDA_INDEX_DONE"
