#!/bin/bash
# Re-run the personalised-genome arms with the library-appropriate counting settings.
#
# Same pipeline as run_phase2_continuation.sh + run_quantify_cleanup.sh, with three deliberate changes:
#   * featureCounts gains  -s 2 --countReadPairs  (reverse-stranded, fragment-level counting).
#     salmon `-l A` auto-detected ISR for all 12 libraries; `-p` alone does not switch subread 2.x
#     from read-level to fragment-level counting.
#   * condition B reverse-conversion uses g2g_bam_convert_v2.py, which gives contigs declared in the
#     VCI header but carrying no records (chrM, chrY, 39 unplaced scaffolds) an identity mapping,
#     instead of dropping every alignment on them.
#   * STAR for the next sample runs while the current sample is converted and counted, so the
#     single-threaded conversion does not idle the machine. Only one STAR process exists at a time
#     (the 27 GB index does not fit twice in 58 GB).
# Condition A needs no coordinate conversion: an SNV-only personalised genome has reference coordinates.
#
# Usage: run_phase2_s2pairs.sh <root> <outdir> <converter_dir> <gtf> <job...>
#   <job> is SAMPLE:COND, e.g. B6135-1:B
set -uo pipefail
ROOT=$1; OUT=$2; CONV_DIR=$3; GTF=$4; shift 4
JOBS=("$@")
SAMTOOLS=${SAMTOOLS:-/usr/local/bin/samtools}
STAR_BIN=${STAR_BIN:-/usr/local/bin/STAR}
FC=${FEATURECOUNTS:-featureCounts}
PY=${PYTHON:-/usr/bin/python3}
cd "$ROOT"
mkdir -p "$OUT"/{star,reverted,counts,logs}

star_one() {                     # $1 sample, $2 cond
  local s=$1 c=$2
  local d="$OUT/star/${s}.cond${c}"
  mkdir -p "$d"
  "$STAR_BIN" --genomeDir "refs/star_cond${c}" \
    --readFilesIn "results/phase0/fastp/${s}_trimmed_R1.fastq.gz" "results/phase0/fastp/${s}_trimmed_R2.fastq.gz" \
    --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --runThreadN 8 \
    --outFileNamePrefix "$d/" > "$d/star.log" 2>&1
  "$SAMTOOLS" index "$d/Aligned.sortedByCoord.out.bam"
}

quantify_one() {                 # $1 sample, $2 cond
  local s=$1 c=$2
  local bam="$OUT/star/${s}.cond${c}/Aligned.sortedByCoord.out.bam"
  local target="$bam"
  if [ ! -s "$bam" ]; then echo "[$(date +%T)] FAIL no BAM: $s cond$c"; return 1; fi
  if [ "$c" = "B" ]; then
    local rev="$OUT/reverted/${s}.condB.reverted.bam"
    local srt="$OUT/reverted/${s}.condB.reverted.sorted.bam"
    echo "[$(date +%T)] convert(v2): $s condB"
    "$PY" "$CONV_DIR/g2g_bam_convert_v2.py" -c results/phase2/g2g/condB/FVB_full.vci.gz \
      -i "$bam" -o "$rev" --reverse > "$OUT/logs/${s}.condB.convert.log" 2>&1
    if [ ! -s "$rev" ] || ! "$SAMTOOLS" quickcheck "$rev" 2>/dev/null; then
      echo "[$(date +%T)] FAIL convert: $s condB (BAM kept)"; return 1; fi
    "$SAMTOOLS" sort -@ 4 -o "$srt" "$rev" && "$SAMTOOLS" index "$srt"
    if [ ! -s "$srt" ]; then echo "[$(date +%T)] FAIL sort: $s condB (BAM kept)"; return 1; fi
    rm -f "$rev"
    target="$srt"
  fi
  echo "[$(date +%T)] featureCounts -s 2 --countReadPairs: $s cond$c"
  "$FC" -T 4 -p --countReadPairs -s 2 -a "$GTF" -t exon -g gene_id \
    -o "$OUT/counts/${s}.cond${c}.counts.txt" "$target" > "$OUT/counts/${s}.cond${c}.featureCounts.log" 2>&1
  if [ ! -s "$OUT/counts/${s}.cond${c}.counts.txt" ]; then
    echo "[$(date +%T)] FAIL featureCounts: $s cond$c (BAMs kept)"; return 1; fi
  rm -rf "$OUT/star/${s}.cond${c}"
  [ "$c" = "B" ] && rm -f "$target" "${target}.bai"
  echo "[$(date +%T)] done: $s cond$c  (disk free $(df -h "$ROOT" | tail -1 | awk '{print $4}'))"
}

"$FC" -v 2>&1 | grep -o 'v[0-9.]*' | head -1 > "$OUT/featureCounts_version.txt" || true
"$STAR_BIN" --version > "$OUT/STAR_version.txt" 2>&1 || true

n=${#JOBS[@]}
echo "[$(date +%T)] $n jobs: ${JOBS[*]}"
prefetch_pid=""
for i in "${!JOBS[@]}"; do
  s="${JOBS[$i]%%:*}"; c="${JOBS[$i]##*:}"
  if [ -s "$OUT/counts/${s}.cond${c}.counts.txt" ]; then
    echo "[$(date +%T)] SKIP (already counted): $s cond$c"; continue
  fi
  if [ -n "$prefetch_pid" ]; then wait "$prefetch_pid"; prefetch_pid=""
  else echo "[$(date +%T)] STAR: $s cond$c"; star_one "$s" "$c"; fi
  # start the next alignment while this sample is converted and counted
  j=$((i+1))
  if [ $j -lt $n ]; then
    ns="${JOBS[$j]%%:*}"; nc="${JOBS[$j]##*:}"
    if [ ! -s "$OUT/counts/${ns}.cond${nc}.counts.txt" ]; then
      echo "[$(date +%T)] STAR (prefetch): $ns cond$nc"
      star_one "$ns" "$nc" & prefetch_pid=$!
    fi
  fi
  quantify_one "$s" "$c" || echo "[$(date +%T)] job failed: $s cond$c"
done
[ -n "$prefetch_pid" ] && wait "$prefetch_pid"
echo "ALL_PHASE2_S2PAIRS_DONE"
