#!/usr/bin/env bash
#
# STEP 4 — Mark PCR/optical duplicates
#
#   ./scripts/04_mark_duplicates.sh
#
# A duplicate is the same original DNA molecule sequenced more than once,
# usually because PCR amplified it. It looks like independent evidence and
# is not. If a duplicated molecule carries a sequencing error, duplicates
# turn one error into "5 reads support this allele" -- a false variant with
# convincing depth.
#
# MarkDuplicates finds them by ALIGNMENT COORDINATES (both ends of the
# fragment), not by sequence identity: two copies of one molecule have
# independent sequencing errors, so their sequences differ slightly, but
# they start and end at exactly the same place.
#
# It MARKS them (flag 0x400) rather than deleting them. Nothing is thrown
# away; HaplotypeCaller simply ignores flagged reads.
#
# Output: results/03_dedup/<sample>.dedup.bam (+ .bai) and *.metrics.txt

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

require gatk samtools
mkdir -p "$DEDUP"

for s in $SAMPLES; do
  need_file "$MAP/${s}.sorted.bam"

  log "[$s] gatk MarkDuplicates"
  run gatk_run MarkDuplicates \
    -I "$MAP/${s}.sorted.bam" \
    -O "$DEDUP/${s}.dedup.bam" \
    -M "$DEDUP/${s}.metrics.txt" \
    --CREATE_INDEX true

  # GATK writes the index as sample.dedup.bai; samtools and IGV also look
  # for sample.dedup.bam.bai. Make both names available.
  [ -f "$DEDUP/${s}.dedup.bai" ] && cp -f "$DEDUP/${s}.dedup.bai" "$DEDUP/${s}.dedup.bam.bai"

  log "[$s] duplicate fraction"
  run bash -c "grep -A2 '^LIBRARY' '$DEDUP/${s}.metrics.txt' | cut -f1-10 | column -t"
done

cat <<'TXT'

--------------------------------------------------------------------------
PERCENT_DUPLICATION in the tables above is the number to report. One sample
is far worse than the other two -- that library was over-amplified.

Two questions to answer in your report:

  1. sample02 has ~30x raw coverage. After marking duplicates, what is its
     EFFECTIVE coverage -- the number of independent molecules? Compute it.
     Compare it to sample03, which has ~12x and almost no duplicates.

  2. We marked duplicates but did not remove them. Give one concrete reason
     why keeping them in the file is better than deleting them.

Confirm the marking landed, using the SAM flag directly:

    samtools view -c -f 1024 results/03_dedup/sample02.dedup.bam
    samtools flagstat results/03_dedup/sample02.dedup.bam
--------------------------------------------------------------------------
TXT
