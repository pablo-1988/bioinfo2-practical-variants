#!/usr/bin/env bash
#
# STEP 3 — Align the reads to the reference with bwa mem
#
#   ./scripts/03_map_bwa.sh
#
# This is the step that turns "a pile of 150 bp strings" into "evidence about
# positions in a genome". Everything afterwards is interpretation of it.
#
# Output: results/02_mapping/<sample>.sorted.bam (+ .bai)

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

require bwa samtools
need_file "$REF.bwt"
mkdir -p "$MAP"

for s in $SAMPLES; do
  need_file "$RAW/${s}_R1.fastq.gz" "$RAW/${s}_R2.fastq.gz"

  log "[$s] bwa mem -> sorted BAM"
  #
  #  -R '@RG\t...'
  #        The read group. This is NOT bookkeeping you can skip: GATK
  #        refuses to run without it, because BQSR models error rates per
  #        read group and every downstream tool identifies the sample by
  #        the SM field. Get SM wrong and your VCF has the wrong name on it.
  #
  #        ID  unique id for this run of this library
  #        SM  SAMPLE NAME -- this is what ends up in the VCF column header
  #        PL  platform, ILLUMINA
  #        LB  library; MarkDuplicates uses it to decide what may duplicate
  #        PU  platform unit, flowcell.lane.barcode -- the real unit of error
  #
  #  -M    flag split alignments as secondary. Picard-era tools expect this.
  #
  # bwa writes SAM to stdout; we never let that hit the disk. samtools sort
  # reads the stream and writes a coordinate-sorted BAM directly.
  #
  run bash -c "bwa mem -t $THREADS -M \
      -R '@RG\tID:${s}\tSM:${s}\tPL:ILLUMINA\tLB:lib1\tPU:FLOWCELLX.1.${s}' \
      '$REF' '$RAW/${s}_R1.fastq.gz' '$RAW/${s}_R2.fastq.gz' \
    | samtools sort -@ $THREADS -o '$MAP/${s}.sorted.bam' -"

  run samtools index -@ "$THREADS" "$MAP/${s}.sorted.bam"

  log "[$s] samtools flagstat"
  run samtools flagstat "$MAP/${s}.sorted.bam"
done

cat <<'TXT'

--------------------------------------------------------------------------
Read the flagstat output above, do not scroll past it.

  - "properly paired" should be high (>98%). If it is not, the insert size
    or the orientation is wrong, and every variant call downstream inherits
    the problem.
  - "mapped" being ~100% is a property of SIMULATED data. Real data has
    contamination, adapters, and unplaced sequence. Never assume it.

Now look at a read by hand:

    samtools view results/02_mapping/sample01.sorted.bam | head -3

Identify, in the columns: the FLAG, the position, the MAPQ, and the CIGAR.
Then find reads the aligner could not place uniquely:

    samtools view -c -q 1 results/02_mapping/sample01.sorted.bam   # MAPQ >= 1
    samtools view -c results/02_mapping/sample01.sorted.bam        # all

The difference is the MAPQ 0 reads. Where are they?

    samtools view results/02_mapping/sample01.sorted.bam \
      | awk '$5 == 0 {print $3"\t"int($4/10000)*10000}' | sort | uniq -c | sort -rn | head

Two coordinate windows dominate. That is not noise -- the genome contains
the same 2 kb segment twice, and bwa cannot know which copy a read came
from, so it honestly reports MAPQ 0. Remember those coordinates.
--------------------------------------------------------------------------
TXT
