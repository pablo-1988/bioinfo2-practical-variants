#!/usr/bin/env bash
#
# STEP 6 — Judge the alignments before you trust any variant from them
#
#   ./scripts/06_alignment_qc.sh
#
# Every variant call is a claim about read evidence. If the coverage is
# uneven, or a region has no uniquely-mapping reads, the caller will still
# emit something -- confidently. Measure the alignments first.
#
# Output: results/05_alignment_qc/*.txt and a MultiQC report

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

require gatk samtools multiqc
mkdir -p "$ALNQC"

for s in $SAMPLES; do
  need_file "$BQSR/${s}.recal.bam"

  log "[$s] samtools flagstat / stats"
  run bash -c "samtools flagstat '$BQSR/${s}.recal.bam' | tee '$ALNQC/${s}.flagstat.txt'"
  run bash -c "samtools stats '$BQSR/${s}.recal.bam' > '$ALNQC/${s}.stats.txt'"

  log "[$s] CollectAlignmentSummaryMetrics"
  run gatk_run CollectAlignmentSummaryMetrics \
    -R "$REF" -I "$BQSR/${s}.recal.bam" \
    -O "$ALNQC/${s}.alignment_metrics.txt"

  log "[$s] CollectInsertSizeMetrics"
  gatk_run CollectInsertSizeMetrics \
    -I "$BQSR/${s}.recal.bam" \
    -O "$ALNQC/${s}.insert_size_metrics.txt" \
    -H "$ALNQC/${s}.insert_size_histogram.pdf" \
    || echo "    (histogram PDF needs R -- the metrics file was still written.)"

  log "[$s] Depth — twice: all reads, and only reads the caller will use"
  #
  # samtools depth -a reports every position, including the zeros. Without
  # -a you would average only over covered positions and flatter yourself.
  #
  # -Q 20 keeps only reads with MAPQ >= 20. That is roughly what
  # HaplotypeCaller uses, and the difference between the two numbers is the
  # part of your coverage that exists but cannot be used.
  #
  # Note also what samtools depth does silently: it skips reads flagged as
  # duplicates. So sample02's depth here is ALREADY its effective depth, not
  # its nominal 30x. Check that against the arithmetic you did in step 4.
  #
  run bash -c "samtools depth -a '$BQSR/${s}.recal.bam' \
      | awk '{sum+=\$3; n++; if(\$3<10) low++} \
             END{printf \"all reads      : mean %.1fx, %.2f%% below 10x\n\", \
                 sum/n, 100*low/n}' | tee '$ALNQC/${s}.depth_summary.txt'"

  run bash -c "samtools depth -a -Q 20 '$BQSR/${s}.recal.bam' \
      | awk '{sum+=\$3; n++; if(\$3<10) low++} \
             END{printf \"MAPQ >= 20 only: mean %.1fx, %.2f%% below 10x\n\", \
                 sum/n, 100*low/n}' | tee -a '$ALNQC/${s}.depth_summary.txt'"

  # Per-window depth, so the unevenness is visible rather than averaged away.
  # Partial windows at the end of a contig are dropped -- they are an artefact
  # of the binning, not of the data.
  for Q in 0 20; do
    run bash -c "samtools depth -a -Q $Q '$BQSR/${s}.recal.bam' \
        | awk '{w=int(\$2/1000); c[\$1\"\t\"w]+=\$3; n[\$1\"\t\"w]++} \
               END{for(k in c) if(n[k]>=500) printf \"%s\t%.1f\n\", k, c[k]/n[k]}' \
        | sort -k1,1 -k2,2n > '$ALNQC/${s}.depth_1kb_q${Q}.tsv'"
  done
done

log "MultiQC — every metric from every sample in one page"
run multiqc --force --outdir "$ALNQC" --filename multiqc_alignment.html "$ALNQC" "$DEDUP"

cat <<'TXT'

--------------------------------------------------------------------------
Open results/05_alignment_qc/multiqc_alignment.html.

Then compare the two depth numbers printed above for each sample. "All
reads" and "MAPQ >= 20 only" are not the same genome. The gap is coverage
that physically exists in your BAM and that the caller will refuse to use.

Find where the gap lives:

    join -j 1 \
      <(awk '{print $1"_"$2"\t"$3}' results/05_alignment_qc/sample01.depth_1kb_q0.tsv  | sort) \
      <(awk '{print $1"_"$2"\t"$3}' results/05_alignment_qc/sample01.depth_1kb_q20.tsv | sort) \
      | awk '{d=$2-$3; if (d > 10) print $0"\tlost="d}' | sort -k4,4

The windows that lose almost all of their coverage are the two copies of the
duplicated segment from step 3: the reads are there, they are ambiguous, and
MAPQ 0 makes them invisible to the caller. In a naive depth plot that region
looks FINE -- even slightly above average, because reads from both copies
pile onto both.

Now look at the other end of the distribution:

    sort -k3,3nr results/05_alignment_qc/sample01.depth_1kb_q20.tsv | head -5

One window on chr2 has roughly TWICE the mean depth, with perfectly normal
MAPQ. Nothing was ambiguous, nothing was filtered -- there is simply twice
as much sequence there as the reference can account for. Write down the
coordinates and form a hypothesis now; you will see the consequence in
step 8.

For your report:

  1. Mean depth of each sample, both ways, and the fraction below 10x.
  2. The insert size distribution: mean and standard deviation. Does it
     match what a 380 bp library should look like?
  3. The two anomalous regions, with coordinates: one that loses coverage to
     MAPQ, and one that has too much of it. Say what you expect each to do
     to the variant calls.
--------------------------------------------------------------------------
TXT
