#!/usr/bin/env bash
#
# STEP 5 — Base Quality Score Recalibration (BQSR)
#
#   ./scripts/05_bqsr.sh
#
# The instrument reports a quality score per base: "I am 99.9% sure this is
# a T". Those numbers are systematically wrong. They drift with the cycle
# number, the sequence context, and the machine. HaplotypeCaller's model is
# built on those probabilities being honest, so we fix them first.
#
# The trick is how you measure the true error rate without knowing the truth:
#
#   1. Align the reads. Any base that mismatches the reference is either a
#      sequencing error or a real variant.
#   2. Mask out every position in the known-sites database (dbSNP). What is
#      left is dominated by errors.
#   3. Tabulate the observed error rate by reported quality, by cycle, and
#      by dinucleotide context. That is the recalibration table.
#   4. Rewrite every base quality with the empirical value.
#
# Note the consequence: real variants NOT in the database get counted as
# errors and quietly push qualities down. That is why the database matters,
# and why BQSR is a bad idea on a species with no variant catalogue.
#
# Output: results/04_bqsr/<sample>.recal.bam, before/after tables and a plot

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

require gatk
need_file "$KNOWN"
mkdir -p "$BQSR"

for s in $SAMPLES; do
  need_file "$DEDUP/${s}.dedup.bam"

  log "[$s] BaseRecalibrator — build the model"
  run gatk_run BaseRecalibrator \
    -R "$REF" \
    -I "$DEDUP/${s}.dedup.bam" \
    --known-sites "$KNOWN" \
    -O "$BQSR/${s}.recal.table"

  log "[$s] ApplyBQSR — rewrite the base qualities"
  run gatk_run ApplyBQSR \
    -R "$REF" \
    -I "$DEDUP/${s}.dedup.bam" \
    --bqsr-recal-file "$BQSR/${s}.recal.table" \
    -O "$BQSR/${s}.recal.bam"

  [ -f "$BQSR/${s}.recal.bai" ] && cp -f "$BQSR/${s}.recal.bai" "$BQSR/${s}.recal.bam.bai"

  log "[$s] BaseRecalibrator again, on the corrected BAM — did it work?"
  run gatk_run BaseRecalibrator \
    -R "$REF" \
    -I "$BQSR/${s}.recal.bam" \
    --known-sites "$KNOWN" \
    -O "$BQSR/${s}.after.table"

  log "[$s] AnalyzeCovariates — before/after plot"
  # This one needs R. If it fails, the pipeline continues: the plot is nice
  # to have, the tables above are the actual evidence.
  gatk_run AnalyzeCovariates \
    -before "$BQSR/${s}.recal.table" \
    -after "$BQSR/${s}.after.table" \
    -plots "$BQSR/${s}.bqsr.pdf" \
    || echo "    (AnalyzeCovariates failed -- needs R + ggplot2. Not fatal.)"
done

cat <<'TXT'

--------------------------------------------------------------------------
Open one of the tables and read it as a table, not as output to scroll past:

    grep -A 12 '#:GATKTable:.*RecalTable1' results/04_bqsr/sample01.recal.table

The columns that matter are QualityScore (what the machine claimed),
EmpiricalQuality (what the data says), and Observations. Find a reported
quality where the empirical value is clearly lower. The machine was
over-confident there.

Then compare before and after:

    grep -A 12 '#:GATKTable:.*RecalTable1' results/04_bqsr/sample01.after.table

After recalibration, reported and empirical should sit close to each other.
That is the entire goal -- not "better" qualities, HONEST ones.

If the PDF was produced, look at results/04_bqsr/sample01.bqsr.pdf, panel
"Quality Score Covariate": the before points sit off the diagonal, the after
points sit on it.

Report: does the size of the correction differ between the three samples?
Which covariate -- reported quality, cycle, or context -- carries most of it
in this dataset?
--------------------------------------------------------------------------
TXT
