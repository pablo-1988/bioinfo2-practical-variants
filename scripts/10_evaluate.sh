#!/usr/bin/env bash
#
# STEP 10 — Compare your calls against the truth
#
#   ./scripts/10_evaluate.sh
#
# This is the step you can only do on simulated data, which is exactly why
# the data is simulated. In a real project nobody hands you the answer, so
# you never learn what your pipeline missed. Here you can.
#
# data/truth/truth.vcf holds the variants that were actually planted, with
# the real genotype of each individual. gatk Concordance compares your calls
# against it and reports, per variant type:
#
#   sensitivity (recall) = of the real variants, how many did you find?
#   precision            = of your calls, how many are real?
#
# Output: results/09_evaluation/

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

require gatk bcftools bgzip tabix
mkdir -p "$EVAL"
need_file "data/truth/truth.vcf" "$FILT/cohort.filtered.vcf.gz"

log "Index the truth set"
run bash -c "bgzip -f -c data/truth/truth.vcf > '$TRUTH'"
run tabix -f -p vcf "$TRUTH"

# --- Normalisation, and why it is not a formality ---------------------
#
# The same indel can be written in several ways. A deletion of "CT" inside
# CTCTCT can be placed at three different positions, all of them describing
# the identical molecule. If your truth set writes it one way and your caller
# writes it another, a naive comparison scores ONE false positive and ONE
# false negative for a variant you got exactly right.
#
# bcftools norm fixes this by left-aligning every indel against the reference
# and splitting multi-allelic records, so both files use the same convention.
# Skip this and your indel numbers are simply wrong -- in this dataset it
# costs about 20% of the indel accuracy, invented out of nothing.
log "Normalise both sides before comparing (left-align indels)"
run bash -c "bcftools norm -f '$REF' -m -any -O z \
    -o '$EVAL/truth.norm.vcf.gz' '$TRUTH' 2>&1 | tail -2"
run bcftools index -f -t "$EVAL/truth.norm.vcf.gz"

run bash -c "bcftools norm -f '$REF' -m -any -O z \
    -o '$EVAL/calls.norm.vcf.gz' '$FILT/cohort.filtered.vcf.gz' 2>&1 | tail -2"
run bcftools index -f -t "$EVAL/calls.norm.vcf.gz"

run bash -c "bcftools norm -f '$REF' -m -any -O z \
    -o '$EVAL/raw.norm.vcf.gz' '$JOINT/cohort.raw.vcf.gz' 2>&1 | tail -2"
run bcftools index -f -t "$EVAL/raw.norm.vcf.gz"

for s in $SAMPLES; do
  log "[$s] extract this sample from the truth and from your calls"
  run gatk_run SelectVariants \
    -R "$REF" -V "$EVAL/truth.norm.vcf.gz" -sn "$s" \
    --exclude-non-variants --remove-unused-alternates \
    -O "$EVAL/${s}.truth.vcf.gz"

  # PASS sites only: filtered sites are calls you already rejected, so
  # counting them against yourself would be dishonest in the other direction.
  run gatk_run SelectVariants \
    -R "$REF" -V "$EVAL/calls.norm.vcf.gz" -sn "$s" \
    --exclude-non-variants --remove-unused-alternates \
    --exclude-filtered \
    -O "$EVAL/${s}.calls.vcf.gz"

  log "[$s] Concordance against the truth"
  run gatk_run Concordance \
    -R "$REF" \
    -eval "$EVAL/${s}.calls.vcf.gz" \
    --truth "$EVAL/${s}.truth.vcf.gz" \
    --summary "$EVAL/${s}.concordance.tsv"

  run bash -c "column -t '$EVAL/${s}.concordance.tsv'"

  log "[$s] genotype concordance (did you get het vs hom right?)"
  run bash -c "bcftools stats -s '$s' '$EVAL/${s}.calls.vcf.gz' \
      | grep -E '^PSC' | cut -f3-8 | column -t"
done

log "Same comparison WITHOUT filtering, to price the filters"
for s in $SAMPLES; do
  gatk_run SelectVariants -R "$REF" -V "$EVAL/raw.norm.vcf.gz" -sn "$s" \
    --exclude-non-variants --remove-unused-alternates \
    -O "$EVAL/${s}.raw_calls.vcf.gz" > /dev/null 2>&1
  gatk_run Concordance -R "$REF" \
    -eval "$EVAL/${s}.raw_calls.vcf.gz" --truth "$EVAL/${s}.truth.vcf.gz" \
    --summary "$EVAL/${s}.concordance_raw.tsv" > /dev/null 2>&1
  echo "--- $s, before filtering ---"
  column -t "$EVAL/${s}.concordance_raw.tsv"
done

cat <<'TXT'

--------------------------------------------------------------------------
Build this table for your report, from the files in results/09_evaluation/:

              SNP sens.  SNP prec.  INDEL sens.  INDEL prec.
  sample01
  sample02
  sample03

and the same table before filtering.

Then answer, with the numbers in front of you:

  1. sample01 and sample02 both have ~30x raw coverage. Do they perform
     equally? Explain the difference using the duplicate metrics from step 4.

  2. sample03 has ~12x. Where does it lose -- sensitivity or precision? Why
     is that the expected direction? Which genotype does it get wrong most
     often, and in which direction (het called as hom-ref, or the reverse)?

  3. Compare filtered against raw. Filtering always trades sensitivity for
     precision. Quantify the trade for SNPs: how many true calls did you
     lose, and how many false ones did you remove? Was it worth it?

  4. Locate your false positives:

         bcftools isec -C results/09_evaluation/sample01.calls.vcf.gz \
             results/09_evaluation/sample01.truth.vcf.gz | head -20

     They are not scattered. Where are they, and is the count similar in all
     three samples? Explain the mechanism -- you have everything you need
     from steps 6, 8 and 9. Note that this class of error does not go away
     with more coverage: sequencing deeper makes it MORE convincing.

  5. Now locate your false negatives:

         bcftools isec -C results/09_evaluation/sample01.truth.vcf.gz \
             results/09_evaluation/sample01.calls.vcf.gz | head -20

     A different region, a different mechanism. Which one, and why is "no
     call" the honest output there rather than a failure?

  6. Normalisation. This script left-aligned both VCFs before comparing.
     Find out what that was worth -- compare against the un-normalised truth:

         gatk SelectVariants -R data/reference/reference.fasta \
           -V data/truth/truth.vcf.gz -sn sample01 --exclude-non-variants \
           -O /tmp/s1.raw_truth.vcf.gz
         gatk Concordance -R data/reference/reference.fasta \
           -eval results/09_evaluation/sample01.calls.vcf.gz \
           --truth /tmp/s1.raw_truth.vcf.gz --summary /tmp/s1.nonorm.tsv
         column -t /tmp/s1.nonorm.tsv

     Report both indel numbers. The difference is entirely representation:
     the same deletion written at a different position. Nothing about the
     biology changed. This is the single most common way published variant
     comparisons are quietly wrong.
--------------------------------------------------------------------------
TXT
