#!/usr/bin/env bash
#
# STEP 9 — Filter the raw calls (GATK hard filtering)
#
#   ./scripts/09_filter_variants.sh
#
# Best practice for large human cohorts is VQSR / VariantRecalibrator, which
# learns the boundary between true and false calls from a training set of
# known variants. It needs tens of thousands of variants to fit its Gaussian
# mixture. We have about a thousand. VQSR would fail, loudly, and GATK's own
# documentation says to use hard filters in exactly this situation -- small
# cohorts, non-model organisms, targeted panels.
#
# So: fixed thresholds on the annotations, applied separately to SNPs and to
# indels because they fail in different ways.
#
# Note VariantFiltration MARKS, it does not delete. Filtered sites stay in
# the file with a reason in the FILTER column. You can always audit them.
#
# Output: results/08_filtered/cohort.filtered.vcf.gz  (all sites, flagged)
#         results/08_filtered/cohort.pass.vcf.gz      (PASS only)

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

require gatk bcftools
mkdir -p "$FILT"
need_file "$JOINT/cohort.raw.vcf.gz"

log "Split SNPs and indels — they need different thresholds"
run gatk_run SelectVariants \
  -R "$REF" -V "$JOINT/cohort.raw.vcf.gz" \
  --select-type-to-include SNP \
  -O "$FILT/raw.snps.vcf.gz"

run gatk_run SelectVariants \
  -R "$REF" -V "$JOINT/cohort.raw.vcf.gz" \
  --select-type-to-include INDEL \
  --select-type-to-include MIXED \
  -O "$FILT/raw.indels.vcf.gz"

log "VariantFiltration on SNPs"
#  QD  < 2    variant confidence divided by depth. A high QUAL that is only
#             high because the site has 500x is not evidence of anything.
#  FS  > 60   strand bias (Fisher). Real alleles appear on both strands.
#  SOR > 3    strand bias again, better behaved at high coverage.
#  MQ  < 40   ROOT MEAN SQUARE MAPPING QUALITY. This is the one that removes
#             your duplicated segment: reads there have MAPQ 0.
#  MQRankSum      < -12.5  do ALT reads have systematically worse MAPQ?
#  ReadPosRankSum < -8     does the ALT only ever appear at read ends?
run gatk_run VariantFiltration \
  -R "$REF" -V "$FILT/raw.snps.vcf.gz" \
  -O "$FILT/filtered.snps.vcf.gz" \
  --filter-name "QD2"      --filter-expression "QD < 2.0" \
  --filter-name "QUAL30"   --filter-expression "QUAL < 30.0" \
  --filter-name "SOR3"     --filter-expression "SOR > 3.0" \
  --filter-name "FS60"     --filter-expression "FS > 60.0" \
  --filter-name "MQ40"     --filter-expression "MQ < 40.0" \
  --filter-name "MQRS-12.5"  --filter-expression "MQRankSum < -12.5" \
  --filter-name "RPRS-8"     --filter-expression "ReadPosRankSum < -8.0"

log "VariantFiltration on indels — looser, because indels are messier"
run gatk_run VariantFiltration \
  -R "$REF" -V "$FILT/raw.indels.vcf.gz" \
  -O "$FILT/filtered.indels.vcf.gz" \
  --filter-name "QD2"      --filter-expression "QD < 2.0" \
  --filter-name "QUAL30"   --filter-expression "QUAL < 30.0" \
  --filter-name "FS200"    --filter-expression "FS > 200.0" \
  --filter-name "RPRS-20"  --filter-expression "ReadPosRankSum < -20.0"

log "Merge the two back into one VCF"
run gatk_run MergeVcfs \
  -I "$FILT/filtered.snps.vcf.gz" \
  -I "$FILT/filtered.indels.vcf.gz" \
  -O "$FILT/cohort.filtered.vcf.gz"

log "And a PASS-only copy, for anyone downstream who wants the short answer"
run bash -c "bcftools view -f PASS -O z -o '$FILT/cohort.pass.vcf.gz' '$FILT/cohort.filtered.vcf.gz'"
run bcftools index -f -t "$FILT/cohort.pass.vcf.gz"

log "What each filter removed"
run bash -c "bcftools query -f '%FILTER\n' '$FILT/cohort.filtered.vcf.gz' \
    | sort | uniq -c | sort -rn"

cat <<'TXT'

--------------------------------------------------------------------------
The table above is the honest summary of this step: how many sites passed,
and which criterion each rejected site failed.

Read it before you feel good about it. The hard filters removed on the order
of a dozen sites -- and the cluster of suspicious calls you found on chr2 in
step 8 is still mostly there:

    bcftools query -r chr2:20000-23000 -f '%POS\t%FILTER\tDP=%DP\tMQ=%MQ\n' \
      results/08_filtered/cohort.filtered.vcf.gz | head

Most of them say PASS. That is not a bug in the filters; it is what they are
for. Every threshold in this step tests whether the READS AT A SITE ARE
CONSISTENT -- balanced strands, good mapping quality, ALT not confined to
read ends. Those reads are perfectly consistent. They are consistent
evidence for the wrong genome, because the sample contains a copy of that
segment that the reference does not have, so its reads have nowhere else to
go and pile onto the single reference copy. Their differences then read as
heterozygous SNPs.

No per-site annotation can see that. What gives it away is DEPTH:

    bcftools query -f '%CHROM\t%POS\t%DP\n' results/08_filtered/cohort.pass.vcf.gz \
      | awk '{s+=$3; n++} END{print "mean cohort DP:", s/n}'

    bcftools query -f '%CHROM\t%POS\t%DP\n' results/08_filtered/cohort.pass.vcf.gz \
      | awk '$3 > 130' | head

EXERCISE. Add a depth filter of your own and see what it buys you:

    gatk VariantFiltration -R data/reference/reference.fasta \
      -V results/08_filtered/cohort.filtered.vcf.gz \
      -O results/08_filtered/cohort.dpfilt.vcf.gz \
      --filter-name "DPexcess" --filter-expression "DP > 130.0"

Then re-run step 10 against it and report the change in precision and in
sensitivity. Both will move. Decide, with the numbers, whether you would
keep the filter.

Also check the filter did what you predicted where it SHOULD work:

    bcftools view -H results/08_filtered/cohort.filtered.vcf.gz chr1:40000-42000 \
      | cut -f1,2,4,5,7

For your report: which single filter removed the most sites, and did it
remove the calls you actually distrust? If not, what would?
--------------------------------------------------------------------------
TXT
