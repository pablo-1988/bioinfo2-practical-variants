#!/usr/bin/env bash
#
# STEP 8 — Joint genotyping across all samples
#
#   ./scripts/08_joint_genotyping.sh
#
# Genotyping the samples together is not a convenience. It changes the
# answer. A site with weak evidence in sample03 (2 reads, low confidence)
# would be dropped if that sample were called alone; when the same allele is
# clean and heterozygous in sample01 and sample02, the prior at that site
# shifts and the marginal call in sample03 becomes defensible. This is where
# the low-coverage sample gets rescued -- partially.
#
# Two commands:
#   GenomicsDBImport   loads the per-sample GVCFs into a columnar store
#   GenotypeGVCFs      does the actual joint calling over that store
#
# Output: results/07_joint_genotyping/cohort.raw.vcf.gz

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

require gatk bcftools
mkdir -p "$JOINT"

DB="$JOINT/genomicsdb"
VARGS=()
for s in $SAMPLES; do
  need_file "$GVCF/${s}.g.vcf.gz"
  VARGS+=(-V "$GVCF/${s}.g.vcf.gz")
done

log "GenomicsDBImport — one store for all samples"
# The workspace must NOT exist beforehand; GATK refuses to overwrite it.
run rm -rf "$DB"
run gatk_run GenomicsDBImport \
  "${VARGS[@]}" \
  --genomicsdb-workspace-path "$DB" \
  $(intervals_args)

log "GenotypeGVCFs — joint calling"
run gatk_run GenotypeGVCFs \
  -R "$REF" \
  -V "gendb://$DB" \
  -O "$JOINT/cohort.raw.vcf.gz"

log "How many raw calls did we get?"
run bash -c "bcftools stats '$JOINT/cohort.raw.vcf.gz' | grep -E '^SN' | cut -f3-"

cat <<'TXT'

--------------------------------------------------------------------------
This VCF is RAW. Every site the caller could justify is in it, including the
ones it should not have justified. Nothing has been filtered yet -- the
FILTER column says '.' on every line.

Read the header, it documents everything you are about to filter on:

    bcftools view -h results/07_joint_genotyping/cohort.raw.vcf.gz | grep '^##INFO'

Now go to the two regions you flagged in step 6, and notice that they fail in
OPPOSITE directions.

1. The duplicated segment, chr1:40,000-42,000:

    bcftools view -H results/07_joint_genotyping/cohort.raw.vcf.gz chr1:40000-42000 | wc -l

Almost nothing. The caller did not make mistakes there -- it made no
statements at all, because MAPQ 0 reads are below its threshold. Silence.
Real variants are in that region, and you will miss every one of them.

2. The over-covered window on chr2, around 20,000-23,000:

    bcftools query -r chr2:20000-23000 \
      -f '%POS\t%REF\t%ALT\tDP=%DP\tMQ=%MQ\tQD=%QD\t[%GT %AD; ]\n' \
      results/07_joint_genotyping/cohort.raw.vcf.gz | head -15

Dozens of calls, densely packed, most of them heterozygous in ALL THREE
samples, with allele depths near 50/50 and MQ of 60. By every annotation in
the VCF they look excellent.

Ask yourself what could produce a dense cluster of confident heterozygous
sites, at double the normal depth, identical in three unrelated individuals.
Real polymorphism does not behave like that.

Reproducibility across samples is not evidence of truth. A systematic
artefact reproduces perfectly.
--------------------------------------------------------------------------
TXT
