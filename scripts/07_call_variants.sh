#!/usr/bin/env bash
#
# STEP 7 — Call variants per sample, in GVCF mode
#
#   ./scripts/07_call_variants.sh
#
# HaplotypeCaller does not compare base by base. In any region that looks
# active it throws away the alignment, reassembles the reads into candidate
# haplotypes with a de Bruijn graph, realigns every read against each
# haplotype (pair-HMM), and only then genotypes. That is why it handles
# indels and clustered variants that a pileup caller mangles.
#
# -ERC GVCF changes the output, not the calling: instead of only variant
# sites, it emits a record for EVERY position, including blocks of reference
# with a confidence that they are non-variant. That distinction matters:
#
#     no variant in the VCF  =  "I found nothing here"
#     GVCF reference block   =  "I looked, I had 32x, and it is reference"
#
# Without it, joint genotyping in step 8 could not tell a homozygous
# reference call from a position nobody sequenced.
#
# Output: results/06_gvcf/<sample>.g.vcf.gz

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

require gatk
mkdir -p "$GVCF"

for s in $SAMPLES; do
  need_file "$BQSR/${s}.recal.bam"

  log "[$s] HaplotypeCaller -ERC GVCF"
  #  -ERC GVCF                     per-sample intermediate, see above
  #  --native-pair-hmm-threads     the pair-HMM is the slow part
  run gatk_run HaplotypeCaller \
    -R "$REF" \
    -I "$BQSR/${s}.recal.bam" \
    -O "$GVCF/${s}.g.vcf.gz" \
    -ERC GVCF \
    --native-pair-hmm-threads "$THREADS"
done

cat <<'TXT'

--------------------------------------------------------------------------
Look at what a GVCF actually contains:

    zcat < results/06_gvcf/sample01.g.vcf.gz | grep -v '^##' | head -20

Find a <NON_REF> reference block and read its END= and its GQ. Then find a
line that is a real candidate variant. Notice that the genotypes here are
provisional -- they get decided in step 8, not now.

Question: HaplotypeCaller is the slowest step in this pipeline even on a
240 kb genome. Name the operation that costs the time, and explain why
reassembling a region is worth that cost compared to reading the pileup.

Optional but strongly recommended, if you have IGV installed: load
data/reference/reference.fasta and results/04_bqsr/sample01.recal.bam, and
go to a site the caller flagged. Seeing the reads under a call once is worth
more than any amount of reading about it.
--------------------------------------------------------------------------
TXT
