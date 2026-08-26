#!/usr/bin/env bash
#
# STEP 1 — Look at the reads before you do anything with them
#
#   ./scripts/01_qc_raw.sh
#
# You already did read QC in the previous practical, so this is short. The
# point here is not to repeat FastQC: it is to decide, on evidence, whether
# these libraries can go into an aligner as they are.
#
# Output: results/01_qc_raw/*.html  and  results/01_qc_raw/multiqc_report.html

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

require fastqc multiqc
mkdir -p "$QC_RAW"

for s in $SAMPLES; do
  need_file "$RAW/${s}_R1.fastq.gz" "$RAW/${s}_R2.fastq.gz"
done

log "FastQC on all raw FASTQ files"
run fastqc --threads "$THREADS" --outdir "$QC_RAW" \
  $(for s in $SAMPLES; do echo "$RAW/${s}_R1.fastq.gz" "$RAW/${s}_R2.fastq.gz"; done)

log "MultiQC — one report for the whole set"
run multiqc --force --outdir "$QC_RAW" --filename multiqc_report.html "$QC_RAW"

cat <<'TXT'

--------------------------------------------------------------------------
Open results/01_qc_raw/multiqc_report.html.

Three things to write down before moving on:

  1. How many reads does each sample have? sample03 has far fewer. Predict
     now, in writing, what that will do to its variant calls.

  2. The per-base quality decays towards the 3' end. Note roughly where.
     You will meet that same decay again in step 5, from the other side.

  3. FastQC reports a duplication level per sample. One of the three is
     clearly worse. Which one -- and does FastQC's estimate, computed from
     sequence identity alone, actually prove those are PCR duplicates?

Note what we are NOT doing: trimming. bwa mem soft-clips adapters and bad
tails by itself, and GATK reads the soft-clip. Trimming before alignment
throws away information the caller could have used. See PIPELINE.md.
--------------------------------------------------------------------------
TXT
