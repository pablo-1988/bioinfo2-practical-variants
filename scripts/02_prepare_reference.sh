#!/usr/bin/env bash
#
# STEP 2 — Index the reference genome
#
#   ./scripts/02_prepare_reference.sh
#
# A reference FASTA is useless to the tools until it is indexed, and each
# tool wants a different index:
#
#   bwa index    -> .amb .ann .bwt .pac .sa   the FM-index bwa searches
#   samtools faidx -> .fai                    random access by coordinate
#   gatk CreateSequenceDictionary -> .dict     contig names and lengths
#
# GATK refuses to start if any of these is missing, and the error message it
# gives is famously unhelpful. Do this once, first.
#
# Output: index files next to data/reference/reference.fasta
#         data/reference/known_sites.vcf.gz (+ .tbi)

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

require bwa samtools gatk bgzip tabix
need_file "$REF" "$REFDIR/known_sites.vcf"

log "bwa index — the aligner's index of the reference"
run bwa index "$REF"

log "samtools faidx — coordinate index"
run samtools faidx "$REF"

log "gatk CreateSequenceDictionary — the .dict GATK insists on"
rm -f "${REF%.fasta}.dict"
run gatk_run CreateSequenceDictionary -R "$REF"

log "Compress and index the known-sites VCF (BQSR needs it indexed)"
run bash -c "bgzip -f -c '$REFDIR/known_sites.vcf' > '$KNOWN'"
run tabix -f -p vcf "$KNOWN"

log "What the reference actually contains"
run cat "$REF.fai"

cat <<'TXT'

--------------------------------------------------------------------------
Two contigs, 240 kb in total. That is a toy, and deliberately so: the exact
same commands on a human genome would take hours per sample instead of
seconds, and you would learn nothing extra while waiting.

known_sites.vcf is our stand-in for dbSNP: a catalogue of positions already
known to be variable in the population. BQSR needs it in step 5, for a
reason worth understanding before you get there -- see PIPELINE.md.
--------------------------------------------------------------------------
TXT
