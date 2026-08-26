#!/usr/bin/env bash
#
# config.sh — shared settings and helpers for every step script.
#
# You never run this file directly. Each numbered script sources it, so all
# of them agree on where the data lives and where the results go.
#
# Override anything from the command line, e.g.:
#     SAMPLES=sample01 THREADS=8 ./scripts/03_map_bwa.sh

set -euo pipefail

# Always work from the project root, no matter where the script is called from
cd "$(dirname "${BASH_SOURCE[0]}")/.."

# --- What to process --------------------------------------------------
SAMPLES="${SAMPLES:-sample01 sample02 sample03}"
THREADS="${THREADS:-4}"

# GATK is a Java program; without this it may try to grab far more memory
# than a laptop has. Two gigabytes is plenty for a 240 kb genome.
GATK_JAVA="${GATK_JAVA:--Xmx2g -XX:ParallelGCThreads=2}"
gatk_run() { gatk --java-options "$GATK_JAVA" "$@"; }

# --- Where the data is ------------------------------------------------
RAW="data/raw"
REFDIR="data/reference"
REF="$REFDIR/reference.fasta"
KNOWN="$REFDIR/known_sites.vcf.gz"
TRUTH="data/truth/truth.vcf.gz"

# --- Where the results go ---------------------------------------------
RESULTS="results"
QC_RAW="$RESULTS/01_qc_raw"
MAP="$RESULTS/02_mapping"
DEDUP="$RESULTS/03_dedup"
BQSR="$RESULTS/04_bqsr"
ALNQC="$RESULTS/05_alignment_qc"
GVCF="$RESULTS/06_gvcf"
JOINT="$RESULTS/07_joint_genotyping"
FILT="$RESULTS/08_filtered"
EVAL="$RESULTS/09_evaluation"

# --- Helpers ----------------------------------------------------------

# Print a step header
log() { printf "\n\033[1;34m==> %s\033[0m\n" "$*"; }

# Echo a command in grey, then run it. This is why you can copy any line
# from the terminal output and run it by hand.
run() { printf "    \033[0;90m$ %s\033[0m\n" "$*"; "$@"; }

# Fail early and clearly if a tool is missing
require() {
  for tool in "$@"; do
    command -v "$tool" >/dev/null 2>&1 || {
      echo "ERROR: '$tool' not found. Activate the environment first:" >&2
      echo "       conda activate bioinfo2-variants" >&2
      exit 1
    }
  done
}

# Fail early and clearly if an input file is missing
need_file() {
  for f in "$@"; do
    [ -f "$f" ] || {
      echo "ERROR: missing input file: $f" >&2
      echo "       Did you run the previous step?" >&2
      exit 1
    }
  done
}

# The contigs of the reference, as "-L chr1 -L chr2 ..."
intervals_args() {
  need_file "$REF.fai"
  awk '{printf "-L %s ", $1}' "$REF.fai"
}
