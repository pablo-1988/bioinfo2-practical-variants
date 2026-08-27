#!/usr/bin/env bash
#
# run_pipeline.sh — runs every step script in order.
#
# This file contains NO analysis commands of its own. Each step lives in its
# own script under scripts/, and you are meant to run those one at a time,
# reading the output before moving on:
#
#     ./scripts/01_qc_raw.sh
#     ./scripts/02_prepare_reference.sh
#     ./scripts/03_map_bwa.sh
#     ./scripts/04_mark_duplicates.sh
#     ./scripts/05_bqsr.sh
#     ./scripts/06_alignment_qc.sh
#     ./scripts/07_call_variants.sh
#     ./scripts/08_joint_genotyping.sh
#     ./scripts/09_filter_variants.sh
#     ./scripts/10_evaluate.sh
#     python3 scripts/11_report.py
#
# Use this wrapper only to run the whole thing unattended, or to catch up if
# you fell behind:
#
#     ./run_pipeline.sh              # every step
#     ./run_pipeline.sh 03 04 05     # only those steps
#     SAMPLES=sample01 ./run_pipeline.sh
#     THREADS=8 ./run_pipeline.sh

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

STEPS=(
  scripts/01_qc_raw.sh
  scripts/02_prepare_reference.sh
  scripts/03_map_bwa.sh
  scripts/04_mark_duplicates.sh
  scripts/05_bqsr.sh
  scripts/06_alignment_qc.sh
  scripts/07_call_variants.sh
  scripts/08_joint_genotyping.sh
  scripts/09_filter_variants.sh
  scripts/10_evaluate.sh
  scripts/11_report.py
)

[ -d data/raw ] || {
  echo "ERROR: data/raw not found. Generate the data first:" >&2
  echo "       python3 scripts/simulate_reads.py --outdir data" >&2
  exit 1
}

# With no arguments, run everything. With arguments, run only the steps whose
# number was given (e.g. "03" matches scripts/03_map_bwa.sh).
if [ $# -eq 0 ]; then
  selected=("${STEPS[@]}")
else
  selected=()
  for want in "$@"; do
    match=""
    for step in "${STEPS[@]}"; do
      case "$(basename "$step")" in
        "$want"_*) match="$step" ;;
      esac
    done
    [ -n "$match" ] || { echo "ERROR: no step numbered '$want'" >&2; exit 1; }
    selected+=("$match")
  done
fi

for step in "${selected[@]}"; do
  printf "\n\033[1;32m######## %s ########\033[0m\n" "$step"
  case "$step" in
    *.py) python3 "$step" ;;
    *)    bash "$step" ;;
  esac
done

printf "\n\033[1;32mAll requested steps finished. Open results/10_report/report.html\033[0m\n"
