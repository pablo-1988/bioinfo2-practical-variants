# Bioinformatics II — Read Mapping and Variant Calling (BWA + GATK)

Hands-on practical: align Illumina paired-end reads to a reference genome,
process the alignments the way GATK Best Practices require, call germline
variants with HaplotypeCaller, filter them, and — because the data is
simulated — measure exactly how right you were.

The full teaching guide is in **[PIPELINE.md](PIPELINE.md)** — read it, the
scripts below are just the executable version of it.

---

# ⚠️ BEFORE THE SESSION — do this at home, not in class

**Installing the software takes 10–30 minutes and downloads a few GB.** If
thirty of us do that simultaneously on the classroom wifi, we lose half the
session to a progress bar. Come to class with steps A–C already working.

Budget an evening, not five minutes. If something breaks, you want it to
break while you still have time to ask.

### A. Install the tools

You need **conda**. If you do not have it, install
[Miniforge](https://github.com/conda-forge/miniforge) first.

```bash
git clone https://github.com/pablo-1988/bioinfo2-practical-variants.git
cd bioinfo2-practical-variants

conda env create -f environment.yml
conda activate bioinfo2-variants
```

Unlike the assembly practical, this one does **not** need the
`CONDA_SUBDIR=osx-64` trick: every tool here (including GATK, which is Java)
has an Apple Silicon build. If your solver still gets stuck, that prefix
remains a valid fallback.

**Windows users: start early.** These tools do not run natively on Windows.
You must install **WSL2** (Windows Subsystem for Linux) with Ubuntu, then
install Miniforge *inside* WSL2 and run everything there. WSL2 requires a
reboot and sometimes a BIOS change, so do not leave this for the night
before. Guide: <https://learn.microsoft.com/windows/wsl/install>

### B. Check that it worked

```bash
conda activate bioinfo2-variants

bwa 2>&1 | head -3        # prints its usage and a version
samtools --version | head -1
bcftools --version | head -1
gatk --version | head -2
fastqc --version && multiqc --version && bgzip --version | head -1
```

**No "command not found", and GATK prints a version.** That is the whole
test. GATK is the one that fails most often — it needs a working Java, and
the conda environment provides it, so run it *inside* the activated
environment.

### C. Generate the data

```bash
python3 scripts/simulate_reads.py --outdir data
```

Takes about a minute and writes ~16 MB. Nothing is downloaded — reads are
generated on your machine from a fixed random seed, so everyone's files are
byte-for-byte identical and every number in your report is comparable to
your classmates'.

### If something fails

Do not show up with it broken and unsaid. Send the **exact error message**
(copy the text, not a photo of the screen) ahead of time.

| What you see | What it means |
|---|---|
| `conda: command not found` | Conda is not installed, or you need to open a new terminal after installing it. |
| `A USER ERROR has occurred: ... Fasta dict file ... does not exist` | You skipped step 2. Run `./scripts/02_prepare_reference.sh`. |
| `SAM/BAM/CRAM file ... has no read groups` | You ran `bwa mem` without `-R`. Redo step 3 — GATK genuinely cannot proceed. |
| `Error: Unable to access jarfile` / Java errors from `gatk` | The environment is not activated, or another Java is shadowing it. `conda activate bioinfo2-variants` and retry. |
| `AnalyzeCovariates` fails with an R error | **Not fatal.** Only the before/after PDF is lost; the recalibration tables it plots are already written and contain the same information. |
| GATK warns `Cannot use GATK jar ... spark` | Harmless. Ignore. |
| `WARN NativeLibraryLoader - Unable to load libgkl_compression.dylib ... incompatible architecture` on Apple Silicon | **Harmless.** GATK ships Intel-only native libraries and falls back to the Java implementations. Slightly slower, identical results. |
| `WARNING: A restricted method in java.lang.System has been called` | Harmless Java 17+ noise from GATK. |

---

## Quick start

```bash
# 1. Install and activate
conda env create -f environment.yml
conda activate bioinfo2-variants

# 2. Generate the simulated dataset (~16 MB, no downloads needed)
python3 scripts/simulate_reads.py --outdir data

# 3. Run the steps one at a time, reading the output of each before moving on
./scripts/01_qc_raw.sh              # FastQC + MultiQC on the raw reads
./scripts/02_prepare_reference.sh   # bwa index, faidx, .dict, index known sites
./scripts/03_map_bwa.sh             # bwa mem  ->  sorted, indexed BAM
./scripts/04_mark_duplicates.sh     # gatk MarkDuplicates
./scripts/05_bqsr.sh                # BaseRecalibrator + ApplyBQSR
./scripts/06_alignment_qc.sh        # flagstat, depth, insert size, MultiQC
./scripts/07_call_variants.sh       # HaplotypeCaller -ERC GVCF   <- the slow one
./scripts/08_joint_genotyping.sh    # GenomicsDBImport + GenotypeGVCFs
./scripts/09_filter_variants.sh     # hard filtering (SNPs and indels apart)
./scripts/10_evaluate.sh            # concordance against the known truth
python3 scripts/11_report.py        # one HTML report with every figure

# ... or run them all unattended (~5 min for three samples, 4 threads)
./run_pipeline.sh
./run_pipeline.sh 07 08 09          # or just some of them, by number
```

Each script prints every command before running it, so you can copy any line
out of the terminal and run it by hand. Each one also ends with what to open
and which questions to answer.

Restrict to one sample or change the thread count:

```bash
SAMPLES=sample01 THREADS=8 ./run_pipeline.sh
```

## Layout

```
.
├── PIPELINE.md              # the practical: theory, commands, questions
├── INSTRUCTOR_KEY.md        # answers to every question — do not hand out
├── environment.yml          # conda environment
├── run_pipeline.sh          # wrapper that runs the step scripts in order
├── scripts/
│   ├── simulate_reads.py         # generates the dataset (stdlib only)
│   ├── config.sh                 # shared paths and helpers, sourced by all steps
│   ├── 01_qc_raw.sh              # ┐
│   ├── 02_prepare_reference.sh   # │
│   ├── 03_map_bwa.sh             # │
│   ├── 04_mark_duplicates.sh     # │
│   ├── 05_bqsr.sh                # │ one step per script,
│   ├── 06_alignment_qc.sh        # │ each runnable on its own
│   ├── 07_call_variants.sh       # │
│   ├── 08_joint_genotyping.sh    # │
│   ├── 09_filter_variants.sh     # │
│   ├── 10_evaluate.sh            # ┘
│   └── 11_report.py              # builds the visual report (stdlib only)
├── data/                    # created by the script (not in git)
│   ├── reference/           # reference.fasta + known_sites.vcf (the "dbSNP")
│   ├── raw/                 # sample0{1,2,3}_R{1,2}.fastq.gz
│   └── truth/               # the real genotypes — DO NOT OPEN before step 10
└── results/                 # created by the pipeline
    ├── 01_qc_raw/  02_mapping/  03_dedup/  04_bqsr/  05_alignment_qc/
    ├── 06_gvcf/  07_joint_genotyping/  08_filtered/  09_evaluation/
    └── 10_report/report.html   # <- the figures, open this in a browser
```

## The dataset

Three fictitious individuals of one species, sharing a reference of **chr1
(180 kb) + chr2 (60 kb)**, into which ~900 SNPs and ~110 indels have been
planted with diploid genotypes per individual.

Two traps are built into that genome on purpose, and they fail in opposite
directions:

- `chr1` carries **the same 2 kb segment at two positions** (40,000 and
  130,000). Reads map equally well to both, so bwa reports MAPQ 0, the caller
  discards them, and the region goes silent — real variants there are simply
  missed. **A false-negative trap.**
- Every sample's genome contains a **3 kb paralogue of chr2:20,000–23,000
  that the reference is missing** (a collapsed segmental duplication — very
  common in real references). Its reads pile onto the single reference copy
  at double depth and look like clean heterozygous SNPs, identically in all
  three samples. **A false-positive trap**, and the standard hard filters do
  not catch it.

Finding both before the scripts point them out is part of the exercise.

| | Coverage | Duplicates | Why it is here |
|---|---|---|---|
| `sample01` | ~30x | ~2% | The control — what "fine" looks like |
| `sample02` | ~30x | ~25% | Same nominal coverage, over-amplified library |
| `sample03` | ~12x | ~2% | Clean, but not enough of it |

Nominal coverage is not information content. Proving that, with numbers from
your own run, is most of the point.

The data is regenerated deterministically from `--seed 42`, so everyone in
the class gets identical files.

## The figures

`scripts/11_report.py` turns everything the pipeline wrote into a single
self-contained HTML page — coverage along the genome, the position of every
wrong call, what filtering bought and cost, and whether BQSR actually made
the base qualities honest. It needs no plotting library (the SVG is written
by hand from the standard library), follows your system's light/dark
setting, and the figures are inline SVG so you can screenshot a panel or copy
the `<svg>` block straight into your report.

```bash
python3 scripts/11_report.py
open results/10_report/report.html      # xdg-open on Linux
```

Use the figures to *find* things, not to replace the numbers: the report also
prints the concordance table, and your write-up needs both.

## What you hand in

A report of at most four pages: the QC table, the MAPQ 0 region with
coordinates, the BQSR before/after evidence, the concordance table (filtered
and unfiltered, all three samples, SNPs and indels apart), and the answers to
the questions printed at the end of each script. Section 3 of
[PIPELINE.md](PIPELINE.md) has the full list.
