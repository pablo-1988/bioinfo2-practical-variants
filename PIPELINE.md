# Read Mapping and Variant Calling with BWA + GATK

**Bioinformatics II — practical session 2**

This is the guide. The scripts under `scripts/` are the executable version of
what is written here; if you only run the scripts you will produce files, and
if you only read this you will understand nothing you can operate. Do both.

The whole run takes about 5 minutes on a recent laptop with 4 threads. It is meant to be run one
step at a time, reading the output of each before starting the next.

---

## 0. What the question is

You have short reads from three individuals of the same species. You want to
know, for each individual and each position of the genome, **which alleles
they carry**.

That question decomposes into three that are answerable:

1. **Where did each read come from?** — alignment (`bwa mem`)
2. **How much do I trust each base and each read?** — MarkDuplicates, BQSR,
   mapping quality
3. **Given all the reads over a position, what genotype best explains
   them?** — `HaplotypeCaller`, then joint genotyping

Almost everything that goes wrong in variant calling goes wrong in step 2,
and it goes wrong silently. A caller that receives bad evidence does not
crash; it produces confident, well-formatted, wrong answers. That is why
this pipeline has more QC steps than calling steps.

---

## 1. The dataset

Generated locally by `scripts/simulate_reads.py`, deterministically from
seed 42, so everyone's files are identical.

**Reference:** two contigs, `chr1` (180 kb) and `chr2` (60 kb).

**Two traps are built in on purpose, and they fail in opposite directions.**
You are expected to find both before the scripts point them out.

- `chr1` contains **the same 2 kb segment at two positions** (40,000 and
  130,000). Reads from it map equally well to both copies, so bwa reports
  MAPQ 0, HaplotypeCaller discards them, and the region goes silent: no
  calls at all, and the real variants there are missed. **A false-negative
  trap.**
- Every sample's real genome contains a **3 kb paralogue of chr2:20,000–23,000
  with ~1.2% divergence that is absent from the reference** — a segmental
  duplication collapsed in the assembly, which is very common in real
  references. Its reads have nowhere else to go, so they pile onto the single
  reference copy at double depth and their divergent bases look like clean,
  well-supported heterozygous SNPs, identically in all three samples.
  **A false-positive trap**, and one the standard hard filters do not catch.

**Variants:** ~1,010 sites (900 SNPs, 110 indels of 1–11 bp) planted into
the genome, with a diploid genotype per individual, mostly shared across
individuals at realistic allele frequencies, ~12% private to one person.
A full run recovers roughly 1,026 raw sites, of which ~28 per sample are
false and all of them are in one place.

**Samples:**

| | Coverage | Duplicates | What it is for |
|---|---|---|---|
| `sample01` | ~30x | ~2% | The control. This is what "fine" looks like. |
| `sample02` | ~30x | ~25% | Same nominal coverage, over-amplified library. |
| `sample03` | ~12x | ~2% | Clean, but not enough of it. |

Nominal coverage is not information content. Proving that is most of the
point of the exercise.

**The truth:** `data/truth/truth.vcf` contains the real genotypes. **Do not
open it before step 10.** If you look at the answer first you will
unconsciously tune the pipeline towards it, which is the single most common
way real benchmarks get corrupted.

---

## 2. Step by step

### Step 1 — QC of the raw reads

```bash
./scripts/01_qc_raw.sh
```

You did this in the previous practical, so it is brief. What is new is the
decision at the end: **we do not trim.**

That may surprise you. The reasoning:

- `bwa mem` **soft-clips** — it aligns the part of the read that matches and
  marks the rest in the CIGAR (`30S120M`). Adapter and low-quality tail get
  excluded from the alignment without being deleted from the file.
- GATK reads that soft-clip and ignores clipped bases when calling.
- Trimming, in contrast, is irreversible. A trimmed base cannot be recovered
  if it turns out the clip was wrong.
- Worse, aggressive quality trimming **biases** the data: it preferentially
  removes bases that disagree with the reference, because errors and real
  alternate alleles both lower local quality. You lose reference-discordant
  evidence — the very thing you are trying to measure.

Trim when there is a concrete reason (heavy adapter contamination, a known
bad cycle range). Not as a reflex.

### Step 2 — Index the reference

```bash
./scripts/02_prepare_reference.sh
```

Three indexes for three consumers: `bwa index` (FM-index for search),
`samtools faidx` (`.fai`, random access by coordinate), `CreateSequenceDictionary`
(`.dict`, the contig list GATK validates every input against).

GATK compares the `.dict` against the BAM header on every run. If they
disagree — a different reference version, a renamed contig, `chr1` vs `1` —
it refuses to start. That check is a feature, and it will save you one day.

### Step 3 — Align with `bwa mem`

```bash
./scripts/03_map_bwa.sh
```

```bash
bwa mem -t 4 -M -R '@RG\tID:sample01\tSM:sample01\tPL:ILLUMINA\tLB:lib1\tPU:FLOWCELLX.1.sample01' \
    reference.fasta sample01_R1.fastq.gz sample01_R2.fastq.gz \
  | samtools sort -o sample01.sorted.bam -
```

**The read group (`-R`) is not optional.** GATK will not run without it, and
for good reasons rather than bureaucratic ones:

- `SM` is the sample name that appears as the column header in your VCF. Get
  it wrong and every result is mislabelled.
- `PU` (flowcell.lane.barcode) is the unit over which sequencing error
  behaves consistently — BQSR builds one error model per read group.
- `LB` tells MarkDuplicates which reads *could* be duplicates of each other.
  Two reads from different libraries at the same coordinates are not
  duplicates; they are independent confirmation.

**Note the pipe.** The SAM never touches the disk. On real data that saves
hundreds of gigabytes of I/O per sample.

#### What MAPQ means, and why it is the most useful number in the BAM

MAPQ is `-10 log10 P(this placement is wrong)`. It is **not** alignment
quality. A read can match a location with zero mismatches and still get MAPQ
0 — if it matches a second location just as well.

That is exactly what happens in the duplicated segment. Find it yourself:

```bash
samtools view results/02_mapping/sample01.sorted.bam \
  | awk '$5 == 0 {print $3"\t"int($4/10000)*10000}' | sort | uniq -c | sort -rn | head
```

Two windows dominate. Write down the coordinates. You will meet this region
three more times, and each time it will look like a different problem.

### Step 4 — Mark duplicates

```bash
./scripts/04_mark_duplicates.sh
```

A PCR duplicate is **one molecule sequenced twice**. It carries no
independent information about the genome, but the caller counts reads, so it
inflates apparent depth. If the original molecule had a sequencing error, its
duplicates carry that same error — and 5 reads supporting a wrong allele
looks exactly like a real heterozygous site.

MarkDuplicates identifies them by **both alignment coordinates** (start of R1
and start of R2), not by sequence: independent copies of one molecule have
independent errors, so their sequences differ, but their fragment boundaries
are identical. This is also why single-end data is much harder to deduplicate
honestly.

It **marks** (SAM flag `0x400`) rather than deletes. HaplotypeCaller skips
flagged reads; you keep the ability to audit.

`sample02` will show ~25% duplication. Its effective coverage is therefore
around 22x, not 30x — and those 22x are less evenly distributed than
sample01's 30x, because the duplicates pile onto whichever molecules PCR
happened to favour.

### Step 5 — BQSR

```bash
./scripts/05_bqsr.sh
```

The instrument's quality scores are systematically wrong — biased by cycle
number, by sequence context, by the machine's state. HaplotypeCaller's
genotype likelihoods are built directly on `P(error)` from those scores, so
a systematic overestimate of quality becomes a systematic overestimate of
confidence in false variants.

BQSR measures the real error rate empirically:

1. Every base that mismatches the reference is either a sequencing error or a
   real variant.
2. **Mask every position in the known-sites database** (here `known_sites.vcf`,
   our fictitious dbSNP). What remains is dominated by errors.
3. Tabulate observed error rate by reported quality × cycle × dinucleotide
   context.
4. Rewrite each base quality with its empirical value.

Two consequences worth stating out loud:

- Real variants **not** in the database are counted as errors. With a good
  catalogue that is a small bias; with no catalogue it is a disaster. **On a
  species without a variant database, do not run BQSR** — or bootstrap one by
  calling variants, taking the high-confidence set as "known", and iterating.
- BQSR needs a lot of bases to fit its model. On a 240 kb toy genome the
  correction is real but modest. On a human exome it is substantial.

In this dataset the simulator deliberately makes the reported qualities
optimistic in the second half of each read. You should be able to see exactly
that in the recalibration table, in the `Cycle` covariate.

### Step 6 — Alignment QC

```bash
./scripts/06_alignment_qc.sh
```

Coverage is never uniform, and the mean hides everything interesting. Note
`samtools depth -a`: without `-a`, positions with zero coverage are omitted
and your "mean depth" is computed only over the parts that worked.

The script measures depth **twice**: over all reads, and over reads with
MAPQ ≥ 20 — roughly what HaplotypeCaller will actually use. The two are not
the same genome.

In the duplicated segment, total depth looks *fine*, even slightly above
average: reads from both copies pile onto both. Callable depth there is
close to zero. A naive depth plot shows nothing wrong. **Coverage that
exists and coverage that counts are different quantities**, and only the
second one calls variants.

At the other end, one window on chr2 sits at roughly **twice** the mean with
completely normal mapping qualities. Nothing was filtered; there is simply
more sequence in the sample than the reference can account for. Note the
coordinates. Step 8 shows what it does.

### Step 7 — HaplotypeCaller in GVCF mode

```bash
./scripts/07_call_variants.sh
```

HaplotypeCaller does not read the pileup column by column. In each region
that looks active it:

1. collects the reads,
2. builds a de Bruijn graph and extracts candidate **haplotypes** — whole
   local sequences, not isolated positions,
3. realigns every read against every candidate haplotype with a pair-HMM,
   producing `P(read | haplotype)`,
4. marginalises to per-allele likelihoods and genotypes the site.

This is why it handles indels, and why two nearby variants on the same
molecule are called consistently instead of as two independent, mutually
confusing observations. It is also why it is the slow step: the pair-HMM is
quadratic in read × haplotype length.

`-ERC GVCF` changes the *output*, not the calling. It emits a record for
every position, including reference blocks with a confidence:

```
no line in a VCF          =  "I found nothing"        (could mean anything)
GVCF <NON_REF> block      =  "I looked, 32x, reference, GQ 99"
```

The difference is what makes step 8 possible.

### Step 8 — Joint genotyping

```bash
./scripts/08_joint_genotyping.sh
```

`GenomicsDBImport` then `GenotypeGVCFs`.

Genotyping the cohort together **changes the calls**, it is not a merge. At a
site where sample01 and sample02 are confidently heterozygous, the prior
probability that the site is polymorphic rises, and marginal evidence in
sample03 (say 3 reads) that would be discarded in isolation becomes a
defensible call. Conversely, an artefact appearing in one sample only, at a
site where the others are cleanly reference, gets less support.

This is the standard architecture for cohorts precisely because it scales:
each new sample costs one GVCF, and the joint step is re-run — you never
recall the old samples.

Now compare the two suspect regions in the raw VCF. The duplicated segment
is nearly **empty** — the caller made no statements, because MAPQ 0 reads sit
below its threshold. The over-covered chr2 window is the opposite: a dense
cluster of confident heterozygous calls, present in all three samples, with
balanced allele depths and MQ 60.

Nothing in that VCF looks wrong. **Reproducibility across samples is not
evidence of truth** — a systematic artefact reproduces perfectly.

### Step 9 — Hard filtering

```bash
./scripts/09_filter_variants.sh
```

Best practice for large human cohorts is **VQSR** (`VariantRecalibrator`),
which learns the boundary between true and false calls from known-variant
training sets. It needs tens of thousands of variants to fit its Gaussian
mixture model. We have ~1,000, so it would fail — and GATK's own
documentation prescribes hard filters in exactly this situation: small
cohorts, non-model organisms, targeted panels.

SNPs and indels are filtered separately, because they fail differently.

| Annotation | Filter | What it detects |
|---|---|---|
| `QD` | `< 2.0` | QUAL normalised by depth. High QUAL that is only high because depth is 500x is not evidence. |
| `QUAL` | `< 30` | Weak overall evidence. |
| `FS` | `> 60` (SNP) | Strand bias, Phred-scaled Fisher test. A real allele appears on both strands. |
| `SOR` | `> 3.0` | Strand bias, better behaved at high depth. |
| `MQ` | `< 40` | **RMS mapping quality.** This is the one that catches the duplicated segment. |
| `MQRankSum` | `< -12.5` | Do ALT reads have systematically worse MAPQ than REF reads? |
| `ReadPosRankSum` | `< -8` | Does the ALT only ever appear near read ends — where errors live? |

Indels get looser thresholds (`FS > 200`, `ReadPosRankSum < -20`) because
indel alignment is inherently noisier.

`VariantFiltration` **marks**; the sites remain in the file with a reason in
the `FILTER` column. `bcftools view -f PASS` gives the clean subset when you
need it. Never delete the rejected calls — the day you doubt your pipeline,
they are the evidence.

### Step 10 — Evaluate against the truth

```bash
./scripts/10_evaluate.sh
```

The step you can only do on simulated data. `gatk Concordance` gives, per
sample and per variant type:

- **sensitivity (recall)** — of the real variants, how many did you find?
- **precision** — of the variants you called, how many are real?

The script also runs the comparison on the *unfiltered* calls, so you can
price the filters: exactly how many true calls did they cost, and how many
false ones did they buy?

**Both files are normalised with `bcftools norm` before comparison, and that
is not a formality.** The same indel can be written at several positions —
a deletion of `CT` inside `CTCTCT` has three equally valid representations.
If truth and calls choose differently, a naive comparison charges you one
false positive *and* one false negative for a variant you got exactly right.
In this dataset skipping normalisation destroys about a fifth of the indel
accuracy, entirely fictitiously. It is one of the most common ways published
variant comparisons are quietly wrong.

You should end up around **0.96–0.99 SNP sensitivity** and **~0.95 SNP
precision**, with indels close to perfect once normalised. Every false
positive will be in the collapsed-duplication window on chr2; nearly every
false negative will be in the duplicated segment on chr1. The errors are not
random, and knowing *where* your pipeline fails is more useful than knowing
*how often*.

### Step 11 — Build the visual report

```bash
python3 scripts/11_report.py
open results/10_report/report.html
```

Reads everything the previous steps wrote and produces one self-contained
HTML page: coverage along the genome with the two suspect regions shaded, the
coordinate of every false positive and false negative, the filtering
trade-off as raw → filtered, and the BQSR before/after calibration curve
(which also replaces the `AnalyzeCovariates` PDF when R refuses to build it).

Standard library only — no plotting package to install. The figures are
inline SVG, so screenshot a panel or copy the `<svg>` block into your report.

**A figure is for finding things, not for proving them.** Use the coverage
panel to locate the anomalies, then go back to the BAM and the VCF and get
the numbers. A report that shows the picture without the coordinates and
counts behind it has not done the work.

---

## 3. What to hand in

A report — text and tables, no more than four pages — containing:

1. **A QC table** for the three samples: read count, duplication rate, mean
   depth, effective (deduplicated) depth, fraction below 10x, insert size
   mean and SD.
2. **Both anomalous regions**, with coordinates, and how you found them: the
   one that loses coverage to MAPQ 0, and the one with excess depth. Explain
   in one paragraph why bwa is *right* to report MAPQ 0 in the first, and
   why nothing in the second looks wrong from inside the VCF.
3. **The BQSR evidence**: one before/after comparison showing reported vs
   empirical quality, and which covariate carried most of the correction.
4. **The concordance table**, filtered and unfiltered, for all three samples,
   SNPs and indels separately.
5. **Answers** to the questions printed at the end of each script.
6. **One paragraph** on this: `sample02` and `sample03` both perform worse
   than `sample01`, for different reasons. If you could only fix one of the
   two libraries, which would you fix, and what would you tell the lab to do
   differently?
7. **One paragraph** on the false positives: they are identical in all three
   samples and would survive deeper sequencing. What would actually remove
   them, and what does that imply about trusting a reference genome?

---

## 4. Things that are true here and false in real life

The dataset is honest about being a toy. Know where it lies:

- **~100% of reads map.** Real libraries carry adapter dimers, contamination,
  and sequence absent from the reference. 85–95% is normal; 100% would mean
  something is wrong with your QC.
- **The reference is the true genome.** Real references are a consensus of
  other individuals, with errors, gaps and missing structural variation.
  Reads from a real sample carry variation the reference simply lacks.
- **The errors are substitutions with a clean quality model.** Real errors
  are context-dependent, cluster in homopolymers, and include chimeric
  fragments and index hopping.
- **240 kb, three samples, five minutes.** A human trio is ~1,000× more data,
  needs a cluster, and would use interval scattering, VQSR, and days.
- **One duplicated segment and one collapsed paralogue.** ~5% of the human
  genome is segmentally duplicated, much of it misrepresented in the
  reference, and short-read calling there is unreliable by construction. Real
  callsets deal with this using blacklists of problematic regions
  (ENCODE/GIAB), not with per-site annotations.

The commands do not change. The scale, the failure modes and the amount of
scepticism required all do.

---

## 5. Reference — the commands, without the prose

```bash
# reference
bwa index ref.fasta
samtools faidx ref.fasta
gatk CreateSequenceDictionary -R ref.fasta

# mapping
bwa mem -t 4 -M -R '@RG\tID:s\tSM:s\tPL:ILLUMINA\tLB:lib1\tPU:fc.1.s' \
    ref.fasta s_R1.fq.gz s_R2.fq.gz | samtools sort -o s.bam -
samtools index s.bam

# duplicates
gatk MarkDuplicates -I s.bam -O s.dedup.bam -M s.metrics.txt --CREATE_INDEX true

# BQSR
gatk BaseRecalibrator -R ref.fasta -I s.dedup.bam --known-sites known.vcf.gz -O s.table
gatk ApplyBQSR       -R ref.fasta -I s.dedup.bam --bqsr-recal-file s.table -O s.recal.bam

# calling
gatk HaplotypeCaller -R ref.fasta -I s.recal.bam -O s.g.vcf.gz -ERC GVCF
gatk GenomicsDBImport -V a.g.vcf.gz -V b.g.vcf.gz --genomicsdb-workspace-path db -L chr1
gatk GenotypeGVCFs   -R ref.fasta -V gendb://db -O cohort.vcf.gz

# filtering
gatk SelectVariants     -V cohort.vcf.gz --select-type-to-include SNP -O snps.vcf.gz
gatk VariantFiltration  -V snps.vcf.gz --filter-name QD2 --filter-expression "QD < 2.0" -O snps.filt.vcf.gz
gatk MergeVcfs -I snps.filt.vcf.gz -I indels.filt.vcf.gz -O cohort.filtered.vcf.gz

# evaluation
gatk Concordance -R ref.fasta -eval calls.vcf.gz --truth truth.vcf.gz --summary out.tsv
```
