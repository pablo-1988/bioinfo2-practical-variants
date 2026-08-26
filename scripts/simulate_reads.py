#!/usr/bin/env python3
"""
simulate_reads.py — Generate the fictitious dataset for the mapping and
variant calling practical.

No external dependencies (Python 3.8+ standard library only).

It creates:

  * a reference genome        chr1 (180 kb) + chr2 (60 kb), haploid FASTA
  * a truth VCF               the variants that were actually planted, with
                              the real genotype of each of the three samples
  * a "known sites" VCF       a population database, like dbSNP, used by BQSR
  * paired-end Illumina reads 2 x 150 bp, for three fictitious individuals

The three samples are not equivalent, and that is the point:

  sample01  ~30x, clean library                 -> the well-behaved control
  sample02  ~30x, but ~25% PCR duplicates       -> MarkDuplicates earns its keep
  sample03  ~12x, clean library                 -> heterozygotes start to be
                                                   missed; low GQ everywhere

Two traps are built into the genome on purpose, and they fail in opposite
directions:

  * chr1 carries two identical copies of a 2 kb segment. Reads from them map
    equally well to both, so bwa reports MAPQ 0, HaplotypeCaller discards
    them, and the region goes SILENT: apparent low coverage, no calls, real
    variants missed. A false-negative trap.

  * every sample's genome contains a 3 kb paralogue of chr2:20,000-23,000
    with ~1.2% divergence, which is ABSENT from the reference -- a segmental
    duplication collapsed in the assembly, which is extremely common in real
    references. Its reads have nowhere else to go, so they pile onto the one
    reference copy at double depth and their divergent bases look like clean
    heterozygous SNPs, in every sample, reproducibly. A false-positive trap.

Learning to distrust both regions is most of the exercise.

Usage:
    python3 scripts/simulate_reads.py --outdir data
"""

import argparse
import gzip
import math
import os
import random

# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

# Real Illumina TruSeq adapters, read through in short fragments
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Where the duplicated 2 kb segment is planted on chr1
SEG_DUP_LEN = 2_000
SEG_DUP_POS = (40_000, 130_000)   # start of each copy, 0-based

# The paralogue that the reference is missing: a copy of this slice of chr2,
# present in every sample's real genome, diverged by ~1.2%
PARALOG_CHROM = "chr2"
PARALOG_START = 20_000            # 0-based
PARALOG_LEN = 3_000
PARALOG_DIVERGENCE = 0.012


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def phred(q):
    """Phred+33 character for an integer quality."""
    return chr(min(max(int(q), 2), 41) + 33)


def qual_to_prob(q):
    return 10 ** (-q / 10.0)


# ----------------------------------------------------------------------
# Reference genome
# ----------------------------------------------------------------------

def random_seq(rng, length, gc=0.41):
    """Random sequence with a given GC content."""
    at, gcb = "AT", "GC"
    return "".join(rng.choice(gcb) if rng.random() < gc else rng.choice(at)
                   for _ in range(length))


def make_chr1(rng, length=180_000):
    """chr1, with the same 2 kb segment planted at two distant positions."""
    seq = list(random_seq(rng, length))
    segment = random_seq(rng, SEG_DUP_LEN, gc=0.52)
    for start in SEG_DUP_POS:
        seq[start:start + SEG_DUP_LEN] = list(segment)
    return "".join(seq)


def write_fasta(path, records, width=60):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(">%s\n" % name)
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


# ----------------------------------------------------------------------
# Variants
# ----------------------------------------------------------------------

class Variant:
    """One planted variant, with the genotype of every sample."""

    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos          # 1-based, VCF convention
        self.ref = ref
        self.alt = alt
        self.genotypes = {}     # sample -> (allele_hap1, allele_hap2)

    @property
    def is_indel(self):
        return len(self.ref) != len(self.alt)


def plant_variants(rng, genome, samples, n_snps=900, n_indels=110,
                   min_gap=60):
    """
    Choose variant sites across the genome and assign each sample a diploid
    genotype, the way a real population would: most variants are shared and
    segregate at some allele frequency, a few are private to one individual.
    """
    lengths = {name: len(seq) for name, seq in genome}
    seqs = dict(genome)
    total = sum(lengths.values())

    taken = {name: [] for name in lengths}       # occupied intervals per chrom
    variants = []

    def free(chrom, start, end):
        for s, e in taken[chrom]:
            if start < e + min_gap and s - min_gap < end:
                return False
        return True

    def pick_site(size):
        for _ in range(200):
            chrom = rng.choices(list(lengths), weights=list(lengths.values()))[0]
            # Keep away from the very ends so reads always flank the variant
            pos0 = rng.randrange(1_000, lengths[chrom] - 1_000 - size)
            if "N" in seqs[chrom][pos0:pos0 + size + 1]:
                continue
            if free(chrom, pos0, pos0 + size):
                taken[chrom].append((pos0, pos0 + size))
                return chrom, pos0
        return None

    for _ in range(n_snps):
        site = pick_site(1)
        if site is None:
            continue
        chrom, pos0 = site
        ref = seqs[chrom][pos0]
        alt = rng.choice(BASES.replace(ref, ""))
        variants.append(Variant(chrom, pos0 + 1, ref, alt))

    for _ in range(n_indels):
        length = rng.choice([1, 1, 2, 2, 3, 4, 5, 6, 8, 11])
        deletion = rng.random() < 0.5
        site = pick_site(length + 1)
        if site is None:
            continue
        chrom, pos0 = site
        anchor = seqs[chrom][pos0]
        if deletion:
            ref = seqs[chrom][pos0:pos0 + length + 1]
            alt = anchor
        else:
            ref = anchor
            alt = anchor + random_seq(rng, length, gc=0.45)
        variants.append(Variant(chrom, pos0 + 1, ref, alt))

    variants.sort(key=lambda v: (v.chrom, v.pos))

    # --- genotypes ----------------------------------------------------
    # An allele frequency per site, then Hardy-Weinberg per individual.
    # Every site is forced to be non-reference in at least one sample,
    # otherwise it is not a variant at all and only confuses the count.
    for v in variants:
        private = rng.random() < 0.12
        if private:
            carrier = rng.choice(samples)
            for s in samples:
                if s == carrier:
                    v.genotypes[s] = (0, 1) if rng.random() < 0.85 else (1, 1)
                else:
                    v.genotypes[s] = (0, 0)
        else:
            af = rng.uniform(0.15, 0.85)
            while True:
                for s in samples:
                    v.genotypes[s] = (int(rng.random() < af),
                                      int(rng.random() < af))
                if any(sum(v.genotypes[s]) > 0 for s in samples):
                    break

    return variants


def make_paralog(rng, genome):
    """
    A diverged copy of a chr2 segment that exists in the samples but not in
    the reference. Reads from it have no correct place to map, so they land
    on the single reference copy and masquerade as heterozygous variation.
    """
    seq = dict(genome)[PARALOG_CHROM][PARALOG_START:PARALOG_START + PARALOG_LEN]
    out = []
    for base in seq:
        if rng.random() < PARALOG_DIVERGENCE:
            out.append(rng.choice(BASES.replace(base, "")))
        else:
            out.append(base)
    return "".join(out)


def build_haplotypes(genome, variants, sample, extra=None):
    """
    Apply the sample's genotype to the reference and return its two
    haplotype sequences. Reads are then drawn from these, not from the
    reference — which is exactly how a real genome differs from the one
    you align against.
    """
    haps = []
    for hap_index in (0, 1):
        out = {}
        for chrom, seq in genome:
            pieces, cursor = [], 0
            for v in variants:
                if v.chrom != chrom:
                    continue
                if v.genotypes[sample][hap_index] == 0:
                    continue
                start = v.pos - 1
                pieces.append(seq[cursor:start])
                pieces.append(v.alt)
                cursor = start + len(v.ref)
            pieces.append(seq[cursor:])
            out[chrom] = "".join(pieces)
        if extra:
            out.update(extra)
        haps.append(out)
    return haps


def write_vcf(path, contigs, variants, samples=None, source="simulate_reads.py"):
    """Write a plain (uncompressed) VCF. samples=None -> sites-only."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=%s\n" % source)
        for name, length in contigs:
            fh.write("##contig=<ID=%s,length=%d>\n" % (name, length))
        if samples:
            fh.write('##FORMAT=<ID=GT,Number=1,Type=String,'
                     'Description="Genotype">\n')
        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        if samples:
            header += ["FORMAT"] + list(samples)
        fh.write("\t".join(header) + "\n")

        for v in variants:
            row = [v.chrom, str(v.pos), ".", v.ref, v.alt, ".", ".", "."]
            if samples:
                row.append("GT")
                for s in samples:
                    a, b = v.genotypes[s]
                    row.append("%d/%d" % (a, b))
            fh.write("\t".join(row) + "\n")


# ----------------------------------------------------------------------
# Illumina simulation
# ----------------------------------------------------------------------

def quality_profile(rng, read_len):
    """
    Per-base REPORTED qualities: high at the start, decaying towards the 3'
    end, as a real run does.

    These are the numbers the instrument claims. The actual error rate used
    below is deliberately worse than claimed in the last third of the read
    (see miscalibration in mutate_by_quality) — that mismatch between claimed
    and observed quality is precisely what BQSR is built to correct.
    """
    drop = rng.uniform(5, 10)
    quals = []
    for i in range(read_len):
        frac = i / float(read_len)
        base_q = 37 - drop * (frac ** 2.0)
        quals.append(max(2, int(rng.gauss(base_q, 2.0))))
    return quals


def mutate_by_quality(rng, seq, quals):
    """
    Introduce substitution errors at the rate implied by the quality score —
    but inflate the true error rate in the second half of the read, so the
    reported qualities there are optimistic.
    """
    out = []
    n = len(seq)
    for i, (base, q) in enumerate(zip(seq, quals)):
        p = qual_to_prob(q)
        if i > n * 0.55:
            p *= 3.0                      # reported Q is ~5 points too generous
        if base != "N" and rng.random() < p:
            out.append(rng.choice(BASES.replace(base, "")))
        else:
            out.append(base)
    return "".join(out)


def make_fragment(rng, haps, weights, names, read_len, frag_mean, frag_sd,
                  short_frag_rate):
    """Draw one fragment from a random haplotype. Returns (r1_seq, r2_seq)."""
    hap = haps[rng.randrange(2)]
    chrom = rng.choices(names, weights=weights, k=1)[0]
    seq = hap[chrom]

    if rng.random() < short_frag_rate:
        frag_len = rng.randint(70, read_len - 10)
    else:
        frag_len = max(read_len + 10, int(rng.gauss(frag_mean, frag_sd)))
    if frag_len >= len(seq):
        frag_len = read_len + 10

    start = rng.randrange(0, len(seq) - frag_len)
    frag = seq[start:start + frag_len]
    if rng.random() < 0.5:
        frag = revcomp(frag)

    r1 = frag[:read_len]
    r2 = revcomp(frag)[:read_len]

    # Adapter read-through when the fragment is shorter than the read
    if len(r1) < read_len:
        r1 = (r1 + ADAPTER_R1 + random_seq(rng, read_len))[:read_len]
        r2 = (r2 + ADAPTER_R2 + random_seq(rng, read_len))[:read_len]
    return r1, r2


def simulate_sample(rng, haps, n_pairs, read_len=150,
                    frag_mean=380, frag_sd=70, short_frag_rate=0.02,
                    dup_rate=0.0):
    """
    Yield (name, r1, q1, r2, q2) for one library.

    dup_rate: fraction of the output that comes from re-sequencing a fragment
    already seen. A PCR duplicate is the SAME molecule read twice: same
    coordinates, same insert size, but its own independent sequencing errors.
    That is what MarkDuplicates finds, and why it uses coordinates rather than
    sequence identity.
    """
    # Note: drawn from the SAMPLE's sequences, which include the paralogue.
    # The reference does not have it; that is the whole point.
    names = list(haps[0])
    weights = [len(haps[0][n]) for n in names]

    pool = []          # recently seen fragments, available for duplication
    i = 0
    while i < n_pairs:
        if pool and rng.random() < dup_rate:
            r1_src, r2_src = rng.choice(pool)
        else:
            r1_src, r2_src = make_fragment(rng, haps, weights, names, read_len,
                                           frag_mean, frag_sd, short_frag_rate)
            pool.append((r1_src, r2_src))
            if len(pool) > 4000:
                pool.pop(rng.randrange(len(pool)))

        q1 = quality_profile(rng, len(r1_src))
        q2 = quality_profile(rng, len(r2_src))
        r1 = mutate_by_quality(rng, r1_src, q1)
        r2 = mutate_by_quality(rng, r2_src, q2)

        # Illumina-style name: instrument:run:flowcell:lane:tile:x:y
        name = "SIM01:1:FLOWCELLX:1:%d:%d:%d" % (
            rng.randrange(1101, 1112), rng.randrange(1000, 30000),
            rng.randrange(1000, 30000))
        yield (name,
               r1, "".join(phred(q) for q in q1),
               r2, "".join(phred(q) for q in q2))
        i += 1


def write_fastq(path, records):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt", compresslevel=6) as fh:
        for name, seq, qual in records:
            fh.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))


# ----------------------------------------------------------------------

SAMPLE_SPECS = [
    # name        coverage  dup_rate  note
    ("sample01",  30,       0.02,     "clean library, ~30x"),
    ("sample02",  30,       0.25,     "~30x, heavy PCR duplication"),
    ("sample03",  12,       0.02,     "clean library, only ~12x"),
]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="data", help="output directory")
    ap.add_argument("--seed", type=int, default=42, help="random seed")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    samples = [s[0] for s in SAMPLE_SPECS]

    print("Building the reference genome ...")
    chr1 = make_chr1(rng, 180_000)
    chr2 = random_seq(rng, 60_000, gc=0.45)
    genome = [("chr1", chr1), ("chr2", chr2)]
    contigs = [(n, len(s)) for n, s in genome]
    genome_size = sum(len(s) for _, s in genome)

    write_fasta(os.path.join(args.outdir, "reference", "reference.fasta"), genome)
    print("  chr1: %d bp (two identical copies of a %d bp segment at %s)"
          % (len(chr1), SEG_DUP_LEN, " and ".join(str(p) for p in SEG_DUP_POS)))
    print("  chr2: %d bp" % len(chr2))

    print("\nPlanting variants ...")
    variants = plant_variants(rng, genome, samples)
    n_snp = sum(1 for v in variants if not v.is_indel)
    n_ind = sum(1 for v in variants if v.is_indel)
    print("  %d sites: %d SNPs, %d indels" % (len(variants), n_snp, n_ind))
    for s in samples:
        het = sum(1 for v in variants if sum(v.genotypes[s]) == 1)
        hom = sum(1 for v in variants if sum(v.genotypes[s]) == 2)
        print("    %s: %d het, %d hom-alt" % (s, het, hom))

    # The truth. Students must not look at this until the evaluation step.
    write_vcf(os.path.join(args.outdir, "truth", "truth.vcf"),
              contigs, variants, samples=samples)

    # The "known sites" database that BQSR uses. It is deliberately
    # INCOMPLETE (85% of the real sites) and contains some sites that are
    # not variable in these three individuals — like any real dbSNP.
    known = [v for v in variants if rng.random() < 0.85]
    decoys = plant_variants(rng, genome, samples, n_snps=200, n_indels=0)
    occupied = {(v.chrom, v.pos) for v in variants}
    decoys = [d for d in decoys
              if not any((d.chrom, d.pos + off) in occupied
                         for off in range(-30, 31))]
    known = sorted(known + decoys, key=lambda v: (v.chrom, v.pos))
    write_vcf(os.path.join(args.outdir, "reference", "known_sites.vcf"),
              contigs, known, samples=None, source="fictitious dbSNP")
    print("  known_sites.vcf: %d sites (85%% of the truth + 200 decoys)"
          % len(known))

    paralog = make_paralog(rng, genome)
    print("  paralogue of %s:%d-%d (%d bp, ~%.1f%% divergent) is present in"
          % (PARALOG_CHROM, PARALOG_START + 1, PARALOG_START + PARALOG_LEN,
             PARALOG_LEN, 100 * PARALOG_DIVERGENCE))
    print("  every sample but ABSENT from the reference.")

    read_len = 150
    for name, coverage, dup_rate, note in SAMPLE_SPECS:
        n_pairs = int(genome_size * coverage / (2.0 * read_len))
        print("\nSimulating %s (%s) ..." % (name, note))
        haps = build_haplotypes(genome, variants, name,
                                extra={"paralog": paralog})
        r1_recs, r2_recs = [], []
        for rname, s1, q1, s2, q2 in simulate_sample(
                rng, haps, n_pairs, read_len=read_len,
                dup_rate=dup_rate):
            r1_recs.append((rname + " 1:N:0:ATCACG", s1, q1))
            r2_recs.append((rname + " 2:N:0:ATCACG", s2, q2))

        # R1 and R2 must stay in the same order or bwa will pair them wrongly
        order = list(range(len(r1_recs)))
        rng.shuffle(order)
        r1_recs = [r1_recs[i] for i in order]
        r2_recs = [r2_recs[i] for i in order]

        ill = os.path.join(args.outdir, "raw")
        write_fastq(os.path.join(ill, "%s_R1.fastq.gz" % name), r1_recs)
        write_fastq(os.path.join(ill, "%s_R2.fastq.gz" % name), r2_recs)
        print("  %d read pairs (~%dx)" % (len(r1_recs), coverage))

    print("\nDone. Files written under %s/" % args.outdir)
    print("Do NOT open data/truth/truth.vcf before step 10.")


if __name__ == "__main__":
    main()
