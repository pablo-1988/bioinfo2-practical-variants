#!/usr/bin/env python3
"""
11_report.py — turn the pipeline's output into one visual HTML report.

    python3 scripts/11_report.py

No external dependencies (Python 3.8+ standard library only) — the figures are
plain SVG written by hand, so this runs in the same environment as everything
else and needs no plotting library.

It reads what the earlier steps already wrote:

    results/03_dedup/<sample>.metrics.txt          duplicate rate
    results/04_bqsr/<sample>.recal.table           BQSR, before
    results/04_bqsr/<sample>.after.table           BQSR, after
    results/05_alignment_qc/<sample>.depth_*.tsv   coverage
    results/09_evaluation/<sample>.concordance*.tsv  accuracy
    results/09_evaluation/<sample>.{truth,calls}.vcf.gz  where the errors are

and writes:

    results/10_report/report.html

Open that file in a browser. Every figure is inline SVG, so you can screenshot
any panel for your report, or open the file and copy the SVG out.
"""

import argparse
import gzip
import math
import os
import re
from collections import defaultdict

# ----------------------------------------------------------------------
# The two regions the dataset was built around. Keep in sync with
# scripts/simulate_reads.py — these are the same constants.
# ----------------------------------------------------------------------

REGIONS = [
    ("chr1", 40_000, 42_000, "segmental duplication (MAPQ 0)"),
    ("chr1", 130_000, 132_000, "segmental duplication (MAPQ 0)"),
    ("chr2", 20_000, 23_000, "paralogue missing from the reference"),
]

# ----------------------------------------------------------------------
# Palette. Two categorical slots per chart, validated for colour-vision
# deficiency against both surfaces; see references in the report footer.
# ----------------------------------------------------------------------

CSS = """
:root {
  --page:#f9f9f7; --surface:#fcfcfb; --ink:#0b0b0b; --ink-2:#52514e;
  --muted:#898781; --grid:#e1e0d9; --axis:#c3c2b7; --border:rgba(11,11,11,.10);
  --s1:#2a78d6; --s2:#eb6834; --s3:#e34948;
}
@media (prefers-color-scheme: dark) {
  :root {
    --page:#0d0d0d; --surface:#1a1a19; --ink:#fff; --ink-2:#c3c2b7;
    --muted:#898781; --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,.10);
    --s1:#3987e5; --s2:#d95926; --s3:#e66767;
  }
}
* { box-sizing:border-box; }
body { margin:0; padding:32px 24px 64px; background:var(--page); color:var(--ink);
  font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif; }
main { max-width:960px; margin:0 auto; }
h1 { font-size:26px; margin:0 0 4px; letter-spacing:-.01em; }
h2 { font-size:18px; margin:0 0 6px; letter-spacing:-.01em; }
.sub { color:var(--ink-2); margin:0 0 28px; }
section { background:var(--surface); border:1px solid var(--border); border-radius:10px;
  padding:20px 22px 18px; margin:0 0 20px; overflow-x:auto; }
.caption { color:var(--ink-2); font-size:13.5px; margin:2px 0 14px; max-width:72ch; }
.note { color:var(--ink-2); font-size:13px; margin:12px 0 0; max-width:72ch; }
.legend { display:flex; gap:18px; flex-wrap:wrap; margin:0 0 10px; font-size:13px;
  color:var(--ink-2); }
.legend span { display:inline-flex; align-items:center; gap:7px; }
.swatch { width:11px; height:11px; border-radius:2px; display:inline-block; }
.tiles { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:2px;
  background:var(--border); border:1px solid var(--border); border-radius:8px;
  overflow:hidden; }
.tile { background:var(--surface); padding:13px 15px; }
.tile .k { font-size:11.5px; text-transform:uppercase; letter-spacing:.05em;
  color:var(--muted); }
.tile .v { font-size:22px; margin-top:3px; }
.tile .u { font-size:13px; color:var(--ink-2); }
table { border-collapse:collapse; font-size:13.5px; width:100%; }
th, td { text-align:right; padding:7px 10px; border-bottom:1px solid var(--border);
  font-variant-numeric:tabular-nums; }
th:first-child, td:first-child { text-align:left; font-variant-numeric:normal; }
thead th { color:var(--ink-2); font-weight:600; border-bottom:1px solid var(--axis); }
tbody tr:last-child td { border-bottom:none; }
svg { display:block; overflow:visible; }
svg text { font:11.5px system-ui,-apple-system,sans-serif; fill:var(--muted); }
svg .lbl { fill:var(--ink-2); }
svg .ttl { fill:var(--ink); font-size:12.5px; }
svg .grid { stroke:var(--grid); stroke-width:1; }
svg .axis { stroke:var(--axis); stroke-width:1; }
svg .band { fill:var(--muted); opacity:.11; }
svg .a1 { fill:var(--s1); opacity:.20; stroke:none; }
svg .l1 { stroke:var(--s1); stroke-width:2; fill:none; stroke-linejoin:round; }
svg .l2 { stroke:var(--s2); stroke-width:2; fill:none; stroke-linejoin:round; }
svg .f1 { fill:var(--s1); } svg .f2 { fill:var(--s2); } svg .f3 { fill:var(--s3); }
svg .k1 { stroke:var(--s1); stroke-width:2; } svg .k3 { stroke:var(--s3); stroke-width:2; }
svg .conn { stroke:var(--axis); stroke-width:2; }
footer { color:var(--muted); font-size:12.5px; max-width:72ch; margin-top:28px; }
code { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:.92em; }
"""

# ----------------------------------------------------------------------
# Tiny SVG helpers
# ----------------------------------------------------------------------


def esc(t):
    return (str(t).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def txt(x, y, s, cls="", anchor="start", dy=0):
    return ('<text x="%.1f" y="%.1f" text-anchor="%s"%s class="%s">%s</text>'
            % (x, y + dy, anchor, "", cls, esc(s)))


def line(x1, y1, x2, y2, cls):
    return ('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" class="%s"/>'
            % (x1, y1, x2, y2, cls))


def rect(x, y, w, h, cls, rx=0):
    return ('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%g" class="%s"/>'
            % (x, y, max(w, 0), max(h, 0), rx, cls))


def poly(points, cls):
    if not points:
        return ""
    d = " ".join("%.1f,%.1f" % p for p in points)
    return '<polyline points="%s" class="%s"/>' % (d, cls)


def area(points, baseline, cls):
    """A filled band from the baseline up to the given points."""
    if not points:
        return ""
    d = " ".join("%.1f,%.1f" % p for p in points)
    return ('<polygon points="%.1f,%.1f %s %.1f,%.1f" class="%s"/>'
            % (points[0][0], baseline, d, points[-1][0], baseline, cls))


def dot(x, y, r, cls):
    return '<circle cx="%.1f" cy="%.1f" r="%g" class="%s"/>' % (x, y, r, cls)


def nice_ticks(lo, hi, n=4):
    """A short list of round numbers spanning [lo, hi]."""
    if hi <= lo:
        hi = lo + 1
    raw = (hi - lo) / float(n)
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            step = m * mag
            break
    start = math.ceil(lo / step) * step
    out, v = [], start
    while v <= hi + step * 1e-9:
        out.append(round(v, 10))
        v += step
    return out


# ----------------------------------------------------------------------
# Readers
# ----------------------------------------------------------------------


def read_fai(path):
    contigs = []
    with open(path) as fh:
        for ln in fh:
            f = ln.split("\t")
            contigs.append((f[0], int(f[1])))
    return contigs


def read_depth_windows(path):
    """<chrom>\t<window index>\t<mean depth> -> {(chrom, win): depth}"""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for ln in fh:
            f = ln.split()
            if len(f) == 3:
                out[(f[0], int(f[1]))] = float(f[2])
    return out


def read_dup_metrics(path):
    """PERCENT_DUPLICATION and READ_PAIRS_EXAMINED out of a Picard metrics file."""
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        lines = fh.read().splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("LIBRARY\t"):
            keys = ln.split("\t")
            vals = lines[i + 1].split("\t")
            row = dict(zip(keys, vals))
            return {"dup": float(row.get("PERCENT_DUPLICATION", 0) or 0),
                    "pairs": int(row.get("READ_PAIRS_EXAMINED", 0) or 0)}
    return None


def read_depth_summary(path):
    """The two lines step 6 wrote: total depth and MAPQ>=20 depth."""
    out = {}
    if not os.path.exists(path):
        return out
    for ln in open(path):
        m = re.search(r"mean ([\d.]+)x, ([\d.]+)% below 10x", ln)
        if not m:
            continue
        key = "callable" if "MAPQ" in ln else "all"
        out[key] = (float(m.group(1)), float(m.group(2)))
    return out


def read_concordance(path):
    """gatk Concordance summary -> {'SNP': {...}, 'INDEL': {...}}"""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        header = fh.readline().split()
        for ln in fh:
            f = ln.split()
            if len(f) < len(header):
                continue
            row = dict(zip(header, f))
            out[row["type"]] = {
                "TP": int(row["TP"]), "FP": int(row["FP"]), "FN": int(row["FN"]),
                "recall": float(row["RECALL"]), "precision": float(row["PRECISION"])}
    return out


def read_vcf_sites(path):
    """{(chrom, pos, ref, alt)} — gzip reads bgzf files fine."""
    sites = set()
    if not os.path.exists(path):
        return sites
    # errors="replace": GATK stamps its headers with a localised date, which on
    # a non-English machine can contain bytes that are not valid UTF-8. We only
    # ever read the data lines, so mangling a header character is harmless.
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            sites.add((f[0], int(f[1]), f[3], f[4]))
    return sites


def recal_by_cycle(path):
    """
    RecalTable2, Cycle covariate -> {cycle: (reported_Q, empirical_Q)}.

    Reported quality is averaged in PROBABILITY space, not in Phred space:
    Phred is a logarithm, and averaging logarithms of error rates would
    understate the common, low-quality bases. Empirical quality comes from the
    error counts GATK observed at that cycle.
    """
    if not os.path.exists(path):
        return {}
    obs = defaultdict(float)
    err = defaultdict(float)
    perr = defaultdict(float)
    in_t2 = False
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#:GATKTable:RecalTable2"):
                in_t2 = True
                continue
            if in_t2 and ln.startswith("#:GATKTable"):
                break
            if not in_t2:
                continue
            f = ln.split()
            if len(f) < 8 or f[3] != "Cycle" or f[4] != "M":
                continue
            q = float(f[1])
            cycle = abs(int(f[2]))
            n = float(f[6])
            e = float(f[7])
            obs[cycle] += n
            err[cycle] += e
            perr[cycle] += n * (10 ** (-q / 10.0))
    out = {}
    for c in sorted(obs):
        if obs[c] < 200:                      # too few bases to mean anything
            continue
        rep = -10 * math.log10(perr[c] / obs[c])
        emp = -10 * math.log10(max(err[c], 0.5) / obs[c])
        out[c] = (rep, emp)
    return out


# ----------------------------------------------------------------------
# Genome coordinate axis, shared by the two whole-genome figures
# ----------------------------------------------------------------------


class GenomeAxis:
    """Lays the contigs end to end on one x axis."""

    def __init__(self, contigs, x0, x1):
        self.contigs = contigs
        self.total = sum(l for _, l in contigs)
        self.offset = {}
        acc = 0
        for name, length in contigs:
            self.offset[name] = acc
            acc += length
        self.x0, self.x1 = x0, x1

    def x(self, chrom, pos):
        return self.x0 + (self.offset[chrom] + pos) / float(self.total) * (self.x1 - self.x0)

    def chrom_labels(self, y):
        out = []
        for name, length in self.contigs:
            mid = self.x(name, length / 2.0)
            out.append(txt(mid, y, name, cls="lbl", anchor="middle"))
        return out

    def dividers(self, y0, y1):
        out = []
        for name, _ in self.contigs[1:]:
            x = self.x(name, 0)
            out.append(line(x, y0, x, y1, "axis"))
        return out

    def bands(self, y0, h):
        out = []
        for chrom, s, e, _ in REGIONS:
            if chrom not in self.offset:
                continue
            x_a, x_b = self.x(chrom, s), self.x(chrom, e)
            out.append(rect(x_a, y0, max(x_b - x_a, 1.5), h, "band"))
        return out


# ----------------------------------------------------------------------
# Figure 1 — coverage along the genome
# ----------------------------------------------------------------------


def fig_depth(samples, depth_all, depth_q20, contigs):
    W, LEFT, RIGHT = 900, 50, 14
    ROW, GAP, TOP = 104, 42, 16
    H = TOP + len(samples) * (ROW + GAP) + 18
    ax = GenomeAxis(contigs, LEFT, W - RIGHT)

    hi = 0
    for s in samples:
        for d in list(depth_all[s].values()) + list(depth_q20[s].values()):
            hi = max(hi, d)
    hi = math.ceil(hi / 10.0) * 10 or 10

    out = ['<svg viewBox="0 0 %d %d" width="100%%" role="img">' % (W, H)]

    for i, s in enumerate(samples):
        y0 = TOP + i * (ROW + GAP)
        y1 = y0 + ROW

        def y(v):
            return y1 - (v / float(hi)) * ROW

        out += ax.bands(y0, ROW)
        for t in nice_ticks(0, hi, 3):
            out.append(line(LEFT, y(t), W - RIGHT, y(t), "grid"))
            out.append(txt(LEFT - 7, y(t) + 3.5, "%g" % t, anchor="end"))
        out += ax.dividers(y0, y1)
        out.append(line(LEFT, y1, W - RIGHT, y1, "axis"))

        # "all reads" as a filled envelope, "callable" as a line on top of it:
        # the two agree almost everywhere, so a second line would simply hide
        # under the first. The shape of the gap is the information.
        for chrom, length in contigs:
            pts = [(ax.x(chrom, w * 1000 + 500), y(min(d, hi)))
                   for (c, w), d in sorted(depth_all[s].items()) if c == chrom]
            out.append(area(pts, y1, "a1"))
            pts = [(ax.x(chrom, w * 1000 + 500), y(min(d, hi)))
                   for (c, w), d in sorted(depth_q20[s].items()) if c == chrom]
            out.append(poly(pts, "l2"))

        out.append(txt(LEFT, y0 - 10, s, cls="ttl"))

    for chrom, a, b, label in REGIONS:
        out.append(txt(ax.x(chrom, (a + b) / 2.0), TOP - 4, "▼", anchor="middle"))
    out += ax.chrom_labels(H - 3)
    out.append("</svg>")
    return "\n".join(out)


# ----------------------------------------------------------------------
# Figure 2 — where the errors are
# ----------------------------------------------------------------------


def fig_errors(samples, fp, fn, contigs):
    W, LEFT, RIGHT = 900, 46, 14
    LANE, GAP, TOP = 56, 40, 16
    H = TOP + len(samples) * (LANE + GAP) + 18
    ax = GenomeAxis(contigs, LEFT, W - RIGHT)

    out = ['<svg viewBox="0 0 %d %d" width="100%%" role="img">' % (W, H)]

    for i, s in enumerate(samples):
        y0 = TOP + i * (LANE + GAP)
        mid = y0 + LANE / 2.0
        out += ax.bands(y0, LANE)
        out += ax.dividers(y0, y0 + LANE)
        out.append(line(LEFT, mid, W - RIGHT, mid, "axis"))
        for chrom, pos in fp.get(s, []):
            x = ax.x(chrom, pos)
            out.append(line(x, mid, x, mid - LANE / 2.0 + 3, "k3"))
        for chrom, pos in fn.get(s, []):
            x = ax.x(chrom, pos)
            out.append(line(x, mid, x, mid + LANE / 2.0 - 3, "k1"))
        out.append(txt(LEFT, y0 - 11, "%s — %d false positives above, %d false "
                       "negatives below"
                       % (s, len(fp.get(s, [])), len(fn.get(s, []))), cls="ttl"))
    out += ax.chrom_labels(H - 2)
    out.append("</svg>")
    return "\n".join(out)


# ----------------------------------------------------------------------
# Figure 3 — what filtering bought and what it cost
# ----------------------------------------------------------------------


def fig_tradeoff(samples, conc, conc_raw):
    panels = [("SNP", "recall", "SNP sensitivity"),
              ("SNP", "precision", "SNP precision"),
              ("INDEL", "recall", "Indel sensitivity"),
              ("INDEL", "precision", "Indel precision")]
    PW, PH, LEFT, TOP, GAP_X, GAP_Y = 348, 104, 78, 22, 84, 54
    W = LEFT + PW * 2 + GAP_X + 26
    H = TOP + PH * 2 + GAP_Y + 10

    out = ['<svg viewBox="0 0 %d %d" width="100%%" role="img">' % (W, H)]

    for p, (vtype, metric, title) in enumerate(panels):
        col, row = p % 2, p // 2
        x0 = LEFT + col * (PW + GAP_X)
        y0 = TOP + row * (PH + GAP_Y)
        vals = []
        for s in samples:
            for src in (conc, conc_raw):
                d = src.get(s, {}).get(vtype)
                if d:
                    vals.append(d[metric])
        if not vals:
            continue
        lo = min(min(vals), 0.95) - 0.01
        hi = max(max(vals), 1.0) + 0.002
        if hi - lo < 1e-6:
            lo, hi = lo - 0.01, hi + 0.01

        def x(v):
            return x0 + (v - lo) / (hi - lo) * PW

        out.append(txt(x0 - 68, y0 - 8, title, cls="ttl"))
        for t in nice_ticks(lo, hi, 3):
            if t > hi - (hi - lo) * 0.02:
                continue
            out.append(line(x(t), y0, x(t), y0 + PH - 14, "grid"))
            out.append(txt(x(t), y0 + PH - 1, ("%.3f" % t).rstrip("0").rstrip("."),
                           anchor="middle"))
        for j, s in enumerate(samples):
            yy = y0 + 16 + j * 24
            a = conc_raw.get(s, {}).get(vtype)
            b = conc.get(s, {}).get(vtype)
            if not a or not b:
                continue
            xa, xb = x(a[metric]), x(b[metric])
            out.append(txt(x0 - 8, yy + 4, s, cls="lbl", anchor="end"))
            out.append(line(xa, yy, xb, yy, "conn"))
            out.append(dot(xa, yy, 4.5, "f2"))
            out.append(dot(xb, yy, 5, "f1"))
    out.append("</svg>")
    return "\n".join(out)


# ----------------------------------------------------------------------
# Figure 4 — did BQSR make the qualities honest?
# ----------------------------------------------------------------------


def fig_bqsr(samples, before, after):
    W, LEFT, RIGHT = 900, 52, 14
    ROW, GAP, TOP = 96, 42, 16
    H = TOP + len(samples) * (ROW + GAP) + 26

    series = []
    for s in samples:
        for src in (before, after):
            for c, (rep, emp) in src.get(s, {}).items():
                series.append(emp - rep)
    if not series:
        return ""
    lo, hi = min(series + [0]), max(series + [0])
    pad = max((hi - lo) * 0.15, 0.5)
    lo, hi = lo - pad, hi + pad
    max_cycle = max([max(before.get(s, {c: 0 for c in [1]}).keys() or [1])
                     for s in samples] + [1])

    out = ['<svg viewBox="0 0 %d %d" width="100%%" role="img">' % (W, H)]

    for i, s in enumerate(samples):
        y0 = TOP + i * (ROW + GAP)
        y1 = y0 + ROW

        def y(v):
            return y1 - (v - lo) / (hi - lo) * ROW

        def x(c):
            return LEFT + (c - 1) / float(max(max_cycle - 1, 1)) * (W - LEFT - RIGHT)

        for t in nice_ticks(lo, hi, 4):
            out.append(line(LEFT, y(t), W - RIGHT, y(t), "grid"))
            out.append(txt(LEFT - 7, y(t) + 3.5, "%+g" % t, anchor="end"))
        out.append(line(LEFT, y(0), W - RIGHT, y(0), "axis"))
        for src, cls in ((before.get(s, {}), "l1"), (after.get(s, {}), "l2")):
            pts = [(x(c), y(max(min(emp - rep, hi), lo)))
                   for c, (rep, emp) in sorted(src.items())]
            out.append(poly(pts, cls))
        out.append(txt(LEFT, y0 - 10, s, cls="ttl"))
        if i == len(samples) - 1:
            for c in (1, 50, 100, 150):
                if c <= max_cycle:
                    out.append(txt(x(c), y1 + 14, str(c), anchor="middle"))
            out.append(txt((LEFT + W - RIGHT) / 2.0, y1 + 28, "cycle (position in read)",
                           cls="lbl", anchor="middle"))
    out.append("</svg>")
    return "\n".join(out)


# ----------------------------------------------------------------------
# Page assembly
# ----------------------------------------------------------------------


def legend(items):
    parts = []
    for color, label in items:
        parts.append('<span><i class="swatch" style="background:var(--%s)"></i>%s</span>'
                     % (color, esc(label)))
    return '<div class="legend">%s</div>' % "".join(parts)


def tiles(rows):
    cells = "".join('<div class="tile"><div class="k">%s</div>'
                    '<div class="v">%s <span class="u">%s</span></div></div>'
                    % (esc(k), esc(v), esc(u)) for k, v, u in rows)
    return '<div class="tiles">%s</div>' % cells


def table(headers, rows):
    head = "".join("<th>%s</th>" % esc(h) for h in headers)
    body = "".join("<tr>%s</tr>" % "".join("<td>%s</td>" % esc(c) for c in r)
                   for r in rows)
    return "<table><thead><tr>%s</tr></thead><tbody>%s</tbody></table>" % (head, body)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="results")
    ap.add_argument("--reference", default="data/reference/reference.fasta")
    ap.add_argument("--samples", default=os.environ.get(
        "SAMPLES", "sample01 sample02 sample03"))
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    R = args.results
    samples = args.samples.split()
    out_path = args.out or os.path.join(R, "10_report", "report.html")

    fai = args.reference + ".fai"
    if not os.path.exists(fai):
        raise SystemExit("ERROR: %s not found. Run ./scripts/02_prepare_reference.sh"
                         % fai)
    contigs = read_fai(fai)

    depth_all, depth_q20, dup, dsum = {}, {}, {}, {}
    conc, conc_raw, fp, fn = {}, {}, {}, {}
    bqsr_before, bqsr_after = {}, {}

    for s in samples:
        depth_all[s] = read_depth_windows("%s/05_alignment_qc/%s.depth_1kb_q0.tsv" % (R, s))
        depth_q20[s] = read_depth_windows("%s/05_alignment_qc/%s.depth_1kb_q20.tsv" % (R, s))
        dup[s] = read_dup_metrics("%s/03_dedup/%s.metrics.txt" % (R, s))
        dsum[s] = read_depth_summary("%s/05_alignment_qc/%s.depth_summary.txt" % (R, s))
        conc[s] = read_concordance("%s/09_evaluation/%s.concordance.tsv" % (R, s))
        conc_raw[s] = read_concordance("%s/09_evaluation/%s.concordance_raw.tsv" % (R, s))
        bqsr_before[s] = recal_by_cycle("%s/04_bqsr/%s.recal.table" % (R, s))
        bqsr_after[s] = recal_by_cycle("%s/04_bqsr/%s.after.table" % (R, s))

        truth = read_vcf_sites("%s/09_evaluation/%s.truth.vcf.gz" % (R, s))
        calls = read_vcf_sites("%s/09_evaluation/%s.calls.vcf.gz" % (R, s))
        fp[s] = sorted((c, p) for c, p, _, _ in calls - truth)
        fn[s] = sorted((c, p) for c, p, _, _ in truth - calls)

    if not any(depth_all.values()):
        raise SystemExit("ERROR: no results found under %s/. Run the pipeline first." % R)

    # --- the page ------------------------------------------------------
    html = ['<!doctype html><html lang="en"><head><meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width,initial-scale=1">',
            '<title>Variant calling report</title>',
            '<style>%s</style></head><body><main>' % CSS,
            '<h1>Mapping and variant calling — visual report</h1>',
            '<p class="sub">Generated by <code>scripts/11_report.py</code> from the '
            'files under <code>%s/</code>. Every number here is one you can also '
            'read out of those files by hand.</p>' % esc(R)]

    # Tiles
    html.append('<section><h2>The three libraries</h2>'
                '<p class="caption">Nominal coverage is not information content. '
                'Compare the duplicate rate against the gap between total and '
                'callable depth — they are different losses, with different '
                'causes.</p>')
    for s in samples:
        rows = []
        if dup.get(s):
            rows.append(("read pairs", "{:,}".format(dup[s]["pairs"]), ""))
            rows.append(("duplicates", "%.1f" % (100 * dup[s]["dup"]), "%"))
        if dsum.get(s).get("all"):
            rows.append(("depth, all reads", "%.1f" % dsum[s]["all"][0], "x"))
        if dsum.get(s).get("callable"):
            rows.append(("depth, MAPQ ≥ 20", "%.1f" % dsum[s]["callable"][0], "x"))
            rows.append(("below 10x", "%.1f" % dsum[s]["callable"][1], "%"))
        html.append('<h2 style="font-size:14px;margin:16px 0 7px">%s</h2>%s'
                    % (esc(s), tiles(rows)))
    html.append('</section>')

    # Depth
    html.append('<section><h2>Coverage along the genome</h2>'
                '<p class="caption">Mean depth (×) per 1 kb window. Grey bands '
                'mark the three regions the dataset was built around. In the '
                'duplicated segments on chr1 the orange line dives into the blue '
                'area: the reads are there, but MAPQ 0 makes them unusable. On '
                'chr2 area and line rise together — that extra coverage is real '
                'and usable, and should not exist at all.</p>')
    html.append(legend([("s1", "all reads (shaded area)"),
                        ("s2", "MAPQ ≥ 20 — what the caller uses (line)")]))
    html.append(fig_depth(samples, depth_all, depth_q20, contigs))
    html.append('<p class="note">Depth is averaged over 1 kb windows and '
                'duplicate-flagged reads are excluded, so sample02 already shows '
                'its effective coverage rather than its nominal one.</p></section>')

    # Errors
    html.append('<section><h2>Where the errors are</h2>'
                '<p class="caption">Every false positive and false negative, at '
                'its coordinate. This is the figure the whole practical is for: '
                'the errors are not scattered, they are two regions, and the two '
                'fail in opposite directions.</p>')
    html.append(legend([("s3", "false positive (called, not real)"),
                        ("s1", "false negative (real, not called)")]))
    html.append(fig_errors(samples, fp, fn, contigs))
    html.append('<p class="note">False positives sit in the collapsed paralogue '
                'and repeat identically in all three samples — sequencing deeper '
                'would make them more convincing, not fewer. False negatives sit '
                'in the segmental duplications, where "no call" is the honest '
                'answer.</p></section>')

    # Trade-off
    html.append('<section><h2>What filtering bought, and what it cost</h2>'
                '<p class="caption">Each line runs from the unfiltered call set to '
                'the filtered one. Precision should move right, sensitivity left. '
                'The length of the line is the price of the filter.</p>')
    html.append(legend([("s2", "unfiltered"), ("s1", "after hard filtering")]))
    html.append(fig_tradeoff(samples, conc, conc_raw))
    html.append('</section>')

    # BQSR
    if any(bqsr_before.values()):
        html.append('<section><h2>Did recalibration make the qualities honest?</h2>'
                    '<p class="caption">Zero means the instrument\'s reported '
                    'quality matches the error rate actually observed. Below zero '
                    'means it claimed more confidence than it earned. Before '
                    'recalibration the curve should sag towards the end of the '
                    'read; afterwards it should sit on the line.</p>')
        html.append(legend([("s1", "before BQSR"), ("s2", "after BQSR")]))
        html.append(fig_bqsr(samples, bqsr_before, bqsr_after))
        html.append('<p class="note">This replaces the AnalyzeCovariates PDF from '
                    'step 5, which needs R and often fails to build its plots. The '
                    'underlying numbers are the same recalibration tables.</p>'
                    '</section>')

    # Table view
    rows = []
    for s in samples:
        for vtype in ("SNP", "INDEL"):
            a, b = conc_raw.get(s, {}).get(vtype), conc.get(s, {}).get(vtype)
            if not b:
                continue
            rows.append([s, vtype, b["TP"], b["FP"], b["FN"],
                         "%.3f" % b["recall"], "%.3f" % b["precision"],
                         "%.3f" % a["recall"] if a else "—",
                         "%.3f" % a["precision"] if a else "—"])
    html.append('<section><h2>The numbers</h2>'
                '<p class="caption">The same data as the figures above, for the '
                'table in your report. TP + FN is the number of real variants that '
                'sample carries.</p>')
    html.append(table(["sample", "type", "TP", "FP", "FN", "sens.", "prec.",
                       "sens. raw", "prec. raw"], rows))
    html.append('</section>')

    html.append('<footer>Figures are inline SVG — screenshot a panel, or open this '
                'file in an editor and copy the <code>&lt;svg&gt;</code> block out. '
                'The page follows your system light/dark setting. Colours were '
                'checked for colour-vision deficiency against both '
                'backgrounds.</footer>')
    html.append('</main></body></html>')

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("\n".join(html))

    print("Wrote %s" % out_path)
    for s in samples:
        print("  %-10s %3d false positives, %3d false negatives"
              % (s, len(fp[s]), len(fn[s])))
    print("\nOpen it in a browser:\n    open %s" % out_path)


if __name__ == "__main__":
    main()
