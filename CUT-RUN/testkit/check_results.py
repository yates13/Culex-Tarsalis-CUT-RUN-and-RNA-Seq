#!/usr/bin/env python3
"""Check pipeline output against the truth files from make_test_data.py.

Each subcommand prints PASS/FAIL lines and exits non-zero if anything fails,
so it can end a Slurm test job. Standard library only.

  names    contig names agree between FASTA, GTF and (optionally) BAM
  fastq    read pairs survive trimming and R1/R2 stay in sync
  aligned  read pairs map back to where they came from (pipe in samtools view)
  peaks    narrow peaks (AC) recover every true peak with few false calls
  broad    broad peaks (ME) cover the true domains and little else
  counts   gene counts match the true fragment counts
  de       differential expression calls match the true DE genes

`names` works on real data too; the rest need --truth testdata/truth.
"""
import argparse
import csv
import gzip
import os
import sys

FAILS = []


def report(ok, msg):
    print(("PASS  " if ok else "FAIL  ") + msg)
    if not ok:
        FAILS.append(msg)


def hint(msg):
    print("      hint: " + msg)


def read_tsv(path):
    with open(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


def read_bed(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), f))
    return rows


def overlaps(a, b):
    return a[0] == b[0] and a[1] < b[2] and b[1] < a[2]


def truth_contigs(truth):
    return {r[1] for r in read_tsv(os.path.join(truth, "genome_facts.tsv")) if r[0] == "contig"}


# ---------------------------------------------------------------- names
def cmd_names(a):
    with open(a.fasta) as fh:
        fa = {line[1:].split()[0] for line in fh if line.startswith(">")}
    with open(a.gtf) as fh:
        gtf = {line.split("\t")[0] for line in fh if line.strip() and not line.startswith("#")}
    report(bool(fa), "FASTA has %d contigs" % len(fa))
    missing = sorted(gtf - fa)
    report(not missing, "all %d GTF contigs exist in the FASTA" % len(gtf))
    if missing:
        hint("in GTF but not FASTA: " + ", ".join(missing[:5]))
    if a.idxstats:
        bam = {r[0] for r in read_tsv(a.idxstats) if r[0] != "*"}
        report(bam == fa, "BAM contigs (%d) equal FASTA contigs (%d)" % (len(bam), len(fa)))
        shared = bam & gtf
        report(bool(shared), "BAM and GTF share %d contig names" % len(shared))
        if not shared:
            hint("no shared names means every gene will count zero")


# ---------------------------------------------------------------- fastq
def count_reads(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as fh:
        return sum(1 for _ in fh) // 4


def cmd_fastq(a):
    truth = {r[0]: int(r[2]) for r in read_tsv(os.path.join(a.truth, "read_pairs.tsv"))[1:]
             if a.assay in ("all", r[1])}
    files = sorted(os.listdir(a.dir))
    for sample, expected in truth.items():
        mates = {}
        for m in ("1", "2"):
            hits = [f for f in files if sample in f and ("_R%s" % m in f or "_%s.f" % m in f)]
            mates[m] = hits
        if len(mates["1"]) != 1 or len(mates["2"]) != 1:
            report(False, "%s: expected one R1 and one R2 file, found %d and %d"
                   % (sample, len(mates["1"]), len(mates["2"])))
            continue
        n1 = count_reads(os.path.join(a.dir, mates["1"][0]))
        n2 = count_reads(os.path.join(a.dir, mates["2"][0]))
        report(n1 == n2, "%s: R1 and R2 in sync (%d, %d)" % (sample, n1, n2))
        ok = a.min_frac * expected <= n1 <= expected
        report(ok, "%s: %d of %d pairs kept (need >= %.0f%%)"
               % (sample, n1, expected, 100 * a.min_frac))


# ---------------------------------------------------------------- aligned
def cmd_aligned(a):
    truth = {r[0]: int(r[2]) for r in read_tsv(os.path.join(a.truth, "read_pairs.tsv"))[1:]}
    seen, good = {}, {}
    fh = sys.stdin if a.sam == "-" else open(a.sam)
    for line in fh:
        if line.startswith("@"):
            continue
        f = line.split("\t", 9)
        flag = int(f[1])
        if not flag & 0x40 or flag & 0x904 or flag & 0x8:
            continue   # first mate only; skip secondary, supplementary, unmapped
        name = f[0].split(":")
        if len(name) != 6:
            continue
        sample, assay, contig = name[0], name[2], name[3]
        seen[sample] = seen.get(sample, 0) + 1
        ok = f[2] == contig
        if ok and assay == "cr":
            start, end = int(name[4]), int(name[5])
            ok = min(int(f[3]), int(f[7])) == start + 1 and abs(int(f[8])) == end - start
        good[sample] = good.get(sample, 0) + ok
    report(bool(seen), "found alignments with test-data read names")
    for sample in sorted(seen):
        rate = seen[sample] / truth.get(sample, seen[sample])
        acc = good[sample] / seen[sample]
        report(rate >= a.min_mapped, "%s: %.1f%% of pairs aligned (need >= %.0f%%)"
               % (sample, 100 * rate, 100 * a.min_mapped))
        report(acc >= a.min_correct, "%s: %.1f%% aligned to their true position (need >= %.0f%%)"
               % (sample, 100 * acc, 100 * a.min_correct))


# ---------------------------------------------------------------- peaks
def cmd_peaks(a):
    truth = read_bed(os.path.join(a.truth, "peaks_AC.bed"))
    calls = read_bed(a.peaks)
    contigs = truth_contigs(a.truth)
    report(bool(calls), "%d peaks called" % len(calls))
    bad = sorted({c[0] for c in calls} - contigs)
    report(not bad, "peak contig names match the test genome")
    found = [t for t in truth if any(overlaps(t, c) for c in calls)]
    report(len(found) == len(truth), "%d of %d true peaks recovered" % (len(found), len(truth)))
    false = [c for c in calls if not any(overlaps(t, c) for t in truth)]
    report(len(false) <= a.max_false, "%d false peaks (allow %d)" % (len(false), a.max_false))
    if len(calls) and len(false) > 5 * len(truth):
        hint("many false peaks: was the input control used, and -g set to the test genome size?")
    dists = []
    for c in calls:
        if len(c[3]) >= 10 and c[3][9].lstrip("-").isdigit() and int(c[3][9]) >= 0:
            summit = c[1] + int(c[3][9])
            for t in truth:
                if overlaps(t, c):
                    dists.append(abs(summit - (t[1] + t[2]) // 2))
    if dists:
        med = sorted(dists)[len(dists) // 2]
        report(med <= a.max_summit, "median summit is %d bp from the true centre (allow %d)"
               % (med, a.max_summit))


# ---------------------------------------------------------------- broad
def cmd_broad(a):
    truth = read_bed(os.path.join(a.truth, "domains_ME.bed"))
    calls = read_bed(a.peaks)
    report(bool(calls), "%d broad peaks called" % len(calls))
    for t in truth:
        covered = set()
        for c in calls:
            if overlaps(t, c):
                covered.update(range(max(t[1], c[1]), min(t[2], c[2])))
        frac = len(covered) / (t[2] - t[1])
        report(frac >= a.min_cover, "%s %s:%d-%d is %.0f%% covered (need >= %.0f%%)"
               % (t[3][3], t[0], t[1], t[2], 100 * frac, 100 * a.min_cover))
    total = sum(c[2] - c[1] for c in calls)
    inside = 0
    for c in calls:
        for t in truth:
            lo, hi = max(c[1], t[1] - a.slop), min(c[2], t[2] + a.slop)
            inside += max(0, hi - lo)
    outside = (total - inside) / total if total else 0
    report(outside <= a.max_outside, "%.0f%% of called bases lie outside true domains (allow %.0f%%)"
           % (100 * outside, 100 * a.max_outside))


# ---------------------------------------------------------------- counts
def load_counts(a, samples):
    """Return {sample: {gene: count}} from a featureCounts table, a plain
    matrix, or a directory of htseq-count files."""
    obs = {}
    if a.htseq_dir:
        for sample in samples:
            hits = [f for f in sorted(os.listdir(a.htseq_dir)) if sample in f]
            if len(hits) != 1:
                continue
            obs[sample] = {r[0]: float(r[1]) for r in read_tsv(os.path.join(a.htseq_dir, hits[0]))
                           if len(r) >= 2 and not r[0].startswith("__")}
        return obs
    rows = [r for r in read_tsv(a.table) if not r[0].startswith("#")]
    header, body = rows[0], rows[1:]
    first = 6 if header[0] == "Geneid" and len(header) > 6 else 1
    for col in range(first, len(header)):
        match = [s for s in samples if s in header[col]]
        if len(match) == 1:
            obs[match[0]] = {r[0]: float(r[col]) for r in body}
    return obs


def cmd_counts(a):
    rows = read_tsv(os.path.join(a.truth, "counts.tsv"))
    samples = rows[0][1:]
    truth = {s: {r[0]: int(r[i + 1]) for r in rows[1:]} for i, s in enumerate(samples)}
    obs = load_counts(a, samples)
    report(len(obs) == len(samples), "found count columns for %d of %d samples" % (len(obs), len(samples)))
    for sample in samples:
        if sample not in obs:
            continue
        genes = truth[sample]
        missing = [g for g in genes if g not in obs[sample]]
        if missing:
            report(False, "%s: %d genes missing from the table (e.g. %s)" % (sample, len(missing), missing[0]))
            continue
        wrong = [g for g in genes
                 if abs(obs[sample][g] - genes[g]) > max(3, a.tol * genes[g])]
        report(not wrong, "%s: %d of %d genes within %.0f%% of the true count"
               % (sample, len(genes) - len(wrong), len(genes), 100 * a.tol))
        if wrong:
            tot_obs, tot_exp = sum(obs[sample][g] for g in genes), sum(genes.values())
            ratio = tot_obs / tot_exp
            if ratio < 0.05:
                hint("near-zero counts: wrong strandedness setting, mismatched contig names, or wrong feature type")
            elif 1.9 < ratio < 2.1:
                hint("about 2x expected: counting reads, not fragments (featureCounts needs -p --countReadPairs)")
            elif 0.4 < ratio < 0.6:
                hint("about half expected: check strandedness and how mates are paired")
            else:
                g = wrong[0]
                hint("e.g. %s: observed %g, expected %d" % (g, obs[sample][g], genes[g]))


# ---------------------------------------------------------------- de
def cmd_de(a):
    truth = {r[0]: (int(r[1]), r[2] == "yes") for r in read_tsv(os.path.join(a.truth, "de_truth.tsv"))[1:]}
    with open(a.table, newline="") as fh:
        dialect = csv.excel_tab if "\t" in fh.readline() else csv.excel
        fh.seek(0)
        rows = list(csv.reader(fh, dialect))
    header = [h.strip('"') for h in rows[0]]
    if len(rows[1]) == len(header) + 1:      # R row names without a header cell
        header = ["gene"] + header
    try:
        li, pi = header.index("log2FoldChange"), header.index("padj")
    except ValueError:
        report(False, "table needs log2FoldChange and padj columns (found: %s)" % ", ".join(header))
        return
    res = {}
    for r in rows[1:]:
        try:
            lfc = float(r[li]) * (-1 if a.flip else 1)
        except ValueError:
            continue
        try:
            padj = float(r[pi])
        except ValueError:
            padj = 1.0   # NA
        res[r[0].strip('"')] = (lfc, padj)
    de = [g for g, (_, is_de) in truth.items() if is_de]
    hit = [g for g in de if g in res and res[g][1] < a.alpha and res[g][0] * truth[g][0] > 0
           and abs(res[g][0] - truth[g][0]) <= 1]
    report(len(hit) == len(de), "%d of %d true DE genes significant with the right direction and size"
           % (len(hit), len(de)))
    if de and not hit:
        signflip = [g for g in de if g in res and res[g][0] * truth[g][0] < 0]
        if len(signflip) == len(de):
            hint("every direction is reversed: contrast is BF vs SF; rerun with --flip or swap the contrast")
    false = [g for g, (_, is_de) in truth.items() if not is_de and g in res and res[g][1] < a.alpha]
    report(len(false) <= a.max_false, "%d non-DE genes called significant (allow %d)" % (len(false), a.max_false))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("names")
    p.add_argument("--fasta", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--idxstats", help="output of samtools idxstats")
    p.set_defaults(func=cmd_names)

    p = sub.add_parser("fastq")
    p.add_argument("--truth", required=True)
    p.add_argument("--dir", required=True, help="directory of trimmed FASTQ files")
    p.add_argument("--assay", choices=["cutrun", "rnaseq", "all"], default="all")
    p.add_argument("--min-frac", type=float, default=0.98)
    p.set_defaults(func=cmd_fastq)

    p = sub.add_parser("aligned")
    p.add_argument("--truth", required=True)
    p.add_argument("sam", help="SAM file, or - for stdin (samtools view x.bam | ...)")
    p.add_argument("--min-mapped", type=float, default=0.95)
    p.add_argument("--min-correct", type=float, default=0.98)
    p.set_defaults(func=cmd_aligned)

    p = sub.add_parser("peaks")
    p.add_argument("--truth", required=True)
    p.add_argument("--peaks", required=True, help="narrowPeak or BED file")
    p.add_argument("--max-false", type=int, default=2)
    p.add_argument("--max-summit", type=int, default=75)
    p.set_defaults(func=cmd_peaks)

    p = sub.add_parser("broad")
    p.add_argument("--truth", required=True)
    p.add_argument("--peaks", required=True, help="broadPeak or BED file")
    p.add_argument("--min-cover", type=float, default=0.6)
    p.add_argument("--max-outside", type=float, default=0.15)
    p.add_argument("--slop", type=int, default=1000)
    p.set_defaults(func=cmd_broad)

    p = sub.add_parser("counts")
    p.add_argument("--truth", required=True)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--table", help="featureCounts output or a gene x sample TSV")
    g.add_argument("--htseq-dir", help="directory of per-sample htseq-count files")
    p.add_argument("--tol", type=float, default=0.05)
    p.set_defaults(func=cmd_counts)

    p = sub.add_parser("de")
    p.add_argument("--truth", required=True)
    p.add_argument("--table", required=True, help="DESeq2 results as CSV or TSV")
    p.add_argument("--flip", action="store_true", help="table is BF vs SF instead of SF vs BF")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--max-false", type=int, default=1)
    p.set_defaults(func=cmd_de)

    a = ap.parse_args()
    a.func(a)
    print("\n%s: %d check(s) failed" % (a.cmd, len(FAILS)) if FAILS else "\n%s: all checks passed" % a.cmd)
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
