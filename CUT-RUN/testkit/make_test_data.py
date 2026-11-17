#!/usr/bin/env python3
"""Make a tiny synthetic dataset with known answers for pipeline testing.

Creates a 316 kb genome, a GTF, paired-end CUT&RUN reads (H3K27ac-style narrow
peaks, H3K9me3-style broad domains, input) and reverse-stranded paired-end
RNA-seq reads, plus truth files that say exactly what a correct pipeline
should find. Standard library only. Same seed -> identical files.

    python3 make_test_data.py --out testdata

Every read name records where the pair came from:
    cutrun: {sample}:{n}:cr:{contig}:{start}:{end}     (0-based, end-exclusive)
    rnaseq: {sample}:{n}:rna:{contig}:{gene_id}:0
"""
import argparse
import gzip
import os
import random

CONTIGS = [("ctgA", 200000), ("ctgB", 100000), ("Mt_test", 16000)]
GROUPS = ["BF", "SF"]
DAY = "D1"
CR_READ_LEN = 50
RNA_READ_LEN = 75

# Truth layout (0-based, end-exclusive)
AC_CENTERS = [("ctgA", 20000 + k * 22000) for k in range(8)] + \
             [("ctgB", 15000 + k * 20000) for k in range(4)]
AC_HALF_WIDTH = 150
ME_DOMAINS = [("ctgA", 30000, 38000), ("ctgA", 120000, 132000), ("ctgB", 50000, 58000)]

DE_UP = {2, 9, 14, 21}      # 4-fold higher in SF than BF
DE_DOWN = {5, 11, 18, 23}   # 4-fold lower in SF than BF
BASE_EXPR = [150, 300, 600, 1000]

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def make_genome(rng):
    genome = {}
    for name, length in CONTIGS:
        genome[name] = rng.choices("ACGT", weights=[30, 20, 20, 30], k=length)
    return genome


def make_genes(genome):
    """24 nuclear genes (1-3 exons, both strands) + 1 mitochondrial gene.
    Plants canonical splice sites so spliced aligners accept the introns."""
    genes = []
    for i in range(24):
        if i < 16:
            contig, start = "ctgA", 5000 + i * 11000
        else:
            contig, start = "ctgB", 5000 + (i - 16) * 11000
        strand = "+" if i % 2 == 0 else "-"
        n_exons = 1 + (i % 3)
        exons, pos = [], start
        for e in range(n_exons):
            exons.append((pos, pos + 600))
            pos += 600
            if e < n_exons - 1:
                donor, acceptor = ("GT", "AG") if strand == "+" else ("CT", "AC")
                genome[contig][pos:pos + 2] = list(donor)
                genome[contig][pos + 498:pos + 500] = list(acceptor)
                pos += 500
        genes.append({"id": "gene%d" % (i + 1), "idx": i, "contig": contig,
                      "strand": strand, "exons": exons})
    genes.append({"id": "MtGene1", "idx": 24, "contig": "Mt_test",
                  "strand": "+", "exons": [(1000, 1900)]})
    return genes


def transcript_seq(gene, genome_str):
    seq = "".join(genome_str[gene["contig"]][s:e] for s, e in gene["exons"])
    return seq if gene["strand"] == "+" else revcomp(seq)


def write_fasta(path, genome_str):
    with open(path, "w") as fh:
        for name, _ in CONTIGS:
            fh.write(">%s\n" % name)
            seq = genome_str[name]
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def write_gtf(path, genes):
    with open(path, "w") as fh:
        for g in genes:
            gs, ge = g["exons"][0][0] + 1, g["exons"][-1][1]   # GTF is 1-based, inclusive
            ga = 'gene_id "%s";' % g["id"]
            ta = ga + ' transcript_id "%s-RA";' % g["id"]
            row = [g["contig"], "maker", "", "", "", ".", g["strand"], "."]
            for feat, s, e, attr in [("gene", gs, ge, ga), ("transcript", gs, ge, ta)]:
                row[2], row[3], row[4] = feat, str(s), str(e)
                fh.write("\t".join(row + [attr]) + "\n")
            for n, (s, e) in enumerate(g["exons"], 1):
                row[2], row[3], row[4] = "exon", str(s + 1), str(e)
                fh.write("\t".join(row + [ta + ' exon_number "%d";' % n]) + "\n")


def write_fastq_pair(prefix, records):
    """records: list of (name, read1, read2)."""
    for mate in (1, 2):
        path = "%s_R%d.fastq.gz" % (prefix, mate)
        with open(path, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                chunks = []
                for rec in records:
                    seq = rec[mate]
                    chunks.append("@%s\n%s\n+\n%s\n" % (rec[0], seq, "I" * len(seq)))
                gz.write("".join(chunks).encode())


def cutrun_fragments(rng, target, strength):
    """Return a list of (contig, start, end) fragments for one sample."""
    frags = []
    total = sum(length for _, length in CONTIGS)

    def background(n):
        for _ in range(n):
            x = rng.randrange(total)
            for name, length in CONTIGS:
                if x < length:
                    break
                x -= length
            flen = min(300, max(100, int(rng.gauss(160, 25))))
            start = min(max(0, x), length - flen)
            frags.append((name, start, start + flen))

    if target == "IN":
        background(6000)
    elif target == "AC":
        background(3000)
        for contig, center in AC_CENTERS:
            for _ in range(int(250 * strength)):
                flen = min(300, max(100, int(rng.gauss(160, 25))))
                start = int(center + rng.gauss(0, 60) - flen / 2)
                frags.append((contig, start, start + flen))
    elif target == "ME":
        background(3000)
        for contig, ds, de in ME_DOMAINS:
            for _ in range(int((de - ds) / 1000 * 95 * strength)):
                flen = min(300, max(100, int(rng.gauss(160, 25))))
                start = rng.randrange(ds, de - flen)
                frags.append((contig, start, start + flen))
    rng.shuffle(frags)
    return frags


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="testdata", help="output directory")
    ap.add_argument("--seed", type=int, default=20261116)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    out = args.out
    for sub in ("reference", "fastq/cutrun", "fastq/rnaseq", "truth"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    genome = make_genome(rng)
    genes = make_genes(genome)
    genome_str = {name: "".join(seq) for name, seq in genome.items()}
    write_fasta(os.path.join(out, "reference", "test_genome.fa"), genome_str)
    write_gtf(os.path.join(out, "reference", "test_genes.gtf"), genes)

    read_pairs = []   # (sample, assay, pairs)

    # ---- CUT&RUN ----
    cr_samples = []
    for group in GROUPS:
        for target, reps in (("AC", 2), ("ME", 2), ("IN", 1)):
            for rep in range(1, reps + 1):
                sample = "%s-%s-%s-R%d" % (group, target, DAY, rep)
                strength = max(0.6, rng.gauss(1.0, 0.1))
                recs = []
                for n, (contig, s, e) in enumerate(cutrun_fragments(rng, target, strength), 1):
                    frag = genome_str[contig][s:e]
                    name = "%s:%d:cr:%s:%d:%d" % (sample, n, contig, s, e)
                    r1, r2 = frag[:CR_READ_LEN], revcomp(frag)[:CR_READ_LEN]
                    if rng.random() < 0.5:
                        r1, r2 = r2, r1
                    recs.append((name, r1, r2))
                write_fastq_pair(os.path.join(out, "fastq", "cutrun", sample), recs)
                cr_samples.append(sample)
                read_pairs.append((sample, "cutrun", len(recs)))

    # ---- RNA-seq (reverse-stranded: read 1 is antisense to the transcript) ----
    tx = {g["id"]: transcript_seq(g, genome_str) for g in genes}
    rna_samples, counts = [], {g["id"]: {} for g in genes}
    for group in GROUPS:
        for rep in range(1, 4):
            sample = "%s-RNA-%s-R%d" % (group, DAY, rep)
            depth = rng.gauss(1.0, 0.08)
            recs = []
            for g in genes:
                fold = 1.0
                if group == "SF" and g["idx"] in DE_UP:
                    fold = 4.0
                if group == "SF" and g["idx"] in DE_DOWN:
                    fold = 0.25
                n_frag = max(1, round(BASE_EXPR[g["idx"] % 4] * fold * depth * rng.gauss(1.0, 0.05)))
                counts[g["id"]][sample] = n_frag
                seq = tx[g["id"]]
                for _ in range(n_frag):
                    flen = max(RNA_READ_LEN, min(len(seq), int(rng.gauss(250, 30))))
                    s = rng.randrange(0, len(seq) - flen + 1)
                    frag = seq[s:s + flen]
                    recs.append([g, revcomp(frag)[:RNA_READ_LEN], frag[:RNA_READ_LEN]])
            rng.shuffle(recs)
            named = [("%s:%d:rna:%s:%s:0" % (sample, n, r[0]["contig"], r[0]["id"]), r[1], r[2])
                     for n, r in enumerate(recs, 1)]
            write_fastq_pair(os.path.join(out, "fastq", "rnaseq", sample), named)
            rna_samples.append(sample)
            read_pairs.append((sample, "rnaseq", len(named)))

    # ---- Truth files ----
    t = os.path.join(out, "truth")
    with open(os.path.join(t, "peaks_AC.bed"), "w") as fh:
        for n, (contig, c) in enumerate(AC_CENTERS, 1):
            fh.write("%s\t%d\t%d\tAC_peak%d\n" % (contig, c - AC_HALF_WIDTH, c + AC_HALF_WIDTH, n))
    with open(os.path.join(t, "domains_ME.bed"), "w") as fh:
        for n, (contig, s, e) in enumerate(ME_DOMAINS, 1):
            fh.write("%s\t%d\t%d\tME_domain%d\n" % (contig, s, e, n))
    with open(os.path.join(t, "read_pairs.tsv"), "w") as fh:
        fh.write("sample\tassay\tread_pairs\n")
        for row in read_pairs:
            fh.write("%s\t%s\t%d\n" % row)
    with open(os.path.join(t, "counts.tsv"), "w") as fh:
        fh.write("gene_id\t" + "\t".join(rna_samples) + "\n")
        for g in genes:
            fh.write(g["id"] + "\t" + "\t".join(str(counts[g["id"]][s]) for s in rna_samples) + "\n")
    with open(os.path.join(t, "de_truth.tsv"), "w") as fh:
        fh.write("gene_id\tlog2fc_SF_vs_BF\tis_de\n")
        for g in genes:
            lfc = 2 if g["idx"] in DE_UP else -2 if g["idx"] in DE_DOWN else 0
            fh.write("%s\t%d\t%s\n" % (g["id"], lfc, "yes" if lfc else "no"))
    with open(os.path.join(t, "genome_facts.tsv"), "w") as fh:
        fh.write("genome_size\t%d\n" % sum(length for _, length in CONTIGS))
        for name, length in CONTIGS:
            fh.write("contig\t%s\t%d\n" % (name, length))
        fh.write("cutrun_read_length\t%d\n" % CR_READ_LEN)
        fh.write("rnaseq_read_length\t%d\n" % RNA_READ_LEN)
        fh.write("rnaseq_strandedness\treverse (HISAT2 RF, featureCounts -s 2, htseq-count -s reverse)\n")
        fh.write("seed\t%d\n" % args.seed)
    with open(os.path.join(out, "samples_cutrun.txt"), "w") as fh:
        fh.write("\n".join(cr_samples) + "\n")
    with open(os.path.join(out, "samples_rnaseq.txt"), "w") as fh:
        fh.write("\n".join(rna_samples) + "\n")

    print("Wrote %s: %d CUT&RUN samples, %d RNA-seq samples, %d genes, %d AC peaks, %d ME domains"
          % (out, len(cr_samples), len(rna_samples), len(genes), len(AC_CENTERS), len(ME_DOMAINS)))


if __name__ == "__main__":
    main()
