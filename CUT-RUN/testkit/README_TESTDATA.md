# Pipeline test kit

A tiny synthetic dataset where the right answers are known in advance, and a
checker that fails loudly when a pipeline gets them wrong. Use it to test any
script change in minutes before spending cluster hours on real data.

Two scripts, Python 3 standard library only (no conda env needed):

- `make_test_data.py` writes the genome, GTF, reads and truth files (about 7 MB).
- `check_results.py` compares pipeline output with the truth.

## What is in the data

| Item | Content |
| --- | --- |
| Genome | 3 contigs, 316,000 bp: `ctgA`, `ctgB`, `Mt_test` |
| Genes | 25 (24 nuclear with 1-3 exons on both strands, 1 mitochondrial) |
| CUT&RUN | 10 samples: `BF` and `SF`, each with `AC` x2, `ME` x2, `IN` x1, day `D1`; paired 2x50 |
| AC truth | 12 narrow peaks (`truth/peaks_AC.bed`) |
| ME truth | 3 broad domains of 8-12 kb (`truth/domains_ME.bed`) |
| RNA-seq | 6 samples: `BF-RNA-D1-R1..3`, `SF-RNA-D1-R1..3`; paired 2x75, **reverse-stranded** |
| RNA truth | exact fragments per gene (`truth/counts.tsv`); 8 DE genes, 4 up and 4 down by 4-fold in SF (`truth/de_truth.tsv`) |

Sample names follow `{group}-{target}-{day}-R{rep}`. Files are
`{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`. Each read name records where
the pair came from, so alignments can be checked position by position.

## Use

Make the data once (same seed gives identical files):

```bash
python3 make_test_data.py --out testdata
```

Run your pipeline scripts with their paths pointed at `testdata/` and a test
output directory. Then check each stage (set `T=testdata/truth` first):

```bash
python3 check_results.py names --fasta testdata/reference/test_genome.fa --gtf testdata/reference/test_genes.gtf --idxstats sample.idxstats.txt
python3 check_results.py fastq --truth $T --dir trimmed_dir --assay cutrun
samtools view BF-AC-D1-R1.bam | python3 check_results.py aligned --truth $T -
python3 check_results.py peaks --truth $T --peaks BF-AC-D1_peaks.narrowPeak
python3 check_results.py broad --truth $T --peaks BF-ME-D1_peaks.broadPeak
python3 check_results.py counts --truth $T --table featureCounts_out.txt
python3 check_results.py counts --truth $T --htseq-dir htseq_out_dir
python3 check_results.py de --truth $T --table deseq2_results.csv
```

Each command prints PASS or FAIL per check and exits non-zero on any failure,
so it can be the last line of a Slurm test job. `names` also works on real data.

## Settings the test data needs

These differ from the real runs. Check each flag against `--help` for your
installed version.

- **HISAT2 index:** build one from `test_genome.fa` (seconds).
- **MACS2 genome size:** `-g 316000`.
- **RNA-seq strandedness:** reverse. HISAT2 `--rna-strandness RF`;
  featureCounts `-s 2`; htseq-count `-s reverse`.
- **featureCounts pairs:** count fragments, not reads (`-p`, plus
  `--countReadPairs` in newer versions).
- **Counting feature:** the GTF has `gene`, `transcript` and `exon` rows with
  `gene_id`; count on `exon` grouped by `gene_id`.
- **DESeq2 contrast:** SF vs BF. If yours is BF vs SF, add `--flip`.
- **Array sizes:** 10 CUT&RUN and 6 RNA-seq samples; sample lists are
  `samples_cutrun.txt` and `samples_rnaseq.txt`.

## What the checks catch

| Check | Catches |
| --- | --- |
| `names` | Contig names that differ between FASTA, GTF and BAM (the zero-count bug) |
| `fastq` | Lost reads, R1/R2 out of sync, missing or duplicated sample files |
| `aligned` | Low mapping rate, wrong index, pairs placed in the wrong position |
| `peaks` | Missed peaks, false peaks, summits off-centre, wrong contig names |
| `broad` | Domains fragmented or missed (for example `--broad` left off), calls outside domains |
| `counts` | Wrong strandedness, reads counted instead of fragments, wrong feature type, dropped Mt gene |
| `de` | Reversed contrast, missed or false DE genes, wrong fold-change size |

## What has and has not been tested

Tested when this kit was built:

- The generator is deterministic, and every CUT&RUN read pair matches the
  genome at the position in its name.
- Perfectly aligned BAMs made from the reads, run through MACS3 with the input
  control and `-g 316000`, recovered 12 of 12 AC peaks with 0 false peaks and
  all 3 ME domains at 97-99% coverage.
- The checker passes correct count and DE tables and fails broken ones (all
  zeros, doubled counts, reversed contrast, renamed contig).

Not tested: an actual run through fastp, HISAT2, featureCounts, htseq-count or
DESeq2. Expect to adjust a flag or tolerance on the first run.

## Limits

- The genome is random sequence with no repeats, so this tests plumbing and
  settings, not multimapper or duplicate handling.
- Peaks are strong on purpose. The test will not tell you whether thresholds
  are well tuned for real data; in testing, peaks still passed with the wrong
  `-g` and with no input control.
- Both groups have the same peaks, so there is no differential-binding truth.
- Reads have no adapters or sequencing errors, so trimming should keep
  essentially everything.
