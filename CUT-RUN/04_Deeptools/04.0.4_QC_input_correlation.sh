#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=QC_input_corr
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

################################################################################
# QC: do the input replicates of each site and day agree?
# Run after Step 2.3.1 (needs the filtered individual input BAMs), before deciding how to
# normalize the bigWigs (Step 4.0): pooled same-day input vs one input each.
#
# Counts reads in 10 kb bins for every input BAM (multiBamSummary), correlates
# all inputs with each other (Spearman and Pearson on log counts), then for each
# site x day trio reports each input's mean correlation with its two siblings.
#
# OUTPUT in Round1_Round2_With_Duplicates/QC/input_correlation/:
#   inputs_spearman_heatmap.png, inputs_pearson_heatmap.png   clustered heatmaps
#   inputs_spearman.tab, inputs_pearson.tab                   correlation matrices
#   input_trio_summary.tsv                                    one row per input
# A row is flagged OUTLIER when its mean Spearman correlation with its siblings
# is below MIN_COR, or more than GAP below the best member of its trio.
# Inputs should cluster by site and day on the heatmap; one that sits apart from
# its trio is a candidate to leave out of the pooled input.
#
# Run from 04_Deeptools (it needs logs/):   sbatch 04.0.4_QC_input_correlation.sh
################################################################################

set -eo pipefail
MIN_COR=0.80
GAP=0.10

SCRATCH=/scratch/alpine/.colostate.edu/c832500103
PROJ=/projects/.colostate.edu/c832500103
PIPE=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates
BAMDIR=${PIPE}/02_Cxt_Alignment/02.3.1_Cxt_BamFiltered
OUT=${PIPE}/QC/input_correlation
mkdir -p "$OUT"

source ${PROJ}/miniforge3/etc/profile.d/conda.sh
conda activate ${PROJ}/conda_envs/deeptools_kernel_v2

BAMS=$(ls ${BAMDIR}/*-IN-*-Cxt.bam | sort -V)
LABELS=$(for b in $BAMS; do basename "$b" -Cxt.bam; done)
echo "Inputs: $(echo $BAMS | wc -w)"

multiBamSummary bins -b $BAMS -l $LABELS --binSize 10000 -p 8 \
    -o "$OUT/inputs_10kb.npz" --outRawCounts "$OUT/inputs_10kb_counts.tab"

plotCorrelation -in "$OUT/inputs_10kb.npz" -c spearman -p heatmap --skipZeros --plotNumbers \
    --plotTitle "Input BAMs, Spearman (10 kb bins)" -o "$OUT/inputs_spearman_heatmap.png" \
    --outFileCorMatrix "$OUT/inputs_spearman.tab"
plotCorrelation -in "$OUT/inputs_10kb.npz" -c pearson --log1p -p heatmap --skipZeros --plotNumbers \
    --plotTitle "Input BAMs, Pearson log1p (10 kb bins)" -o "$OUT/inputs_pearson_heatmap.png" \
    --outFileCorMatrix "$OUT/inputs_pearson.tab"

python - "$OUT" "$MIN_COR" "$GAP" <<'EOF'
import sys
out, min_cor, gap = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
def read(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path) if not l.startswith("#")]
    labels = [x.strip("'") for x in rows[0][1:]]
    m = {r[0].strip("'"): dict(zip(labels, map(float, r[1:]))) for r in rows[1:]}
    return labels, m
labels, sp = read(f"{out}/inputs_spearman.tab")
_, pe = read(f"{out}/inputs_pearson.tab")
grp = lambda s: "-".join(s.split("-")[:3]).replace("-IN-", "-")   # BF-IN-D1-R1 -> BF-D1
site = lambda s: s.split("-")[0]
rows = []
for s in labels:
    sib = [t for t in labels if t != s and grp(t) == grp(s)]
    other_day = [t for t in labels if site(t) == site(s) and grp(t) != grp(s)]
    other_site = [t for t in labels if site(t) != site(s)]
    mean = lambda d, ts: sum(d[s][t] for t in ts) / len(ts) if ts else float("nan")
    rows.append([s, grp(s), mean(sp, sib), mean(pe, sib), mean(sp, other_day), mean(sp, other_site)])
best = {}
for r in rows: best[r[1]] = max(best.get(r[1], -1), r[2])
with open(f"{out}/input_trio_summary.tsv", "w") as f:
    f.write("input\tsite_day\tspearman_vs_siblings\tpearson_vs_siblings\tspearman_vs_other_day\tspearman_vs_other_sites\tflag\n")
    for r in rows:
        flag = "OUTLIER" if (r[2] < min_cor or best[r[1]] - r[2] > gap) else "ok"
        f.write("%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n" % (*r, flag))
        print("%-14s %-6s sib_sp=%.3f sib_pe=%.3f other_day=%.3f other_site=%.3f %s" % (*r, flag))
EOF
echo "Done. Summary: $OUT/input_trio_summary.tsv   Heatmaps: $OUT/inputs_*_heatmap.png"
