#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=artifact_candidates_check
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=01:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

################################################################################
# ONE-OFF CHECK (not a pipeline step; nothing in the pipeline is changed)
#
# Question: do the top Greenscreen-overlap genes meet the gene5911 rule,
#   i.e. extreme signal in ALL targets (AC, ME, SR), both days and all sites?
#
# For every sample in the combined matrix: mean signal per gene over the
# +/-2 kb window (nan skipped), then for each candidate gene:
#   fold  = gene mean / median gene mean in that sample
#   pct   = percentile rank of the gene in that sample (100 = highest)
#
# INPUT:   04_Deeptools/matrix_combined_gene_only_bed_file.gz (Step 4.1.2, 11/23)
# OUTPUT:  QC/oneoff_checks/artifact_candidates/  (NEW folder; refuses to run if it exists)
#            per_sample.tsv   gene, sample, target, day, site, mean, fold, pct
#            summary.txt      per gene and target x day: samples, median fold,
#                             samples in the top 1%
# ENV:     deeptools_kernel_v2 (numpy)
# NOTE:    Run "mkdir -p logs" in the folder you submit from.
################################################################################

set -eo pipefail

PROJ=/projects/c832500103@colostate.edu
source ${PROJ}/miniforge3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found" >&2; exit 1; }
conda activate ${PROJ}/conda_envs/deeptools_kernel_v2 || { echo "ERROR: conda activate deeptools_kernel_v2 failed" >&2; exit 1; }

RUN=/scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates
MATRIX=${RUN}/04_Deeptools/matrix_combined_gene_only_bed_file.gz
OUT=${RUN}/QC/oneoff_checks/artifact_candidates
GENES="gene12689 gene2811 gene223 gene225 gene248 gene226 gene224 gene247 gene2805 gene222"

[[ -s "$MATRIX" ]] || { echo "ERROR: matrix not found: $MATRIX" >&2; exit 1; }
if [[ -e "$OUT" ]]; then
    echo "ERROR: $OUT already exists; rename it first. Nothing was overwritten." >&2
    exit 1
fi
mkdir -p "$OUT"

python3 - "$MATRIX" "$OUT" $GENES <<'EOF'
import sys, gzip, json
import numpy as np

matrix, out = sys.argv[1], sys.argv[2]
cands = sys.argv[3:]

with gzip.open(matrix, 'rt') as fh:
    header = json.loads(fh.readline()[1:])
    labels = header['sample_labels']
    bounds = header['sample_boundaries']
    names, means = [], []
    for line in fh:
        p = line.rstrip('\n').split('\t')
        v = np.array(p[6:], dtype=float)
        names.append(p[3])
        with np.errstate(all='ignore'):
            means.append([np.nanmean(v[bounds[i]:bounds[i+1]]) for i in range(len(labels))])

means = np.array(means)              # genes x samples
print(f"rows {means.shape[0]}, samples {means.shape[1]}")
missing = [g for g in cands if g not in names]
if missing:
    print("WARNING: not in matrix:", " ".join(missing))

med = np.nanmedian(means, axis=0)
rows = []
for g in cands:
    if g not in names:
        continue
    r = names.index(g)
    for s, lab in enumerate(labels):
        col = means[:, s]
        x = means[r, s]
        ok = ~np.isnan(col)
        pct = 100.0 * np.sum(col[ok] <= x) / ok.sum() if not np.isnan(x) else float('nan')
        site, target, day = (lab.split('-') + ['?', '?', '?'])[:3]
        fold = x / med[s] if med[s] > 0 else float('nan')
        rows.append((g, lab, target, day, site, x, fold, pct))

with open(f"{out}/per_sample.tsv", 'w') as f:
    f.write("gene\tsample\ttarget\tday\tsite\tmean\tfold_vs_median\tpct_rank\n")
    for t in rows:
        f.write("%s\t%s\t%s\t%s\t%s\t%.4f\t%.1f\t%.2f\n" % t)

lines = ["Per gene, by target x day: samples | median fold vs a typical gene | samples in the top 1% | sites in the top 1%", ""]
for g in cands:
    sub = [t for t in rows if t[0] == g]
    if not sub:
        continue
    lines.append(g)
    tot_top = 0
    for target in ('AC', 'ME', 'SR'):
        for day in ('D1', 'D3'):
            ss = [t for t in sub if t[2] == target and t[3] == day]
            if not ss:
                continue
            top = [t for t in ss if t[7] >= 99]
            tot_top += len(top)
            sites = sorted({t[4] for t in top})
            lines.append("  %s %s: n=%d  median fold %.1f  top1%% %d/%d  sites %s"
                         % (target, day, len(ss), np.nanmedian([t[6] for t in ss]),
                            len(top), len(ss), ",".join(sites) if sites else "-"))
    lines.append("  TOTAL top1%%: %d/%d samples" % (tot_top, len(sub)))
    lines.append("")
txt = "\n".join(lines)
print(txt)
open(f"{out}/summary.txt", 'w').write(txt + "\n")
EOF

echo "Output: $OUT"
