#!/usr/bin/env bash

#SBATCH --partition=acpu
#SBATCH --job-name=plotHeatMap_By_Target_Day
#SBATCH --output=logs/%x.%j.log
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=4:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu

################################################################################
# STEP 4.3: Heatmaps of ALL replicates, one heatmap per target and day
#
# RUN THIS AFTER Step 4.1.2 has made matrix_combined_gene_only_bed_file.gz.
#
# USAGE (a title is required):
#   sbatch 04.3_HeatMap_Individual_Bams_By_Treatment_updated.sh "Base Title"
#   Each heatmap is titled "<Base Title> - <TARGET> <DAY>".
#   The script exits immediately with an error if no title is given.
#
# INPUT:   matrix_combined_gene_only_bed_file.gz (from Step 4.1.2, in 04_Deeptools)
#            The sample list is read from the matrix's own header, so no other
#            input files are needed.
# OUTPUT:  6 heatmaps, one per TARGET (AC, ME, SR) x DAY (D1, D3):
#            04_Deeptools/Figures/matrix_<TARGET>_<DAY>_all_pops_heatmap.png
#          6 subset matrices in 04_Deeptools:
#            matrix_<TARGET>_<DAY>_all_pops.gz
# ENV:     deeptools_kernel_v2 (deepTools computeMatrixOperations, plotHeatmap)
# RESOURCES: 8 CPUs, 32 GB, 4 h
#
# WHAT EACH HEATMAP SHOWS
#   For one target and day, EVERY replicate of every group (BF, MP, SF, ZH)
#   appears as its own panel around the TSS, in the order the samples appear in
#   the combined matrix. TARGET: AC = H3K27ac, ME = H3K9me, SR = SREBP.
#   Compare with 04.3_Heat_Maps.sh, which shows one representative replicate
#   per group.
#
# DIRECTORIES
#   Run from 04_Deeptools. Job logs go to logs/ in the folder
#   you submit from. Everything else uses full paths.
#   DEEPTOOLS_DIR = <path>/04_Deeptools
#   MATRIX        = DEEPTOOLS_DIR/matrix_combined_gene_only_bed_file.gz
#   OUTDIR        = DEEPTOOLS_DIR/Figures (created automatically)
#
# BEFORE SUBMITTING (run from 04_Deeptools)
#   1. mkdir -p logs
#        (SLURM opens the log files before the script starts, so these folders
#         must exist ahead of time or the job fails with no log.)
#   2. Confirm Step 4.1.2 finished:
#        ls -lh <DEEPTOOLS_DIR>/matrix_combined_gene_only_bed_file.gz
#   3. sbatch 04.3_HeatMap_Individual_Bams_By_Treatment_updated.sh "Your Base Title"
#
# NOT AN ARRAY JOB
#   There is no array size to change. One job loops over the target/day pairs
#   listed in TREATMENTS. New samples in the combined matrix are picked up
#   automatically, because the sample list comes from the matrix header.
#
# HOW SAMPLES ARE CHOSEN
#   Sample names in the matrix look like <GROUP>-<TARGET>-<DAY>-R<n> (e.g.
#   BF-AC-D1-R1). For each target/day pair, every sample whose name contains
#   -<TARGET>-<DAY>- is included. A pair with no matching samples is skipped
#   and the job exits with an error at the end.
#
# THINGS YOU MAY WANT TO EDIT
#   TREATMENTS   The target/day pairs to plot. Add or remove pairs here.
#                (The name is historical; each entry is a target and a day.)
#   COLORMAP, Z_MIN, Z_MAX, REF_POINT_LABEL, WHAT_TO_SHOW   Heatmap appearance.
#
# CHECKING THE RESULTS
#   ls <OUTDIR>/*_all_pops_heatmap.png     # 6 heatmaps
#   grep "samples:" logs/*.log          # which samples went into each plot
#   grep "skipping" logs/*.log          # pairs that were left out
#   The log ends with a list of the output files and their sizes.
#
# NOTES
#   - The script stops on any error (set -e), for example if plotHeatmap fails.
#   - Paths to check before running: SCRATCH, DEEPTOOLS_DIR, and the email above.
################################################################################

set -e

# ----------------------------
# Environment setup
# ----------------------------
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh
conda activate ${SCRATCH}/conda_envs/deeptools_kernel_v2

# ----------------------------
# Input/output paths
# ----------------------------
DEEPTOOLS_DIR="${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/04_Deeptools"
MATRIX="${DEEPTOOLS_DIR}/matrix_combined_gene_only_bed_file.gz"
OUTDIR="${DEEPTOOLS_DIR}/Figures"
mkdir -p "$OUTDIR"

# ----------------------------
# Target x Day pairs to process
# ----------------------------
TREATMENTS=(
    "AC D1"
    "AC D3"
    "ME D1"
    "ME D3"
    "SR D1"
    "SR D3"
)

# ----------------------------
# Checks: tools, input matrix, title argument
# ----------------------------
if ! command -v plotHeatmap &> /dev/null; then
    echo "ERROR: plotHeatmap not found. Install with: conda install deeptools"
    exit 1
fi

if ! command -v computeMatrixOperations &> /dev/null; then
    echo "ERROR: computeMatrixOperations not found."
    exit 1
fi

if [[ ! -f "$MATRIX" ]]; then
    echo "ERROR: Combined matrix not found: $MATRIX"
    exit 1
fi

if [[ $# -lt 1 ]]; then
    echo "ERROR: Missing title argument"
    echo "USAGE: sbatch 04.3_HeatMap_Individual_Bams_By_Treatment_updated.sh 'Base Title'"
    exit 1
fi
BASE_TITLE=$1

# ============================================================
# HEATMAP SETTINGS
# ============================================================
COLORMAP="RdBu"
Z_MIN="auto"
Z_MAX="auto"
REF_POINT_LABEL="TSS"
WHAT_TO_SHOW="heatmap and colorbar"

# ============================================================
# EXTRACT ALL SAMPLE LABELS FROM THE COMBINED MATRIX'S HEADER
# ============================================================
get_sample_labels() {
    local matrix="$1"
    zcat "$matrix" | head -n 1 | python3 -c '
import sys, json
line = sys.stdin.readline().strip()
if line.startswith("@"):
    line = line[1:]
data = json.loads(line)
labels = data.get("sample_labels", [])
for l in labels:
    print(l)
'
}

mapfile -t ALL_SAMPLES < <(get_sample_labels "$MATRIX")

if [[ ${#ALL_SAMPLES[@]} -eq 0 ]]; then
    echo "ERROR: Could not parse sample labels from $MATRIX"
    exit 1
fi

echo "Found ${#ALL_SAMPLES[@]} samples total in combined matrix."

# ============================================================
# MAIN LOOP - ONE HEATMAP PER TARGET x DAY
# ============================================================
FAILED=0

for pair in "${TREATMENTS[@]}"; do
    TARGET=$(echo "$pair" | awk '{print $1}')
    DAY=$(echo "$pair" | awk '{print $2}')

    echo "============================================================"
    echo "TARGET: ${TARGET}   DAY: ${DAY}"
    echo "============================================================"

    # Filter to samples matching this target AND this day
    # Sample names look like: BF-AC-D1-R1  (site-target-day-rep)
    GROUP_SAMPLES=()
    for s in "${ALL_SAMPLES[@]}"; do
        [[ "$s" == *-${TARGET}-${DAY}-* ]] && GROUP_SAMPLES+=("$s")
    done

    if [[ ${#GROUP_SAMPLES[@]} -eq 0 ]]; then
        echo "  No samples found for ${TARGET} ${DAY} — skipping"
        FAILED=1
        echo ""
        continue
    fi

    echo "  ${#GROUP_SAMPLES[@]} samples: ${GROUP_SAMPLES[*]}"

    GROUP_MATRIX="${DEEPTOOLS_DIR}/matrix_${TARGET}_${DAY}_all_pops.gz"
    TITLE="${BASE_TITLE} - ${TARGET} ${DAY}"
    filename="matrix_${TARGET}_${DAY}_all_pops"

    # Subset the combined matrix to this target/day
    echo "  Subsetting matrix to ${TARGET} ${DAY} samples..."
    computeMatrixOperations subset \
        -m "$MATRIX" \
        --samples "${GROUP_SAMPLES[@]}" \
        -o "$GROUP_MATRIX"

    if [[ ! -f "$GROUP_MATRIX" ]]; then
        echo "  ERROR: Subset failed for ${TARGET} ${DAY}, skipping"
        FAILED=1
        echo ""
        continue
    fi

    # Plot the heatmap
    echo "  Creating heatmap..."
    plotHeatmap \
        -m "$GROUP_MATRIX" \
        -out "${OUTDIR}/${filename}_heatmap.png" \
        --colorMap "$COLORMAP" \
        --whatToShow "$WHAT_TO_SHOW" \
        --refPointLabel "$REF_POINT_LABEL" \
        --samplesLabel "${GROUP_SAMPLES[@]}" \
        -T "$TITLE" \
        --zMin "$Z_MIN" \
        --zMax "$Z_MAX" \
        --dpi 300 \
        --regionsLabel " "

    if [[ $? -eq 0 ]]; then
        echo "  ✅ Done: ${OUTDIR}/${filename}_heatmap.png"
        ls -lh "${OUTDIR}/${filename}_heatmap.png"
    else
        echo "  🚨 Heatmap creation FAILED for ${TARGET} ${DAY}"
        FAILED=1
    fi
    echo ""
done

# ============================================================
# SUMMARY
# ============================================================
echo "============================================================"
echo "VISUALIZATION COMPLETE"
echo "============================================================"
echo ""
echo "Output files in ${OUTDIR}/:"
ls -lh "${OUTDIR}/"*_all_pops_heatmap.png 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
echo ""

if [[ "$FAILED" -eq 1 ]]; then
    echo "One or more target/day combinations failed or were skipped - see warnings above."
    exit 1
fi
