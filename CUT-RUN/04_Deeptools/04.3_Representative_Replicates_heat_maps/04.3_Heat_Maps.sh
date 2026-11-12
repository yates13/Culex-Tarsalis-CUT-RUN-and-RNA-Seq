#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=All_Marks_By_Day_Heatmap
#SBATCH --output=out_log/%x.%A-%a.log
#SBATCH --error=err_log/%x.%A-%a.err
#SBATCH --time=4:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shawn.yates@colostate.edu

################################################################################
# STEP 4.3: Heatmaps of all four groups together, one heatmap per target and day
#           (representative replicates only)
#
# RUN THIS AFTER Step 4.1.2 has made matrix_combined_gene_only_bed_file.gz.
#
# USAGE (a title is required):
#   sbatch 04.3_Heat_Maps.sh "Base Title"
#   Each heatmap is titled "<Base Title> - <TARGET> <DAY>".
#   The script exits immediately with an error if no title is given.
#
# INPUT:   matrix_combined_gene_only_bed_file.gz (from Step 4.1.2, in 04_Deeptools)
#          BigWigs in 00_BigWig_Files_Individual_Bams (used only to look up the
#            sample names, which must match the matrix sample labels)
# OUTPUT:  6 heatmaps, one per TARGET (AC, ME, SR) x DAY (D1, D3):
#            04.3_Representative_Replicates_heat_maps/Figures/All_Marks_By_Day_Heatmap/<DAY>/<TARGET>_<DAY>_AllMarks_heatmap.png
#          6 subset matrices in 04_Deeptools:
#            matrix_<TARGET>_<DAY>_AllMarks_Heatmap.gz
# ENV:     deeptools_kernel_v2 (deepTools computeMatrixOperations, plotHeatmap)
# RESOURCES: 8 CPUs, 64 GB, 4 h
#
# WHAT EACH HEATMAP SHOWS
#   For one target and day, the representative replicate of each group
#   (SF, BF, MP, ZH) is shown as its own side-by-side panel around the TSS,
#   in that order. TARGET: AC = H3K27ac, ME = H3K9me3, SR = SREBP.
#
# DIRECTORIES
#   Run from 04_Deeptools/04.3_Representative_Replicates_heat_maps. Job logs go
#   to out_log/ and err_log/ in the folder you submit from. Everything else uses
#   full paths.
#   DEEPTOOLS_DIR = <path>/04_Deeptools
#   BW_DIR        = DEEPTOOLS_DIR/00_BigWig_Files_Individual_Bams
#   MATRIX        = DEEPTOOLS_DIR/matrix_combined_gene_only_bed_file.gz
#   OUTDIR        = DEEPTOOLS_DIR/04.3_Representative_Replicates_heat_maps/Figures/All_Marks_By_Day_Heatmap
#                   (created automatically, with D1/ and D3/ subfolders)
#
# BEFORE SUBMITTING (run from 04.3_Representative_Replicates_heat_maps)
#   1. mkdir -p out_log err_log
#        (SLURM opens the log files before the script starts, so these folders
#         must exist ahead of time or the job fails with no log.)
#   2. Confirm Step 4.1.2 finished:
#        ls -lh <DEEPTOOLS_DIR>/matrix_combined_gene_only_bed_file.gz
#   3. sbatch 04.3_Heat_Maps.sh "Your Base Title"
#
# NOT AN ARRAY JOB
#   There is no array size to change. One job loops over all target/day
#   combinations.
#
# THINGS YOU MAY WANT TO EDIT
#   REP_MAP     Which replicate represents each group/target/day, in the form
#               [GROUP_TARGET_DAY]="R<n>". KEEP THIS IDENTICAL TO REP_MAP in
#               04.2_All_Marks_By_Day_Profile_HardCoded_Yaxis.sh so the heatmaps
#               and profile plots show the same samples. If you change one,
#               change the other. A missing entry, or a replicate that isn't
#               found, leaves that group off the heatmap (a message is printed).
#   COLORMAP, Z_MIN, Z_MAX, REF_POINT_LABEL, WHAT_TO_SHOW   Heatmap appearance.
#   Panel order  Set by the "for MARK in SF BF MP ZH" loop. Change the order there.
#   Adding or removing a group: update REP_MAP and that loop.
#   Adding a target or day: update the "for TARGET" / "for DAY" loops and give
#               REP_MAP entries for the new combinations.
#
# CHECKING THE RESULTS
#   ls <OUTDIR>/D1 <OUTDIR>/D3          # 3 heatmaps in each folder
#   grep "skipping" out_log/*.log       # groups or plots that were left out
#   grep "using" out_log/*.log          # which sample was chosen for each group
#   The job exits with an error at the end if any target/day was skipped or failed.
#
# NOTES
#   - The script stops on any error (set -e), for example if plotHeatmap fails.
#   - Sample names look like <GROUP>-<TARGET>-<DAY>-R<n>; replicates are matched
#     by the R<n> part ("R1", "R2", ...), not "Rep1".
#   - Paths to check before running: SCRATCH, DEEPTOOLS_DIR, and the email above.
################################################################################

set -e

# ----------------------------
# Environment setup
# ----------------------------
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
PROJ=/projects/c832500103@colostate.edu
source ${PROJ}/miniforge3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${PROJ}/miniforge3" >&2; exit 1; }
conda activate ${PROJ}/conda_envs/deeptools_kernel_v2 || { echo "ERROR: conda activate deeptools_kernel_v2 failed" >&2; exit 1; }

# ----------------------------
# Input/output paths
# ----------------------------
DEEPTOOLS_DIR="${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/04_Deeptools"
BW_DIR="${DEEPTOOLS_DIR}/00_BigWig_Files_Individual_Bams"
MATRIX="${DEEPTOOLS_DIR}/matrix_combined_gene_only_bed_file.gz"
OUTDIR="${DEEPTOOLS_DIR}/04.3_Representative_Replicates_heat_maps/Figures/All_Marks_By_Day_Heatmap"
mkdir -p "$OUTDIR" out_log err_log

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
    echo "USAGE: sbatch 04.7_All_Marks_By_Day_Heatmap.sh 'Base Title'"
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
# Representative replicate for every MARK_TARGET_DAY combination
# (identical to 04.6 — keep these two scripts in sync if you
# update representative reps in one, update the other)
# ============================================================
declare -A REP_MAP=(
    [BF_AC_D1]="R2"  [BF_AC_D3]="R4"  [BF_ME_D1]="R2"  [BF_ME_D3]="R3"  [BF_SR_D1]="R3"  [BF_SR_D3]="R2"
    [MP_AC_D1]="R2"  [MP_AC_D3]="R1"  [MP_ME_D1]="R1"  [MP_ME_D3]="R1"  [MP_SR_D1]="R3"  [MP_SR_D3]="R3"
    [ZH_AC_D1]="R3"  [ZH_AC_D3]="R3"  [ZH_ME_D1]="R3"  [ZH_ME_D3]="R3"  [ZH_SR_D1]="R1"  [ZH_SR_D3]="R3"
    [SF_AC_D1]="R2"  [SF_AC_D3]="R3"  [SF_ME_D1]="R2"  [SF_ME_D3]="R3"  [SF_SR_D1]="R3"  [SF_SR_D3]="R2"
)

# ============================================================
# For every TARGET + DAY, gather one representative sample per
# MARK and plot all four marks as panels in a single heatmap
# ============================================================
FAILED=0

for TARGET in AC ME SR; do
    for DAY in D1 D3; do
        echo "============================================================"
        echo "TARGET: $TARGET   DAY: $DAY"
        echo "============================================================"

        SAMPLES=()
        LABELS=()

        for MARK in SF BF MP ZH; do
            REP="${REP_MAP[${MARK}_${TARGET}_${DAY}]:-}"

            if [[ -z "$REP" ]]; then
                echo "  No representative defined for ${MARK} ${TARGET} ${DAY} — skipping mark"
                continue
            fi

            mapfile -t MARK_SAMPLES < <(
                ls "${BW_DIR}"/${MARK}-${TARGET}-${DAY}-*_inputsubtracted.bw 2>/dev/null \
                    | xargs -n1 basename \
                    | sed 's/_inputsubtracted\.bw$//' \
                    | sort
            )

            SAMPLE=""
            for s in "${MARK_SAMPLES[@]}"; do
                if [[ "$s" == *"-${REP}" || "$s" == *"-${REP}-"* ]]; then
                    SAMPLE="$s"
                    break
                fi
            done

            if [[ -z "$SAMPLE" ]]; then
                echo "  Representative '${REP}' not found among ${MARK} samples for ${TARGET} ${DAY} — skipping mark"
                continue
            fi

            echo "  ${MARK}: using ${SAMPLE}"
            SAMPLES+=("$SAMPLE")
            LABELS+=("$MARK")
        done

        if [[ ${#SAMPLES[@]} -eq 0 ]]; then
            echo "  No marks found for ${TARGET} ${DAY} — skipping"
            FAILED=1
            continue
        fi

        DAY_DIR="${OUTDIR}/${DAY}"
        mkdir -p "$DAY_DIR"

        MATRIX_TD="${DEEPTOOLS_DIR}/matrix_${TARGET}_${DAY}_AllMarks_Heatmap.gz"
        TITLE="${BASE_TITLE} - ${TARGET} ${DAY}"

        echo "  Subsetting matrix to representative samples..."
        computeMatrixOperations subset \
            -m "$MATRIX" \
            --samples "${SAMPLES[@]}" \
            -o "$MATRIX_TD"

        if [[ ! -f "$MATRIX_TD" ]]; then
            echo "  ERROR: Subset failed for ${TARGET} ${DAY}, skipping"
            FAILED=1
            continue
        fi

        echo "  Creating heatmap..."
        plotHeatmap \
            -m "$MATRIX_TD" \
            -out "${DAY_DIR}/${TARGET}_${DAY}_AllMarks_heatmap.png" \
            --colorMap "$COLORMAP" \
            --whatToShow "$WHAT_TO_SHOW" \
            --refPointLabel "$REF_POINT_LABEL" \
            --samplesLabel "${LABELS[@]}" \
            -T "$TITLE" \
            --zMin "$Z_MIN" \
            --zMax "$Z_MAX" \
            --dpi 300 \
            --regionsLabel " "

        if [[ $? -eq 0 ]]; then
            echo " Created: ${DAY_DIR}/${TARGET}_${DAY}_AllMarks_heatmap.png"
        else
            echo " Heatmap creation FAILED for ${TARGET} ${DAY}"
            FAILED=1
        fi
        echo ""
    done
done

echo "============================================================"
echo "All target/day AllMarks heatmaps complete."
echo "============================================================"

if [[ "$FAILED" -eq 1 ]]; then
    echo "One or more target/day combinations failed or were skipped - see warnings above."
    exit 1
fi
