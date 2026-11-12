#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=All_Marks_By_Day_Profile
#SBATCH --output=out_log/%x.%A-%a.log
#SBATCH --error=err_log/%x.%A-%a.err
#SBATCH --time=1:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shawn.yates@colostate.edu

################################################################################
# STEP 4.2: Profile plots of all four groups together, one plot per target and day
#           (representative replicates only)
#
# RUN THIS AFTER Step 4.1.2 has made matrix_combined_gene_only_bed_file.gz.
#
# INPUT:   matrix_combined_gene_only_bed_file.gz (from Step 4.1.2, in 04_Deeptools)
#          BigWigs in 00_BigWig_Files_Individual_Bams (used only to look up the
#            sample names, which must match the matrix sample labels)
# OUTPUT:  6 plots, one per TARGET (AC, ME, SR) x DAY (D1, D3):
#            04.2_Representative_Replicates/Figures/All_Marks_By_Day/<DAY>/<TARGET>_<DAY>_AllMarks.png
#          6 subset matrices in 04_Deeptools (kept for the heatmap step):
#            matrix_<TARGET>_<DAY>_AllMarks.gz
# ENV:     deeptools_kernel_v2 (deepTools computeMatrixOperations, plotProfile)
# RESOURCES: 2 CPUs, 64 GB, 1 h
#
# WHAT EACH PLOT SHOWS
#   For one target and day, the representative replicate of each group
#   (BF, MP, ZH, SF) is plotted on one graph, one colored line per group,
#   all on the same fixed y-axis so plots can be compared directly.
#   TARGET: AC = H3K27ac, ME = H3K9me, SR = SREBP.
#
# DIRECTORIES
#   Run from 04_Deeptools/04.2_Representative_Replicates. Job logs go to out_log/
#   and err_log/ in the folder you submit from. Everything else uses full paths.
#   DEEPTOOLS_DIR = <path>/04_Deeptools
#   BW_DIR        = DEEPTOOLS_DIR/00_BigWig_Files_Individual_Bams
#   OUTDIR        = DEEPTOOLS_DIR/04.2_Representative_Replicates/Figures/All_Marks_By_Day
#                   (created automatically, with D1/ and D3/ subfolders)
#
# BEFORE SUBMITTING (run from 04.2_Representative_Replicates)
#   1. mkdir -p out_log err_log
#        (SLURM opens the log files before the script starts, so these folders
#         must exist ahead of time or the job fails with no log.)
#   2. Confirm Step 4.1.2 finished:
#        ls -lh <DEEPTOOLS_DIR>/matrix_combined_gene_only_bed_file.gz
#   3. sbatch 04.2_All_Marks_By_Day_Profile_HardCoded_Yaxis.sh
#
# NOT AN ARRAY JOB
#   There is no array size to change. One job loops over all target/day
#   combinations.
#
# THINGS YOU MAY WANT TO EDIT
#   REP_MAP       Which replicate represents each group/target/day, in the form
#                 [GROUP_TARGET_DAY]="R<n>". Every combination needs an entry. A
#                 missing entry, or a replicate that isn't in the matrix, causes
#                 that group to be left off the plot (a message is printed and
#                 the job still finishes).
#   YMIN, YMAX    Shared y-axis range for every plot (currently -0.10 to 0.55,
#                 chosen by looking at earlier plots). Change these if new data
#                 falls outside the range.
#   MARK_COLORS   One color per group, the same on every plot.
#   Adding or removing a group: update MARK_COLORS, REP_MAP, the
#                 "for MARK in BF MP ZH SF" loop, and the plot title (-T).
#   Adding a target or day: update the "for TARGET" / "for DAY" loops and give
#                 REP_MAP entries for the new combinations.
#
# CHECKING THE RESULTS
#   ls <OUTDIR>/D1 <OUTDIR>/D3          # 3 plots in each folder
#   grep "skipping" out_log/*.log       # groups or plots that were left out
#   grep "using" out_log/*.log          # which sample was chosen for each group
#
# NOTES
#   - The script stops on any error (set -e), for example if a sample name is
#     not found in the combined matrix.
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
#MATRIX="${DEEPTOOLS_DIR}/bigwiglist_Individual_Bams_gene_only_bed_file.gz"
MATRIX="${DEEPTOOLS_DIR}/matrix_combined_gene_only_bed_file.gz"
OUTDIR="${DEEPTOOLS_DIR}/04.2_Representative_Replicates/Figures/All_Marks_By_Day"
mkdir -p "$OUTDIR" out_log err_log

# Hard-coded shared y-axis range for every graph (from visual inspection of prior plots)
YMIN="-0.10"
YMAX="0.55"

# ============================================================
# One color per mark (BF, MP, ZH, SF) — consistent across every plot
# ============================================================
declare -A MARK_COLORS=(
    [BF]="darkorange"
    [MP]="forestgreen"
    [ZH]="crimson"
    [SF]="purple"
)

# ============================================================
# Representative replicate for every MARK_TARGET_DAY combination
# Note: actual sample names use "R1", "R2", "R3"... not "Rep1", "Rep2"...
# ============================================================
declare -A REP_MAP=(
    [BF_AC_D1]="R2"  [BF_AC_D3]="R4"  [BF_ME_D1]="R2"  [BF_ME_D3]="R3"  [BF_SR_D1]="R3"  [BF_SR_D3]="R2"
    [MP_AC_D1]="R2"  [MP_AC_D3]="R1"  [MP_ME_D1]="R1"  [MP_ME_D3]="R1"  [MP_SR_D1]="R3"  [MP_SR_D3]="R3"
    [ZH_AC_D1]="R3"  [ZH_AC_D3]="R3"  [ZH_ME_D1]="R3"  [ZH_ME_D3]="R3"  [ZH_SR_D1]="R1"  [ZH_SR_D3]="R3"
    [SF_AC_D1]="R2"  [SF_AC_D3]="R3"  [SF_ME_D1]="R2"  [SF_ME_D3]="R3"  [SF_SR_D1]="R3"  [SF_SR_D3]="R2"
)

# ============================================================
# For every TARGET + DAY, gather one representative sample per MARK
# and plot all four marks together in a single graph, sharing a
# fixed y-axis range across every graph
# ============================================================
for TARGET in AC ME SR; do
    for DAY in D1 D3; do
        echo "============================================================"
        echo "TARGET: $TARGET   DAY: $DAY"
        echo "============================================================"

        SAMPLES=()
        LABELS=()
        COLORS=()

        for MARK in BF MP ZH SF; do
            REP="${REP_MAP[${MARK}_${TARGET}_${DAY}]:-}"

            if [[ -z "$REP" ]]; then
                echo "  No representative defined for ${MARK} ${TARGET} ${DAY} — skipping mark"
                continue
            fi

            # All sample names for this group/target/day (from the bigWig file names)
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
            COLORS+=("${MARK_COLORS[$MARK]}")
        done

        if [[ ${#SAMPLES[@]} -eq 0 ]]; then
            echo "  No marks found for ${TARGET} ${DAY} — skipping"
            continue
        fi

        DAY_DIR="${OUTDIR}/${DAY}"
        mkdir -p "$DAY_DIR"

        # Subset the combined matrix to the chosen samples
        MATRIX_TD="${DEEPTOOLS_DIR}/matrix_${TARGET}_${DAY}_AllMarks.gz"
        computeMatrixOperations subset \
            -m "$MATRIX" \
            --samples "${SAMPLES[@]}" \
            -o "$MATRIX_TD"

        # Plot all groups together on the fixed y-axis
        plotProfile \
            -m "$MATRIX_TD" \
            -out "${DAY_DIR}/${TARGET}_${DAY}_AllMarks.png" \
            --perGroup \
            --colors "${COLORS[@]}" \
            --samplesLabel "${LABELS[@]}" \
            --regionsLabel "" \
            -T "${TARGET} ${DAY} — All Marks (BF, MP, ZH, SF)" \
            -y "Reads per Genomic Content" \
            --yMin "$YMIN" \
            --yMax "$YMAX" \
            --plotHeight 8 \
            --plotWidth 12 \
            --dpi 300

        echo "  Created: ${DAY_DIR}/${TARGET}_${DAY}_AllMarks.png"
    done
done

echo "All target/day AllMarks plots complete with fixed y-axis (${YMIN} to ${YMAX})."
