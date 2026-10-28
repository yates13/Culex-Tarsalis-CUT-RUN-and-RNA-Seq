#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=MP_Targets_Individual_Bams_dynamic_naming
#SBATCH --output=logs/%x.%j.log
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=1:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu

################################################################################
# STEP 4.2: Profile plots of ALL replicates for the MP group, one plot per
#           target and day
#
# RUN THIS AFTER Step 4.1.2 has made matrix_combined_gene_only_bed_file.gz.
#
# INPUT:   matrix_combined_gene_only_bed_file.gz (from Step 4.1.2, in 04_Deeptools)
#          BigWigs in 00_BigWig_Files_Individual_Bams (used only to look up which
#            samples exist; names must match the matrix sample labels)
# OUTPUT:  6 plots, one per TARGET (AC, ME, SR) x DAY (D1, D3):
#            04_Deeptools/Figures/<TARGET>_<DAY>_MP_profile.png
#          12 subset matrices in 04_Deeptools (2 per plot):
#            matrix_<TARGET>_<DAY>_ALL.gz   (all groups for that target/day)
#            matrix_<TARGET>_<DAY>_MP.gz    (MP samples only; used for the plot)
# ENV:     deeptools_kernel_v2 (deepTools computeMatrixOperations, plotProfile)
# RESOURCES: 2 CPUs, 16 GB, 1 h
#
# WHAT EACH PLOT SHOWS
#   For one target and day, every MP replicate is drawn as its own colored line
#   around the TSS (-2 kb to +2 kb). TARGET: AC = H3K27ac, ME = H3K9me,
#   SR = SREBP. The y-axis is fixed (-0.3 to 0.8) so plots can be compared.
#
# ONE OF FOUR SIBLING SCRIPTS
#   04.2_BF_..., 04.2_MP_..., 04.2_SF_... and 04.2_ZH_Plot_profile_Individual_Bams_
#   Dynamic_sample_naming.sh are the same script for different groups. To make
#   another, copy this one and change MP to the new group in: the job name, the
#   "MP-*" filter, the MP_SAMPLES/MATRIX_MP names, the output file name, and the
#   plot title. The palette (COLOR_PALETTE) is the only per-group setting.
#
# DIRECTORIES
#   Run from 04_Deeptools. Job logs go to out_log/ and err_log/ in the folder
#   you submit from. Everything else uses full paths.
#   DEEPTOOLS_DIR = <path>/04_Deeptools
#   BW_DIR        = DEEPTOOLS_DIR/00_BigWig_Files_Individual_Bams
#   MATRIX        = DEEPTOOLS_DIR/matrix_combined_gene_only_bed_file.gz
#   OUTDIR        = DEEPTOOLS_DIR/Figures (created automatically)
#
# BEFORE SUBMITTING (run from 04_Deeptools)
#   1. mkdir -p out_log err_log
#        (SLURM opens the log files before the script starts, so these folders
#         must exist ahead of time or the job fails with no log.)
#   2. Confirm Step 4.1.2 finished:
#        ls -lh <DEEPTOOLS_DIR>/matrix_combined_gene_only_bed_file.gz
#   3. sbatch 04.2_MP_Plot_profile_Individual_Bams_Dynamic_sample_naming.sh
#
# NOT AN ARRAY JOB
#   There is no array size to change. One job loops over every target/day. The
#   sample lists are built from the bigWig files on disk, so new MP replicates
#   are picked up automatically.
#
# HOW SAMPLES ARE CHOSEN
#   Sample names look like <GROUP>-<TARGET>-<DAY>-R<n> (e.g. MP-AC-D1-R1). For
#   each target/day, every bigWig named *-<TARGET>-<DAY>-*_inputsubtracted.bw is
#   found, then those starting with "MP-" are kept. A target/day with no MP
#   samples is skipped with a message and the job still finishes.
#
# LIMIT: 8 REPLICATES PER PLOT
#   COLOR_PALETTE has 8 colors, and each MP sample needs one. If a target/day
#   has more than 8 MP samples, plotProfile gets fewer colors than samples and
#   will likely fail. Add colors to COLOR_PALETTE if that happens.
#
# THINGS YOU MAY WANT TO EDIT
#   COLOR_PALETTE   Line colors, used in order, one per replicate.
#   --yMin/--yMax   Shared y-axis range (-0.3 to 0.8). Change if new data falls
#                   outside it. The representative-replicate profile script
#                   (04.2_All_Marks_By_Day_...) uses -0.10 to 0.55 instead.
#   Targets/days    The "for TARGET" and "for DAY" loops.
#
# CHECKING THE RESULTS
#   ls <OUTDIR>/*_MP_profile.png     # 6 plots
#   grep "MP samples:" out_log/*.log # which samples went into each plot
#   grep "skipping" out_log/*.log    # plots that were left out
#
# NOTES
#   - The script stops on any error (set -e), for example if a sample name is
#     not found in the combined matrix.
#   - Paths to check before running: SCRATCH, DEEPTOOLS_DIR, and the email above.
################################################################################

set -e
shopt -s nullglob

SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh
conda activate ${SCRATCH}/conda_envs/deeptools_kernel_v2

# ----------------------------
# Environment setup
# ----------------------------
DEEPTOOLS_DIR="${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/04_Deeptools"
BW_DIR="${DEEPTOOLS_DIR}/00_BigWig_Files_Individual_Bams"
#MATRIX="${DEEPTOOLS_DIR}/bigwiglist_Individual_Bams_gene_only_bed_file.gz"
MATRIX="${DEEPTOOLS_DIR}/matrix_combined_gene_only_bed_file.gz"
OUTDIR="${DEEPTOOLS_DIR}/Figures"
mkdir -p "$OUTDIR"

# ----------------------------
# Input/output paths
# ----------------------------
COLOR_PALETTE=(orange darkorange peru chocolate goldenrod sandybrown tan orangered)

# ============================================================
# Loop through targets and days, building sample lists live
# from whatever bigwigs actually exist on disk
# ============================================================
for TARGET in AC ME SR; do
    for DAY in D1 D3; do
        echo "============================================================"
        echo "TARGET: $TARGET   DAY: $DAY"
        echo "============================================================"

        # All sites, this target + day
        mapfile -t ALL_SAMPLES < <(
            ls "${BW_DIR}"/*-${TARGET}-${DAY}-*_inputsubtracted.bw 2>/dev/null \
                | xargs -n1 basename \
                | sed 's/_inputsubtracted\.bw$//' \
                | sort
        )
        if [[ ${#ALL_SAMPLES[@]} -eq 0 ]]; then
            echo "  No samples found for ${TARGET} ${DAY} — skipping"
            continue
        fi
        echo "  Found ${#ALL_SAMPLES[@]} samples: ${ALL_SAMPLES[*]}"

        # MP-only subset of the above
        MP_SAMPLES=()
        for s in "${ALL_SAMPLES[@]}"; do
            [[ "$s" == MP-* ]] && MP_SAMPLES+=("$s")
        done
        if [[ ${#MP_SAMPLES[@]} -eq 0 ]]; then
            echo "  No MP samples for ${TARGET} ${DAY} — skipping MP plot"
            continue
        fi
        echo "  MP samples: ${MP_SAMPLES[*]}"

        # Step 1 — subset all-samples matrix down to this target+day
        MATRIX_TARGET="${DEEPTOOLS_DIR}/matrix_${TARGET}_${DAY}_ALL.gz"
        computeMatrixOperations subset \
            -m "$MATRIX" \
            --samples "${ALL_SAMPLES[@]}" \
            -o "$MATRIX_TARGET"

        # Step 2 — subset down further to MP only
        MATRIX_MP="${DEEPTOOLS_DIR}/matrix_${TARGET}_${DAY}_MP.gz"
        computeMatrixOperations subset \
            -m "$MATRIX_TARGET" \
            --samples "${MP_SAMPLES[@]}" \
            -o "$MATRIX_MP"

        # Enough colors for however many MP reps exist this round
        n=${#MP_SAMPLES[@]}
        PLOT_COLORS=("${COLOR_PALETTE[@]:0:$n}")

        plotProfile \
            -m "$MATRIX_MP" \
            -out "${OUTDIR}/${TARGET}_${DAY}_MP_profile.png" \
            --perGroup \
            --colors "${PLOT_COLORS[@]}" \
            --samplesLabel "${MP_SAMPLES[@]}" \
            --regionsLabel "" \
            -T "${TARGET} ${DAY} — MP Only" \
            -y "Input-subtracted signal (BPM)" \
            --yMin -0.3 \
            --yMax 0.8 \
            --plotHeight 8 \
            --plotWidth 12 \
            --dpi 300
        echo "  ✔ Created: ${OUTDIR}/${TARGET}_${DAY}_MP_profile.png"
    done
done
echo "All MP-only target/day plots complete."
