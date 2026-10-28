#!/usr/bin/env bash
#SBATCH --qos=cpu-normal
#SBATCH --job-name=computeMatrix_array
#SBATCH --output=logs/%x.%A_%a.log
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --partition=acpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu
#SBATCH --array=0-3		# 0-(number of groups in POPS - 1); 4 groups = 0-3 (see below)


################################################################################
# STEP 4.1: Build a TSS signal matrix for each group (one group per array task)
#
# INPUT:   bw_chunks/chunk_<GROUP>.txt = the bigWig files for one group (BF, MP,
#            SF, or ZH), from the input-subtracted bigWigs made in Step 4.0
#            (00_BigWig_Files_Individual_Bams). Each entry is a full path ending
#            in <sample>_inputsubtracted.bw, separated by whitespace or newlines.
#          CtarK1_TSS.bed = transcription start sites (reference regions)
# OUTPUT:  in the folder you submit from, per group:
#            matrix_<GROUP>.gz     (deepTools matrix, TSS -2 kb to +2 kb)
#            region_<GROUP>.bed    (regions kept, in matrix order)
#          Sample labels in the matrix = bigWig names without _inputsubtracted.bw.
# ENV:     deeptools_kernel_v2 (deepTools computeMatrix)
# RESOURCES per task: 4 CPUs, 16 GB, 6 h
#
# DIRECTORIES
#   Run from 04_Deeptools. The chunk files, matrices and region files all use
#   paths relative to the folder you submit from.
#   bw_chunks/       Input. Must contain one chunk_<GROUP>.txt per group.
#   logs/            Job logs. Must exist before submitting.
#   REFERENCE_BED    <path>/00_Cxt_Genome/Genome_gff/CtarK1_TSS.bed (outside the repo)
#
# BEFORE SUBMITTING (run from 04_Deeptools)
#   1. mkdir -p logs
#        (SLURM opens the log files before the script starts, so this folder
#         must exist ahead of time or the job fails with no log.)
#   2. Make sure bw_chunks/chunk_<GROUP>.txt exists for every group in POPS.
#        ls bw_chunks/
#        wc -w bw_chunks/chunk_*.txt       # bigWigs per group
#   3. sbatch 04.1_ComputeMatrix_Individual_array.sh
#
# CHANGING THE ARRAY SIZE
#   Here the array counts GROUPS, not samples. Task 0 = first entry of POPS,
#   task 1 = second, and so on:
#        POPS=(BF MP SF ZH)      # 0=BF  1=MP  2=SF  3=ZH
#   To add or remove a group, edit BOTH the POPS list and the array range, so
#   the array is 0-(number of groups - 1). Examples:
#        5 groups: --array=0-4
#        3 groups: --array=0-2
#     - Array smaller than POPS: the extra groups are silently skipped.
#     - Array larger than POPS: extra tasks get an empty group name, find no
#       chunk file, and fail.
#   Adding samples to an existing group does NOT change the array. Add them to
#   that group's chunk file instead.
#   Other useful forms:
#        sbatch --array=1 script.sh       # rerun only MP
#        sbatch --array=0,3 script.sh     # rerun BF and ZH
#
# MATRIX SETTINGS
#   reference-point at TSS, 2000 bp before and after (-b 2000 -a 2000),
#   --skipZeros (drops regions with no signal in any sample), -p 4.
#
# CHECKING THE RESULTS
#   ls matrix_*.gz region_*.bed        # one of each per group
#   The log prints "Population: <GROUP> | Files: <n>" -- check n matches the
#   number of bigWigs you expect in that group.
#
# NOTES
#   - The script stops on any error (set -e), including a missing chunk file.
#   - Paths to check before running: SCRATCH, REFERENCE_BED, and the email above.
#   - The follow-up script 04.1.2_ComputeMatrix_Individual_array_RUN_THIS_AFTER.sh
#     is named to run after this one.
################################################################################

set -e

# ----------------------------
# Environment setup
# ----------------------------
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh
conda activate ${SCRATCH}/conda_envs/deeptools_kernel_v2

# ----------------------------
# Select this task's group
# ----------------------------
# Task 0 = BF, 1 = MP, 2 = SF, 3 = ZH. Edit this list AND the --array range
# together (see CHANGING THE ARRAY SIZE above).
POPS=(BF MP SF ZH)
POP=${POPS[$SLURM_ARRAY_TASK_ID]}
CHUNK="bw_chunks/chunk_${POP}.txt"
BIGWIGS=$(cat "$CHUNK")

# ----------------------------
# Reference regions
# ----------------------------
REFERENCE_BED="${SCRATCH}/Cxt_Cut_Run_Pipeline/00_Cxt_Genome/Genome_gff/CtarK1_TSS.bed"

# ----------------------------
# Sample labels (bigWig name without _inputsubtracted.bw)
# ----------------------------
LABELS=()
for bw in $BIGWIGS; do
    LABELS+=("$(basename "$bw" _inputsubtracted.bw)")
done

echo "Population: $POP  |  Files: $(echo $BIGWIGS | wc -w)"

# ----------------------------
# Compute the matrix
# ----------------------------
computeMatrix reference-point \
    --referencePoint TSS \
    -b 2000 -a 2000 \
    -R "$REFERENCE_BED" \
    -S $BIGWIGS \
    --skipZeros \
    -o "matrix_${POP}.gz" \
    -p 4 \
    --outFileSortedRegions "region_${POP}.bed" \
    --samplesLabel "${LABELS[@]}"

echo "Done: $POP"
