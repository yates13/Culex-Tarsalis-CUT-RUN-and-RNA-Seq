#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=DiffBind_contrasts
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=1-24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shawn.yates@colostate.edu

################################################################################
# STEP 5.0: DiffBind, one array task per contrast (24 = 2 days x 3 marks x 4 pairs)
#
# RUNS:    05.0_DiffBind_Contrasts.R <task number>   (same folder as this script)
#          "Rscript 05.0_DiffBind_Contrasts.R list" prints which number is which.
# NEEDS:   Step 3.1 pooled peaks (done 11/20), filtered BAMs (2.3.1), pooled inputs (2.4)
# OUTPUT:  Round1_Round2_With_Duplicates/05_DiffBind_2026-11-24/{d1,d3}/  (NEW)
#          (the output folder is OUT_BASE in 05.0_DiffBind_Contrasts.R, not set here)
# ENV:     diffbind_env (R + DiffBind + edgeR)
# RESOURCES per task: 1 CPU, 32 GB, 12 h
#
# BEFORE SUBMITTING (from the folder you submit from):
#   1. mkdir -p logs
#   2. Rscript 05.0_DiffBind_Contrasts.R list        # check the task numbers
#   3. sbatch 05.0_DiffBind_array.sh                 # all tasks
#      sbatch --array=4,8,12,16,20,24 05.0_DiffBind_array.sh   # only some tasks
#   A task whose output CSV already exists stops without overwriting anything,
#   so resubmitting the whole array is safe.
#
# USER SETTINGS
#   #SBATCH lines (not tagged, so Slurm reads them exactly; edit the values above):
#     --array=1-24        must equal the number of rows printed by "list"
#                         (it changes when DAYS, MARKS or PAIRS change in the R script)
#     --time=12:00:00     raise if a task times out
#     --mem=32G           raise if a task is killed for memory
#     --mail-user         your email
#   Code lines tagged [EDIT] (list them: grep -n "\[EDIT\]" 05.0_DiffBind_array.sh):
#     PROJ                your /projects folder (miniforge3 + conda_envs)
#     SCRIPT_DIR          folder holding 05.0_DiffBind_Contrasts.R (the repo clone)
#
# Analysis settings (contrasts, peaks, normalization, paths to data) are NOT
#   here: they are in the SETTINGS and PATHS blocks of 05.0_DiffBind_Contrasts.R.
################################################################################

set -eo pipefail

PROJ=/projects/c832500103@colostate.edu   # [EDIT]
source ${PROJ}/miniforge3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found" >&2; exit 1; }
conda activate ${PROJ}/conda_envs/diffbind_env || { echo "ERROR: conda activate diffbind_env failed" >&2; exit 1; }

SCRIPT_DIR=/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/05_DiffBind   # [EDIT]
Rscript "${SCRIPT_DIR}/05.0_DiffBind_Contrasts.R" "${SLURM_ARRAY_TASK_ID}"
