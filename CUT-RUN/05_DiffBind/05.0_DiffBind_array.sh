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
# ENV:     diffbind_env (R + DiffBind + edgeR)
# NOTES:   Run "mkdir -p logs" in the folder you submit from.
################################################################################

set -eo pipefail

PROJ=/projects/c832500103@colostate.edu
source ${PROJ}/miniforge3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found" >&2; exit 1; }
conda activate ${PROJ}/conda_envs/diffbind_env || { echo "ERROR: conda activate diffbind_env failed" >&2; exit 1; }

SCRIPT_DIR=/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/05_DiffBind
Rscript "${SCRIPT_DIR}/05.0_DiffBind_Contrasts.R" "${SLURM_ARRAY_TASK_ID}"
