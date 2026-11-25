#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=getGeneDistances
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=06:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shawn.yates@colostate.edu

################################################################################
# STEP 6.1: peak-to-gene distances for the 18 DiffBind contrasts (one job)
#
# RUNS:    06.1_getGeneDistances_Culex_Alpine.R (same folder as this script)
# NEEDS:   Step 6.0 (CtarK1_MAKER_genes.gtf) and Step 5.0 (DiffBind CSVs)
# OUTPUT:  Round1_Round2_With_Duplicates/06_Gene_Distances_2026-11-24/  (NEW)
#          contrasts whose output already exists are skipped, not overwritten
# ENV:     diffbind_env (R; readr is used if installed, base R otherwise)
# NOTES:   Run "mkdir -p logs" in the folder you submit from.
################################################################################

set -eo pipefail

PROJ=/projects/c832500103@colostate.edu
source ${PROJ}/miniforge3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found" >&2; exit 1; }
conda activate ${PROJ}/conda_envs/diffbind_env || { echo "ERROR: conda activate diffbind_env failed" >&2; exit 1; }

SCRIPT_DIR=/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/06_Gene_Distances
Rscript "${SCRIPT_DIR}/06.1_getGeneDistances_Culex_Alpine.R"
