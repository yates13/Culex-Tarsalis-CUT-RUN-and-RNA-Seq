#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=Gviz_07.1
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=02:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shawn.yates@colostate.edu

################################################################################
# STEP 7.1: Gviz genome-track figures for one DiffBind contrast (one job)
#
# RUNS:    07.1_Gviz_Culex_Alpine.R (same folder as this script)
# NEEDS:   Step 7.0 (CtarK1_Gviz_Meta_Alpine.csv), Step 5.0 (DiffBind CSVs),
#          Step 4 bigWigs, Step 3.1 pooled peaks
# OUTPUT:  Round1_Round2_With_Duplicates/07_Gviz/<contrast>_FDR_<fdr>_<direction>_<timestamp>/  (NEW)
# USAGE:   mkdir -p logs   (in the folder you submit from, e.g. $RUN/07_Gviz)
#          sbatch 07.1_Gviz_Culex.sh                        # USER SETTINGS in the R script
#          sbatch 07.1_Gviz_Culex.sh ZH BF SR D1 0.1 up 5    # targetPop refPop mark day FDR direction topN
# ENV:     gviz_env
################################################################################

set -eo pipefail

PROJ=/projects/c832500103@colostate.edu
source ${PROJ}/miniforge3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found" >&2; exit 1; }
conda activate ${PROJ}/conda_envs/gviz_env || { echo "ERROR: conda activate gviz_env failed" >&2; exit 1; }

SCRIPT_DIR=/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/07_Gviz
Rscript "${SCRIPT_DIR}/07.1_Gviz_Culex_Alpine.R" "$@"
