#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=MultiQC_fastp
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

################################################################################
# STEP 1.2: Summarize the fastp reports for all samples with MultiQC
#
# INPUT:   the *_fastp.json files in 01.1_Cxt_fastp_out (from 01.1_fastp_Cxt.sh)
# OUTPUT:  MultiQC report in 01.2_FastP_MultiQC/MultiQC_out
# ENV:     rnaPseudo_clean (contains multiqc)
#
# Submit from the folder that holds a logs/ directory:
#     mkdir -p logs
#     sbatch 01.2_multiqc_fastp.sh
# (It is light enough to run directly with `bash 01.2_multiqc_fastp.sh` too.)
################################################################################

SCRATCH=/scratch/alpine/.colostate.edu/c832500103
PIPE=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates

FASTP_OUT=${PIPE}/01_Cxt_FastP/01.1_Cxt_fastp_out
MQC_OUT=${PIPE}/01_Cxt_FastP/01.2_FastP_MultiQC/MultiQC_out

source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${SCRATCH}/miniconda3" >&2; exit 1; }
conda activate ${SCRATCH}/conda_envs/rnaPseudo_clean || { echo "ERROR: conda activate rnaPseudo_clean failed" >&2; exit 1; }

set -eo pipefail

[[ -d "$FASTP_OUT" ]] || { echo "ERROR: fastp output folder not found: $FASTP_OUT" >&2; exit 1; }
mkdir -p "$MQC_OUT"

multiqc "$FASTP_OUT" -o "$MQC_OUT"
