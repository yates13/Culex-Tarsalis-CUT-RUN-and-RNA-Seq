#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=Ctark1Build
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=4:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c832500103@colostate.edu  # CHANGE to your own address

################################################################################
# STEP 2.1: Build the HISAT2 index for the C. tarsalis CtarK1 genome
#
# INPUT:   Cxt_Genome/CtarK1.fa   (relative to the folder you submit from)
# OUTPUT:  CtarK1_index.*.ht2     (8 files, written one folder up from Cxt_Genome)
# ENV:     rnaPseudo (must contain hisat2)
#
# Submit from the folder that CONTAINS Cxt_Genome/:
#     mkdir -p logs
#     sbatch 02.1_runBuildCtark1Genome_apline.sh
################################################################################

SCRATCH=/scratch/alpine/.colostate.edu/c832500103
module purge

# Load conda
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${SCRATCH}/miniconda3" >&2; exit 1; }
conda activate rnaPseudo || { echo "ERROR: conda activate rnaPseudo failed" >&2; exit 1; }

set -eo pipefail

# Move into genome directory
cd Cxt_Genome

# Build HISAT2 index
hisat2-build -p 8 CtarK1.fa ../CtarK1_index
