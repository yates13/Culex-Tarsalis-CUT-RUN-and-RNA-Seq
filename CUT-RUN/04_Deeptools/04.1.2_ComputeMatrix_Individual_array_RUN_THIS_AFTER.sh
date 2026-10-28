#!/usr/bin/env bash
#SBATCH --qos=cpu-normal
#SBATCH --job-name=computeMatrix_merge
#SBATCH --output=logs/%x.%j.log
#SBATCH --error=logs/%x.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --partition=acpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu

################################################################################
# STEP 4.1.2: Combine the per-group matrices into one matrix (single job)
#
# RUN THIS AFTER 04.1_ComputeMatrix_Individual_array.sh has finished for every group.
#
# INPUT:   matrix_BF.gz, matrix_MP.gz, matrix_SF.gz, matrix_ZH.gz
#            (made by Step 4.1, in the folder you submit from)
# OUTPUT:  matrix_combined_gene_only_bed_file.gz
#            One matrix with the samples from all four groups side by side.
# ENV:     deeptools_kernel_v2 (deepTools computeMatrixOperations)
# RESOURCES: 4 CPUs, 8 GB, 6 h (this step is light; it does not use the extra CPUs)
#
# DIRECTORIES
#   Run from 04_Deeptools, the same folder used for Step 4.1. Input and output
#   paths are relative to the folder you submit from.
#   logs/  Job logs. Must exist before submitting.
#
# BEFORE SUBMITTING (run from 04_Deeptools)
#   1. mkdir -p logs
#        (SLURM opens the log files before the script starts, so this folder
#         must exist ahead of time or the job fails with no log.)
#   2. Confirm all four input matrices exist:
#        ls matrix_BF.gz matrix_MP.gz matrix_SF.gz matrix_ZH.gz
#   3. sbatch 04.1.2_ComputeMatrix_Individual_array_RUN_THIS_AFTER.sh
#
# NOT AN ARRAY JOB
#   There is no array size to change. To add or remove a group, edit the list of
#   matrix files in the cbind command below.
#
# HOW IT WORKS
#   cbind puts the samples side by side (columns) and requires every input matrix
#   to have the same regions in the same order. The commented-out rbind line
#   below does the opposite: it stacks regions (rows) and requires the same
#   samples in every matrix. Use cbind here, since each group has different
#   samples over the same regions.
#
# CHECKING THE RESULTS
#   ls -lh matrix_combined_gene_only_bed_file.gz
#   computeMatrixOperations info -m matrix_combined_gene_only_bed_file.gz
#     (lists the sample labels and the number of regions)
#
# NOTES
#   - The script stops on any error (set -e).
#   - Paths to check before running: SCRATCH and the email above.
################################################################################

set -e

# ----------------------------
# Environment setup
# ----------------------------
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh
conda activate ${SCRATCH}/conda_envs/deeptools_kernel_v2

# ----------------------------
# Combine the group matrices (samples side by side)
# ----------------------------
computeMatrixOperations cbind -m matrix_BF.gz matrix_MP.gz matrix_SF.gz matrix_ZH.gz -o matrix_combined_gene_only_bed_file.gz
#computeMatrixOperations rbind -m matrix_BF.gz matrix_MP.gz matrix_SF.gz matrix_ZH.gz -o matrix_combined_gene_only_bed_file.gz
