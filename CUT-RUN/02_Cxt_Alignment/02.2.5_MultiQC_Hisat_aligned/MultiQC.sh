#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --job-name=MultiQC
#SBATCH --mail-user=shawn.yates@colostate.edu
#SBATCH --mail-type=all
#SBATCH --output=logs/%x.%A.log
#SBATCH --error=logs/%x.%A.err

################################################################################
# STEP 2.2.5: Combined QC report for trimming and alignment (single job)
#
# RUN THIS AFTER Step 1 (fastp) and Step 2.2 (HISAT2) have finished for every
# sample.
#
# INPUT:   01_Cxt_FastP/01.1_Cxt_fastp_out            (fastp .json reports)
#          02_Cxt_Alignment/02.2_Cxt_HisatAligned_Cleaned   (.hisat2.log files)
#          MultiQC scans both folders and picks up the files it recognizes.
# OUTPUT:  in 02_Cxt_Alignment/02.2.5_MultiQC_Hisat_aligned:
#            Cxt-CUTRUN-Pipeline-QC-fastp--HISAT2_multiqc_report.html
#            Cxt-CUTRUN-Pipeline-QC-fastp--HISAT2_multiqc_report_data/
#          Open the .html file in a browser (download it from Alpine first).
# ENV:     rnaPseudo_clean (MultiQC)
# RESOURCES: 1 CPU, 16 GB, 1 h
#
# DIRECTORIES
#   Run from 02_Cxt_Alignment. Job logs go to logs/ in the folder you
#   submit from. The MultiQC input and output folders are full paths in the
#   multiqc command below.
#
# BEFORE SUBMITTING (run from 02_Cxt_Alignment)
#   1. mkdir -p logs
#        (SLURM opens the log files before the script starts, so this folder
#         must exist ahead of time or the job fails with no log.)
#   2. Confirm both input folders are complete:
#        ls <fastp_out_dir>/*.json | wc -l           # one per sample
#        ls <HisatAligned_dir>/*.hisat2.log | wc -l   # one per sample
#   3. sbatch MultiQC.sh
#
# NOT AN ARRAY JOB
#   There is no array size to change. New samples are picked up automatically,
#   since MultiQC reads whatever is in the two input folders.
#
# HOW IT WORKS
#   - -f overwrites any existing report in the output folder, so rerunning
#     replaces the old report.
#   - --ignore "*.err" and --ignore "logs/" keep job logs out of the scan.
#   - --title sets the report title and the report file name.
#   - An optional MultiQC config file (multiqc_config.yaml) can be added with
#     -c to control sample naming. It is commented out at the end of the command.
#
# CHECKING THE RESULTS
#   ls <output_dir>/*_multiqc_report.html
#   Open the report and check the General Statistics table: the number of rows
#   should match your sample count, and the alignment rates should look sensible.
#   A sample missing from the table means its fastp or HISAT2 file is missing.
#
# NOTES
#   - Paths to check before running: the input and output paths in the multiqc
#     command, and the email above.
#   - The SCRATCH and PIPELINE variables below are defined but not used by the
#     multiqc command, which uses full paths.
################################################################################

# ----------------------------
# Environment setup
# ----------------------------
source /scratch/alpine/.colostate.edu/c832500103/miniconda3/etc/profile.d/conda.sh
conda activate /scratch/alpine/.colostate.edu/c832500103/conda_envs/rnaPseudo_clean

# ----------------------------
# Paths (not used by the command below)
# ----------------------------
#SCRATCH=/scratch/alpine/.colostate.edu/c832500103
#PIPELINE=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment

# ----------------------------
# Run MultiQC on the fastp and HISAT2 outputs
# ----------------------------
multiqc \
    /scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/01_Cxt_FastP/01.1_Cxt_fastp_out \
    /scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/02.2_Cxt_HisatAligned_Cleaned \
    --ignore "*.err" \
    --ignore "logs/" \
    -o /scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/02.2.5_MultiQC_Hisat_aligned \
    -f \
    --title "Cxt CUT&RUN Pipeline QC - fastp + HISAT2"
#-c multiqc_config.yaml
