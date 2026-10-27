#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=FastP_Cxt_array
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --array=1-96		# MUST equal the number of lines in r1_file_list.txt (see below)
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c832500103@colostate.edu  # CHANGE to your own address


################################################################################
# STEP 1: Trim paired-end reads with fastp (one sample per array task)
#
# INPUT:   <sample>_R1_001.fastq.gz and <sample>_R2_001.fastq.gz in INDIR
# OUTPUT:  in OUTDIR, per sample:
#            <sample>_R1_trimmed.fastq.gz, <sample>_R2_trimmed.fastq.gz
#            <sample>_fastp.html, <sample>_fastp.json  (QC reports)
# ENV:     fastp_env
# RESOURCES per task: 8 CPUs, 8 GB, 2 h (raise --time if large samples time out)
#
# BEFORE SUBMITTING (run from the folder you submit from):
#   1. mkdir -p logs
#   2. Make the R1 file list. Task N processes line N of this file:
#        cd /scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Raw_Files
#        ls *R1_001.fastq.gz | sort > r1_file_list.txt
#        wc -l r1_file_list.txt        # <- this number is your array size
#
# CHANGING THE ARRAY SIZE TO MATCH YOUR SAMPLE COUNT
#   The "#SBATCH --array=1-96" line above must equal the line count of
#   r1_file_list.txt.
#     - Array smaller than the list: extra samples are silently skipped.
#     - Array larger than the list: extra tasks fail with
#       "No R1 file found for array index N".
#   Either edit that line, or override it at submission (no editing needed):
#        N=$(wc -l < r1_file_list.txt)
#        sbatch --array=1-${N} 01_fastp_Cxt.sh
#   Other useful forms:
#        sbatch --array=1-${N}%20 01.1_fastp_Cxt.sh   # at most 20 tasks at once
#        sbatch --array=5,17 01.1_fastp_Cxt.sh        # rerun only tasks 5 and 17
#
# NOTES
#   - If a sample has no matching R2 file, the task prints a message and exits
#     without error. If the output count is low, check the logs for
#     "No R2 found".
#   - Paths to check before running: INDIR, OUTDIR, and the email address above.
################################################################################

SCRATCH=/scratch/alpine/.colostate.edu/c832500103
PIPE=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates

# Directory containing your FASTQ files
INDIR="${SCRATCH}/Cxt_Cut_Run_Raw_Files"
OUTDIR="${PIPE}/01_Cxt_FastP/01.1_Cxt_fastp_out"
R1_LIST="${INDIR}/r1_file_list.txt"

source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${SCRATCH}/miniconda3"; exit 1; }
conda activate ${SCRATCH}/conda_envs/fastp_env || { echo "ERROR: conda activate rnaPseudo failed"; exit 1; }

# Stop on errors (after conda activation, which does not play well with set -e)
set -eo pipefail

mkdir -p $OUTDIR

# Move into the FASTQ directory
cd "$INDIR"

######################## Loop through each file (takes to long)
# Loop through all R1 files
#for R1 in *R1_001.fastq.gz; do

 # Skip md5 files
 #   [[ "$R1" == *.md5 ]] && continue

    # Derive the matching R2 filename
  #  R2=${R1/_R1_/_R2_}

   # if [[ ! -f "$R2" ]]; then
    #    echo "No R2 found for $R1 — skipping."
     #   continue
   # fi
########################

######################## Run as an array
# Get this task's R1 file (line number = array task ID)
R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$R1_LIST")

if [[ -z "$R1" ]]; then
    echo "ERROR: No R1 file found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Matching R2 file: swap _R1_ for _R2_
R2=${R1/_R1_/_R2_}
if [[ ! -f "$INDIR/$R2" ]]; then
    echo "No R2 found for $R1 - skipping."
    exit 1
fi

    # Extract sample name prefix (everything before _R1)
    SAMPLE=${R1%%_R1_001.fastq.gz}
    echo "Task $SLURM_ARRAY_TASK_ID: Processing sample: $SAMPLE"

    fastp \
        -i "$INDIR/$R1" \
        -I "$INDIR/$R2" \
        -o "$OUTDIR/${SAMPLE}_R1_trimmed.fastq.gz" \
        -O "$OUTDIR/${SAMPLE}_R2_trimmed.fastq.gz" \
        --html "$OUTDIR/${SAMPLE}_fastp.html" \
        --json "$OUTDIR/${SAMPLE}_fastp.json" \
        --thread 8

echo "Task $SLURM_ARRAY_TASK_ID: done ($SAMPLE)"

