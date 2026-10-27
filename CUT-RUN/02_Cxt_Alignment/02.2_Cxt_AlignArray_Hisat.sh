#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --job-name=HisatArrayCxt
#SBATCH --mail-user=shawn.yates@colostate.edu  # Change to your email
#SBATCH --mail-type=all
#SBATCH --output=logs/%x.%A-%a.log
#SBATCH --error=logs/%x.%A-%a.err
#SBATCH --array=0-95%20

################################################################################
# STEP 2.2: Align trimmed paired-end reads with HISAT2 (one sample per array task)
#
# INPUT:   trimmed FASTQ pairs from Step 1 (01.1_Cxt_fastp_out):
#            <sample>_R1_trimmed.fastq.gz and <sample>_R2_trimmed.fastq.gz
#          FastqTrimmed.txt = list of the R1 trimmed file names (one per line)
# OUTPUT:  in 02.2_Cxt_HisatAligned_Cleaned, per sample:
#            <sample>.sam            (alignments)
#            <sample>.hisat2.log     (HISAT2 alignment summary)
# INDEX:   00_Cxt_Indexed_Genome/CtarK1_index (built by 02.1_runBuildCtark1Genome_apline.sh)
# ENV:     rnaPseudo_clean
# RESOURCES per task: 4 CPUs, 16 GB, 4 h
# HISAT2 flags: --no-spliced-alignment (CUT&RUN is genomic, not RNA),
#               --no-mixed (only pairs where both mates align), --mm (memory-mapped index)
#
# BEFORE SUBMITTING (run from 02_Cxt_Alignment):
#   1. mkdir -p logs
#        (SLURM opens the log files before the script starts, so this folder
#         must exist ahead of time or the job fails with no log.)
#   2. Build the file list. Task N processes line N+1 of this file:
#        ls <fastp_out_dir> | grep '_R1_' | grep 'trimmed' > FastqTrimmed.txt
#        wc -l FastqTrimmed.txt        # <- this number is N
#
# CHANGING THE ARRAY SIZE TO MATCH YOUR SAMPLE COUNT
#   This script counts tasks from 0, not 1. If FastqTrimmed.txt has N lines,
#   the array must be 0-(N-1). For 96 samples that is 0-95.
#     - Array smaller than the list: extra samples are silently skipped.
#     - Array larger than the list: extra tasks find no line and quietly
#       do nothing, so a successful job does NOT prove every sample ran.
#   Either edit the "#SBATCH --array" line, or override it at submission:
#        N=$(wc -l < FastqTrimmed.txt)
#        sbatch --array=0-$((N-1))%20 02.2_Cxt_AlignArray_Hisat.sh
#   Other useful forms:
#        sbatch --array=0-$((N-1)) script.sh    # no limit on simultaneous tasks
#        sbatch --array=5,17 script.sh          # rerun only tasks 5 and 17
#
# CHECKING THE RESULTS
#   ls <OUT_DIR>/*.sam | wc -l      # should equal N
#   Each .hisat2.log ends with the overall alignment rate.
#
# NOTES
#   - Paths to check before running: SCRATCH, FASTQ_DIR, OUT_DIR, INDEX, filename.
#   - If a sample's R1 or R2 file is missing, that task exits with an error.
################################################################################

echo "[$0] Job: $SLURM_JOB_NAME | Task ID: $SLURM_ARRAY_TASK_ID"

# ----------------------------
# Environment setup
# ----------------------------
module purge

SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${SCRATCH}/miniconda3"; exit 1; }
conda activate ${SCRATCH}/conda_envs/rnaPseudo_clean || { echo "ERROR: conda activate rnaPseudo_clean failed"; exit 1; }

# Stop on errors (after conda activation, which does not play well with set -e)
set -eo pipefail

date

# ----------------------------
# Input/output paths
# ----------------------------
FASTQ_DIR=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/01_Cxt_FastP/01.1_Cxt_fastp_out
OUT_DIR=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/02.2_Cxt_HisatAligned_Cleaned
INDEX=${SCRATCH}/Cxt_Cut_Run_Pipeline/00_Cxt_Indexed_Genome/CtarK1_index
filename=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/FastqTrimmed.txt
mkdir -p $OUT_DIR
mkdir -p ${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/slurm_logs

# ----------------------------
# Check the file list exists
# ----------------------------
# Build it once with:
#   ls ${FASTQ_DIR} | grep "_R1_" | grep "trimmed" > ${filename}
if [[ ! -f "$filename" ]]; then
    echo "ERROR: $filename not found — build it first with:"
    echo "  ls ${FASTQ_DIR} | grep '_R1_' | grep 'trimmed' > ${filename}"
    exit 1
fi


# Debug (VERY helpful)
echo "Total files found: $(wc -l < ${filename})"

# =========================
# Select file for this task
# =========================
linenum=0
while read -r line
do
  if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
  then
    R1=${FASTQ_DIR}/${line}
    R2=${R1/_R1_/_R2_}
    BASE=$(basename $R1)
    SAMPLE=${BASE%%_R1_trimmed*}

    # Safety check
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
      echo "ERROR: Missing FASTQ files"
      exit 1
    fi

    # =========================
    # Run HISAT2
    # =========================
    hisat2 \
      --phred33 \
      --mm \
      --no-mixed \
      --no-spliced-alignment \
      -p $SLURM_CPUS_PER_TASK \
      -x $INDEX \
      -1 $R1 \
      -2 $R2 \
      -S $OUT_DIR/${SAMPLE}.sam \
      2> $OUT_DIR/${SAMPLE}.hisat2.log
    echo "Completed: $SAMPLE"
  fi

  linenum=$((linenum + 1))
done < ${filename}

date
echo "Job complete"
