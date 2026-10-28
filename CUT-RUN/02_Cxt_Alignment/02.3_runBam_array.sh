#!/usr/bin/env bashy
#SBATCH --partition=acpu
#SBATCH --job-name=BamConvert_array
#SBATCH --output=logs/%x.%A_%a.log
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --time=3:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shawn.yates@colostate.edu
#SBATCH --array=1-96             # MUST equal the number of lines in sam_file_list.txt (see below)

################################################################################
# STEP 2.3: Convert SAM to sorted, indexed BAM (one sample per array task)
#
# INPUT:   <sample>.sam files from Step 2.2 (02.2_Cxt_HisatAligned_Cleaned)
#          sam_file_list.txt = full paths of the SAM files, one per line
# OUTPUT:  in 02.3_Cxt_BamConverted, per sample:
#            <sample>.bam            (sorted)
#            <sample>.bam.bai        (index)
#            <sample>_flagstat.txt   (alignment stats from samtools flagstat)
# ENV:     rnaPseudo_clean (samtools)
# RESOURCES per task: 8 CPUs, 32 GB, 3 h
#
# BEFORE SUBMITTING (run from 02_Cxt_Alignment):
#   1. mkdir -p out_log err_log
#        (SLURM opens the log files before the script starts, so these folders
#         must exist ahead of time or the job fails with no log.)
#   2. Build the SAM list. Task N processes line N of this file:
#        ls <samdir>/*.sam > <samdir>/sam_file_list.txt
#        wc -l <samdir>/sam_file_list.txt      # <- this number is your array size
#
# CHANGING THE ARRAY SIZE TO MATCH YOUR SAMPLE COUNT
#   The "#SBATCH --array=1-96" line must equal the line count of
#   sam_file_list.txt (task N processes line N, counting from 1).
#     - Array smaller than the list: extra samples are silently skipped.
#     - Array larger than the list: extra tasks fail with
#       "No SAM file found for array index N".
#   Either edit that line, or override it at submission (no editing needed):
#        N=$(wc -l < <samdir>/sam_file_list.txt)
#        sbatch --array=1-${N} 02.3_runBamm_array.sh
#   Other useful forms:
#        sbatch --array=1-${N}%20 02.3_runBamm_array.sh     # at most 20 tasks at once
#        sbatch --array=5,17 02.3_runBamm_array.sh          # rerun only tasks 5 and 17
#
# RERUNNING
#   Safe to resubmit. A sample is skipped if its .bam.bai already exists and is
#   newer than its SAM file. To force a redo, delete that sample's .bam.bai.
#
# CHECKING THE RESULTS
#   ls <outdir>/*.bam | wc -l       # should equal N
#   cat <outdir>/<sample>_flagstat.txt   # mapped read counts and pairing stats
#
# NOTES
#   - Paths to check before running: SCRATCH, samdir, outdir, and the email above.
#   - SAM files are large. Once every BAM is confirmed good, the SAMs can be
#     deleted to free up scratch space.
################################################################################

# ----------------------------
# Environment setup
# ----------------------------
module purge
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${SCRATCH}/miniconda3" >&2; exit 1; }
conda activate ${SCRATCH}/conda_envs/rnaPseudo_clean || { echo "ERROR: conda activate rnaPseudo_clean failed" >&2; exit 1; }
date
nthreads=$SLURM_CPUS_PER_TASK

# Stop on errors, including a failure in the middle of a pipe
# (after conda activation, which does not play well with set -e)
set -eo pipefail


# ----------------------------
# Input/output paths
# ----------------------------
samdir=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/02.2_Cxt_HisatAligned_Cleaned
outdir=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/02.3_Cxt_BamConverted
SAM_LIST=${samdir}/sam_file_list.txt

mkdir -p "$outdir"

if [[ ! -f "$SAM_LIST" ]]; then
    echo "ERROR: $SAM_LIST not found — build it first with:"
    echo "  ls ${samdir}/*.sam > ${SAM_LIST}"
    exit 1
fi

# ----------------------------
# Select this task's SAM file
# ----------------------------
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAM_LIST")

if [[ -z "$file" ]]; then
    echo "ERROR: No SAM file found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

base=$(basename "$file" .sam)
echo "Task $SLURM_ARRAY_TASK_ID: Processing $base"

# ----------------------------
# Skip if already done
# ----------------------------
if [[ -f "${outdir}/${base}.bam.bai" && "${outdir}/${base}.bam.bai" -nt "$file" ]]; then
    echo "EXISTS: ${base} — skipping"
    exit 0
fi

# ----------------------------
# Bam creation, sorting, and indexing
# ----------------------------
samtools view -bS -@$nthreads "$file" | \
samtools sort -@$nthreads -o ${outdir}/${base}.bam
samtools flagstat ${outdir}/${base}.bam > ${outdir}/${base}_flagstat.txt
samtools index -@$nthreads ${outdir}/${base}.bam

echo "Completed: $base"
date
