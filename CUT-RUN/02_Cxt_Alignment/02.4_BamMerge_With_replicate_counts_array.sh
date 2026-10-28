#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --job-name=BamMerge_with_replicate_count
#SBATCH --mail-user=shawn.yates@colostate.edu
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x.%A_%a.log
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --array=1-32		 # MUST equal the number of CONDITIONS in conditions.txt (not samples)

################################################################################
# STEP 2.4: Merge replicate BAMs for each condition (one condition per array task)
#
# INPUT:   <condition>-R<n>-Cxt.bam files from Step 2.3 (02.3_Cxt_BamConverted)
#            e.g. <condition>-R1-Cxt.bam, <condition>-R2-Cxt.bam
#          conditions.txt = one condition name per line (built below)
# OUTPUT:  in 02.4_Cxt_BamMerged_with_replicate_count, per condition:
#            <condition>-rep<digits>-merged-sorted.bam (+ .bai index)
#            The rep tag lists which replicates went in: rep123 = R1, R2, R3;
#            rep12 = R1 and R2 only.
#            <condition>-rep<digits>-merged.bam (unsorted intermediate, safe to delete)
# ENV:     rnaPseudo_clean (samtools)
# RESOURCES per task: 8 CPUs, 16 GB, 4 h
# If a condition has only one replicate, its BAM is copied instead of merged.
#
# BEFORE SUBMITTING
#   SLURM needs the array size at submit time, so this script can't set its own
#   range. Two steps:
#
#   1. Build the conditions list (one line per condition, replicate suffix removed):
#        bamdir=<path>/02.3_Cxt_BamConverted
#        mergedir=<path>/02.4_Cxt_BamMerged_with_replicate_count
#        mkdir -p $mergedir
#        ls ${bamdir}/*.bam | xargs -n1 basename | sed 's/-R[0-9]-Cxt\.bam//' \
#            | sort -u > ${mergedir}/conditions.txt
#
#   2. Submit with the array size taken from that file:
#        sbatch --array=1-$(wc -l < ${mergedir}/conditions.txt) \
#            02.4_BamMerge_With_replicate_counts_array.sh
#
# CHANGING THE ARRAY SIZE
#   The array size is the number of CONDITIONS, not the number of samples.
#   With 96 BAMs at 3 replicates each, that is 32. Task N merges line N of
#   conditions.txt (counting from 1).
#     - Array smaller than the list: extra conditions are silently skipped.
#     - Array larger than the list: extra tasks fail with
#       "no condition found for array index N".
#   Using the wc -l form above avoids editing the "#SBATCH --array" line.
#   Other useful forms:
#        sbatch --array=1-32%10 02.4_BamMerge_With_replicate_counts_array.sh   # at most 10 tasks at once
#        sbatch --array=5,17 02.4_BamMerge_With_replicate_counts_array.sh      # rerun only tasks 5 and 17
#
# CHECKING THE RESULTS
#   ls <mergedir>/*-merged-sorted.bam | wc -l    # should equal the number of conditions
#   Each task's log ends with a line like:
#     "Condition X: 3 replicate(s) (rep123) were merged."
#
# NOTES
#   - Requires BAM names to follow <condition>-R<n>-Cxt.bam. Other names break the
#     conditions list and the replicate detection.
#   - Logs are written to the folder you submit from, one .log and .err per task.
#   - Paths to check before running: SCRATCH, bamdir, mergedir, and the email above.
#   - The script stops on any error (set -euo pipefail).
################################################################################

# -----------------------------
# Load environment
# -----------------------------
#module purge
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh
conda activate ${SCRATCH}/conda_envs/rnaPseudo_clean
date
nthreads=$SLURM_CPUS_PER_TASK

set -euo pipefail

# -----------------------------
# Define directories
# -----------------------------
samdir=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/02.2_Cxt_HisatAligned_Cleaned
bamdir=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/02.3_Cxt_BamConverted
mergedir=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/02.4_Cxt_BamMerged_with_replicate_count
mkdir -p $mergedir

CONDITIONS_LIST=${mergedir}/conditions.txt
if [ ! -s "$CONDITIONS_LIST" ]; then
    echo "ERROR: $CONDITIONS_LIST missing or empty. Run the prep step (see header comment) before submitting this array job." >&2
    exit 1
fi

# -----------------------------
# Pick this task's condition from the pre-built list
# -----------------------------
cond=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CONDITIONS_LIST")
if [ -z "$cond" ]; then
    echo "ERROR: no condition found for array index ${SLURM_ARRAY_TASK_ID} in $CONDITIONS_LIST" >&2
    exit 1
fi

echo "----------------------------------------"
echo "Array task ${SLURM_ARRAY_TASK_ID}: merging condition: $cond"

# Find BAMs for this condition
bams=$(ls ${bamdir}/${cond}-R*-Cxt.bam)

# Extract replicate numbers (R1, R2, R3 -> 1 2 3)
reps=$(echo "$bams" | sed -E 's/.*-R([0-9]+)-Cxt\.bam/\1/' | tr '\n' ' ')
repstring="rep$(echo $reps | tr -d ' ')"

# Count replicates
count=$(echo "$reps" | wc -w)
echo "Found $count replicate(s): $repstring"

# Output filename includes repstring
merged_base="${mergedir}/${cond}-${repstring}-merged.bam"

if [ "$count" -gt 1 ]; then
    echo "Merging $count BAM files..."
    samtools merge -f "$merged_base" $bams
    merge_status="merged"
else
    echo "Only one BAM found - copying instead of merging."
    cp $bams "$merged_base"
    merge_status="copied"
fi

# Sort + index
sorted="${mergedir}/${cond}-${repstring}-merged-sorted.bam"
samtools sort -@$nthreads -o "$sorted" "$merged_base"
samtools index "$sorted"

echo "Completed: $sorted"
echo "Condition $cond: $count replicate(s) ($repstring) were $merge_status."
date
