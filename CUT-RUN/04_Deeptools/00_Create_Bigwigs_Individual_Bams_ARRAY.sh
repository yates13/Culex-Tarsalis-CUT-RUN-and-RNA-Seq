#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=BamToBigwig_Individual_bams
#SBATCH --output=logs/%x.%A_%a.log
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --time=5:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu
#SBATCH --array=1-72		# MUST equal the number of lines in treatment_bam_list.txt (see below)

################################################################################
# STEP 4.0: Make input-subtracted bigWigs from individual BAMs (one BAM per array task)
#
# INPUT:   <TREATMENT>-<TARGET>-<DAY>-R<n>-Cxt.bam files from Step 2.3
#            (02.3_Cxt_BamConverted), e.g. ZH-AC-D1-R1-Cxt.bam
#          treatment_bam_list.txt = full paths of the treatment BAMs, one per line
#            (input controls, names containing "-IN-", are NOT in this list)
# OUTPUT:  in 00_BigWig_Files_Individual_Bams, per sample:
#            <sample>_inputsubtracted.bw   (treatment minus input, BPM-normalized)
#          These feed the matrix and heatmap steps (04.1 onward).
# ENV:     deeptools_kernel_v2 (deepTools bamCompare)
# RESOURCES per task: 8 CPUs, 32 GB, 5 h
#
# DIRECTORIES
#   BAM_INPUT_DIR = <path>/02_Cxt_Alignment/02.3_Cxt_BamConverted
#                   Input. Must contain the sorted, indexed BAMs (.bam + .bai),
#                   including the input controls (-IN-).
#   BW_OUTPUT_DIR = <path>/04_Deeptools/00_BigWig_Files_Individual_Bams
#                   Output. Created automatically.
#   TREATMENT_LIST = <path>/04_Deeptools/treatment_bam_list.txt
#                   Must be built before submitting (see below).
#   Job logs go to logs/ in the folder you submit from.
#
# BEFORE SUBMITTING (run from 04_Deeptools)
#   1. mkdir -p logs
#        (SLURM opens the log files before the script starts, so this folder
#         must exist ahead of time or the job fails with no log.)
#   2. Build the treatment list (all non-input BAMs, full paths). See also
#      How_to_make_Treatment_list.md:
#        ls <BAM_INPUT_DIR>/*-Cxt.bam | grep -v -- '-IN-' > treatment_bam_list.txt
#        wc -l treatment_bam_list.txt      # <- this number is your array size
#   3. Submit with the array size taken from that file:
#        sbatch --array=1-$(wc -l < treatment_bam_list.txt) 00_Create_Bigwigs_Individual_Bams_ARRAY.sh
#
# CHANGING THE ARRAY SIZE TO MATCH YOUR SAMPLE COUNT
#   The "#SBATCH --array=1-72" line must equal the line count of
#   treatment_bam_list.txt (task N processes line N, counting from 1). 72 =
#   the 96 BAMs minus the input controls.
#     - Array smaller than the list: extra samples are silently skipped.
#     - Array larger than the list: extra tasks fail with
#       "No treatment BAM found for array index N".
#   Either edit that line, or use the wc -l form above so you never have to.
#   Put --array BEFORE the script name. Anything after the script name is passed
#   to the script, not to sbatch, and is ignored.
#   Other useful forms:
#        sbatch --array=1-72%10 00_Create_Bigwigs_Individual_Bams_ARRAY.sh   # at most 10 tasks at once
#        sbatch --array=5,17 00_Create_Bigwigs_Individual_Bams_ARRAY.sh      # rerun only tasks 5 and 17
#
# HOW THE INPUT CONTROL IS CHOSEN (TEMPORARY)
#   Each treatment BAM is compared against ANY Day 1 input BAM for the same
#   TREATMENT (first match of <TREATMENT>-IN-D1-*-Cxt.bam), NOT the input with
#   the same replicate number, because not every treatment replicate has a
#   same-numbered input yet. Once all input replicates exist, switch to
#   rep-matched pairing:
#        rep=$(echo "$basename" | cut -d'-' -f4)
#        input_bam="${BAM_INPUT_DIR}/${treatment}-IN-D1-${rep}-Cxt.bam"
#   All days (D1 and D3) currently use a Day 1 input.
#
# RERUNNING
#   Safe to resubmit. A sample is skipped if its .bw already exists. To redo a
#   sample, delete its .bw first. Samples with no input control are skipped
#   and the task still reports success.
#
# CHECKING THE RESULTS
#   ls <BW_OUTPUT_DIR>/*.bw | wc -l    # should equal the array size
#   grep -l "SKIP" logs/*.log          # tasks skipped for no input control
#
# BIGWIG SETTINGS
#   --normalizeUsing BPM, --scaleFactorsMethod None, --operation subtract,
#   --extendReads 150 (fragment length for paired-end reads), -p 8.
#
# NOTES
#   - The script stops on any error (set -e).
#   - Paths to check before running: SCRATCH, BAM_INPUT_DIR, BW_OUTPUT_DIR,
#     TREATMENT_LIST, and the email above.
################################################################################

set -e

# ----------------------------
# Environment setup
# ----------------------------
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh
conda activate ${SCRATCH}/conda_envs/deeptools_kernel_v2

# ----------------------------
# Input/output paths
# ----------------------------
BAM_INPUT_DIR="${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/02.3_Cxt_BamConverted"
BW_OUTPUT_DIR="${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/04_Deeptools/00_BigWig_Files_Individual_Bams"
TREATMENT_LIST="${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/04_Deeptools/treatment_bam_list.txt"

mkdir -p "$BW_OUTPUT_DIR"
mkdir -p logs

# ----------------------------
# Select this task's treatment BAM (task ID 1 = first line)
# ----------------------------
treatment_bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TREATMENT_LIST")

if [[ -z "$treatment_bam" ]]; then
    echo "ERROR: No treatment BAM found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# --- OLD SERIAL VERSION (built + looped through ALL treatments in one job) ---
# declare -a TREATMENTS
# for bam_file in "${BAM_INPUT_DIR}"/*-Cxt.bam; do
#     basename=$(basename "$bam_file")
#     if [[ "$basename" == *"-IN-"* ]]; then
#         continue
#     fi
#     TREATMENTS+=("$bam_file")
# done
# total=${#TREATMENTS[@]}
# echo "Found $total treatment samples"
# echo ""
# count=0
# for treatment_bam in "${TREATMENTS[@]}"; do
#     count=$((count + 1))
# -------------------------------------------------------------------------

basename=$(basename "$treatment_bam" -Cxt.bam)
treatment=$(echo "$basename" | cut -d'-' -f1)

#-------------------------------------------------------------------------------------------------------
#These are for merged bams
#
#    input_bam="${BAM_INPUT_DIR}/${treatment}-IN-D1-rep"*"-merged-sorted.bam"
#    input_bam=$(ls "${BAM_INPUT_DIR}/${treatment}-IN-D1-rep"*"-merged-sorted.bam" 2>/dev/null | head -1)
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# TEMPORARY: pairing to ANY available input replicate for this treatment site,
# not matched by replicate number - not every treatment rep has a same-numbered
# input yet. Revisit once all input replicates exist and switch to exact
# rep-matched pairing (see commented-out version below).
input_bam=$(ls "${BAM_INPUT_DIR}/${treatment}-IN-D1-"*"-Cxt.bam" 2>/dev/null | head -1)
# --- exact rep-matched version to swap back in later ---
# rep=$(echo "$basename" | cut -d'-' -f4)
# input_bam="${BAM_INPUT_DIR}/${treatment}-IN-D1-${rep}-Cxt.bam"
#---------------------------------------------------------------------------------------------------------

echo "Task $SLURM_ARRAY_TASK_ID: $basename"
echo "  treatment: $treatment_bam"
echo "  input:     $input_bam"

if [[ -z "$input_bam" || ! -f "$input_bam" ]]; then
    echo "SKIP: $basename (no input control)"
    exit 0
fi

# ----------------------------
# Skip if already done
# ----------------------------
output_bw="${BW_OUTPUT_DIR}/${basename}_inputsubtracted.bw"
if [[ -f "$output_bw" ]]; then
    echo "EXISTS: $basename"
    exit 0
fi

# ----------------------------
# Create the input-subtracted bigWig
# ----------------------------
echo "Creating: $basename (treatment=$treatment_bam, input=$input_bam)"
bamCompare \
    -b1 "$treatment_bam" \
    -b2 "$input_bam" \
    -o "$output_bw" \
    --normalizeUsing BPM \
    --scaleFactorsMethod None \
    --operation subtract \
    -p 8 \
    --extendReads 150
echo "Done: $basename"

# --- OLD SERIAL VERSION continued (closing loop + final summary) ---
# done
# echo ""
# echo "BigWigs created in: $BW_OUTPUT_DIR"
# ls -lh "$BW_OUTPUT_DIR"/*.bw
# -------------------------------------------------------------------------
