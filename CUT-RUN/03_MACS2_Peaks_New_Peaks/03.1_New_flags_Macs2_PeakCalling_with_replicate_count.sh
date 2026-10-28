#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=MACS2_NewFlags_with_replicate_counts
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=24:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu

################################################################################
# STEP 3.1: Call peaks with MACS2 on every merged BAM (one job, loops over samples)
#
# INPUT:   *-merged-sorted.bam files from Step 2.4
#            (02.4_Cxt_BamMerged_with_replicate_count)
#          Names must look like <TREATMENT>-<TARGET>-<DAY>-rep<digits>-merged-sorted.bam
#            TREATMENT = field 1 (e.g. ZH, MP)
#            TARGET    = field 2 (AC = H3K27ac, ME = H3K9me, SR = SREBP, IN = input)
# OUTPUT:  in 03.1_MacsPeaks_with_replicate_count/<sample>/ : MACS2 peak files
#          macs2.log (all samples, appended) in 03.1_MacsPeaks_with_replicate_count/
# ENV:     macs2_env
# RESOURCES: 8 CPUs, 32 GB, 24 h for the whole loop
#
# HOW IT WORKS
#   - Input controls (names containing "-IN-") are skipped as samples.
#   - ME samples (H3K9me, broad domains) are called with --broad.
#     AC and SR samples use MACS2's default narrow peak calling.
#   - Each sample is compared against the input control matching its TREATMENT:
#       <TREATMENT>-IN-D1*-merged-sorted.bam  (first match is used)
#     If no input is found, the sample is skipped with a message.
#   - Samples with no "rep<digits>" in the name are reported as an error and skipped.
#
# NOT AN ARRAY JOB
#   There is no array size to change. The loop uses every *-merged-sorted.bam in
#   the folder, so new samples are picked up automatically. Rerunning redoes ALL
#   samples and overwrites their outputs. To run only some samples, move the
#   other BAMs out of the folder or add a filter to the loop.
#
# SETTINGS TO CHECK
#   GENOME_SIZE=800000000   effective genome size used by MACS2 (-g)
#   FORMAT=BAMPE            paired-end BAMs; MACS2 uses real fragment sizes
#   KEEP_DUP=all            duplicates kept (this is the "With_Duplicates" run)
#   --slocal 2372           local background window
#
# CHECKING THE RESULTS
#   ls <PEAK_OUT> | wc -l    # one folder per non-input sample
#   The job output (.out) prints "[n/total] sample" per sample, with SKIP or
#   ERROR lines for anything not run. macs2.log has full MACS2 messages.
#
# NOTES
#   - Job logs go to logs/ in the folder you submit from. Run "mkdir -p logs" before sbatch.
#   - Paths to check before running: SCRATCH, BAMDIR, PEAK_OUT, and the email above.
################################################################################

set -eo pipefail

# ----------------------------
# Environment setup
# ----------------------------
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh
conda activate ${SCRATCH}/conda_envs/macs2_env


# ----------------------------
# Input/output paths
# ----------------------------
BAMDIR=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment/02.4_Cxt_BamMerged_with_replicate_count
PEAK_OUT=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/03_MACS2_Peaks_New_Peaks/03.1_MacsPeaks_with_replicate_count

mkdir -p logs
mkdir -p ${PEAK_OUT}

# ----------------------------
# MACS2 settings
# ----------------------------
GENOME_SIZE=800000000
FORMAT="BAMPE"
KEEP_DUP="all"

# ----------------------------
# Count the merged BAMs (used for the [n/total] progress lines)
# ----------------------------
total=$(ls ${BAMDIR}/*-merged-sorted.bam 2>/dev/null | wc -l)
count=0

# ----------------------------
# Loop over every merged BAM
# ----------------------------
for BAM in ${BAMDIR}/*-merged-sorted.bam; do
    count=$((count + 1))
    SAMPLE=$(basename "$BAM" | sed 's/-merged-sorted\.bam$//')

    # Skip input controls
    if [[ "$SAMPLE" == *"-IN-"* ]]; then
        echo "[$count/$total] SKIP (input): $SAMPLE"
        continue
    fi

    # Extract target (AC, ME, or SR)
    TARGET=$(echo $SAMPLE | cut -d'-' -f2)
    TREATMENT=$(echo $SAMPLE | cut -d'-' -f1)

    # Extract replicate count (rep123, rep23, rep1, etc.)
    REPS=$(echo "$SAMPLE" | grep -o 'rep[0-9]\+' || true)

    if [[ -z "$REPS" ]]; then
        echo "[$count/$total] ERROR: No replicate count found in $SAMPLE"
        continue
    fi

    # Find matching input control
    INPUT_BAM=$(ls ${BAMDIR}/${TREATMENT}-IN-D1*-merged-sorted.bam 2>/dev/null | head -1 || true)

    if [[ -z "$INPUT_BAM" ]]; then
        echo "[$count/$total] SKIP: $SAMPLE (no input control for $TREATMENT)"
        continue
    fi

    # Get the input name for logging
    INPUT_NAME=$(basename "$INPUT_BAM")

    # One output folder per sample
    PEAK_DIR="${PEAK_OUT}/${SAMPLE}"
    mkdir -p ${PEAK_DIR}

    echo "[$count/$total] $SAMPLE (Target: $TARGET)"

    # ----------------------------
    # Choose peak type by target
    # ----------------------------
    if [[ "$TARGET" == "ME" ]]; then
        # H3K9me3 (methylated) - broad peaks
        MACS_FLAGS="--broad"
        echo "        Using --broad for H3K9me3 methylation marks"
    else
        # AC (H3K27ac) and SR (SREBP) - narrow peaks
        MACS_FLAGS=""
        echo "        Using narrow peaks for $TARGET"
    fi

    # ----------------------------
    # Run MACS2
    # ----------------------------
    macs2 callpeak \
        -t ${BAM} \
        -c ${INPUT_BAM} \
        -f ${FORMAT} \
        -g ${GENOME_SIZE} \
        -n ${SAMPLE} \
        --keep-dup ${KEEP_DUP} \
        --extsize 250 \
        --slocal 2372 \
        ${MACS_FLAGS} \
        --outdir ${PEAK_DIR} \
        2>&1 | tee -a ${PEAK_OUT}/macs2.log

    echo "     ✓ Done"
    echo ""
done

echo "Peak calling complete."
echo "Output: ${PEAK_OUT}"

