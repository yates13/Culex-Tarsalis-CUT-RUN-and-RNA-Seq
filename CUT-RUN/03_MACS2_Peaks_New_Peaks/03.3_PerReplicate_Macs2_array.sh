#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=MACS2_PerReplicate
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=1-72
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shawn.yates@colostate.edu

################################################################################
# STEP 3.3: Per-replicate MACS2 peaks (for the peak reproducibility check only)
#
# PURPOSE: Call peaks on each of the 72 target replicates separately, against the
#          same pooled same-day input used in Step 3.1. Step 3.4 then counts how
#          many pooled (Step 3.1) peaks are found in >= 2 of 3 replicates.
#          REPORT ONLY: nothing downstream (DiffBind, figures) uses these peaks.
#
# USAGE:   mkdir -p logs && sbatch 03.3_PerReplicate_Macs2_array.sh
#          (array 1-72, one task per target replicate; inputs are skipped)
#
# INPUT:   filtered BAMs  02.3.1_Cxt_BamFiltered/<SITE>-<MARK>-<DAY>-R<n>-Cxt.bam
#          pooled inputs  02.4_Cxt_BamMerged_with_replicate_count/<SITE>-IN-<DAY>-*-merged-sorted.bam
# OUTPUT:  03_MACS2_Peaks_New_Peaks/03.3_PerReplicate_Peaks_2026-11-27/<sample>/
#          (refuses to run if the sample folder already exists)
# ENV:     macs2_env (MACS2 2.2.9.1)
#
# SETTINGS: identical to Step 3.1 so per-replicate and pooled peaks are comparable:
#   -f BAMPE, -g 800000000, --keep-dup all, --extsize 250 (ignored with BAMPE),
#   --slocal 2372, default q 0.05, --broad for ME (H3K9me3).
#   Lines tagged [EDIT] are paths or settings you may change; [FIXED] match Step 3.1.
################################################################################

set -eo pipefail

# ---------------- USER SETTINGS ----------------
SCRATCH=/scratch/alpine/.colostate.edu/c832500103                                   # [EDIT] @-free alias of the scratch path
PROJ=/projects/c832500103@colostate.edu                                             # [EDIT]
RUN=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates                   # [EDIT]
BAMDIR=${RUN}/02_Cxt_Alignment/02.3.1_Cxt_BamFiltered                                # [EDIT]
INDIR=${RUN}/02_Cxt_Alignment/02.4_Cxt_BamMerged_with_replicate_count               # [EDIT]
OUTDIR=${RUN}/03_MACS2_Peaks_New_Peaks/03.3_PerReplicate_Peaks_2026-11-27            # [EDIT]
EXPECTED=72                                                                         # [EDIT] number of target replicates
GENOME_SIZE=800000000                                                               # [FIXED]
FORMAT="BAMPE"                                                                      # [FIXED]
KEEP_DUP="all"                                                                      # [FIXED]
# ------------------------------------------------

source ${PROJ}/miniforge3/etc/profile.d/conda.sh
conda activate ${PROJ}/conda_envs/macs2_env

# Target replicate BAMs only (inputs excluded), in a fixed order
mapfile -t BAMS < <(ls ${BAMDIR}/*-Cxt.bam | grep -v -- '-IN-' | sort)
if [[ ${#BAMS[@]} -ne ${EXPECTED} ]]; then
    echo "ERROR: found ${#BAMS[@]} target BAMs, expected ${EXPECTED}" >&2; exit 1
fi

BAM=${BAMS[$((SLURM_ARRAY_TASK_ID - 1))]}
SAMPLE=$(basename "$BAM" | sed 's/-Cxt\.bam$//')
SITE=$(echo "$SAMPLE" | cut -d- -f1)
MARK=$(echo "$SAMPLE" | cut -d- -f2)
DAY=$(echo "$SAMPLE" | cut -d- -f3)

# Same pooled same-day input as Step 3.1; exactly one must match
mapfile -t INS < <(ls ${INDIR}/${SITE}-IN-${DAY}-*-merged-sorted.bam 2>/dev/null)
if [[ ${#INS[@]} -ne 1 ]]; then
    echo "ERROR: ${SAMPLE}: found ${#INS[@]} pooled inputs for ${SITE} ${DAY}, expected 1" >&2; exit 1
fi
INPUT_BAM=${INS[0]}

PEAK_DIR=${OUTDIR}/${SAMPLE}
if [[ -e ${PEAK_DIR} ]]; then
    echo "ERROR: ${PEAK_DIR} already exists; refusing to overwrite" >&2; exit 1
fi
mkdir -p ${PEAK_DIR}

if [[ "$MARK" == "ME" ]]; then MACS_FLAGS="--broad"; else MACS_FLAGS=""; fi       # [FIXED]

echo "Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE} | input $(basename ${INPUT_BAM}) | flags '${MACS_FLAGS}'"

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
    > ${PEAK_DIR}/${SAMPLE}.macs2.log 2>&1

echo "Done: $(ls ${PEAK_DIR}/*_peaks.*Peak | xargs -n1 basename) ($(cat ${PEAK_DIR}/*_peaks.*Peak | wc -l) peaks)"
