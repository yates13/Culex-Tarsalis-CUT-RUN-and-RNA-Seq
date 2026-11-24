#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=Greenscreen_Inputs
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=24:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shawn.yates@colostate.edu

################################################################################
# STEP 3.0: Build a greenscreen (artifact-region list) from the input controls
#
# METHOD:  Klasfeld, Roule & Wagner 2022, Plant Cell 34:4795 (doi 10.1093/plcell/koac282)
#          github.com/sklasfeld/GreenscreenProject
#            1. MACS2 broad peaks on EACH input separately, no local lambda
#            2. keep input peaks with -log10(q) >= 10 (broadPeak column 9)
#            3. pool all inputs' peaks, merge regions within MERGE_DIST bp
#            4. keep merged regions found in >= MIN_INPUTS distinct inputs
#
# INPUT:   the 24 individual filtered input BAMs from Step 2.3.1
#            02.3.1_Cxt_BamFiltered/<TREATMENT>-IN-<DAY>-R<n>-Cxt.bam
#            (-q 10 -f 2 -F 2308: the same reads used for the bigWigs)
# OUTPUT:  in 03.0_Greenscreen/ (NEW folder; nothing existing is changed):
#            macs2_inputs/<sample>/        MACS2 output per input
#            input_peaks_q10/<sample>.bed  peaks passing the q filter
#            greenscreen_all_merged.bed    every merged region + input count (col 4)
#            Cxt_greenscreen_<N>of24_<D>bp.bed   FINAL list (regions in >= N inputs)
#            greenscreen_summary.txt       settings, genome size, counts, checks
#            macs2.log                     all MACS2 messages
# ENV:     macs2_env (macs2 2.2.9.1); bedtools taken from another env (see below)
# RESOURCES: 2 CPUs, 32 GB, 24 h (24 MACS2 runs in a loop; MACS2 is single-threaded)
#
# SAFETY
#   Refuses to run if 03.0_Greenscreen/ already exists, so earlier results are
#   never overwritten. To rerun with other settings, rename the old folder or
#   change OUT_NAME below.
#   This step changes no existing peaks, matrices or figures. Applying the list
#   (03.1 peaks, 4.1 --blackListFileName) is done later, into new output folders.
#
# SETTINGS TO CHECK
#   QVAL_MIN=10       -log10 q cutoff for input peaks (paper: 10)
#   MERGE_DIST=10000  merge distance in bp (paper: 5 kb Arabidopsis, 10 kb rice,
#                     20 kb human; 10 kb chosen for the ~800 Mb CtarK1 genome)
#   MIN_INPUTS=12     at least half of the 24 inputs (paper: half)
#   GENOME_SIZE       computed below as the non-N bases of CtarK1.fa
#   FORMAT=BAMPE      paired-end; real fragments are used, so no --extsize
#
# NOTES
#   - Run "mkdir -p logs" in the folder you submit from before sbatch.
################################################################################

set -eo pipefail

# ----------------------------
# Environment setup
# ----------------------------
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
PROJ=/projects/c832500103@colostate.edu
source ${PROJ}/miniforge3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found at ${PROJ}/miniforge3" >&2; exit 1; }
conda activate ${PROJ}/conda_envs/macs2_env || { echo "ERROR: conda activate macs2_env failed" >&2; exit 1; }

# bedtools is not in macs2_env; use the first conda env that has it
BEDTOOLS=$(ls ${PROJ}/conda_envs/*/bin/bedtools 2>/dev/null | head -1 || true)
if [[ -z "$BEDTOOLS" || ! -x "$BEDTOOLS" ]]; then
    echo "ERROR: no bedtools found in ${PROJ}/conda_envs/*/bin/" >&2
    exit 1
fi

# ----------------------------
# Input/output paths
# ----------------------------
RUN=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates
BAMDIR=${RUN}/02_Cxt_Alignment/02.3.1_Cxt_BamFiltered
GENOME_FA=${SCRATCH}/Cxt_Cut_Run_Pipeline/00_Cxt_Indexed_Genome/CtarK1.fa
OUT_NAME=03.0_Greenscreen
OUT=${RUN}/03_MACS2_Peaks_New_Peaks/${OUT_NAME}

# ----------------------------
# Greenscreen settings
# ----------------------------
FORMAT="BAMPE"
QVAL_MIN=10
MERGE_DIST=10000
MIN_INPUTS=12
EXPECTED_INPUTS=24
EXCLUDE_CONTIG="Mt"

# ----------------------------
# Checks before anything is written
# ----------------------------
if [[ -e "$OUT" ]]; then
    echo "ERROR: $OUT already exists. Rename it or change OUT_NAME; nothing was overwritten." >&2
    exit 1
fi
[[ -s "$GENOME_FA" ]] || { echo "ERROR: genome not found: $GENOME_FA" >&2; exit 1; }
command -v macs2 &> /dev/null || { echo "ERROR: macs2 not found" >&2; exit 1; }

INPUT_BAMS=( $(ls ${BAMDIR}/*-IN-*-Cxt.bam 2>/dev/null || true) )
if [[ ${#INPUT_BAMS[@]} -ne $EXPECTED_INPUTS ]]; then
    echo "ERROR: found ${#INPUT_BAMS[@]} input BAMs in $BAMDIR, expected $EXPECTED_INPUTS" >&2
    exit 1
fi

mkdir -p logs
mkdir -p ${OUT}/macs2_inputs ${OUT}/input_peaks_q10
SUMMARY=${OUT}/greenscreen_summary.txt

# ----------------------------
# Effective genome size = non-N bases in CtarK1.fa
# ----------------------------
GENOME_SIZE=$(grep -v '^>' "$GENOME_FA" | tr -cd 'ACGTacgt' | wc -c)
TOTAL_BP=$(grep -v '^>' "$GENOME_FA" | tr -d '\n' | wc -c)
echo "Genome: $GENOME_FA"
echo "  total bases: $TOTAL_BP   non-N bases (MACS2 -g): $GENOME_SIZE"

{
    echo "Greenscreen built $(date)"
    echo "macs2: $(macs2 --version 2>&1)"
    echo "bedtools: $BEDTOOLS ($($BEDTOOLS --version))"
    echo "BAMs: $BAMDIR (${#INPUT_BAMS[@]} inputs)"
    echo "Genome: $GENOME_FA  total=$TOTAL_BP  non-N=$GENOME_SIZE"
    echo "Settings: FORMAT=$FORMAT QVAL_MIN=$QVAL_MIN MERGE_DIST=$MERGE_DIST MIN_INPUTS=$MIN_INPUTS EXCLUDE_CONTIG=$EXCLUDE_CONTIG"
    echo ""
    echo -e "input\tbroad_peaks\tpeaks_q>=${QVAL_MIN}"
} > "$SUMMARY"

# ----------------------------
# 1-2. MACS2 on each input, then the q-value filter
# ----------------------------
count=0
for BAM in "${INPUT_BAMS[@]}"; do
    count=$((count + 1))
    SAMPLE=$(basename "$BAM" | sed 's/-Cxt\.bam$//')
    PEAK_DIR=${OUT}/macs2_inputs/${SAMPLE}
    mkdir -p "$PEAK_DIR"
    echo "[$count/${#INPUT_BAMS[@]}] $SAMPLE"

    # --nomodel/--extsize are not used: with BAMPE MACS2 uses the real fragments
    macs2 callpeak \
        -t "$BAM" \
        -f ${FORMAT} \
        -g ${GENOME_SIZE} \
        -n ${SAMPLE} \
        --keep-dup auto \
        --broad \
        --nolambda \
        --outdir "$PEAK_DIR" \
        2>&1 | tee -a ${OUT}/macs2.log

    BROAD=${PEAK_DIR}/${SAMPLE}_peaks.broadPeak
    [[ -f "$BROAD" ]] || { echo "ERROR: MACS2 made no broadPeak for $SAMPLE" >&2; exit 1; }

    # keep -log10(q) >= QVAL_MIN, drop the mitochondrial contig, label by sample
    awk -v q=${QVAL_MIN} -v s=${SAMPLE} -v mt=${EXCLUDE_CONTIG} 'BEGIN{OFS="\t"} $9>=q && $1!=mt {print $1,$2,$3,s}' \
        "$BROAD" > ${OUT}/input_peaks_q10/${SAMPLE}.bed

    n_all=$(wc -l < "$BROAD")
    n_q=$(wc -l < ${OUT}/input_peaks_q10/${SAMPLE}.bed)
    echo -e "${SAMPLE}\t${n_all}\t${n_q}" >> "$SUMMARY"
    echo "     ✓ $n_all broad peaks, $n_q with -log10(q) >= $QVAL_MIN"
done

# ----------------------------
# 3. Pool, sort, merge within MERGE_DIST; col 4 = number of distinct inputs
# ----------------------------
cat ${OUT}/input_peaks_q10/*.bed \
    | sort -k1,1 -k2,2n \
    | "$BEDTOOLS" merge -i stdin -d ${MERGE_DIST} -c 4 -o count_distinct \
    > ${OUT}/greenscreen_all_merged.bed

# ----------------------------
# 4. Keep regions found in >= MIN_INPUTS inputs
# ----------------------------
FINAL=${OUT}/Cxt_greenscreen_${MIN_INPUTS}of${EXPECTED_INPUTS}_${MERGE_DIST}bp.bed
awk -v m=${MIN_INPUTS} 'BEGIN{OFS="\t"} $4>=m {print $1,$2,$3,"greenscreen",$4}' \
    ${OUT}/greenscreen_all_merged.bed > "$FINAL"

n_merged=$(wc -l < ${OUT}/greenscreen_all_merged.bed)
n_final=$(wc -l < "$FINAL")
bp_final=$(awk '{s+=$3-$2} END{print s+0}' "$FINAL")
pct=$(awk -v b=$bp_final -v g=$GENOME_SIZE 'BEGIN{printf "%.4f", 100*b/g}')

# Sanity check (report only): does the list cover the known gene5911 artifact
# at the start of contig tig00005036 (TSS at 116; +2 kb window)?
hit5911=$(awk '$1=="tig00005036" && $2<2200' "$FINAL" | wc -l)

{
    echo ""
    echo "Merged regions (any support): $n_merged"
    echo "Final regions (>= $MIN_INPUTS of $EXPECTED_INPUTS inputs): $n_final"
    echo "Final bp: $bp_final  (${pct}% of non-N genome)"
    echo "gene5911 region (tig00005036:0-2200) covered: $([[ $hit5911 -gt 0 ]] && echo YES || echo NO)"
    echo "Final file: $FINAL"
} | tee -a "$SUMMARY"

echo "Greenscreen complete."
echo "Output: ${OUT}"
