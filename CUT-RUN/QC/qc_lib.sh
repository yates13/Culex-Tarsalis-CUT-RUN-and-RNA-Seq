#!/usr/bin/env bash
################################################################################
# qc_lib.sh: shared paths and helpers for the per-step QC scripts.
# Not run on its own. Each QC script sets STEP and then sources this file:
#   STEP=align; source <repo>/CUT-RUN/QC/qc_lib.sh
# Paths are defined here once. Change them here if the layout changes.
#
# Each QC script writes, in Round1_Round2_With_Duplicates/QC/:
#   <step>_qc.tsv        one row per sample or condition
#   <step>_summary.txt   PASS / WARN / FAIL lines
# FAIL = missing, broken, or paired with the wrong day: the job exits 1, so
#        anything queued behind it with --dependency=afterok is cancelled.
# WARN = outlier (below WARN_FRAC x the median of its group); job exits 0.
################################################################################

set -uo pipefail

SCRATCH=/scratch/alpine/.colostate.edu/c832500103
PROJ=/projects/.colostate.edu/c832500103          # no "@" (see GOTCHAS.md)
PIPE=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates
RAW=${SCRATCH}/Cxt_Cut_Run_Raw_Files/00_Merged_Raw
FASTP=${PIPE}/01_Cxt_FastP/01.1_Cxt_fastp_out
SAMDIR=${PIPE}/02_Cxt_Alignment/02.2_Cxt_HisatAligned_Cleaned
BAMDIR=${PIPE}/02_Cxt_Alignment/02.3_Cxt_BamConverted
FILTDIR=${PIPE}/02_Cxt_Alignment/02.3.1_Cxt_BamFiltered   # Step 2.3.1 output (MAPQ>=10, primary, proper pairs)
MERGED=${PIPE}/02_Cxt_Alignment/02.4_Cxt_BamMerged_with_replicate_count
LIST=${PIPE}/02_Cxt_Alignment/FastqTrimmed.txt
D4=${PIPE}/04_Deeptools
BWDIR=${D4}/00_BigWig_Files_Individual_Bams
PEAKS=${PIPE}/03_MACS2_Peaks_New_Peaks/03.1_MacsPeaks_with_replicate_count
TSS=${SCRATCH}/Cxt_Cut_Run_Pipeline/00_Cxt_Genome/Genome_gff/CtarK1_TSS_excl_gene5911.bed
QCDIR=${PIPE}/QC
WARN_FRAC=0.25      # WARN if below 25% of the group median
MIN_ALIGN=20        # WARN if HISAT2 overall alignment rate is below this %

mkdir -p "$QCDIR"
TSV="$QCDIR/${STEP}_qc.tsv"
SUM="$QCDIR/${STEP}_summary.txt"
: > "$SUM"
nfail=0; nwarn=0
pass() { echo "PASS  $*" | tee -a "$SUM"; }
warn() { echo "WARN  $*" | tee -a "$SUM"; nwarn=$((nwarn+1)); }
fail() { echo "FAIL  $*" | tee -a "$SUM"; nfail=$((nfail+1)); }
need() { [[ -e "$1" ]] || { fail "missing: $1"; finish; }; }

# Samples expected, from the Step 2.2 list (one R1 trimmed FASTQ per sample)
samples() { sed 's/_R1_trimmed.fastq.gz$//' "$LIST"; }

# Flag rows whose column $2 is below WARN_FRAC x median of their group (column $1).
# Input: TSV with a header. Prints "name group value median" for each outlier.
low_outliers() {   # $1 group column, $2 value column, $3 file
    local med="$QCDIR/.medians.$$"
    awk -F'\t' -v g="$1" -v c="$2" 'NR>1 {print $g "\t" $c}' "$3" | sort -k1,1 -k2,2g |
        awk -F'\t' '{a[$1]=a[$1] " " $2; n[$1]++}
            END {for (k in a) {split(substr(a[k],2), x, " ");
                 m = (n[k]%2 ? x[(n[k]+1)/2] : (x[n[k]/2]+x[n[k]/2+1])/2); print k "\t" m}}' > "$med"
    awk -F'\t' -v g="$1" -v c="$2" -v f="$WARN_FRAC" 'NR==FNR {m[$1]=$2; next}
        FNR>1 && $c < f*m[$g] {print $1, $g, $c, m[$g]}' "$med" "$3"
    rm -f "$med"
}

load_samtools() {
    set +u
    source ${PROJ}/miniforge3/etc/profile.d/conda.sh && conda activate ${PROJ}/conda_envs/hisat2_env_fixed
    set -u
    command -v samtools >/dev/null || { fail "samtools not found in hisat2_env_fixed"; finish; }
}

finish() {
    echo "---- ${STEP}: $nfail FAIL, $nwarn WARN  ($(date))" | tee -a "$SUM"
    [[ -f "$TSV" ]] && echo "table:   $TSV"
    echo "summary: $SUM"
    [[ $nfail -eq 0 ]] && exit 0 || exit 1
}

need "$LIST"
N=$(samples | wc -l)
