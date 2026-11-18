#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=CheckMD5_array
#SBATCH --time=1:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu  # CHANGE to your own address

################################################################################
# STEP 1.0.2: Verify the checksum of every raw FASTQ (one file per array task)
#
# INPUT:   INDIR (00_Merged_Raw from step 1.0.1): <file>.fastq.gz and
#          <file>.fastq.gz.md5 side by side (symlinks are fine)
# OUTPUT:  INDIR/md5_file_list.txt     line N = array task N
#          INDIR/md5_results.tsv       status<TAB>file   (OK, FAIL, MISSING_MD5)
#          INDIR/logs/CheckMD5_array.<jobid>_<task>.out/.err
#
# USAGE (login node; builds the list and submits itself with the right
# array size, so no editing --array):
#   bash 01.0.2_Md5sum_Check.sh
#
# CHECK RESULTS (after the job finishes):
#   cd $INDIR
#   cut -f1 md5_results.tsv | sort | uniq -c    # expect only OK
#   grep -v '^OK' md5_results.tsv               # anything that failed
#   wc -l < md5_file_list.txt; wc -l < md5_results.tsv   # must match
# Resubmitting clears md5_results.tsv first.
################################################################################

set -euo pipefail

INDIR="${INDIR:-/scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Raw_Files/00_Merged_Raw}"
FILE_LIST="${INDIR}/md5_file_list.txt"
RESULTS="${INDIR}/md5_results.tsv"

SCRIPT=$(readlink -f "$0")   # resolve before cd, or a relative $0 breaks
cd "$INDIR"

# ---------------------------------------------------------------- submit mode
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    mkdir -p logs
    ls *.fastq.gz | sort > "$FILE_LIST"
    N=$(wc -l < "$FILE_LIST")
    : > "$RESULTS"
    echo "${N} FASTQ files listed in ${FILE_LIST}"
    sbatch --array=1-"${N}" \
        --output="logs/%x.%A_%a.out" \
        --error="logs/%x.%A_%a.err" \
        "$SCRIPT"
    exit 0
fi

# ------------------------------------------------------------------ task mode
gz=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")
if [[ -z "$gz" ]]; then
    echo "ERROR: no file for array index $SLURM_ARRAY_TASK_ID in $FILE_LIST"
    exit 1
fi

md5file="${gz}.md5"
echo "Task $SLURM_ARRAY_TASK_ID: checking $gz -> $(readlink -f "$gz")"

if [[ ! -f "$md5file" ]]; then
    printf 'MISSING_MD5\t%s\n' "$gz" >> "$RESULTS"
    echo "MISSING checksum file: $md5file"
    exit 1
fi

expected=$(awk 'NF{print $1; exit}' "$md5file")
actual=$(md5sum "$gz" | awk '{print $1}')

if [[ "$expected" == "$actual" ]]; then
    printf 'OK\t%s\n' "$gz" >> "$RESULTS"
    echo "OK: $gz"
else
    printf 'FAIL\t%s\n' "$gz" >> "$RESULTS"
    echo "FAIL: $gz (expected $expected, got $actual)"
    exit 1
fi
