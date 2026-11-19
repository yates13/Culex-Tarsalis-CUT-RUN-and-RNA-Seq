#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=QC_fastp
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=0:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

################################################################################
# QC for Step 1.1 (fastp). Run after 01.1_fastp_Cxt.sh and 01.1_fix_sample_names.sh.
# Checks: trimmed R1/R2 and fastp report for every sample. WARN: reads after filtering < 25% of the target median.
# Output: QC/fastp_qc.tsv and QC/fastp_summary.txt in Round1_Round2_With_Duplicates.
# FAIL exits 1 (cancels anything queued behind it with afterok); WARN exits 0.
#
# Run from 01_Cxt_FastP (it needs a logs/ folder):
#   bash   01.1.5_QC_fastp.sh                       # quick; sbatch also works
#   sbatch --dependency=afterok:<jobid of the step> 01.1.5_QC_fastp.sh
################################################################################

STEP=fastp
source /scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/QC/qc_lib.sh

printf 'sample\ttarget\treads_after_filter\tdup_pct\n' > "$TSV"
nmiss=0
for s in $(samples); do
    j="$FASTP/${s}_fastp.json"; [[ -f $j ]] || j="$FASTP/${s/-Cxt/_Cxt}_fastp.json"
    [[ -f "$FASTP/${s}_R1_trimmed.fastq.gz" && -f "$FASTP/${s}_R2_trimmed.fastq.gz" ]] || { fail "$s: trimmed R1/R2 missing"; continue; }
    [[ -f $j ]] || { fail "$s: fastp json missing"; nmiss=$((nmiss+1)); continue; }
    reads=$(grep -m2 '"total_reads"' "$j" | tail -1 | tr -dc '0-9')
    dup=$(grep -A1 '"duplication"' "$j" | grep '"rate"' | head -1 | sed 's/.*: *//; s/,//')
    printf '%s\t%s\t%s\t%.2f\n' "$s" "$(echo $s | cut -d- -f2)" "$reads" "$(echo "$dup*100" | bc -l)" >> "$TSV"
done
[[ $nmiss -eq 0 ]] && pass "$N samples: trimmed R1/R2 and fastp report present"
while read -r s g v m; do warn "$s: $v reads after filtering (median for $g: $m)"; done < <(low_outliers 2 3 "$TSV")

finish
