#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=QC_align
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

################################################################################
# QC for Steps 2.2 + 2.3 (HISAT2, BAMs). Run after 02.3_runBam_array.sh has finished.
# Checks: HISAT2 log complete, BAM passes quickcheck, .bai present, flagstat. WARN: alignment rate < 20%.
# Output: QC/align_qc.tsv and QC/align_summary.txt in Round1_Round2_With_Duplicates.
# FAIL exits 1 (cancels anything queued behind it with afterok); WARN exits 0.
#
# Run from 02_Cxt_Alignment (it needs a logs/ folder):
#   sbatch 02.3.5_QC_align.sh                       # uses samtools, so submit it
#   sbatch --dependency=afterok:<jobid of the step> 02.3.5_QC_align.sh
################################################################################

STEP=align
source /scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/QC/qc_lib.sh

load_samtools
printf 'sample\ttarget\thisat2_rate\tmapped_pct\tproperly_paired_pct\n' > "$TSV"
for s in $(samples); do
    log="$SAMDIR/${s}.hisat2.log"; bam="$BAMDIR/${s}.bam"; fs="$BAMDIR/${s}_flagstat.txt"
    rate=$(grep 'overall alignment rate' "$log" 2>/dev/null | sed 's/%.*//')
    [[ -n "$rate" ]] || { fail "$s: no overall alignment rate in $log (missing or truncated)"; continue; }
    [[ -f "$bam" && -f "$bam.bai" ]] || { fail "$s: BAM or .bai missing"; continue; }
    samtools quickcheck "$bam" || { fail "$s: BAM fails samtools quickcheck"; continue; }
    [[ -f "$fs" ]] || samtools flagstat "$bam" > "$fs"
    mp=$(grep -m1 'mapped (' "$fs" | sed 's/.*(\([0-9.]*\)%.*/\1/')
    pp=$(grep 'properly paired' "$fs" | sed 's/.*(\([0-9.]*\)%.*/\1/')
    printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$(echo $s | cut -d- -f2)" "$rate" "$mp" "$pp" >> "$TSV"
    awk -v r="$rate" -v m="$MIN_ALIGN" 'BEGIN{exit !(r<m)}' && warn "$s: HISAT2 alignment rate $rate% (below $MIN_ALIGN%)"
done
n=$(($(wc -l < "$TSV")-1)); [[ $n -eq $N ]] && pass "$N samples: HISAT2 log complete, BAM intact and indexed"
awk -F'\t' 'NR>1{print $3}' "$TSV" | sort -n | awk '{a[NR]=$1} END{if(NR) printf "INFO  HISAT2 rate min %s, median %s, max %s\n", a[1], a[int((NR+1)/2)], a[NR]}' | tee -a "$SUM"

finish
