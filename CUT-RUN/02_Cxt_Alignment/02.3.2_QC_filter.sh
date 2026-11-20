#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=QC_filter
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

################################################################################
# QC for Step 2.3.1 (filtered BAMs). Run after 02.3.1_Filter_Bams_array.sh, before 2.4.
# Checks: a filtered BAM + .bai per sample that passes quickcheck; reads kept > 0;
#   WARN if a sample kept far fewer reads than others of its target (low_outliers).
# Output: QC/filter_qc.tsv and QC/filter_summary.txt in Round1_Round2_With_Duplicates.
# FAIL exits 1 (cancels anything queued behind it with afterok); WARN exits 0.
#
# Run from 02_Cxt_Alignment (it needs a logs/ folder):
#   sbatch --dependency=afterok:<jobid of 2.3.1> 02.3.2_QC_filter.sh
################################################################################

STEP=filter
source /scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/QC/qc_lib.sh

load_samtools
printf 'sample\ttarget\tprimary_before\tkept\tpct_kept\n' > "$TSV"
for s in $(samples); do
    b="$FILTDIR/${s}.bam"; c="$FILTDIR/${s}_filter_counts.txt"
    [[ -f "$b" && -f "$b.bai" && -f "$c" ]] || { fail "$s: filtered BAM, .bai or counts file missing"; continue; }
    samtools quickcheck "$b" || { fail "$s: filtered BAM fails quickcheck"; continue; }
    read -r _ before kept pct < <(sed -n 2p "$c")
    [[ "$kept" -gt 0 ]] || fail "$s: no reads kept"
    printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$(echo $s | cut -d- -f2)" "$before" "$kept" "$pct" >> "$TSV"
done
low_outliers 2 4 "$TSV" | while read -r s g v m; do warn "$s: $v reads kept, below ${WARN_FRAC} x the $g median ($m)"; done
nwarn=$(grep -c '^WARN' "$SUM" || true)
n=$(($(wc -l < "$TSV")-1))
[[ $n -eq $N ]] && pass "$n filtered BAMs present; median % kept by target:" && \
    awk -F'\t' 'NR>1{print $2"\t"$5}' "$TSV" | sort -k1,1 -k2,2g | awk -F'\t' '{a[$1]=a[$1]" "$2;n[$1]++} END{for(k in a){split(substr(a[k],2),x," "); print "      " k ": " x[int((n[k]+1)/2)] "%"}}' | tee -a "$SUM"

finish

