#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=QC_merge
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

################################################################################
# QC for Step 2.4 (merged BAMs). Run after 02.4_BamMerge_With_replicate_counts_array.sh has finished, before Step 3.
# Checks: one merged BAM per condition; rep tag = replicates on disk; merged reads = sum of replicates.
# Output: QC/merge_qc.tsv and QC/merge_summary.txt in Round1_Round2_With_Duplicates.
# FAIL exits 1 (cancels anything queued behind it with afterok); WARN exits 0.
#
# Run from 02_Cxt_Alignment (it needs a logs/ folder):
#   sbatch 02.4.5_QC_merge.sh                       # uses samtools, so submit it
#   sbatch --dependency=afterok:<jobid of the step> 02.4.5_QC_merge.sh
################################################################################

STEP=merge
source /scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/QC/qc_lib.sh

load_samtools
printf 'condition\treps\tmerged_reads\tsum_of_rep_reads\n' > "$TSV"
conds=$(samples | sed 's/-R[0-9]*-Cxt$//' | sort -u)
for c in $conds; do
    m=$(ls "$MERGED/${c}"-rep*-merged-sorted.bam 2>/dev/null | head -1)
    [[ -n "$m" && -f "$m.bai" ]] || { fail "$c: merged BAM or .bai missing"; continue; }
    samtools quickcheck "$m" || { fail "$c: merged BAM fails quickcheck"; continue; }
    tag=$(basename "$m" | grep -o 'rep[0-9]*')
    want="rep$(samples | grep "^${c}-R" | sed -E 's/.*-R([0-9]+)-Cxt$/\1/' | sort -n | tr -d '\n')"
    [[ "$tag" == "$want" ]] || fail "$c: merged as $tag, but the replicates on disk are $want"
    mr=$(samtools idxstats "$m" | awk '{s+=$3}END{print s+0}')
    rr=0; for b in "$BAMDIR/${c}"-R*-Cxt.bam; do rr=$((rr + $(samtools idxstats "$b" | awk '{s+=$3}END{print s+0}'))); done
    printf '%s\t%s\t%s\t%s\n' "$c" "$tag" "$mr" "$rr" >> "$TSV"
    [[ $mr -eq $rr ]] || fail "$c: merged BAM has $mr mapped reads, replicates add up to $rr"
done
n=$(($(wc -l < "$TSV")-1)); nc=$(echo "$conds" | wc -l)
[[ $n -eq $nc && $nfail -eq 0 ]] && pass "$nc conditions merged; replicate tags and read totals match"

finish
