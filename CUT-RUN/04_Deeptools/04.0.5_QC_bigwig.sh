#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=QC_bigwig
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=0:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

################################################################################
# QC for Step 4.0 (bigWigs). Run after 00_Create_Bigwigs_Individual_Bams_ARRAY.sh, before 4.1.
# Checks: a non-empty bigWig per treatment BAM; no input bigWigs; the input used (from the newest run logs) is the pooled input from the SAME DAY.
# Output: QC/bigwig_qc.tsv and QC/bigwig_summary.txt in Round1_Round2_With_Duplicates.
# FAIL exits 1 (cancels anything queued behind it with afterok); WARN exits 0.
#
# Run from 04_Deeptools (it needs a logs/ folder):
#   bash   04.0.5_QC_bigwig.sh                       # quick; sbatch also works
#   sbatch --dependency=afterok:<jobid of the step> 04.0.5_QC_bigwig.sh
################################################################################

STEP=bigwig
source /scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/QC/qc_lib.sh

printf 'sample\tday\tinput_used\tinput_day\tbw_MB\n' > "$TSV"
# Input used: from the newest run's logs ("  input:" line of each task)
jid=$(ls "$D4"/logs/BamToBigwig_Individual_bams.*_*.log 2>/dev/null | sed -E 's/.*\.([0-9]+)_[0-9]+\.log$/\1/' | sort -n | tail -1)
[[ -n "$jid" ]] || warn "no bigWig logs in $D4/logs; cannot check which inputs were used"
for s in $(samples | grep -v -- '-IN-'); do
    b=${s%-Cxt}; bw="$BWDIR/${b}_inputsubtracted.bw"
    [[ -s "$bw" ]] || { fail "$b: bigWig missing or empty"; continue; }
    mb=$(du -m "$bw" | cut -f1); day=$(echo $b | cut -d- -f3)
    inp=""; [[ -n "$jid" ]] && inp=$(grep -l "treatment: .*/${s}.bam\$" "$D4"/logs/BamToBigwig_Individual_bams.${jid}_*.log 2>/dev/null | head -1 | xargs -r grep -m1 '  input:' | awk '{print $2}' | xargs -r basename)
    iday=$(echo "$inp" | cut -d- -f3)
    printf '%s\t%s\t%s\t%s\t%s\n' "$b" "$day" "${inp:-NA}" "${iday:-NA}" "$mb" >> "$TSV"
    [[ -z "$inp" || "$iday" == "$day" ]] || fail "$b ($day) was normalized with $inp ($iday)"
    [[ -z "$inp" || "$inp" == *-merged-sorted.bam ]] || fail "$b: input $inp is one replicate, not the pooled same-day input"
done
n=$(($(wc -l < "$TSV")-1)); [[ $n -eq $(samples | grep -vc -- '-IN-') ]] && pass "$n bigWigs present and non-empty"
[[ -n "$jid" ]] && [[ $(awk -F'\t' 'NR>1 && $2!=$4' "$TSV" | wc -l) -eq 0 ]] && pass "every bigWig used an input from the same day (job $jid)"
[[ $(ls "$BWDIR" | grep -c -- '-IN-') -eq 0 ]] || fail "input samples found among the bigWigs"

finish
