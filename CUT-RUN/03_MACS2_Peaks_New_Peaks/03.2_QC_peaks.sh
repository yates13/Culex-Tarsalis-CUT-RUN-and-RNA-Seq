#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=QC_peaks
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=0:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

################################################################################
# QC for Step 3 (MACS2 peaks). Run after 03.1_New_flags_Macs2_PeakCalling_with_replicate_count.sh.
# Checks: a peak file per condition, > 0 peaks, MACS2 control from the SAME DAY. WARN: peaks < 25% of the target median.
# Output: QC/peaks_qc.tsv and QC/peaks_summary.txt in Round1_Round2_With_Duplicates.
# FAIL exits 1 (cancels anything queued behind it with afterok); WARN exits 0.
#
# Run from 03_MACS2_Peaks_New_Peaks (it needs a logs/ folder):
#   bash   03.2_QC_peaks.sh                       # quick; sbatch also works
#   sbatch --dependency=afterok:<jobid of the step> 03.2_QC_peaks.sh
################################################################################

STEP=peaks
source /scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/QC/qc_lib.sh

printf 'sample\ttarget\tday\tpeaks\tcontrol\tcontrol_day\n' > "$TSV"
conds=$(samples | grep -v -- '-IN-' | sed 's/-R[0-9]*-Cxt$//' | sort -u)
for c in $conds; do
    d=$(ls -d "$PEAKS/${c}"-rep* 2>/dev/null | head -1)
    [[ -n "$d" ]] || { fail "$c: no MACS2 output folder"; continue; }
    s=$(basename "$d"); t=$(echo $c | cut -d- -f2); day=$(echo $c | cut -d- -f3)
    pk=$(ls "$d"/*Peak 2>/dev/null | grep -v gapped | head -1)
    [[ -n "$pk" ]] || { fail "$s: no narrowPeak/broadPeak file"; continue; }
    np=$(wc -l < "$pk")
    ctl=$(grep -m1 '^# control file' "$d"/*_peaks.xls 2>/dev/null | sed 's/.*= *//' | xargs -r basename)
    cday=$(echo "$ctl" | cut -d- -f3)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$t" "$day" "$np" "${ctl:-NA}" "${cday:-NA}" >> "$TSV"
    [[ "$cday" == "$day" ]] || fail "$s ($day) used control ${ctl:-unknown} (${cday:-?})"
    [[ $np -gt 0 ]] || fail "$s: 0 peaks"
done
n=$(($(wc -l < "$TSV")-1)); [[ $n -eq $(echo "$conds" | wc -l) ]] && pass "$n conditions have peaks"
[[ $(awk -F'\t' 'NR>1 && $3!=$6' "$TSV" | wc -l) -eq 0 && $n -gt 0 ]] && pass "every condition used a same-day control"
while read -r s g v m; do warn "$s: $v peaks (median for $g: $m)"; done < <(low_outliers 2 4 "$TSV")

finish
