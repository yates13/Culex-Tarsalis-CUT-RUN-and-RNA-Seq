#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=QC_raw
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=0:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

################################################################################
# QC for Step 1.0 (merge + checksums). Run after 01.0.2_Md5sum_Check.sh has finished.
# Checks: FASTQ count = 2 x samples; every checksum OK. WARN: files with no .md5.
# Output: QC/raw_qc.tsv and QC/raw_summary.txt in Round1_Round2_With_Duplicates.
# FAIL exits 1 (cancels anything queued behind it with afterok); WARN exits 0.
#
# Run from 01_Cxt_FastP (it needs a logs/ folder):
#   bash   01.0.3_QC_raw.sh                       # quick; sbatch also works
#   sbatch --dependency=afterok:<jobid of the step> 01.0.3_QC_raw.sh
################################################################################

STEP=raw
source /scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/QC/qc_lib.sh

need "$RAW/md5_results.tsv"
nfq=$(ls "$RAW"/*.fastq.gz 2>/dev/null | wc -l)
[[ $nfq -eq $((2*N)) ]] && pass "$nfq FASTQs in 00_Merged_Raw (2 x $N samples)" || fail "$nfq FASTQs in 00_Merged_Raw, expected $((2*N))"
cp "$RAW/md5_results.tsv" "$TSV"
ok=$(grep -c '^OK' "$TSV"); bad=$(grep -c '^FAIL' "$TSV"); mis=$(grep -c '^MISSING' "$TSV")
[[ $((ok+bad+mis)) -eq $nfq ]] || fail "md5_results.tsv has $((ok+bad+mis)) rows for $nfq FASTQs (some tasks did not finish)"
[[ $bad -eq 0 ]] && pass "no checksum failures ($ok OK)" || fail "$bad checksum FAIL: $(grep '^FAIL' "$TSV" | cut -f2 | tr '\n' ' ')"
[[ $mis -eq 0 ]] || warn "$mis file(s) with no .md5 to check: $(grep '^MISSING' "$TSV" | cut -f2 | tr '\n' ' ')"

finish
