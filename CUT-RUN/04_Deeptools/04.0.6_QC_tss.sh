#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=QC_tss
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=0:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

################################################################################
# QC for the TSS BED used by Step 4.1. Run after building CtarK1_TSS.bed (00_Cxt_Genome/*instructions.txt), before 4.1.
# Checks: TSS BED rows are 1 bp with a strand; every TSS contig name is in the BAM header. WARN: gene count far from ~20k.
# Output: QC/tss_qc.tsv and QC/tss_summary.txt in Round1_Round2_With_Duplicates.
# FAIL exits 1 (cancels anything queued behind it with afterok); WARN exits 0.
#
# Run from 04_Deeptools (it needs a logs/ folder):
#   sbatch 04.0.6_QC_tss.sh                       # uses samtools, so submit it
#   sbatch --dependency=afterok:<jobid of the step> 04.0.6_QC_tss.sh
################################################################################

STEP=tss
source /scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN/QC/qc_lib.sh

load_samtools
need "$TSS"
n=$(wc -l < "$TSS"); bad=$(awk -F'\t' 'NF<6 || $3!=$2+1 || ($6!="+" && $6!="-")' "$TSS" | wc -l)
[[ $bad -eq 0 ]] && pass "$n TSS rows, all 1 bp with a strand" || fail "$bad malformed TSS rows (need 6 columns, end = start + 1, strand + or -)"
[[ $n -gt 5000 && $n -lt 50000 ]] || warn "$n TSS rows; expected roughly 15-25k genes"
samtools view -H "$(ls "$BAMDIR"/*.bam | head -1)" | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); print $2}' | sort -u > "$QCDIR/.bam_contigs"
cut -f1 "$TSS" | sort -u > "$QCDIR/.tss_contigs"
miss=$(comm -13 "$QCDIR/.bam_contigs" "$QCDIR/.tss_contigs" | wc -l)
{ echo -e "tss_contigs\t$(wc -l < "$QCDIR/.tss_contigs")"; echo -e "bam_contigs\t$(wc -l < "$QCDIR/.bam_contigs")"; echo -e "tss_contigs_not_in_bam\t$miss"; } > "$TSV"
[[ $miss -eq 0 ]] && pass "every TSS contig is in the BAM header" || fail "$miss TSS contig name(s) not in the BAMs, e.g. $(comm -13 "$QCDIR/.bam_contigs" "$QCDIR/.tss_contigs" | head -3 | tr '\n' ' ')"

finish
