#!/usr/bin/env bash
################################################################################
# Rebuild everything downstream of Step 2.3 as one chain of dependent jobs:
#   2.3.1 filter BAMs -> 2.3.2 QC -> 2.4 merge -> 2.4.5 QC -> 3.1 MACS2 -> 3.2 QC
#                                 -> 4.0.4 input-correlation QC
#                                 2.4.5 QC -> 4.0 bigWigs (pooled same-day input)
#                                          -> 4.0.5 QC + 4.0.6 QC -> 4.1 -> 4.1.2 -> 4.3 heatmaps
# Every job waits for the one before it (afterok); if any step or QC fails, the
# jobs behind it are cancelled (--kill-on-invalid-dep). Representative-replicate
# figures (4.2/4.3 representatives) are left out: run them once REP_MAP is chosen.
#
# Run on a login node, NOT with sbatch (it only submits jobs):
#   bash <repo>/CUT-RUN/Rebuild_From_Step2.3.1.sh
# Start from clean outputs: move the old merged BAMs, peaks, bigWigs, matrices
# and figures aside first. bigWigs are skipped when they already exist.
# Job IDs are printed and saved to QC/rebuild_jobs_<date>.txt.
################################################################################

set -euo pipefail
REPO=/scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/CUT-RUN
P=/scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates
A=$P/02_Cxt_Alignment; M=$P/03_MACS2_Peaks_New_Peaks; D=$P/04_Deeptools
FILT=$A/02.3.1_Cxt_BamFiltered; MERGED=$A/02.4_Cxt_BamMerged_with_replicate_count
mkdir -p $A/logs $M/logs $D/logs $FILT $MERGED $P/QC

# Lists each array job reads (all names are known before the jobs run)
ls $A/02.3_Cxt_BamConverted/*-Cxt.bam > $FILT/bam_file_list.txt
xargs -n1 basename < $FILT/bam_file_list.txt | sed -E 's/-R[0-9]+-Cxt\.bam$//' | sort -u > $MERGED/conditions.txt
sed "s|/02.3_Cxt_BamConverted/|/02.3.1_Cxt_BamFiltered/|" $FILT/bam_file_list.txt | grep -v -- '-IN-' > $D/treatment_bam_list.txt
nb=$(wc -l < $FILT/bam_file_list.txt); nc=$(wc -l < $MERGED/conditions.txt); nt=$(wc -l < $D/treatment_bam_list.txt)
echo "BAMs: $nb   conditions: $nc   treatment BAMs: $nt"

sub() { local dir=$1 dep=$2; shift 2; sbatch --parsable -D "$dir" ${dep:+--dependency=afterok:$dep --kill-on-invalid-dep=yes} "$@"; }

j_filt=$(sub $A ""        --array=1-$nb $REPO/02_Cxt_Alignment/02.3.1_Filter_Bams_array.sh)
q_filt=$(sub $A $j_filt   $REPO/02_Cxt_Alignment/02.3.2_QC_filter.sh)
q_icor=$(sub $D $q_filt   $REPO/04_Deeptools/04.0.4_QC_input_correlation.sh)
j_merg=$(sub $A $q_filt   --array=1-$nc $REPO/02_Cxt_Alignment/02.4_BamMerge_With_replicate_counts_array.sh)
q_merg=$(sub $A $j_merg   $REPO/02_Cxt_Alignment/02.4.5_QC_merge.sh)
j_macs=$(sub $M $q_merg   $REPO/03_MACS2_Peaks_New_Peaks/03.1_New_flags_Macs2_PeakCalling_with_replicate_count.sh)
q_peak=$(sub $M $j_macs   $REPO/03_MACS2_Peaks_New_Peaks/03.2_QC_peaks.sh)
j_bw=$(sub $D $q_merg     --array=1-$nt $REPO/04_Deeptools/00_Create_Bigwigs_Individual_Bams_ARRAY.sh)
q_bw=$(sub $D $j_bw       $REPO/04_Deeptools/04.0.5_QC_bigwig.sh)
q_tss=$(sub $D $q_bw      $REPO/04_Deeptools/04.0.6_QC_tss.sh)
j_mat=$(sub $D $q_tss     $REPO/04_Deeptools/04.1_ComputeMatrix_Individual_array.sh)
j_cb=$(sub $D $j_mat      $REPO/04_Deeptools/04.1.2_ComputeMatrix_Individual_array_RUN_THIS_AFTER.sh)
j_hm=$(sub $D $j_cb       $REPO/04_Deeptools/04.3_HeatMap_Individual_Bams_By_Treatment_updated.sh 'Input-subtracted signal at TSS')

LOG=$P/QC/rebuild_jobs_$(date +%F).txt
cat <<EOF | tee "$LOG"
filter      $j_filt   QC_filter  $q_filt   QC_input_corr $q_icor
merge       $j_merg   QC_merge   $q_merg
MACS2       $j_macs   QC_peaks   $q_peak
bigWigs     $j_bw     QC_bigwig  $q_bw     QC_tss $q_tss
matrix      $j_mat    cbind      $j_cb     heatmaps $j_hm
EOF
echo "Watch: squeue -u \$USER     Saved: $LOG"
