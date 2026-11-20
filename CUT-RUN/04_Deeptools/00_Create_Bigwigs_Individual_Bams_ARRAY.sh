#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=BamToBigwig_Individual_bams
#SBATCH --output=logs/%x.%A_%a.log
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --time=5:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shawn.yates@colostate.edu
#SBATCH --array=1-72		# MUST equal the number of lines in treatment_bam_list.txt (see below)

################################################################################
# STEP 4.0: Make input-subtracted bigWigs from individual BAMs (one BAM per array task)
#
# INPUT:   <TREATMENT>-<TARGET>-<DAY>-R<n>-Cxt.bam files from Step 2.3.1
#            (02.3.1_Cxt_BamFiltered), e.g. ZH-AC-D1-R1-Cxt.bam
#          treatment_bam_list.txt = full paths of the treatment BAMs, one per line
#            (input controls, names containing "-IN-", are NOT in this list)
# OUTPUT:  in 00_BigWig_Files_Individual_Bams, per sample:
#            <sample>_inputsubtracted.bw   (treatment minus input, BPM-normalized)
#          These feed the matrix and heatmap steps (04.1 onward).
# ENV:     deeptools_kernel_v2 (deepTools bamCompare)
# RESOURCES per task: 8 CPUs, 32 GB, 5 h
#
# DIRECTORIES
#   BAM_INPUT_DIR = <path>/02_Cxt_Alignment/02.3.1_Cxt_BamFiltered
#                   Input. Must contain the sorted, indexed BAMs (.bam + .bai),
#                   including the input controls (-IN-).
#   BW_OUTPUT_DIR = <path>/04_Deeptools/00_BigWig_Files_Individual_Bams
#                   Output. Created automatically.
#   TREATMENT_LIST = <path>/04_Deeptools/treatment_bam_list.txt
#                   Must be built before submitting (see below).
#   Job logs go to logs/ in the folder you submit from.
#
# BEFORE SUBMITTING (run from 04_Deeptools)
#   1. mkdir -p logs
#        (SLURM opens the log files before the script starts, so this folder
#         must exist ahead of time or the job fails with no log.)
#   2. Build the treatment list (all non-input BAMs, full paths). See also
#      How_to_make_Treatment_list.md:
#        ls <BAM_INPUT_DIR>/*-Cxt.bam | grep -v -- '-IN-' > treatment_bam_list.txt
#        wc -l treatment_bam_list.txt      # <- this number is your array size
#   3. Submit with the array size taken from that file:
#        sbatch --array=1-$(wc -l < treatment_bam_list.txt) 00_Create_Bigwigs_Individual_Bams_ARRAY.sh
#
# CHANGING THE ARRAY SIZE TO MATCH YOUR SAMPLE COUNT
#   The "#SBATCH --array=1-72" line must equal the line count of
#   treatment_bam_list.txt (task N processes line N, counting from 1). 72 =
#   the 96 BAMs minus the input controls.
#     - Array smaller than the list: extra samples are silently skipped.
#     - Array larger than the list: extra tasks fail with
#       "No treatment BAM found for array index N".
#   Either edit that line, or use the wc -l form above so you never have to.
#   Put --array BEFORE the script name. Anything after the script name is passed
#   to the script, not to sbatch, and is ignored.
#   Other useful forms:
#        sbatch --array=1-72%10 00_Create_Bigwigs_Individual_Bams_ARRAY.sh   # at most 10 tasks at once
#        sbatch --array=5,17 00_Create_Bigwigs_Individual_Bams_ARRAY.sh      # rerun only tasks 5 and 17
#
# Input = all input replicates of this site and day, pooled (merged input BAM from Step 2.4)
day=$(echo "$basename" | cut -d- -f3)
input_bam=$(ls "${MERGED_DIR}/${treatment}-IN-${day}-"rep*-merged-sorted.bam 2>/dev/null | head -1)
