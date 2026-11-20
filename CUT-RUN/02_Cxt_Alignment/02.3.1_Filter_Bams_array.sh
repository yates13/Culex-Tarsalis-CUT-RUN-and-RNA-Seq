#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=BamFilter_array
#SBATCH --output=logs/%x.%A_%a.log
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shawn.yates@colostate.edu
#SBATCH --array=1-96             # MUST equal the number of lines in bam_file_list.txt

################################################################################
# STEP 2.3.1: Keep only uniquely placed, properly paired primary reads
#
# WHY: HISAT2 reports up to 5 alignments per read. The unfiltered BAMs hold
#   secondary alignments (about half as many as primary reads in BF-ME-D1-R1)
#   and multi-mapped primary reads with MAPQ < 10 (about 18%). bamCompare counts
#   secondary alignments, so reads in repeats were counted up to 5 times in the
#   bigWigs, and MACS2 counted the low-MAPQ reads. Every step after this one
#   reads the filtered BAMs, so peaks and bigWigs are built from the same reads.
#
# FILTER:  samtools view -q 10 -f 2 -F 2308
#            -q 10   MAPQ >= 10 (drops reads that map equally well elsewhere)
#            -f 2    properly paired
#            -F 2308 drop unmapped (4), secondary (256), supplementary (2048)
#          The input is coordinate-sorted, so the output stays sorted.
#
# INPUT:   <sample>.bam from Step 2.3 (02.3_Cxt_BamConverted)
#          bam_file_list.txt in 02.3.1_Cxt_BamFiltered = full paths, one per line
# OUTPUT:  in 02.3.1_Cxt_BamFiltered, per sample (same file names as Step 2.3):
#            <sample>.bam, <sample>.bam.bai, <sample>_flagstat.txt
#            <sample>_filter_counts.txt   primary reads before, kept, % kept
# ENV:     hisat2_env_fixed (samtools)
# QC:      02.3.2_QC_filter.sh
#
# BEFORE SUBMITTING (run from 02_Cxt_Alignment):
#   mkdir -p logs 02.3.1_Cxt_BamFiltered   (pipeline folder, not the repo)
#   ls <path>/02.3_Cxt_BamConverted/*-Cxt.bam > <path>/02.3.1_Cxt_BamFiltered/bam_file_list.txt
#   sbatch --array=1-$(wc -l < <path>/02.3.1_Cxt_BamFiltered/bam_file_list.txt) 02.3.1_Filter_Bams_array.sh
#
# RERUNNING: a sample is skipped if its filtered .bam.bai exists and is newer
#   than the unfiltered BAM. Delete the .bai to force a redo.
################################################################################

module purge
SCRATCH=/scratch/alpine/.colostate.edu/c832500103
PROJ=/projects/.colostate.edu/c832500103   # no "@" (see GOTCHAS.md)
source ${PROJ}/miniforge3/etc/profile.d/conda.sh || { echo "ERROR: conda.sh not found" >&2; exit 1; }
conda activate ${PROJ}/conda_envs/hisat2_env_fixed || { echo "ERROR: conda activate hisat2_env_fixed failed" >&2; exit 1; }
set -eo pipefail
nthreads=$SLURM_CPUS_PER_TASK

ALN=${SCRATCH}/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/02_Cxt_Alignment
indir=${ALN}/02.3_Cxt_BamConverted
outdir=${ALN}/02.3.1_Cxt_BamFiltered
BAM_LIST=${outdir}/bam_file_list.txt
MAPQ=10
mkdir -p "$outdir"

[[ -f "$BAM_LIST" ]] || { echo "ERROR: $BAM_LIST not found. Build it: ls ${indir}/*-Cxt.bam > ${BAM_LIST}" >&2; exit 1; }
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_LIST")
[[ -n "$file" ]] || { echo "ERROR: no BAM for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
base=$(basename "$file" .bam)
out=${outdir}/${base}.bam
echo "Task $SLURM_ARRAY_TASK_ID: $base"; date

if [[ -f "${out}.bai" && "${out}.bai" -nt "$file" ]]; then
    echo "EXISTS: $base, skipping"; exit 0
fi

samtools view -b -@$nthreads -q $MAPQ -f 2 -F 2308 -o "$out" "$file"
samtools index -@$nthreads "$out"
samtools flagstat -@$nthreads "$out" > ${outdir}/${base}_flagstat.txt

before=$(samtools view -c -@$nthreads -F 2308 "$file")
kept=$(samtools view -c -@$nthreads "$out")
pct=$(awk -v k=$kept -v b=$before 'BEGIN{printf "%.2f", b ? 100*k/b : 0}')
printf 'sample\tprimary_before\tkept\tpct_kept\n%s\t%s\t%s\t%s\n' "$base" "$before" "$kept" "$pct" > ${outdir}/${base}_filter_counts.txt
echo "Kept $kept of $before primary reads ($pct%)"
echo "Completed: $base"; date
