#!/usr/bin/env bash
#SBATCH --partition=acpu
#SBATCH --job-name=gs_tss_signal_check
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --time=01:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

################################################################################
# ONE-OFF CHECK (not a pipeline step; nothing in the pipeline is changed)
#
# Question: do the TSS windows that overlap Greenscreen regions have extreme
# signal in the figures' matrix (like gene5911 did), or are they ordinary?
#
# INPUT:   04_Deeptools/matrix_combined_gene_only_bed_file.gz (Step 4.1.2, 11/23)
#          Cxt_reference/CtarK1_TSS_excl_gene5911.bed (14,712 TSS)
#          03.0_Greenscreen/Cxt_greenscreen_12of24_10000bp.bed (468 regions)
# OUTPUT:  QC/oneoff_checks/greenscreen_tss_overlap/   (NEW folder; refuses to
#          run if it exists)
#            overlap_genes.txt        the overlapping gene names
#            row_stats.tsv            per gene: mean and max over all 72 samples
#                                     and bins, percentile rank, overlap flag
#            summary.txt              counts, percentiles, effect on the mean
# NOTE:    Run "mkdir -p logs" in the folder you submit from.
################################################################################

set -eo pipefail

RUN=/scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates
MATRIX=${RUN}/04_Deeptools/matrix_combined_gene_only_bed_file.gz
TSS=/projects/c832500103@colostate.edu/Cxt_reference/CtarK1_TSS_excl_gene5911.bed
GSB=${RUN}/03_MACS2_Peaks_New_Peaks/03.0_Greenscreen/Cxt_greenscreen_12of24_10000bp.bed
BT=/projects/c832500103@colostate.edu/conda_envs/homer_env/bin/bedtools
OUT=${RUN}/QC/oneoff_checks/greenscreen_tss_overlap

for f in "$MATRIX" "$TSS" "$GSB" "$BT"; do
    [[ -e "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done
if [[ -e "$OUT" ]]; then
    echo "ERROR: $OUT already exists; rename it first. Nothing was overwritten." >&2
    exit 1
fi
mkdir -p "$OUT"

# 1. Genes whose TSS +/-2 kb window overlaps a Greenscreen region
awk 'BEGIN{OFS="\t"}{s=$2-2000; if(s<0)s=0; print $1,s,$3+2000,$4}' "$TSS" \
    | "$BT" intersect -u -a stdin -b "$GSB" | cut -f4 | sort -u > "$OUT/overlap_genes.txt"
n_ov=$(wc -l < "$OUT/overlap_genes.txt")
echo "Overlapping genes: $n_ov"

# 2. Per-row mean and max over all samples and bins (nan skipped).
#    Matrix rows: chrom start end name score strand value1 ... valueN; header starts with @
zcat "$MATRIX" | awk -v OFS="\t" '
    NR==FNR {ov[$1]=1; next}
    /^@/ {next}
    {
        s=0; n=0; mx="";
        for (i=7; i<=NF; i++) {
            if ($i=="nan" || $i=="NaN") continue;
            v=$i+0; s+=v; n++;
            if (mx=="" || v>mx) mx=v;
        }
        m = (n>0) ? s/n : "NA";
        print $4, m, (mx=="" ? "NA" : mx), ($4 in ov ? 1 : 0)
    }' "$OUT/overlap_genes.txt" - > "$OUT/row_stats.unsorted.tsv"

# 3. Percentile rank by row mean (100 = highest)
{
    echo -e "gene\trow_mean\trow_max\tgreenscreen\tpct_rank_mean"
    grep -v $'\tNA\t' "$OUT/row_stats.unsorted.tsv" | sort -t$'\t' -k2,2g \
        | awk -v OFS="\t" '{a[NR]=$0} END{for(i=1;i<=NR;i++) print a[i], sprintf("%.2f",100*i/NR)}'
} > "$OUT/row_stats.tsv"
rm "$OUT/row_stats.unsorted.tsv"

# 4. Summary
awk -F'\t' -v nov=$n_ov '
    NR==1 {next}
    {
        all_s+=$2; all_n++
        if ($4==1) { ov_n++; ov_s+=$2; pr[ov_n]=$5; if ($5>=99) top1++; if ($5>=95) top5++ }
        else { rest_s+=$2; rest_n++ }
    }
    END{
        n=ov_n   # rows are sorted by mean, so pr[] is already in ascending order
        printf "Matrix rows with data: %d\n", all_n
        printf "Greenscreen-overlap genes in the matrix: %d (of %d listed)\n", ov_n, nov
        printf "Their percentile rank by row mean: min %.1f, median %.1f, max %.1f\n", pr[1], pr[int((n+1)/2)], pr[n]
        printf "Of them, in the top 5%%: %d   in the top 1%%: %d   (expected by chance: %.1f and %.1f)\n", top5, top1, ov_n*0.05, ov_n*0.01
        printf "Mean of row means, all genes: %.5f\n", all_s/all_n
        printf "Mean of row means, without the overlap genes: %.5f\n", rest_s/rest_n
        printf "Change from removing them: %.2f%%\n", 100*((rest_s/rest_n)-(all_s/all_n))/(all_s/all_n)
    }' "$OUT/row_stats.tsv" | tee "$OUT/summary.txt"

echo ""
echo "Top 10 overlap genes by row mean:"
awk -F'\t' 'NR>1 && $4==1' "$OUT/row_stats.tsv" | sort -t$'\t' -k2,2gr | head -10 | tee -a "$OUT/summary.txt"
echo "Output: $OUT"
