#!/usr/bin/env bash
################################################################################
# STEP 3.4: Peak reproducibility across replicates (report only)
#
# PURPOSE: For each of the 24 pooled peak sets (Step 3.1), count how many pooled
#          peaks overlap a peak in at least 1, at least 2, and all 3 replicates
#          (per-replicate peaks from Step 3.3). A pooled peak counts as
#          reproducible if it overlaps peaks in >= 2 of 3 replicates.
#
# USAGE:   on a compute node (acompile), after all 72 Step 3.3 tasks finish:
#            conda activate /projects/c832500103@colostate.edu/conda_envs/homer_env
#            bash 03.4_Peak_Reproducibility.sh
#
# INPUT:   03.1_MacsPeaks_with_replicate_count/<SITE>-<MARK>-<DAY>-rep*/..._peaks.<narrow|broad>Peak
#          03.3_PerReplicate_Peaks_2026-11-27/<SITE>-<MARK>-<DAY>-R<n>/..._peaks.<narrow|broad>Peak
# OUTPUT:  QC/peak_reproducibility_2026-11-27/peak_reproducibility.tsv
#          (refuses to overwrite)
# ENV:     homer_env (bedtools 2.31.1)
################################################################################

set -eo pipefail

# ---------------- USER SETTINGS ----------------
RUN=/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates   # [EDIT]
POOLED=${RUN}/03_MACS2_Peaks_New_Peaks/03.1_MacsPeaks_with_replicate_count                       # [EDIT]
PERREP=${RUN}/03_MACS2_Peaks_New_Peaks/03.3_PerReplicate_Peaks_2026-11-27                        # [EDIT]
OUTDIR=${RUN}/QC/peak_reproducibility_2026-11-27                                                 # [EDIT]
# ------------------------------------------------

command -v bedtools >/dev/null || { echo "ERROR: bedtools not found; activate homer_env" >&2; exit 1; }
OUT=${OUTDIR}/peak_reproducibility.tsv
if [[ -e ${OUT} ]]; then echo "ERROR: ${OUT} exists; refusing to overwrite" >&2; exit 1; fi
mkdir -p ${OUTDIR}

printf "condition\tpooled_peaks\treplicate_peaks\tin_ge1\tin_ge2\tpct_ge2\tin_all\n" > ${OUT}

for D in ${POOLED}/*-rep*/; do
    COND=$(basename "$D")                       # e.g. BF-AC-D1-rep123
    SMD=$(echo "$COND" | cut -d- -f1-3)         # e.g. BF-AC-D1
    MARK=$(echo "$COND" | cut -d- -f2)
    if [[ "$MARK" == "ME" ]]; then EXT=broadPeak; else EXT=narrowPeak; fi

    P=$(ls ${D}*_peaks.${EXT})
    mapfile -t REPS < <(ls ${PERREP}/${SMD}-R*/*_peaks.${EXT} 2>/dev/null)
    if [[ ${#REPS[@]} -ne 3 ]]; then
        echo "ERROR: ${SMD}: found ${#REPS[@]} replicate peak files, expected 3" >&2; exit 1
    fi

    NP=$(wc -l < "$P")
    RN=$(for r in "${REPS[@]}"; do wc -l < "$r"; done | paste -sd, -)

    # For each replicate, the pooled peaks it overlaps (each pooled peak once per replicate)
    COUNTS=$(for r in "${REPS[@]}"; do bedtools intersect -u -a "$P" -b "$r" | cut -f1-3; done \
             | sort | uniq -c | awk '{print $1}')
    GE1=$(echo "$COUNTS" | grep -c . || true)
    GE2=$(echo "$COUNTS" | awk '$1>=2' | wc -l)
    ALL=$(echo "$COUNTS" | awk '$1>=3' | wc -l)
    PCT=$(awk -v a=$GE2 -v b=$NP 'BEGIN{printf "%.1f", 100*a/b}')

    printf "%s\t%d\t%s\t%d\t%d\t%s\t%d\n" "$COND" "$NP" "$RN" "$GE1" "$GE2" "$PCT" "$ALL" >> ${OUT}
done

column -t ${OUT}
echo "Written: ${OUT}"

