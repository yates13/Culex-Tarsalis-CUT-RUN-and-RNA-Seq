#!/usr/bin/env bash
################################################################################
# STEP 7.0: Make the Gviz sample sheet for 07.1_Gviz_Culex_Alpine.R
#           (run ONCE, with bash, on a login node; takes seconds)
#
# WHY:     The old CtarK1_Meta_File.csv (NAS) pointed at pre-rerun BAMs/peaks.
#          This builds a new sheet from the files that actually exist after the
#          2026-11 rerun, so nothing is typed by hand.
#
# INPUT:   filtered BAMs   02_Cxt_Alignment/02.3.1_Cxt_BamFiltered/<sample>-Cxt.bam (+ .bam.bai)
#          pooled peaks    03_MACS2_Peaks_New_Peaks/03.1_MacsPeaks_with_replicate_count/
#                          <SITE>-<MARK>-<DAY>-rep*/<...>_peaks.narrowPeak (AC, SR) / .broadPeak (ME)
#          bigWigs         04_Deeptools/00_BigWig_Files_Individual_Bams/<sample>_inputsubtracted.bw
#          REP_MAP         04.2_All_Marks_By_Day_Profile_HardCoded_Yaxis.sh (representative replicates)
# OUTPUT:  CUT-RUN/07_Gviz/CtarK1_Gviz_Meta_Alpine.csv   (NEW; refuses to overwrite)
#          one row per ChIP library (expected 72; inputs are not listed)
#          columns: Population,Factor,Day,Replicate,Sample,Peaks,bamReads,bigWig,Representative
#          (Population/Factor/Day/Peaks/bamReads = same columns as the old CtarK1_Meta_File.csv)
#
# USAGE:   bash 07.0_Make_Gviz_Meta.sh
################################################################################

set -eo pipefail

RUN=/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates        # [EDIT]
REPO=/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq   # [EDIT]

BAM_DIR=${RUN}/02_Cxt_Alignment/02.3.1_Cxt_BamFiltered
PEAK_DIR=${RUN}/03_MACS2_Peaks_New_Peaks/03.1_MacsPeaks_with_replicate_count
BW_DIR=${RUN}/04_Deeptools/00_BigWig_Files_Individual_Bams
REP_SCRIPT=${REPO}/CUT-RUN/04_Deeptools/04.2_Representative_Replicates/04.2_All_Marks_By_Day_Profile_HardCoded_Yaxis.sh
OUT=${REPO}/CUT-RUN/07_Gviz/CtarK1_Gviz_Meta_Alpine.csv
EXPECTED_ROWS=72
EXPECTED_REPS=24

[[ -e "$OUT" ]] && { echo "ERROR: $OUT already exists; nothing was overwritten." >&2; exit 1; }
for d in "$BAM_DIR" "$PEAK_DIR" "$BW_DIR"; do
    [[ -d "$d" ]] || { echo "ERROR: folder not found: $d" >&2; exit 1; }
done
[[ -s "$REP_SCRIPT" ]] || { echo "ERROR: REP_MAP script not found: $REP_SCRIPT" >&2; exit 1; }

# Representative replicate per SITE_MARK_DAY, read from the 04.2 REP_MAP (single source of truth)
declare -A REP
while read -r key rep; do
    REP[$key]=$rep
done < <(grep -oE '\[[A-Z]{2}_[A-Z]{2}_D[13]\]="R[0-9]+"' "$REP_SCRIPT" | sed -E 's/^\[([A-Z0-9_]+)\]="(R[0-9]+)"$/\1 \2/')
[[ ${#REP[@]} -eq $EXPECTED_REPS ]] || { echo "ERROR: expected $EXPECTED_REPS REP_MAP entries, found ${#REP[@]}" >&2; exit 1; }

tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT
echo "Population,Factor,Day,Replicate,Sample,Peaks,bamReads,bigWig,Representative" > "$tmp"

n=0; nrep=0
for s in BF MP SF ZH; do
  for m in AC ME SR; do
    for d in D1 D3; do
      ext=narrowPeak; [[ $m == ME ]] && ext=broadPeak       # ME = broad peaks (H3K9me3)
      pk=( "${PEAK_DIR}/${s}-${m}-${d}-rep"*/"${s}-${m}-${d}-rep"*"_peaks.${ext}" )
      if [[ ${#pk[@]} -ne 1 || ! -s "${pk[0]}" ]]; then
          echo "ERROR: expected exactly one peak file for ${s}-${m}-${d}: ${pk[*]}" >&2; exit 1
      fi
      for bam in "${BAM_DIR}/${s}-${m}-${d}-R"*"-Cxt.bam"; do
        samp=$(basename "$bam" -Cxt.bam)                    # e.g. ZH-SR-D1-R1
        r=${samp##*-}                                       # e.g. R1
        bw=${BW_DIR}/${samp}_inputsubtracted.bw
        for f in "$bam" "${bam}.bai" "$bw"; do
            [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
        done
        isrep=FALSE
        if [[ "${REP[${s}_${m}_${d}]}" == "$r" ]]; then isrep=TRUE; nrep=$((nrep + 1)); fi
        echo "${s},${m},${d},${r},${samp},${pk[0]},${bam},${bw},${isrep}" >> "$tmp"
        n=$((n + 1))
      done
    done
  done
done

if [[ $n -ne $EXPECTED_ROWS || $nrep -ne $EXPECTED_REPS ]]; then
    echo "ERROR: got $n rows ($nrep representative); expected $EXPECTED_ROWS ($EXPECTED_REPS). Nothing written." >&2
    exit 1
fi

cp "$tmp" "$OUT"
echo "Wrote $OUT: $n libraries, $nrep representative replicates"
head -3 "$OUT"
