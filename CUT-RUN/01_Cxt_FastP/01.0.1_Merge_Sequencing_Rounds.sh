#!/usr/bin/env bash
################################################################################
# STEP 1.0.1: Merge the sequencing rounds into one flat folder of raw reads
#
# Raw reads arrived in three deliveries, each nested differently:
#   First_Round_Sequencing/CutRun_DNA_Rd1_30-1326595428/00_fastq/*.fastq.gz
#   Second_Round_Sequencing/CutRun_DNA_Rd2_30-1352505226/*.fastq.gz
#   Third_Round_Sequencing/CutRun_MP12_30-1360459436/...
# Some samples were sequenced again in a later round. The LATEST round wins:
# a sample in rounds 2 and 3 uses only the round 3 copy. A sample is taken
# as a unit (R1 + R2 + both .md5 files from the same round), so mates are
# never mixed across rounds.
#
# OUTPUT:  RAW_DIR/00_Merged_Raw with symlinks, per sample, to
#            <sample>_R1_001.fastq.gz, <sample>_R2_001.fastq.gz and their .md5
#          00_Merged_Raw/merge_manifest.tsv
#            sample  round_used  superseded_rounds  r1_md5  r2_md5
# Symlinks leave the round folders untouched and use no extra space.
#
# USAGE (login node, takes seconds, no sbatch):
#   bash 01.0.1_Merge_Sequencing_Rounds.sh
# Rerunning is safe: old symlinks are replaced.
# Exits 1 if any sample is missing a mate or a .md5 (still links the rest).
# Next step: 01.0.2_Md5sum_Check.sh
################################################################################

set -euo pipefail

RAW_DIR="${RAW_DIR:-/scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Raw_Files}"
MERGED_DIR="${MERGED_DIR:-${RAW_DIR}/00_Merged_Raw}"
ROUNDS=(First_Round_Sequencing Second_Round_Sequencing Third_Round_Sequencing)  # oldest -> newest

cd "$RAW_DIR"
for r in "${ROUNDS[@]}"; do
    [[ -d "$r" ]] || { echo "ERROR: ${RAW_DIR}/${r} not found"; exit 1; }
done

mkdir -p "$MERGED_DIR"
if find "$MERGED_DIR" -mindepth 1 -maxdepth 1 -name '*.fastq.gz*' ! -type l | grep -q .; then
    echo "ERROR: $MERGED_DIR contains real files, not symlinks. Not touching it."
    exit 1
fi
find "$MERGED_DIR" -mindepth 1 -maxdepth 1 -type l -delete

declare -A SRC_DIR SRC_ROUND OLDER
problems=0

# Later rounds overwrite earlier ones, so the newest copy is kept.
for r in "${ROUNDS[@]}"; do
    declare -A seen_this_round=()
    while IFS= read -r f; do
        s=$(basename "$f"); s=${s%_R1_001.fastq.gz}
        if [[ -n "${seen_this_round[$s]:-}" ]]; then
            echo "PROBLEM: $s appears twice inside $r ($f and ${seen_this_round[$s]})"
            problems=$((problems + 1)); continue
        fi
        seen_this_round[$s]=$f
        if [[ -n "${SRC_ROUND[$s]:-}" ]]; then
            OLDER[$s]="${OLDER[$s]:+${OLDER[$s]},}${SRC_ROUND[$s]}"
        fi
        SRC_DIR[$s]=$(dirname "$f")
        SRC_ROUND[$s]=$r
    done < <(find "$r" -name '*_R1_001.fastq.gz' \( -type f -o -type l \) | sort)

    # An R2 with no R1 in the same folder would be silently dropped; flag it.
    while IFS= read -r f; do
        [[ -e "${f/_R2_001.fastq.gz/_R1_001.fastq.gz}" ]] || {
            echo "PROBLEM: R2 without R1: $f"; problems=$((problems + 1)); }
    done < <(find "$r" -name '*_R2_001.fastq.gz' \( -type f -o -type l \) | sort)
    unset seen_this_round
done

MANIFEST="${MERGED_DIR}/merge_manifest.tsv"
printf 'sample\tround_used\tsuperseded_rounds\tr1_md5\tr2_md5\n' > "$MANIFEST"

linked=0
for s in $(printf '%s\n' "${!SRC_ROUND[@]}" | sort); do
    d="${RAW_DIR}/${SRC_DIR[$s]}"
    r1="${d}/${s}_R1_001.fastq.gz"
    r2="${d}/${s}_R2_001.fastq.gz"
    if [[ ! -e "$r2" ]]; then
        echo "PROBLEM: no R2 for $s in ${SRC_DIR[$s]} - sample NOT linked"
        problems=$((problems + 1)); continue
    fi
    ln -s "$r1" "${MERGED_DIR}/"
    ln -s "$r2" "${MERGED_DIR}/"
    md5_status=()
    for fq in "$r1" "$r2"; do
        if [[ -f "${fq}.md5" ]]; then
            ln -s "${fq}.md5" "${MERGED_DIR}/"; md5_status+=(yes)
        else
            echo "PROBLEM: missing checksum ${SRC_DIR[$s]}/$(basename "$fq").md5"
            md5_status+=(MISSING); problems=$((problems + 1))
        fi
    done
    printf '%s\t%s\t%s\t%s\t%s\n' "$s" "${SRC_ROUND[$s]}" "${OLDER[$s]:--}" \
        "${md5_status[0]}" "${md5_status[1]}" >> "$MANIFEST"
    linked=$((linked + 1))
done

echo
echo "== Samples linked into $MERGED_DIR: $linked ($((linked * 2)) FASTQs)"
echo "== Samples per round used:"
tail -n +2 "$MANIFEST" | cut -f2 | sort | uniq -c
echo "== Samples taken from a later round (older copy ignored):"
awk -F'\t' 'NR>1 && $3!="-" {print "  "$1"  used "$2", ignored "$3}' "$MANIFEST"
echo "== Problems: $problems"
(( problems == 0 )) || exit 1
