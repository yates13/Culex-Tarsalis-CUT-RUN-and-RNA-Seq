#!/usr/bin/env bash
################################################################################
# STEP 6.0: Make a gene-only GTF of the 14,713 MAKER genes for the R scripts
#           (run ONCE, with bash, on a login node; takes seconds)
#
# WHY:     getGeneDistances (Step 6.1) reads GTF attributes (gene_id "gene123";
#          Name "..."; note "...";). The MAKER genes come as GFF3 (ID=gene123;Name=...).
#          This writes the SAME 14,713 genes as the TSS figures (Step 4) in GTF form,
#          so every analysis uses one gene set.
#
# INPUT:   the original OSF annotation (Main et al. 2021), md5-checked:
#            /projects/.../Cxt_reference/Culex-tarsalis_knwr_BASEFEATURES_CtarK1.gff3.gz
#          filter: column 2 == "maker" and column 3 == "gene" (Mt genes excluded,
#          same filter as genes_only.gff3 and the TSS BED)
# OUTPUT:  /projects/.../Cxt_reference/CtarK1_MAKER_genes.gtf   (NEW; refuses to overwrite)
#          + a copy in Cxt_Cut_Run_Pipeline/00_Cxt_Genome/Genome_gff/
#          GTF column 9: gene_id "<ID>"; Name "<Name>"; note "<Note>";  (Name/note only if present)
#
# USAGE:   bash 06.0_Make_MAKER_gene_GTF.sh
################################################################################

set -eo pipefail

REF=/projects/c832500103@colostate.edu/Cxt_reference
SRC=${REF}/Culex-tarsalis_knwr_BASEFEATURES_CtarK1.gff3.gz
MD5=8a252032e154660d3109c2db815b0ebb
OUT=${REF}/CtarK1_MAKER_genes.gtf
COPY_DIR=/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/00_Cxt_Genome/Genome_gff
EXPECTED=14713

[[ -s "$SRC" ]] || { echo "ERROR: annotation not found: $SRC" >&2; exit 1; }
[[ -e "$OUT" ]] && { echo "ERROR: $OUT already exists; nothing was overwritten." >&2; exit 1; }

got=$(md5sum "$SRC" | cut -d' ' -f1)
[[ "$got" == "$MD5" ]] || { echo "ERROR: md5 mismatch for $SRC ($got, expected $MD5)" >&2; exit 1; }
echo "md5 OK: $SRC"

zcat "$SRC" | awk -F'\t' -v OFS='\t' '
  $2 == "maker" && $3 == "gene" && $1 != "Mt" {
    id = ""; nm = ""; nt = ""
    n = split($9, kv, ";")
    for (i = 1; i <= n; i++) {
      p = index(kv[i], "="); if (!p) continue
      k = substr(kv[i], 1, p - 1); v = substr(kv[i], p + 1)
      gsub(/"/, "\047", v)                       # no double quotes inside GTF values
      if (k == "ID") id = v; else if (k == "Name") nm = v; else if (k == "Note") nt = v
    }
    if (id == "") next
    a = "gene_id \"" id "\";"
    if (nm != "") a = a " Name \"" nm "\";"
    if (nt != "") a = a " note \"" nt "\";"
    print $1, $2, $3, $4, $5, $6, $7, $8, a
  }' > "$OUT"

n=$(wc -l < "$OUT")
nid=$(grep -c -E 'gene_id "gene[0-9]+";' "$OUT" || true)
echo "Genes written: $n   (gene_id like gene123: $nid)"
if [[ "$n" -ne "$EXPECTED" || "$nid" -ne "$EXPECTED" ]]; then
    echo "ERROR: expected $EXPECTED MAKER genes with gene<N> IDs; removing $OUT" >&2
    rm -f "$OUT"; exit 1
fi

cp -n "$OUT" "$COPY_DIR/"
echo "Wrote $OUT (copy in $COPY_DIR)"
head -2 "$OUT"
