################################################################################
# STEP 7.1 (ALPINE): Gviz genome-track figures for DiffBind example genes
#
# Alpine version of Gvis_Analysis_SY.R (NAS). Same functions and the same flow
# (DiffBind results -> genes near significant sites -> pick genes -> plot).
# Lines changed from the original end in  # [CHANGED] <reason>.
# Lines you may edit end in  # [EDIT].  Lines that match the paper end in  # [FIXED].
#
# USAGE:   sbatch 07.1_Gviz_Culex.sh                            (USER SETTINGS below)
#          sbatch 07.1_Gviz_Culex.sh ZH BF SR D1 0.1 up 5        (override on the command line)
#                 args: targetPop refPop mark day FDR direction topN
#          or, in gviz_env:  Rscript 07.1_Gviz_Culex_Alpine.R [same args]
#
# INPUT:   DiffBind results  05_DiffBind_2026-11-24/<d1|d3>/<TARGET>_<MARK>_<day>_vs_<REF>-OUT.csv
#                            (all sites, th = 1; Fold > 0 = higher in TARGET)
#          sample sheet      CUT-RUN/07_Gviz/CtarK1_Gviz_Meta_Alpine.csv (made by 07.0)
#          structural GFF3   Cxt_reference/Culex-tarsalis_knwr_BASEFEATURES_CtarK1.gff3.gz
#                            (MAKER genes, ID=gene<N>; same IDs/coords as CtarK1_MAKER_genes.gtf)
#          notes GTF         Cxt_reference/basefeatures_updated_high_confidence.with_notes.gtf
#                            (OPTIONAL: without it, titles say "function unknown")
# OUTPUT:  <run>/07_Gviz/<TARGET>_<MARK>_<day>_vs_<REF>_FDR_<fdr>_<direction>_<timestamp>/   (NEW folder)
#            NN_<gene>.pdf / .png   one genome-track figure per candidate gene (topN)
#            candidates.csv         all candidate genes, ranked (FDR, then |Fold|)
#            run_settings.txt, sessionInfo.txt
#          <run>/07_Gviz/reference/<structural>_noMt.gff3   (made once, then reused)
#          Nothing is ever overwritten.
# ENV:     gviz_env (Gviz 1.54, txdbmaker 1.6, R 4.5.3; envs/locked/gviz_env.full.yml)
#
# TRACKS:  genome axis, MAKER gene bodies, MAKER transcripts, input-subtracted
#          coverage (BPM, 50 bp) of the REPRESENTATIVE replicate of each site
#          (04.2 REP_MAP) on a SHARED y-axis, pooled MACS2 peaks of each site.
# GENES:   a gene is a candidate when a significant DiffBind site (FDR <= cutoff,
#          chosen direction) overlaps its TSS +/- 2 kb (as Fig. 3B) AND that site
#          lies inside a pooled peak of the site where it is higher.
#
# CAVEAT (2026-11-27): significant sites are mostly peaks called in only one group
#          (peak-origin bias, Lun & Smyth 2014); a group-blind peak set is planned.
#          Figures from these DiffBind results are for checking, not final.
################################################################################

suppressPackageStartupMessages({
  library(Gviz)
  library(GenomicRanges)
  library(rtracklayer)
  library(IRanges)
  library(stringr)
  library(biomartr)          # read_gff()
  library(GenomicFeatures)
  library(txdbmaker)         # [CHANGED] makeTxDbFromGFF() moved here from GenomicFeatures
  library(AnnotationDbi)
})
# [CHANGED] chipseq, biomaRt, urltools no longer loaded: nothing in this script uses them

# Must be run for non-model organisms
options(ucscChromosomeNames = FALSE)   # [FIXED]

# =============================================================================
# USER SETTINGS
# =============================================================================
RUN  <- "/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates"        # [EDIT]
REPO <- "/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq"   # [EDIT]
REF  <- "/projects/c832500103@colostate.edu/Cxt_reference"                                                   # [EDIT]

diffbindDir   <- file.path(RUN, "05_DiffBind_2026-11-24")                                    # [EDIT] DiffBind results (d1/, d3/)
metaFile      <- file.path(REPO, "CUT-RUN/07_Gviz/CtarK1_Gviz_Meta_Alpine.csv")             # [EDIT] made by 07.0
structuralGZ  <- file.path(REF, "Culex-tarsalis_knwr_BASEFEATURES_CtarK1.gff3.gz")          # [EDIT]
notesGTF      <- file.path(REF, "basefeatures_updated_high_confidence.with_notes.gtf")      # [EDIT] optional
outBase       <- file.path(RUN, "07_Gviz")                                                   # [EDIT]

targetPop <- "ZH"    # [EDIT] first site in the DiffBind file name (BF, MP, ZH)
refPop    <- "BF"    # [EDIT] second site (BF, SF, MP)
markName  <- "SR"    # [EDIT] AC, ME or SR
dayName   <- "D1"    # [EDIT] D1 or D3
fdrCut    <- 0.1     # [EDIT] 0.1 = main results (paper: FDR <= 0.10); 0.3 = exploratory
direction <- "up"    # [EDIT] "up" = higher in targetPop (Fold > 0); "down" = higher in refPop (Fold < 0)
topN      <- 5       # [EDIT] how many candidate genes to render

promoterFlank <- 2000                                         # [FIXED] +/- 2 kb TSS window (Fig. 3B)
excludeGenes  <- c("gene5911", "gene2811", "gene12689")       # artifact TSS loci: never used as examples
trackHeights  <- c(1, 1, 2, 2, 0.6, 0.6)                      # [EDIT] genes, transcripts, coverage x2, peaks x2
trackColors   <- c("firebrick3", "grey30")                    # [EDIT] coverage: targetPop, refPop
mtContigName  <- "Mt"                                         # Mt contig name in the structural GFF3 (checked 2026-11-27)

# Command-line overrides (optional): targetPop refPop mark day FDR direction topN
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) { targetPop <- args[1]; refPop <- args[2]; markName <- args[3]; dayName <- args[4] }
if (length(args) >= 5) fdrCut    <- as.numeric(args[5])
if (length(args) >= 6) direction <- args[6]
if (length(args) >= 7) topN      <- as.integer(args[7])
if (!direction %in% c("up", "down") || is.na(fdrCut) || is.na(topN)) {
  stop("Bad settings: direction must be up/down, FDR and topN numeric")
}

# ---- output folder: NEW every run, never overwritten -------------------------
contrastID <- paste0(targetPop, "_", markName, "_", tolower(dayName), "_vs_", refPop)
outDir <- file.path(outBase, paste0(contrastID, "_FDR_", fdrCut, "_", direction, "_",
                                    format(Sys.time(), "%Y%m%d_%H%M%S")))
if (dir.exists(outDir)) stop("Output folder already exists; nothing overwritten: ", outDir)
dir.create(outDir, recursive = TRUE)
refDir <- file.path(outBase, "reference")
dir.create(refDir, recursive = TRUE, showWarnings = FALSE)
cat("Contrast:", contrastID, " FDR <=", fdrCut, " direction:", direction, " topN:", topN,
    "\nOutput:  ", outDir, "\n\n")

# ---------------------------------------------------------------------------
# TWO FILES, TWO JOBS
#
#   structuralGFF -> full gene/transcript/exon/CDS structure (original file).
#                    Used for makeTxDbFromGFF(), makeGViz(), plotAligns().
#
#   notesGTF      -> gene-only file with SwissProt functional annotations
#                    merged in (the "note" attribute). Used ONLY to look up a
#                    gene's description, matched by genomic COORDINATES.
# ---------------------------------------------------------------------------
structuralNoMt <- file.path(refDir, sub("\\.gff3\\.gz$", "_noMt.gff3", basename(structuralGZ)))   # [CHANGED] input is .gz on Alpine; stripped copy kept in 07_Gviz/reference

# Defined ONCE here; referenced everywhere else (makeGViz, importPeaks, etc.)
macsPeakCols <- c(signalValue = "numeric", pValue = "numeric",
                  qValue = "numeric", peak = "integer")

# [CHANGED] MACS2 broadPeak files (ME) have no 10th "peak" column; narrowPeak (AC, SR) do
importPeaks <- function(f) {
  cols <- if (grepl("broadPeak$", f)) macsPeakCols[c("signalValue", "pValue", "qValue")] else macsPeakCols
  rtracklayer::import(f, format = "BED", extraCols = cols)
}

# ---------------------------------------------------------------------------
# Strip the Mt contig from the structural file
# (Mt gene entries are missing proper ID attributes and break makeTxDbFromGFF)
# Made once; reused afterwards (never overwritten).
# ---------------------------------------------------------------------------
if (!file.exists(structuralNoMt)) {
  rawLines  <- readLines(gzfile(structuralGZ))                                               # [CHANGED] read the .gz directly
  keepLines <- grepl("^#", rawLines) | !grepl(paste0("^", mtContigName, "\t"), rawLines)
  writeLines(rawLines[keepLines], structuralNoMt)
  cat("Wrote", structuralNoMt, "(", sum(!keepLines), "Mt lines removed )\n")
}

# Quick sanity check - abort early with a clear message if repair failed
checkLines <- readLines(structuralNoMt, n = 200)
if (any(grepl("%3d|%3b|attribute=", checkLines, ignore.case = TRUE))) {
  stop("GFF repair incomplete - corrupted attribute syntax still present in ",
       structuralNoMt)
}

# Build the TxDb ONCE here and reuse everywhere
testtx <- txdbmaker::makeTxDbFromGFF(structuralNoMt, format = "gff3")                      # [CHANGED] txdbmaker::

# Gene-body GRanges, used for coordinate-based ID lookups and peak overlaps
structGFF   <- rtracklayer::import(structuralNoMt)                                         # [CHANGED] was the raw file (still had Mt)
structGenes <- structGFF[structGFF$type == "gene"]
cat("Structural genes (no Mt):", length(structGenes), "\n")


# ---------------------------------------------------------------------------
# makeGViz: build a list of Gviz tracks (gene bodies, transcripts, coverage, peaks)
# Uses the STRUCTURAL file only.
#
#   colorBySignal = TRUE  -> peak boxes shaded by MACS signalValue (QC/exploration)
#   colorBySignal = FALSE -> flat blue3 peak boxes (cleaner, for figures/posters)
#   showPeakIDs   = TRUE  -> draws each peak's BED "id" text inside the box
#   showPeakIDs   = FALSE -> no in-box text (default; avoids "unknown.####" labels)
#   reads         -> BAM or bigWig paths (bigWig = input-subtracted coverage)
#   txdb          -> reuse an existing TxDb instead of rebuilding it
# ---------------------------------------------------------------------------
makeGViz <- function(annotationFile = "", reads = c(), peaks = c(),
                     geneBodies = TRUE, transcripts = FALSE, chromosome = "",
                     orgName = "Culex tarsalis",
                     peakNames = c(), alignNames = c(),
                     colorBySignal = FALSE, showPeakIDs = FALSE,
                     txdb = NULL) {                                                         # [CHANGED] txdb argument

  peakLists <- c()

  if (geneBodies) {
    tempgff  <- read_gff(annotationFile)
    tempgff  <- tempgff[tempgff$seqid != mtContigName, ]
    geneRows <- tempgff[tempgff$type == "gene", ]

    geneIDS <- str_extract(geneRows$attribute, pattern = "gene[0-9]+")

    geneGR <- GRanges(seqnames = geneRows$seqid,
                      ranges = IRanges(start = geneRows$start, end = geneRows$end),
                      strand = geneRows$strand)

    geneBodyTrack <- AnnotationTrack(range = geneGR, group = geneIDS,
                                     fontcolor.item = "black", just.group = "above",
                                     name = "Culex tarsalis Genes", cex.title = 1)
    peakLists <- c(peakLists, geneBodyTrack)
  }

  if (transcripts) {
    tx <- if (is.null(txdb)) txdbmaker::makeTxDbFromGFF(annotationFile, format = "gff3") else txdb   # [CHANGED] reuse the TxDb built once
    transcriptTrack <- GeneRegionTrack(range = tx,
                                       genome = "Culex tarsalis",
                                       name = "Culex tarsalis Transcripts",
                                       cex.title = 1)
    peakLists <- c(peakLists, transcriptTrack)
  }

  if (length(reads) > 0) {
    if (length(alignNames) == 0) {
      alignNames <- paste0("Alignment_", seq_along(reads))
    }
    for (i in seq_along(reads)) {
      if (grepl("\\.bam$", reads[i])) {
        alTrack <- DataTrack(range = reads[i], type = "histogram",
                             name = alignNames[i], isPaired = TRUE, cex.title = 1)
      } else {                                                                              # [CHANGED] bigWig coverage (no isPaired)
        alTrack <- DataTrack(range = reads[i], type = "histogram",
                             name = alignNames[i], cex.title = 1)
      }
      peakLists <- c(peakLists, alTrack)
    }
  }

  if (length(peaks) > 0) {
    if (length(peakNames) == 0) {
      peakNames <- paste0("Peaks_", seq_along(peaks))
    }
    for (i in seq_along(peaks)) {
      thisPeakGR <- importPeaks(peaks[i])                                                   # [CHANGED] handles broadPeak (ME)

      # color peak boxes by signal strength; fall back to a flat color
      # if there's only one distinct value (cut() fails on a single value)
      if (colorBySignal && length(unique(thisPeakGR$signalValue)) > 1) {
        peakColors <- colorRampPalette(c("lightblue", "darkblue"))(100)[
          as.integer(cut(thisPeakGR$signalValue, 100))
        ]
      } else {
        peakColors <- "blue3"
      }

      tempPeak <- AnnotationTrack(thisPeakGR,
                                  name = peakNames[i],
                                  shape = 'box',
                                  fill = peakColors,
                                  size = 20, cex.title = 1,
                                  featureAnnotation = if (showPeakIDs) "id" else NULL)
      peakLists <- c(peakLists, tempPeak)
    }
  }

  return(peakLists)
}


# [CHANGED] shared y-axis: one ylim from all coverage bigWigs over the plotted region
sharedYlim <- function(bwFiles, chrom, from, to) {
  rng <- GRanges(chrom, IRanges(from, to))
  v <- unlist(lapply(bwFiles, function(f) {
    tryCatch(rtracklayer::import(f, format = "BigWig", which = rng)$score,
             error = function(e) numeric(0))
  }))
  if (length(v) == 0) return(c(0, 1))
  c(min(0, v), max(0.1, v))
}


# ---------------------------------------------------------------------------
# plotAligns: plot a set of tracks over a region, optionally centered on a gene
# Default gff is the Mt-stripped file and default txdb reuses the one built above.
#   bwFiles -> if given, all DataTracks get the same y-axis for this region
# ---------------------------------------------------------------------------
plotAligns <- function(tracks, chromosome, start, end, gene = "", heights = c(),
                       title = "", gff = structuralNoMt, txdb = testtx,
                       showAxis = TRUE, baseSize = 1,
                       background.title = "black",
                       autoTitle = FALSE, geneLookupTable = NULL,
                       titleCex = 1.1, bwFiles = NULL) {                                    # [CHANGED] bwFiles

  if (nchar(gene) > 0) {
    myannot <- AnnotationDbi::select(txdb, keys = c(gene),
                                     columns = c("TXCHROM", "TXSTART", "TXEND", "TXNAME"),
                                     keytype = "GENEID")
    if (nrow(myannot) == 0) {
      stop("Gene '", gene, "' not found in TxDb — check the ID matches structGenes$Name")
    }
    chromosome <- myannot$TXCHROM
    start      <- myannot$TXSTART
    end        <- myannot$TXEND

    if (autoTitle && !is.null(geneLookupTable)) {
      note  <- geneLookupTable$note[geneLookupTable$Name == gene]
      title <- paste0(gene, " - ", if (length(note) > 0) note[1] else "unknown")
    }
  }

  if (length(bwFiles) > 0) {                                                                # [CHANGED] shared y-axis
    yl <- sharedYlim(bwFiles, chromosome[1], min(start), max(end))
    tracks <- lapply(tracks, function(tr) {
      if (is(tr, "DataTrack")) displayPars(tr) <- list(ylim = yl)
      tr
    })
  }

  if (showAxis) {
    tracks <- c(GenomeAxisTrack(), tracks)
    if (length(heights) == 0) {
      heights <- c(0.5, rep(1, length(tracks) - 1))
    } else {
      heights <- c(0.5, heights)
    }
  } else if (length(heights) == 0) {
    heights <- rep(1, length(tracks))
  }


  plotTracks(tracks, sizes = heights, main = title,
             from = start, to = end, transcriptAnnotation = "symbol",
             groupAnnotation = "group", just.group = "above",
             chromosome = chromosome,
             title.width = 1.5,
             cex.main = titleCex,
             background.title = "grey75",
             col.title = "black")
}


# ---------------------------------------------------------------------------
# findGeneID: bridge the notes file's gene_id (e.g. "gene146") to the
# structural file's Name, by matching genomic coordinates.
# ---------------------------------------------------------------------------
findGeneID <- function(seqname, s, e) {
  m <- structGenes[as.character(seqnames(structGenes)) == as.character(seqname) &
                     start(structGenes) == s &
                     end(structGenes) == e]
  if (length(m) == 0) {
    warning("No structural gene found for ", seqname, ":", s, "-", e)
    return(NA_character_)
  }
  if (length(m) > 1) {
    warning("Multiple structural genes matched ", seqname, ":", s, "-", e,
            " - using first")
  }
  mcols(m)$Name[1]
}


# ---------------------------------------------------------------------------
# Gene descriptions from the NOTES file (gene-only, SwissProt notes)
# [CHANGED] optional: if the file is not on Alpine yet, titles say "function unknown"
# ---------------------------------------------------------------------------
if (file.exists(notesGTF)) {
  notesImported <- rtracklayer::import(notesGTF, format = "gtf")
  genesOnly     <- notesImported[notesImported$type == "gene"]
  geneLookup    <- as.data.frame(genesOnly)
  if (!"note" %in% names(geneLookup)) geneLookup$note <- NA_character_
  geneLookup    <- geneLookup[, intersect(c("seqnames", "start", "end", "gene_id", "ID", "Name", "note"),
                                          names(geneLookup))]
} else {
  message("Notes GTF not found; titles will say 'function unknown': ", notesGTF)
  geneLookup <- data.frame(seqnames = character(), start = integer(), end = integer(),
                           note = character())
}

# Example search - replace "kinase" with your keyword/protein of interest
# hits <- geneLookup[grepl("kinase", geneLookup$note, ignore.case = TRUE), ]
# hits


# ---------------------------------------------------------------------------
# Sample metadata (07.0 sample sheet: filtered BAMs, 03.1 pooled peaks, bigWigs)
# ---------------------------------------------------------------------------
allSamps <- read.csv(metaFile, stringsAsFactors = FALSE)
# [CHANGED] removed gsub("03.1_MacsPeaks_with_replicate_count/", ...): that folder holds the current pooled peaks


# ---------------------------------------------------------------------------
# getTopPeakGenes: ad hoc ranking by MACS signal (not a differential test).
# Kept from the original for quick looks; not used below.
# ---------------------------------------------------------------------------
# getTopPeakGenes <- function(population, factorName, day, topN = 50) {
#   samps <- allSamps[allSamps$Population == population &
#                       allSamps$Factor == factorName &
#                       allSamps$Day == day, ]
#   peakGR <- importPeaks(samps$Peaks[1])
#   peakGR <- peakGR[order(-peakGR$signalValue)]
#   topPeaks <- peakGR[seq_len(min(topN, length(peakGR)))]
#   hitGenes <- structGenes[unique(queryHits(findOverlaps(structGenes, topPeaks)))]
#   geneIDs  <- mcols(hitGenes)$ID
#   list(samps = samps,
#        annotated = geneLookup[geneLookup$gene_id %in% geneIDs, ])
# }


# ===========================================================================
# DiffBind results — load, filter, annotate, and find plot-ready examples
#
# The contrast (targetPop vs refPop, mark, day) is set ONCE in USER SETTINGS.
# Everything below (DiffBind file, tracks, peak-overlap check, plots) follows it.
# ===========================================================================

# ---- 1. Load the matching DiffBind comparison file -------------------------
diffFile <- file.path(diffbindDir, tolower(dayName), paste0(contrastID, "-OUT.csv"))     # [CHANGED] new folder and names (AC, never H3K27AC)
if (!file.exists(diffFile)) {
  stop("DiffBind file not found: ", diffFile, "\nAvailable: ",
       paste(list.files(file.path(diffbindDir, tolower(dayName)), pattern = "-OUT\\.csv$"), collapse = ", "))
}
diffResults <- read.csv(diffFile)
diffResults <- diffResults[order(diffResults$FDR, -abs(diffResults$Fold)), ]
sigResults  <- diffResults[!is.na(diffResults$FDR) & diffResults$FDR <= fdrCut, ]       # [CHANGED] <= (paper: FDR <= 0.10)

# Fold sign checked 2026-11-27 for all 24 files: Fold > 0 = higher in the FIRST site (targetPop)
popEnriched <- if (direction == "up") sigResults[sigResults$Fold > 0, ] else sigResults[sigResults$Fold < 0, ]
enrichedPop <- if (direction == "up") targetPop else refPop
cat(sprintf("DiffBind sites: %d | FDR <= %g: %d | higher in %s: %d\n",
            nrow(diffResults), fdrCut, nrow(sigResults), enrichedPop, nrow(popEnriched)))


# ---- 2. Representative replicates of BOTH sites -> tracks ------------------
pickRep <- function(pop) {                                                                  # [CHANGED] both sites, 04.2 representative replicate
  s <- allSamps[allSamps$Population == pop & allSamps$Factor == markName &
                  allSamps$Day == dayName & allSamps$Representative, ]
  if (nrow(s) != 1) stop("Expected 1 representative replicate for ", pop, "-", markName, "-", dayName,
                         ", found ", nrow(s))
  s
}
tgtSamp  <- pickRep(targetPop)
refSamp  <- pickRep(refPop)
popSamps <- if (direction == "up") tgtSamp else refSamp     # site where the chosen sites are higher
bwFiles  <- c(tgtSamp$bigWig, refSamp$bigWig)

myTracksPop <- makeGViz(annotationFile = structuralNoMt,
                        peaks = c(tgtSamp$Peaks, refSamp$Peaks),
                        geneBodies = TRUE,
                        transcripts = TRUE,
                        peakNames = c(paste(targetPop, "peaks"), paste(refPop, "peaks")),
                        reads = bwFiles,                                                    # [CHANGED] input-subtracted bigWigs (paper used input-subtracted bigWigs)
                        alignNames = c(tgtSamp$Sample, refSamp$Sample),
                        colorBySignal = FALSE,
                        showPeakIDs = FALSE,
                        txdb = testtx)

dataIdx <- which(vapply(myTracksPop, is, logical(1), "DataTrack"))                         # [CHANGED] one colour per site
for (k in seq_along(dataIdx)) {
  displayPars(myTracksPop[[dataIdx[k]]]) <- list(fill.histogram = trackColors[k], col.histogram = trackColors[k])
}


# ---- 3. Overlap significant sites with each gene's TSS +/- 2 kb -------------
diffGR <- GRanges(seqnames = popEnriched$chr,                                               # [CHANGED] significant sites in the chosen direction (was ALL sites)
                  ranges = IRanges(start = popEnriched$start, end = popEnriched$end),
                  Fold = popEnriched$Fold,
                  FDR = popEnriched$FDR)

# [CHANGED] window = TSS +/- promoterFlank (strand-aware), as Fig. 3B; was the whole gene body +/- 2 kb
structGenesFlanked <- suppressWarnings(promoters(structGenes, upstream = promoterFlank,
                                                 downstream = promoterFlank))
start(structGenesFlanked) <- pmax(1L, start(structGenesFlanked))

ov <- suppressWarnings(findOverlaps(structGenesFlanked, diffGR))

diffAnnotated <- data.frame(
  gene_ID    = mcols(structGenes)$ID[queryHits(ov)],
  gene_Name  = mcols(structGenes)$Name[queryHits(ov)],
  seqnames   = as.character(seqnames(structGenes))[queryHits(ov)],
  gene_start = start(structGenes)[queryHits(ov)],   # unflanked coords for plotting
  gene_end   = end(structGenes)[queryHits(ov)],
  strand     = as.character(strand(structGenes))[queryHits(ov)],
  site_start = start(diffGR)[subjectHits(ov)],                                              # [CHANGED] keep the site itself (for the peak check)
  site_end   = end(diffGR)[subjectHits(ov)],
  Fold       = mcols(diffGR)$Fold[subjectHits(ov)],
  FDR        = mcols(diffGR)$FDR[subjectHits(ov)],
  stringsAsFactors = FALSE
)

# Bring in SwissProt functional notes by COORDINATE match
diffAnnotated$note <- geneLookup$note[
  match(paste(diffAnnotated$seqnames, diffAnnotated$gene_start, diffAnnotated$gene_end),
        paste(geneLookup$seqnames, geneLookup$start, geneLookup$end))
]

diffAnnotated <- diffAnnotated[order(diffAnnotated$FDR, -abs(diffAnnotated$Fold)), ]
cat("Gene-site pairs with the site in TSS +/- ", promoterFlank, " bp: ", nrow(diffAnnotated), "\n", sep = "")


# ---- 4. Keep sites that lie inside a peak of the site where they are higher -
# [CHANGED] tests the significant SITE itself (was: any peak within gene +/- 2 kb)
popPeakGR <- importPeaks(popSamps$Peaks[1])

diffAnnotated$hasPeak <- suppressWarnings(overlapsAny(
  GRanges(diffAnnotated$seqnames,
          IRanges(diffAnnotated$site_start, diffAnnotated$site_end)),
  popPeakGR
))

diffWithPeaks <- diffAnnotated[diffAnnotated$hasPeak & !(diffAnnotated$gene_ID %in% excludeGenes), ]
diffWithPeaks <- diffWithPeaks[order(diffWithPeaks$FDR, -abs(diffWithPeaks$Fold)), ]
diffWithPeaks <- diffWithPeaks[!duplicated(diffWithPeaks$gene_ID), ]                       # [CHANGED] one row per gene (its best site)

cat("Candidate genes (site in a", enrichedPop, "peak):", nrow(diffWithPeaks), "\n")
write.csv(diffWithPeaks, file.path(outDir, "candidates.csv"), row.names = FALSE)
print(head(diffWithPeaks[, c("gene_ID", "seqnames", "gene_start", "gene_end", "Fold", "FDR", "note")], 10))


# ---- 5/6/7. Render the top candidates to PDF + PNG ---------------------------
# Option B (specific gene): diffWithPeaks <- diffWithPeaks[diffWithPeaks$gene_ID == "gene10137", ]
# Option C (keyword):       diffWithPeaks <- diffWithPeaks[grepl("kinase", diffWithPeaks$note, ignore.case = TRUE), ]
savePlot <- function(filename, width = 10, height = 6, dpi = 300, plotExpr) {               # [CHANGED] PDF + PNG, never overwrite
  for (ext in c("pdf", "png")) {
    f <- paste0(filename, ".", ext)
    if (file.exists(f)) { warning("Exists, not overwritten: ", f); next }
    if (ext == "pdf") pdf(f, width = width, height = height) else png(f, width = width, height = height, units = "in", res = dpi, type = "cairo")
    eval(plotExpr, envir = parent.frame())
    dev.off()
  }
}

nPlot <- min(topN, nrow(diffWithPeaks))
for (i in seq_len(nPlot)) {
  hit <- diffWithPeaks[i, ]   # [CHANGED] starts at rank 1 (original plotted row 2)

  note_i  <- if (is.na(hit$note)) "function unknown" else sub("\\s*\\[.*$", "", hit$note)
  title_i <- substr(paste0(hit$gene_ID, " (", note_i, ") - ", enrichedPop, "-enriched, ",
                           targetPop, " vs ", refPop, " ", markName, " ", dayName,
                           ", Fold=", round(hit$Fold, 2),
                           ", FDR=", format(hit$FDR, digits = 2)), 1, 100)

  savePlot(file.path(outDir, sprintf("%02d_%s", i, hit$gene_ID)),
           plotExpr = quote(
             plotAligns(myTracksPop,
                        chromosome = hit$seqnames,
                        start = max(1, hit$gene_start - promoterFlank),
                        end   = hit$gene_end + promoterFlank,
                        title = title_i,
                        heights = trackHeights,
                        titleCex = 0.9,
                        bwFiles = bwFiles)
           ))
  cat("Plotted", i, "of", nPlot, ":", hit$gene_ID, "\n")
}
if (nPlot == 0) cat("No candidate genes for this contrast/cutoff/direction; nothing plotted.\n")

writeLines(c(paste("contrast:", contrastID), paste("FDR <=", fdrCut), paste("direction:", direction),
             paste("topN:", topN), paste("DiffBind file:", diffFile),
             paste("coverage tracks:", paste(bwFiles, collapse = " ")),
             paste("peak tracks:", tgtSamp$Peaks, refSamp$Peaks),
             paste("notes GTF used:", file.exists(notesGTF))),
           file.path(outDir, "run_settings.txt"))
writeLines(capture.output(sessionInfo()), file.path(outDir, "sessionInfo.txt"))
cat("\nDone. Output in", outDir, "\n")
