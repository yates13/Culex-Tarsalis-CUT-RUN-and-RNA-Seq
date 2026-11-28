library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(chipseq)
library(biomaRt)
library(stringr)
library(biomartr)
library(GenomicFeatures)
library(AnnotationDbi)
library(urltools)       # for repairing the corrupted attribute column

# Must be run for non-model organisms
options(ucscChromosomeNames = FALSE)

# ==== working directory ====
setwd("Z:/NAS_Server_Shared/08-Cxt_RVFV/scripts/R Scripts/gVis")

gVisDir     <- "Z:/NAS_Server_Shared/08-Cxt_RVFV/scripts/R Scripts/gVis"
diffbindDir <- "Z:/NAS_Server_Shared/08-Cxt_RVFV/diffbind_d1/Diff_out"

# structuralGFF <- file.path(gVisDir, "Culex-tarsalis_knwr_BASEFEATURES_CtarK1.gff3")
# notesGTF      <- file.path(gVisDir, "basefeatures_updated_high_confidence.with_notes.gtf")

setwd(gVisDir)  
# ---------------------------------------------------------------------------
# TWO FILES, TWO JOBS
#
#   structuralGFF -> full gene/transcript/exon/CDS structure (original file).
#                    Used for makeTxDbFromGFF(), makeGViz(), plotAligns().
#
#   notesGTF      -> gene-only file with SwissProt functional annotations
#                    merged in (the "note" attribute). Used ONLY to search by
#                    keyword/function and find a gene's coordinates.
#
#   IMPORTANT: the two files do NOT share a gene-ID scheme in the TxDb sense.
#   notesGTF gene_id (e.g. "gene146") must be matched to the structural file's
#   real GENEID (e.g. "maker-tig00000029-exonerate_protein2genome-gene-4.5")
#   by genomic COORDINATES, using findGeneID() below — not by ID string.
# ---------------------------------------------------------------------------
structuralGFF <- "Culex-tarsalis_knwr_BASEFEATURES_CtarK1.gff3"
notesGTF      <- "basefeatures_updated_high_confidence.with_notes.gtf"

structuralNoMt <- sub("\\.gff3$", "_noMt.gff3", structuralGFF)

# Mitochondrial contig name in the structural GFF3 file
mtContigName <- "Mt"

# Defined ONCE here; referenced everywhere else (makeGViz, getTopPeakGenes, etc.)
macsPeakCols <- c(signalValue = "numeric", pValue = "numeric",
                  qValue = "numeric", peak = "integer")


# ---------------------------------------------------------------------------
# Strip the Mt contig from the structural file
# (Mt gene entries are missing proper ID attributes and break makeTxDbFromGFF)
# Only re-runs if the stripped file doesn't exist yet or the source is newer.
# ---------------------------------------------------------------------------
if (!file.exists(structuralNoMt) || file.mtime(structuralGFF) > file.mtime(structuralNoMt)) {
  rawLines  <- readLines(structuralGFF)
  keepLines <- grepl("^#", rawLines) | !grepl(paste0("^", mtContigName, "\t"), rawLines)
  writeLines(rawLines[keepLines], structuralNoMt)
}

# Quick sanity check - abort early with a clear message if repair failed
checkLines <- readLines(structuralNoMt, n = 200)
if (any(grepl("%3d|%3b|attribute=", checkLines, ignore.case = TRUE))) {
  stop("GFF repair incomplete - corrupted attribute syntax still present in ",
       structuralNoMt)
}

# Build the TxDb ONCE here and reuse everywhere (rebuilding per-plot is slow
# and was the source of several earlier bugs when the raw file was used by
# accident instead of the Mt-stripped one).
testtx <- GenomicFeatures::makeTxDbFromGFF(structuralNoMt, format = "gff3")

# Gene-body GRanges, used for coordinate-based ID lookups and peak overlaps
structGFF   <- rtracklayer::import(structuralGFF)
structGenes <- structGFF[structGFF$type == "gene"]


# ---------------------------------------------------------------------------
# makeGViz: build a list of Gviz tracks (gene bodies, transcripts, alignments, peaks)
# Uses the STRUCTURAL file only.
#
#   colorBySignal = TRUE  -> peak boxes shaded by MACS signalValue (QC/exploration)
#   colorBySignal = FALSE -> flat blue3 peak boxes (cleaner, for figures/posters)
#   showPeakIDs   = TRUE  -> draws each peak's BED "id" text inside the box
#   showPeakIDs   = FALSE -> no in-box text (default; avoids "unknown.####" labels)
# ---------------------------------------------------------------------------
makeGViz <- function(annotationFile = "", reads = c(), peaks = c(),
                     geneBodies = TRUE, transcripts = FALSE, chromosome = "",
                     orgName = "Culex tarsalis",
                     peakNames = c(), alignNames = c(),
                     colorBySignal = FALSE, showPeakIDs = FALSE) {
  
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
    noMtFile <- sub("\\.gff3$", "_noMt.gff3", annotationFile)
    tx       <- GenomicFeatures::makeTxDbFromGFF(noMtFile, format = "gff3")
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
      alTrack <- DataTrack(range = reads[i], type = "histogram",
                           name = alignNames[i], isPaired = TRUE, cex.title = 1)
      peakLists <- c(peakLists, alTrack)
    }
  }
  
  if (length(peaks) > 0) {
    if (length(peakNames) == 0) {
      peakNames <- paste0("Peaks_", seq_along(peaks))
    }
    for (i in seq_along(peaks)) {
      thisPeakGR <- import(peaks[i], format = "BED", extraCols = macsPeakCols)
      
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


# ---------------------------------------------------------------------------
# plotAligns: plot a set of tracks over a region, optionally centered on a gene
# Default gff is the Mt-stripped file and default txdb reuses the one built
# above, so this no longer rebuilds/breaks on every call.
# ---------------------------------------------------------------------------

plotAligns <- function(tracks, chromosome, start, end, gene = "", heights = c(),
                       title = "", gff = structuralNoMt, txdb = testtx,
                       showAxis = TRUE, baseSize = 1,
                       background.title = "black",
                       autoTitle = FALSE, geneLookupTable = NULL,
                       titleCex = 1.1) {          # <-- new param
  
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
# structural file's real GENEID, by matching genomic coordinates.
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
# Search genes by function using the NOTES file (gene-only, SwissProt notes)
# ---------------------------------------------------------------------------
notesImported <- rtracklayer::import(notesGTF, format = "gtf")
genesOnly     <- notesImported[notesImported$type == "gene"]

geneLookup <- as.data.frame(genesOnly)[, c("seqnames", "start", "end",
                                           "gene_id", "ID", "Name", "note")]

# Example search - replace "kinase" with your keyword/protein of interest
# hits <- geneLookup[grepl("kinase", geneLookup$note, ignore.case = TRUE), ]
# hits


# ---------------------------------------------------------------------------
# Sample metadata
# ---------------------------------------------------------------------------
allSamps <- read.csv("CtarK1_Meta_File.csv")

# Fix stale subfolder in Peaks paths (03.1_MacsPeaks_with_replicate_count no longer exists)
allSamps$Peaks <- gsub("03.1_MacsPeaks_with_replicate_count/", "", allSamps$Peaks)


# ---------------------------------------------------------------------------
# getTopPeakGenes: find genes with called peaks for a given Population/Factor/Day,
# ranked by MACS signal strength. (Only looks at replicate 1 — see message below
# if you need to merge replicates instead.)
# ---------------------------------------------------------------------------
# getTopPeakGenes <- function(population, factorName, day, topN = 50) {
#   samps <- allSamps[allSamps$Population == population &
#                       allSamps$Factor == factorName &
#                       allSamps$Day == day, ]
#   if (nrow(samps) > 1) {
#     message("Multiple replicates found for ", population, "/", factorName, "/", day,
#             " - using only the first (", samps$Peaks[1], ")")
#   }
#   peakGR <- rtracklayer::import(samps$Peaks[1], format = "BED", extraCols = macsPeakCols)
#   
#   peakGR <- peakGR[order(-peakGR$signalValue)]
#   topPeaks <- peakGR[seq_len(min(topN, length(peakGR)))]
#   
#   hitGenes <- structGenes[unique(queryHits(findOverlaps(structGenes, topPeaks)))]
#   geneIDs  <- mcols(hitGenes)$ID
#   
#   list(samps = samps,
#        annotated = geneLookup[geneLookup$gene_id %in% geneIDs, ])
# }


# ---------------------------------------------------------------------------
# Example: quick single-gene plot using getTopPeakGenes() (ad hoc ranking,
# not tied to a formal differential-binding test). Useful for a fast look,
# but for a poster-quality, statistically supported example use the DiffBind
# section further down instead.
# ---------------------------------------------------------------------------
# mpResult <- getTopPeakGenes("MP", "AC", "D1")
# zhResult <- getTopPeakGenes("ZH", "AC", "D1")
# 
# mpResult$annotated[, c("gene_id", "note")]
# zhResult$annotated[, c("gene_id", "note")]
# 
# mp_target  <- geneLookup[geneLookup$gene_id == "gene10766", ]
# mp_realID  <- findGeneID(mp_target$seqnames, mp_target$start, mp_target$end)
# 
# myTracksMP <- makeGViz(annotationFile = structuralGFF,
#                        peaks = mpResult$samps$Peaks,
#                        geneBodies = TRUE,
#                        transcripts = TRUE,
#                        peakNames = c("MP_AC_D1_Peaks"),
#                        reads = mpResult$samps$bamReads,
#                        alignNames = c("MP_AC_D1"),
#                        colorBySignal = FALSE,
#                        showPeakIDs = FALSE)
# 
# plotAligns(myTracksMP, gene = mp_realID,
#            title = "gene10766 (MP) - Protein pangolin, isoforms A/H/I/S")


# ===========================================================================
# DiffBind results — load, filter, annotate, and find plot-ready examples
#
# Pick ONE population to feature with `targetPop` below. Everything after
# that line (DiffBind file, tracks, peak-overlap check, and the final plot)
# automatically stays consistent with that choice, so you never end up
# plotting a ZH-enriched gene on MP tracks (or vice versa) by accident.
# ===========================================================================
setwd("Z:/NAS_Server_Shared/08-Cxt_RVFV/diffbind_d1/Diff_out")   #ZH directory
#setwd("Z:/NAS_Server_Shared/08-Cxt_RVFV/diffbind_d1")   #MP directory 

targetPop <- "ZH"   # <- set to "ZH" or "MP"
#targetPop <- "MP" 

# ---- 1. Load the matching DiffBind comparison file -------------------------
diffFile <- switch(targetPop,
                   "ZH" = "ZH_H3K27AC_d1_vs_BF-OUT.csv",
                   "MP" = "MP_AC_d1_vs_BF-OUT.csv")   # adjust filename to match yours

setwd(diffbindDir)  #ZH directory
diffResults <- read.csv(diffFile)
diffResults <- diffResults[order(diffResults$FDR, -abs(diffResults$Fold)), ]
sigResults  <- diffResults[diffResults$FDR < 0.1, ]   # tune cutoff as needed

# NOTE: confirm the Fold sign convention for whichever file you load —
# check a row where one Conc column is 0 and verify Fold = Conc_target - Conc_BF.
# For the ZH file this was confirmed: positive Fold = ZH-enriched.
popEnriched <- sigResults[sigResults$Fold > 0, ]


# ---- 2. Pull matching sample metadata + build tracks for that population --

popSamps <- allSamps[allSamps$Population == targetPop &
                       allSamps$Factor == "AC" &
                       allSamps$Day == "D1", ]

setwd(gVisDir)
myTracksPop <- makeGViz(annotationFile = structuralGFF,
                        peaks = popSamps$Peaks,
                        geneBodies = TRUE,
                        transcripts = TRUE,
                        peakNames = "Peaks",                 # was paste0(targetPop, "_AC_D1_Peaks")
                        reads = popSamps$bamReads,
                        alignNames = "Coverage",     # was paste0(targetPop, "_AC_D1")
                        colorBySignal = FALSE,
                        showPeakIDs = FALSE)

# ---- 3. Overlap DiffBind hits with structural genes ------------------------
# promoterFlank widens gene boundaries before overlap so promoter-proximal
# peaks (sitting just outside the gene body) aren't missed.
diffGR <- GRanges(seqnames = diffResults$chr,
                  ranges = IRanges(start = diffResults$start, end = diffResults$end),
                  Fold = diffResults$Fold,
                  FDR = diffResults$FDR)

promoterFlank <- 2000   # bp of flanking sequence to include on each side

structGenesFlanked <- structGenes
start(structGenesFlanked) <- pmax(1, start(structGenes) - promoterFlank)
end(structGenesFlanked)   <- end(structGenes) + promoterFlank

ov <- findOverlaps(structGenesFlanked, diffGR)

diffAnnotated <- data.frame(
  gene_ID    = mcols(structGenes)$ID[queryHits(ov)],
  gene_Name  = mcols(structGenes)$Name[queryHits(ov)],
  seqnames   = as.character(seqnames(structGenes))[queryHits(ov)],
  gene_start = start(structGenes)[queryHits(ov)],   # unflanked coords for plotting
  gene_end   = end(structGenes)[queryHits(ov)],
  Fold       = mcols(diffGR)$Fold[subjectHits(ov)],
  FDR        = mcols(diffGR)$FDR[subjectHits(ov)]
)

# Bring in SwissProt functional notes by COORDINATE match (geneLookup$gene_id
# and structGenes$ID use different ID schemes and can't be string-matched).
diffAnnotated$note <- geneLookup$note[
  match(paste(diffAnnotated$seqnames, diffAnnotated$gene_start, diffAnnotated$gene_end),
        paste(geneLookup$seqnames, geneLookup$start, geneLookup$end))
]

diffAnnotated <- diffAnnotated[order(diffAnnotated$FDR, -abs(diffAnnotated$Fold)), ]
head(diffAnnotated, 20)


# ---- 4. Restrict to genes that actually have a called peak nearby ---------
# hasPeak is computed as a column ON diffAnnotated directly (not a separate
# vector) to avoid row misalignment if diffAnnotated gets reordered later.
# Uses THIS population's own peak file, and the same promoter flank as above,
# and is pre-filtered to this population's enrichment direction (Fold > 0).
popPeakGR <- rtracklayer::import(popSamps$Peaks[1], format = "BED",
                                 extraCols = macsPeakCols)

diffAnnotated$hasPeak <- overlapsAny(
  GRanges(diffAnnotated$seqnames,
          IRanges(pmax(1, diffAnnotated$gene_start - promoterFlank),
                  diffAnnotated$gene_end + promoterFlank)),
  popPeakGR
)

diffWithPeaks <- diffAnnotated[diffAnnotated$hasPeak & diffAnnotated$Fold > 0, ]
diffWithPeaks <- diffWithPeaks[order(diffWithPeaks$FDR, -abs(diffWithPeaks$Fold)), ]

nrow(diffWithPeaks)   # sanity check before plotting
head(diffWithPeaks[, c("gene_ID", "note", "Fold", "FDR")], 10)


# ---- 5. Pick a gene to plot -------------------------------------------------
# Option A: by rank (row number in the sorted table) — default: top hit
topHit <- diffWithPeaks[2, ]
# topHit <- diffWithPeaks[5, ]      # e.g. the 5th-ranked gene instead

# Option B: by specific gene ID
#topHit <- diffWithPeaks[diffWithPeaks$gene_ID == "gene10137", ]

# Option C: search notes for a keyword (e.g. "kinase", "transport", "immune")
# keywordHits <- diffWithPeaks[grepl("kinase", diffWithPeaks$note, ignore.case = TRUE), ]
# head(keywordHits[, c("gene_ID", "note", "Fold", "FDR")])
# topHit <- keywordHits[1, ]


# ---- 6. Build the title and plot -------------------------------------------
noteShort <- if (is.na(topHit$note)) {
  "function unknown"
} else {
  sub("\\s*\\[.*$", "", topHit$note)   # strips everything from "[" onward
}

# titleText <- paste0(topHit$gene_ID, " (", noteShort, ") - ",
#                     targetPop, "-enriched, Fold=", round(topHit$Fold, 2),
#                     ", FDR=", format(topHit$FDR, digits = 2))
# titleText <- substr(titleText, 1, 80)   # hard cap so titles never overflow the plot
sampleID <- paste0(targetPop, "_AC_D1")

titleText <- paste0(topHit$gene_ID, " (", noteShort, ") - ", sampleID)
titleText <- substr(titleText, 1, 60) 

plotAligns(myTracksPop,
           chromosome = topHit$seqnames,
           start = topHit$gene_start - promoterFlank,
           end   = topHit$gene_end + promoterFlank,
           title = titleText)


# ---- 7. (Optional) batch-render several candidates to compare -------------
# Uncomment this whole section to render the top N candidate genes to PNG
# files in your working directory instead of just plotting one at a time.

# savePlot <- function(filename, width = 10, height = 6, dpi = 300, plotExpr) {
#   png(filename, width = width, height = height, units = "in", res = dpi)
#   eval(plotExpr)
#   dev.off()
# }
#
# topN <- 5   # how many candidate genes to render
#
# for (i in seq_len(min(topN, nrow(diffWithPeaks)))) {
#   hit <- diffWithPeaks[i, ]
#
#   note_i  <- if (is.na(hit$note)) "function unknown" else sub("\\s*\\[.*$", "", hit$note)
#   title_i <- substr(paste0(hit$gene_ID, " (", note_i, ") - ",
#                            targetPop, "-enriched, Fold=", round(hit$Fold, 2),
#                            ", FDR=", format(hit$FDR, digits = 2)), 1, 80)
#
#   savePlot(paste0("candidate_", hit$gene_ID, ".png"),
#            plotExpr = quote(
#              plotAligns(myTracksPop,
#                        chromosome = hit$seqnames,
#                        start = hit$gene_start - promoterFlank,
#                        end   = hit$gene_end + promoterFlank,
#                        title = title_i)
#            ))
# }

