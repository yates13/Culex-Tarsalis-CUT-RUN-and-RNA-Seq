library("tools")
if (requireNamespace("readr", quietly = TRUE)) options(readr.show_col_types = FALSE)

# =============================================================================
# 06.1_getGeneDistances_Culex_NAS.R  (RStudio on the NAS)
# Same as getGeneDistances_Culex_SY_Main_GTF_v3_offset.R, updated for the
# 2026-11-24 rerun. The Alpine version (06.1_getGeneDistances_Culex_Alpine.R)
# has the same functions; only the config at the bottom differs.
#
# QUICK START: press Ctrl+F and search for  >>>  to jump to every line you
# need to change. Those marked lines are the ONLY ones to edit.
# =============================================================================

# =============================================================================
# WHAT CHANGED FROM THE v3_offset SCRIPT (and why)
# =============================================================================
#   1. GTF = CtarK1_MAKER_genes.gtf: the 14,713 MAKER genes used for the TSS
#      figures (made by 06.0_Make_MAKER_gene_GTF.sh on Alpine). Same IDs
#      (gene1 ... gene14720), so geneIDPattern "^gene[0-9]+$" still works.
#   2. Peak files = the new DiffBind results (filtered BAMs, 2026-11-24).
#      File names now always use AC (old names mixed AC and H3K27AC), so the
#      5 H3K27AC names in `runs` / `runs_d3` were changed to AC.
#   3. exclude = "^Mt$": the mitochondrial contig in CtarK1.fa is named "Mt";
#      the old "^Mt_CtarK1" matched nothing (no effect on results: the MAKER
#      GTF has no Mt genes).
#   4. GTF files are read through .readGTF() (readr if installed, base R if not).
#   Everything else (FDR <= 0.1, 2000 bp, inGene = TRUE, output columns and
#   names) is unchanged.
#
# ---- HOW TO RUN (quick reference) -------------------------------------------
#   A. Copy from Alpine to the NAS first (folder names below are suggestions):
#        Cxt_reference/CtarK1_MAKER_genes.gtf
#            -> Z:/NAS_Server_Shared/08-Cxt_RVFV/RNASeq/BaseFeatures/
#        Round1_Round2_With_Duplicates/05_DiffBind_2026-11-24/d1 and d3
#            -> Z:/NAS_Server_Shared/08-Cxt_RVFV/diffbind_2026-11-24/d1 and d3
#   B. Set the output folder (setwd), the GTF, and the two peakDir lines (>>>).
#   C. Run the whole script. Outputs go to the output folder:
#        <sampName>_DistanceDF.csv and <sampName>_DistanceDF.RDS
#        sampName = <run>_<day>_vs_<contrast>_FDR_<fdr>_<gtfTag>
#        e.g. ZH_AC_d1_vs_BF_FDR_0.1_CtarK1_MAKER_genes
#      Re-running with the same GTF/FDR/files OVERWRITES those files; use a new
#      output folder (setwd) to keep both versions.
#   D. Peak files need a header with: chr, start, end, Fold, p.value, FDR
#      (DiffBind dba.report output). Missing files are skipped with a warning.
#   E. "Near a gene" = a peak end within `distance` bp of either gene END, or a
#      peak overlapping the gene body (inGene = TRUE). This is not strand-aware;
#      the 5 Prime Dist column is.
#   F. Check it worked: the data.frame() at the end of each section prints rows
#      per run; NA = the run was skipped or failed (see warnings).
# =============================================================================

setwd("Z:/NAS_Server_Shared/08-Cxt_RVFV/scripts/R Scripts/GeneDistance/2026-11-24_MAKER")   # >>> CHANGE: output folder (create it first)

# =============================================================================
# FUNCTIONS (same as getGeneDistances_Culex_SY_Main_GTF_v3_offset.R)
# The only change: GTF files are read through .readGTF(), which uses readr
# when it is installed and base R otherwise (readr may not be in the Alpine
# R env). Both give the same 9 columns; the code below uses column numbers.
# You never need to edit getGeneDistances() or .processChromosome().
# =============================================================================

.readGTF = function(gtfFile) {
  if (requireNamespace("readr", quietly = TRUE)) {
    as.data.frame(readr::read_tsv(gtfFile, col_names = FALSE, comment = "#"))
  } else {
    read.delim(gtfFile, header = FALSE, comment.char = "#", quote = "",
               stringsAsFactors = FALSE)
  }
}

#' Quick diagnostic: see what's in the "feature" column (col 3) of a GTF
#' before deciding whether to filter to feature == "gene".
checkFeatureTypes = function(gtfFile, excludeContig = NULL) {
  gtfDF = .readGTF(gtfFile)
  if (!is.null(excludeContig)) {
    gtfDF = gtfDF[!grepl(excludeContig, gtfDF[[1]]), , drop = FALSE]
  }
  print(table(gtfDF[[3]]))
  invisible(gtfDF)
}

#' Pull a single attribute value by key name out of a GTF-style attribute
#' string like: gene_id "gene123"; Name "maker-..."; note "...";
#' Returns NA if the key is not present in that row's attribute string.
extractAttr = function(attrString, key) {
  pattern = paste0(key, '\\s+"([^"]*)"')
  m = regmatches(attrString, regexpr(pattern, attrString))
  if (length(m) == 0) return(NA_character_)
  sub(pattern, "\\1", m)
}


#' Process a single chromosome's genes against that chromosome's peaks.
#' Returns a list of per-pair result vectors (only for pairs that satisfy
#' nearEnd/inGene), or NULL if nothing on this chromosome qualifies.
.processChromosome = function(rnaIdx, rnaStart, rnaEnd, orientVec,
                               dnaIdxLocal, dnaStart, dnaEnd,
                               distance, inGene, nearEnd, verbose) {
  nGenes = length(rnaIdx)
  nPeaks = length(dnaIdxLocal)
  if (nGenes == 0 || nPeaks == 0) return(NULL)

  pairCount = nGenes * nPeaks   # small: genes on THIS chrom * peaks on THIS chrom

  withinDist  = rep(FALSE, pairCount)
  dnaDist     = rep(0, pairCount)
  dna5Prime   = rep(0, pairCount)
  dnaSamp     = rep(0, pairCount)
  insideGene  = rep(FALSE, pairCount)
  rnaSamp     = rep(0, pairCount)
  over5prime  = rep(FALSE, pairCount)
  over3prime  = rep(FALSE, pairCount)
  dna3Prime   = rep(0, pairCount)

  for (ii in seq_len(nGenes)) {
    i = rnaIdx[ii]
    orientation = orientVec[ii]

    for (jj in seq_len(nPeaks)) {
      j = dnaIdxLocal[jj]
      tempNum = (ii - 1) * nPeaks + jj

      dist = Inf
      temp = abs(as.integer(rnaStart[i]) - as.integer(dnaStart[j])); if (temp < dist) dist = temp
      temp = abs(as.integer(rnaStart[i]) - as.integer(dnaEnd[j]));   if (temp < dist) dist = temp
      temp = abs(as.integer(rnaEnd[i])   - as.integer(dnaStart[j])); if (temp < dist) dist = temp
      temp = abs(as.integer(rnaEnd[i])   - as.integer(dnaEnd[j]));   if (temp < dist) dist = temp

      setCondition = FALSE

      if (nearEnd && (dist < distance)) {
        if (verbose) {
          print("RNA"); print(i); print("Distance is"); print(dist)
        }
        setCondition = TRUE
      }

      if (inGene && (
            (as.integer(rnaEnd[i])   < as.integer(dnaEnd[j])   && as.integer(rnaStart[i]) > as.integer(dnaStart[j])) ||
            (as.integer(rnaStart[i]) < as.integer(dnaStart[j]) && as.integer(rnaEnd[i])   > as.integer(dnaEnd[j])) ||
            (as.integer(rnaEnd[i])   > as.integer(dnaStart[j]) && as.integer(rnaEnd[i])   < as.integer(dnaEnd[j])) ||
            (as.integer(rnaStart[i]) > as.integer(dnaStart[j]) && as.integer(rnaStart[i]) < as.integer(dnaEnd[j]))
          )) {
        if (verbose) print("Inside the gene!")

        if ((as.integer(rnaEnd[i]) > as.integer(dnaStart[j])) && (as.integer(rnaEnd[i]) < as.integer(dnaEnd[j])) && orientation == '-') {
          over5prime[tempNum] = TRUE
        }
        if ((as.integer(rnaStart[i]) > as.integer(dnaStart[j])) && (as.integer(rnaStart[i]) < as.integer(dnaEnd[j])) && orientation == '+') {
          over5prime[tempNum] = TRUE
        }
        if ((as.integer(rnaStart[i]) > as.integer(dnaStart[j])) && (as.integer(rnaStart[i]) < as.integer(dnaEnd[j])) && orientation == '-') {
          over3prime[tempNum] = TRUE
        }
        if ((as.integer(rnaEnd[i]) > as.integer(dnaStart[j])) && (as.integer(rnaEnd[i]) < as.integer(dnaEnd[j])) && orientation == '+') {
          over3prime[tempNum] = TRUE
        }

        insideGene[tempNum] = TRUE
        setCondition = TRUE
      }

      if (setCondition == TRUE) {
        withinDist[tempNum] = TRUE

        if (orientation == '+') {
          dna5Prime[tempNum] = min(c(abs(as.integer(rnaStart[i]) - as.integer(dnaEnd[j])), abs(as.integer(rnaStart[i]) - as.integer(dnaStart[j]))))
          dna3Prime[tempNum] = min(c(abs(as.integer(rnaEnd[i]) - as.integer(dnaEnd[j])), abs(as.integer(rnaEnd[i]) - as.integer(dnaStart[j]))))
        }
        if (orientation == '-') {
          dna5Prime[tempNum] = min(c(abs(as.integer(rnaEnd[i]) - as.integer(dnaEnd[j])), abs(as.integer(rnaEnd[i]) - as.integer(dnaStart[j]))))
          dna3Prime[tempNum] = min(c(abs(as.integer(rnaStart[i]) - as.integer(dnaEnd[j])), abs(as.integer(rnaStart[i]) - as.integer(dnaStart[j]))))
        }

        dnaSamp[tempNum] = j          # global peak index, preserved
        dnaDist[tempNum] = dist
        rnaSamp[tempNum] = i          # global gene index, preserved
      }
    }
  }

  keep = withinDist
  if (!any(keep)) return(NULL)

  list(
    rnaSamp    = rnaSamp[keep],
    dnaSamp    = dnaSamp[keep],
    dna5Prime  = dna5Prime[keep],
    dna3Prime  = dna3Prime[keep],
    insideGene = insideGene[keep],
    over5prime = over5prime[keep],
    over3prime = over3prime[keep]
  )
}


#' gtfFile: 9-column GTF with GTF-style key "value"; attributes (col 9).
#' diffPeakFile: DiffBind differential peak output (CSV with header, or a
#'   headerless tsv/txt fallback, matched by column NAME).
#' excludeContig: optional regex matched against the chromosome column of
#'   BOTH the GTF and the peak file; matching rows are dropped first.
#' featureFilter: optional string; only GTF rows where column 3 equals it are kept.
getGeneDistances = function(gtfFile, diffPeakFile, topN = -1, verbose = FALSE,
                            sampName = "Culex_DistanceDF",
                            distance = 10000, saveResults = TRUE,
                            inGene = TRUE, nearEnd = TRUE, bedFormat = FALSE,
                            excludeContig = NULL, featureFilter = NULL,
                            fdrFilter = NULL,
                            geneIDPattern = "^gene[0-9]+$") {

  orientCol = -1
  startCol  = -1
  endCol    = -1
  nameCol   = -1
  chromCol  = -1
  featureCol = 3

  if (bedFormat) {
    orientCol = 6
    startCol  = 2
    endCol    = 3
    chromCol  = 1
    nameCol   = 10
  } else {
    # GTF/GFF
    orientCol = 7
    startCol  = 4
    endCol    = 5
    chromCol  = 1
    nameCol   = 9
  }

  rnaDF = .readGTF(gtfFile)
  if (!is.null(excludeContig)) {
    keep = !grepl(excludeContig, rnaDF[[chromCol]])
    rnaDF = rnaDF[keep, , drop = FALSE]
  }

  if (!is.null(featureFilter)) {
    rnaDF = rnaDF[rnaDF[[featureCol]] == featureFilter, , drop = FALSE]
  }
  # Keep only rows whose gene_id matches geneIDPattern (MAKER: gene1 ... gene14720)
  if (!is.null(geneIDPattern)) {
    extractedIDs = sapply(rnaDF[[nameCol]], extractAttr, key = "gene_id")
    keepGenes = grepl(geneIDPattern, extractedIDs)
    if (verbose) {
      cat(sprintf("geneIDPattern filter: kept %d of %d GTF rows\n",
                  sum(keepGenes), length(keepGenes)))
    }
    rnaDF = rnaDF[keepGenes, , drop = FALSE]
  }

  if (nrow(rnaDF) == 0) {
    print("No GTF rows remain after filtering - check excludeContig/featureFilter")
    return()
  }

  rnaNames = as.list(rnaDF[[nameCol]])

  diffExt = file_ext(diffPeakFile)
  sigResDF = ""

  if (diffExt == "csv") {
    sigResDF = read.csv(diffPeakFile)
    dnaStart = sigResDF$start
    dnaEnd   = sigResDF$end
    dnaChrom = as.character(sigResDF$chr)
    dnaFold  = sigResDF$Fold
    pvalue   = sigResDF$p.value
    fdr      = sigResDF$FDR
  } else if (diffExt %in% c("tsv", "txt")) {
    sigResDF = as.data.frame(readr::read_tsv(diffPeakFile, col_names = TRUE))
    dnaStart = sigResDF$start
    dnaEnd   = sigResDF$end
    dnaChrom = as.character(sigResDF$chr)
    dnaFold  = sigResDF$Fold
    pvalue   = sigResDF$p.value
    fdr      = sigResDF$FDR
  } else {
    stop(paste("Unrecognized diffPeakFile extension:", diffExt,
               "- expected csv, tsv, or txt with a header row"))
  }

  if (!is.null(excludeContig)) {
    keepPeaks = !grepl(excludeContig, dnaChrom)
    dnaStart = dnaStart[keepPeaks]
    dnaEnd   = dnaEnd[keepPeaks]
    dnaChrom = dnaChrom[keepPeaks]
    dnaFold  = dnaFold[keepPeaks]
    pvalue   = pvalue[keepPeaks]
    fdr      = fdr[keepPeaks]
  }

  rnaStart = rnaDF[[startCol]]
  rnaEnd   = rnaDF[[endCol]]
  rnaChrom = rnaDF[[chromCol]]
  rnaOrient = rnaDF[[orientCol]]

  rnaLen = length(rnaChrom)
  dnaLen = length(dnaStart)

  # Split gene and peak indices by chromosome and only compare pairs on the
  # same chromosome (avoids the genes x peaks integer overflow).
  commonChroms = intersect(unique(rnaChrom), unique(dnaChrom))

  if (length(commonChroms) == 0) {
    print("No chromosome names in common between GTF and peak file - check naming/excludeContig")
    return()
  }

  allResults = list()

  for (chrom in commonChroms) {
    rnaIdx = which(rnaChrom == chrom)
    dnaIdxLocal = which(dnaChrom == chrom)

    if (verbose) {
      cat(sprintf("Chromosome %s: %d genes x %d peaks\n",
                   chrom, length(rnaIdx), length(dnaIdxLocal)))
    }

    res = .processChromosome(
      rnaIdx = rnaIdx, rnaStart = rnaStart, rnaEnd = rnaEnd,
      orientVec = rnaOrient[rnaIdx],
      dnaIdxLocal = dnaIdxLocal, dnaStart = dnaStart, dnaEnd = dnaEnd,
      distance = distance, inGene = inGene, nearEnd = nearEnd, verbose = verbose
    )

    if (!is.null(res)) allResults[[length(allResults) + 1]] = res
  }

  if (length(allResults) == 0) {
    print("No matching regions found for the given parameters")
    return()
  }

  rnaSamp    = unlist(lapply(allResults, `[[`, "rnaSamp"))
  dnaSamp    = unlist(lapply(allResults, `[[`, "dnaSamp"))
  dna5Prime  = unlist(lapply(allResults, `[[`, "dna5Prime"))
  dna3Prime  = unlist(lapply(allResults, `[[`, "dna3Prime"))
  insideGene = unlist(lapply(allResults, `[[`, "insideGene"))
  over5prime = unlist(lapply(allResults, `[[`, "over5prime"))
  over3prime = unlist(lapply(allResults, `[[`, "over3prime"))

  geneChrom  = rep(0, length(rnaSamp))
  geneOrient = rep(0, length(rnaSamp))
  geneStart  = rep(0, length(rnaSamp))
  geneEnd    = rep(0, length(rnaSamp))
  rnaTitle   = rep(0, length(rnaSamp))

  for (i in seq_along(rnaSamp)) {
    k = rnaSamp[i]
    geneChrom[i]  = rnaDF[k, chromCol]
    geneOrient[i] = rnaDF[k, orientCol]
    geneStart[i]  = rnaDF[k, startCol]
    geneEnd[i]    = rnaDF[k, endCol]
    rnaTitle[i]   = rnaDF[k, nameCol]
  }

  reducedDF = data.frame(geneChrom, geneStart, geneEnd, geneOrient)
  colnames(reducedDF) = c("Chromosome", "Gene Start", "Gene End", "Gene Orientation")

  reducedDF["Peak Number"] = dnaSamp
  reducedDF["5 Prime Dist"] = dna5Prime
  reducedDF["3 Prime Dist"] = dna3Prime
  reducedDF["Peak Inside Gene"] = insideGene

  # ---- Attribute extraction by key name ----
  geneID   = sapply(rnaTitle, extractAttr, key = "gene_id")
  geneName = sapply(rnaTitle, extractAttr, key = "Name")
  geneDes  = sapply(rnaTitle, extractAttr, key = "note")
  spanFlag = sapply(rnaTitle, extractAttr, key = "span_flag")
  supportTranscripts = sapply(rnaTitle, extractAttr, key = "support_transcripts")

  reducedDF["Gene ID"] = geneID
  reducedDF["Gene Name"] = geneName
  reducedDF["Gene Description"] = geneDes
  reducedDF["Span Flag"] = spanFlag
  reducedDF["Support Transcripts"] = supportTranscripts

  foldchange = rep(0, length(dnaSamp))
  for (i in seq_along(dnaSamp)) {
    foldchange[i] = dnaFold[dnaSamp[i]]
  }
  reducedDF["Fold Change"] = foldchange

  peak_pvalues = rep(0, length(dnaSamp))
  peak_fdrs    = rep(0, length(dnaSamp))
  peakStart    = rep(0, length(dnaSamp))
  peakEnd      = rep(0, length(dnaSamp))

  for (i in seq_along(dnaSamp)) {
    peakStart[i]    = as.integer(dnaStart[dnaSamp[i]])
    peakEnd[i]      = as.integer(dnaEnd[dnaSamp[i]])
    peak_pvalues[i] = as.double(pvalue[dnaSamp[i]])
    peak_fdrs[i]    = as.double(fdr[dnaSamp[i]])
  }

  reducedDF["Peak Start"] = peakStart
  reducedDF["Peak End"]   = peakEnd
  reducedDF["Over5Prime"] = over5prime
  reducedDF["Over3Prime"] = over3prime

  reducedDF["Peak Width"] = reducedDF["Peak End"] - reducedDF["Peak Start"]
  reducedDF["Gene Width"] = reducedDF["Gene End"] - reducedDF["Gene Start"]

  ori    <- reducedDF[["Gene Orientation"]]
  pS     <- reducedDF[["Peak Start"]]
  pE     <- reducedDF[["Peak End"]]
  gS     <- reducedDF[["Gene Start"]]
  gE     <- reducedDF[["Gene End"]]

  peakStartOffset <- rep(0, nrow(reducedDF))
  peakEndOffset   <- rep(0, nrow(reducedDF))

  plus  <- ori == "+"
  minus <- ori == "-"

  peakStartOffset[plus]  <- pS[plus] - gS[plus]
  peakEndOffset[plus]    <- pE[plus] - gS[plus]

  peakStartOffset[minus] <- gE[minus] - pE[minus]
  peakEndOffset[minus]   <- gE[minus] - pS[minus]

  reducedDF["Peak Start Offset"] <- peakStartOffset
  reducedDF["Peak End Offset"]   <- peakEndOffset
  reducedDF["p-value"] = peak_pvalues
  reducedDF["FDR"]     = peak_fdrs

  if (topN != -1) {
    reducedDF = reducedDF[reducedDF["Peak Number"] <= topN, ]
  }

  if (!is.null(fdrFilter)) {
    nBefore = nrow(reducedDF)
    reducedDF = reducedDF[reducedDF[["FDR"]] <= fdrFilter, , drop = FALSE]
    if (verbose) {
      cat(sprintf("fdrFilter = %.3g: kept %d of %d rows\n", fdrFilter, nrow(reducedDF), nBefore))
    }
  }

  if (saveResults) {
    saveRDS(reducedDF, paste0(sampName, "_DistanceDF.RDS"))
    write.csv(reducedDF, paste0(sampName, "_DistanceDF.csv"), row.names = FALSE)
  }

  return(reducedDF)
}


# =============================================================================
# USAGE
# =============================================================================

# checkFeatureTypes(sharedGtfFile, excludeContig = "^Mt$")   # expect: gene 14713

##------------------------DAY 1-------------------------------------------------
sharedGtfFile <- "Z:/NAS_Server_Shared/08-Cxt_RVFV/RNASeq/BaseFeatures/CtarK1_MAKER_genes.gtf"   # >>> CHANGE: GTF for day 1 AND day 3

gtfTag <- tools::file_path_sans_ext(basename(sharedGtfFile))

cfg <- list(
  gtf      = sharedGtfFile,
  peakDir  = "Z:/NAS_Server_Shared/08-Cxt_RVFV/diffbind_2026-11-24/d1",   # >>> CHANGE: folder with day 1 DiffBind files
  fdr      = 0.1,                    # >>> CHANGE: FDR cutoff (same as the gVis script)
  distance = 2000,                   # >>> CHANGE: bp window
  exclude  = "^Mt$"                  # >>> CHANGE: contigs to drop
)

# >>> CHANGE (day 1): right side = exact file name inside cfg$peakDir; left side = label
runs <- c(
  ZH_AC = "ZH_AC_d1_vs_BF-OUT.csv",
  ZH_ME = "ZH_ME_d1_vs_BF-OUT.csv",
  ZH_SR = "ZH_SR_d1_vs_BF-OUT.csv",
  MP_AC = "MP_AC_d1_vs_BF-OUT.csv",
  MP_ME = "MP_ME_d1_vs_BF-OUT.csv",
  MP_SR = "MP_SR_d1_vs_BF-OUT.csv",
  BF_AC = "BF_AC_d1_vs_SF-OUT.csv",
  BF_ME = "BF_ME_d1_vs_SF-OUT.csv",
  BF_SR = "BF_SR_d1_vs_SF-OUT.csv"
)

# confirm everything exists before running
print(data.frame(run = names(runs),
                 exists = file.exists(file.path(cfg$peakDir, runs))))

results <- list()
for (run in names(runs)) {
  peakFile <- file.path(cfg$peakDir, runs[[run]])
  contrast <- sub(".*_vs_([A-Za-z]+)-OUT\\.csv$", "\\1", runs[[run]])   # "BF" or "SF"

  if (!file.exists(peakFile)) {
    warning("Missing peak file, skipping ", run, ": ", peakFile)
    next
  }

  message("Running ", run, " ...")
  results[[run]] <- tryCatch(
    getGeneDistances(
      gtfFile       = cfg$gtf,
      diffPeakFile  = peakFile,
      bedFormat     = FALSE,
      distance      = cfg$distance,
      inGene        = TRUE,
      excludeContig = cfg$exclude,
      featureFilter = "gene",
      fdrFilter     = cfg$fdr,
      sampName = paste0(run, "_d1_vs_", contrast, "_FDR_", cfg$fdr, "_", gtfTag),
      verbose       = FALSE
    ),
    error = function(e) { warning(run, " failed: ", conditionMessage(e)); NULL }
  )
}

# quick sanity table
data.frame(
  run  = names(runs),
  rows = sapply(names(runs), function(r) if (is.null(results[[r]])) NA else nrow(results[[r]]))
)
##---------------------------------------------------------------------------------

##-------------------------------DAY 3---------------------------------------------
cfg_d3 <- list(
  gtf      = sharedGtfFile,
  peakDir  = "Z:/NAS_Server_Shared/08-Cxt_RVFV/diffbind_2026-11-24/d3",   # >>> CHANGE: folder with day 3 DiffBind files
  fdr      = 0.1,                             # >>> CHANGE: day 3 FDR cutoff
  distance = 2000,                            # >>> CHANGE: day 3 bp window
  exclude  = "^Mt$"                           # >>> CHANGE: contigs to drop
)

# >>> CHANGE (day 3): right side = exact file name inside cfg_d3$peakDir; left side = label
runs_d3 <- c(
  ZH_AC = "ZH_AC_d3_vs_BF-OUT.csv",
  ZH_ME = "ZH_ME_d3_vs_BF-OUT.csv",
  ZH_SR = "ZH_SR_d3_vs_BF-OUT.csv",
  MP_AC = "MP_AC_d3_vs_BF-OUT.csv",
  MP_ME = "MP_ME_d3_vs_BF-OUT.csv",
  MP_SR = "MP_SR_d3_vs_BF-OUT.csv",
  BF_AC = "BF_AC_d3_vs_SF-OUT.csv",
  BF_ME = "BF_ME_d3_vs_SF-OUT.csv",
  BF_SR = "BF_SR_d3_vs_SF-OUT.csv"
)

# confirm everything exists before running
print(data.frame(run = names(runs_d3),
                 exists = file.exists(file.path(cfg_d3$peakDir, runs_d3))))

results_d3 <- list()
for (run in names(runs_d3)) {
  peakFile <- file.path(cfg_d3$peakDir, runs_d3[[run]])
  contrast <- sub(".*_vs_([A-Za-z]+)-OUT\\.csv$", "\\1", runs_d3[[run]])

  if (!file.exists(peakFile)) {
    warning("Missing peak file, skipping ", run, ": ", peakFile)
    next
  }

  message("Running d3 ", run, " ...")
  results_d3[[run]] <- tryCatch(
    getGeneDistances(
      gtfFile       = cfg_d3$gtf,
      diffPeakFile  = peakFile,
      bedFormat     = FALSE,
      distance      = cfg_d3$distance,
      inGene        = TRUE,
      excludeContig = cfg_d3$exclude,
      featureFilter = "gene",
      fdrFilter     = cfg_d3$fdr,
      sampName = paste0(run, "_d3_vs_", contrast, "_FDR_", cfg_d3$fdr, "_", gtfTag),
      verbose       = FALSE
    ),
    error = function(e) { warning(run, " failed: ", conditionMessage(e)); NULL }
  )
}

data.frame(
  run  = names(runs_d3),
  rows = sapply(names(runs_d3), function(r) if (is.null(results_d3[[r]])) NA else nrow(results_d3[[r]]))
)
