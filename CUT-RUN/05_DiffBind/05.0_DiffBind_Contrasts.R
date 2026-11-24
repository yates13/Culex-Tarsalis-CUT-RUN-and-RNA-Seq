################################################################################
# STEP 5.0: DiffBind differential peaks (edgeR), one contrast per run
#
# Usage:   Rscript 05.0_DiffBind_Contrasts.R <contrast number 1-18>
#          Rscript 05.0_DiffBind_Contrasts.R list      # print the 18 contrasts
#          (on Alpine this is run by 05.0_DiffBind_array.sh, one array task each)
#
# Contrasts (same 9 per day as the getGeneDistances runs):
#   for day D1, D3; for mark AC, ME, SR:  ZH vs BF,  MP vs BF,  BF vs SF
#   Fold = log2(group1) - log2(group2): positive Fold = higher in the FIRST site.
#
# Inputs (built automatically from the file names; nothing typed by hand):
#   bamReads   filtered BAMs, Step 2.3.1 (MAPQ >= 10, proper pairs, primary)
#   Peaks      PEAK_MODE = "pooled" (default): the Step 3.1 pooled-replicate peaks
#              (one file per site/mark/day, used for all of that condition's
#              replicates), chosen to stay close to the Aedes paper. "per_replicate" uses 3.3
#              (not run yet) if it is ever needed.
#   bamControl not used (USE_CONTROLS = FALSE), as in the original sample sheets
#   The contrast is chosen by the Condition column (site), not by row number.
#
# Output (NEW folder per day; a contrast stops if its CSV already exists):
#   <OUT_BASE>/d1/ and /d3/
#     <SITE1>_<MARK>_<day>_vs_<SITE2>-OUT.csv   all sites (th = 1); filter FDR later
#     <...>_samplesheet.csv                     the exact sheet used
#     <...>_analysis.RDS                        the analyzed DBA object
#     <...>_plots.pdf                           PCA (all, contrast), MA, volcano
#     <...>_sessionInfo.txt                     R and package versions
#
# Settings (see SETTINGS below): edgeR; library-size normalization (DiffBind
#   default); no blacklist/greylist (no blacklist exists for Culex; matches the
#   Aedes paper); no input controls in the sheet (as in the original scripts);
#   summits: FALSE for all marks (full peak widths, as in the original scripts);
#   ME full peak width (broad domains).
################################################################################

suppressPackageStartupMessages({
  library(BiocParallel)
  library(DiffBind)
  library(edgeR)
})
register(SerialParam())

# ==============================================================================
# PATHS: ALPINE (active)
# ==============================================================================
RUN        <- "/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates"
BAM_DIR    <- file.path(RUN, "02_Cxt_Alignment/02.3.1_Cxt_BamFiltered")
INPUT_DIR  <- file.path(RUN, "02_Cxt_Alignment/02.4_Cxt_BamMerged_with_replicate_count")
POOLED_PEAK_DIR <- file.path(RUN, "03_MACS2_Peaks_New_Peaks/03.1_MacsPeaks_with_replicate_count")
REP_PEAK_DIR    <- file.path(RUN, "03_MACS2_Peaks_New_Peaks/03.3_MacsPeaks_per_replicate")
OUT_BASE   <- file.path(RUN, "05_DiffBind_2026-11-24")

# ==============================================================================
# PATHS: NAS (for later; commented out)
# To run on the NAS: copy the SAME filtered BAMs (+ .bai), the pooled inputs
# (+ .bai) and the Step 3.1 peak folders there, keep the same folder layout,
# then comment out the Alpine block above and uncomment this one.
# Don't point these at the old "with duplicates" BAMs: those are unfiltered.
# ==============================================================================
# RUN        <- "Z:/NAS_Server_Shared/08-Cxt_RVFV"
# BAM_DIR    <- file.path(RUN, "bam/02.3.1_Cxt_BamFiltered")
# INPUT_DIR  <- file.path(RUN, "bam/02.4_Cxt_BamMerged_with_replicate_count")
# POOLED_PEAK_DIR <- file.path(RUN, "peaks/03.1_MacsPeaks_with_replicate_count")
# REP_PEAK_DIR    <- file.path(RUN, "peaks/03.3_MacsPeaks_per_replicate")
# OUT_BASE   <- file.path(RUN, "diffbind_2026-11-24")

# ==============================================================================
# SETTINGS
# ==============================================================================
DAYS     <- c("D1", "D3")
MARKS    <- c("AC", "ME", "SR")
PAIRS    <- list(c("ZH", "BF"), c("MP", "BF"), c("BF", "SF"))
SUMMITS  <- list(AC = FALSE, ME = FALSE, SR = FALSE)   # full peak widths, as in the original DiffBind scripts
MIN_REPS <- 2          # stop if a site has fewer replicates than this
PEAK_MODE <- "pooled"  # "pooled" = Step 3.1 (default); "per_replicate" = Step 3.3 (not run)
USE_CONTROLS <- FALSE  # original sheets had no bamControl; TRUE adds the same-day pooled input

contrasts <- do.call(rbind, lapply(DAYS, function(d)
  do.call(rbind, lapply(MARKS, function(m)
    do.call(rbind, lapply(PAIRS, function(p)
      data.frame(day = d, mark = m, site1 = p[1], site2 = p[2], stringsAsFactors = FALSE)))))))
contrasts$name <- sprintf("%s_%s_%s_vs_%s-OUT", contrasts$site1, contrasts$mark,
                          tolower(contrasts$day), contrasts$site2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript 05.0_DiffBind_Contrasts.R <1-", nrow(contrasts), " | list>")
if (args[1] == "list") {
  print(cbind(n = seq_len(nrow(contrasts)), contrasts))
  quit(save = "no", status = 0)
}
k <- as.integer(args[1])
if (is.na(k) || k < 1 || k > nrow(contrasts)) stop("Contrast number must be 1-", nrow(contrasts))
ct <- contrasts[k, ]
cat(sprintf("Contrast %d: %s  (%s %s: %s vs %s)\n", k, ct$name, ct$mark, ct$day, ct$site1, ct$site2))

outDir <- file.path(OUT_BASE, tolower(ct$day))
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
outCSV <- file.path(outDir, paste0(ct$name, ".csv"))
if (file.exists(outCSV)) stop("Output already exists, nothing overwritten: ", outCSV)

# ------------------------------------------------------------------------------
# Build the sample sheet from the files that exist, and check every path
# ------------------------------------------------------------------------------
peakExt    <- if (ct$mark == "ME") "broadPeak" else "narrowPeak"
peakCaller <- "macs"   # as in the original sample sheets

rowsFor <- function(site) {
  pat  <- sprintf("^%s-%s-%s-R[0-9]+-Cxt\\.bam$", site, ct$mark, ct$day)
  bams <- sort(list.files(BAM_DIR, pattern = pat, full.names = TRUE))
  if (length(bams) < MIN_REPS) stop("Only ", length(bams), " BAMs for ", site, "-", ct$mark, "-", ct$day)
  if (length(bams) != 3) warning(length(bams), " replicates for ", site, "-", ct$mark, "-", ct$day)
  samp  <- sub("-Cxt\\.bam$", "", basename(bams))
  ctrl  <- list.files(INPUT_DIR, pattern = sprintf("^%s-IN-%s-rep[0-9]+-merged-sorted\\.bam$", site, ct$day),
                      full.names = TRUE)
  if (length(ctrl) != 1) stop("Expected 1 pooled ", site, "-IN-", ct$day, " input, found ", length(ctrl))
  if (PEAK_MODE == "pooled") {
    # one Step 3.1 folder per condition, e.g. ZH-AC-D1-rep123 (rep23 etc. if a replicate was dropped)
    cond <- list.files(POOLED_PEAK_DIR, pattern = sprintf("^%s-%s-%s-rep[0-9]+$", site, ct$mark, ct$day))
    if (length(cond) != 1) stop("Expected 1 pooled peak folder for ", site, "-", ct$mark, "-", ct$day,
                                ", found ", length(cond), ": ", paste(cond, collapse = " "))
    peaks <- rep(file.path(POOLED_PEAK_DIR, cond, paste0(cond, "_peaks.", peakExt)), length(samp))
  } else if (PEAK_MODE == "per_replicate") {
    peaks <- file.path(REP_PEAK_DIR, samp, paste0(samp, "_peaks.", peakExt))
  } else stop("PEAK_MODE must be 'pooled' or 'per_replicate'")
  data.frame(SampleID   = samp,
             Factor     = ct$mark,
             Condition  = site,
             Treatment  = ct$day,
             Replicate  = sub(".*-R([0-9]+)$", "\\1", samp),
             bamReads   = bams,
             ControlID  = paste0(site, "-IN-", ct$day),
             bamControl = ctrl,
             Peaks      = peaks,
             PeakCaller = peakCaller,
             stringsAsFactors = FALSE)
}
sheet <- rbind(rowsFor(ct$site1), rowsFor(ct$site2))

need <- c(sheet$bamReads, paste0(sheet$bamReads, ".bai"),
          unique(sheet$bamControl), paste0(unique(sheet$bamControl), ".bai"), unique(sheet$Peaks))
missing <- need[!file.exists(need)]
if (length(missing) > 0) stop("Missing files:\n  ", paste(missing, collapse = "\n  "))
# Same-day check: every control must be from the sample's own day
stopifnot(all(grepl(paste0("-IN-", ct$day, "-"), basename(sheet$bamControl))))
if (!USE_CONTROLS) sheet <- sheet[, setdiff(names(sheet), c("ControlID", "bamControl"))]   # original sheets had no controls

write.csv(sheet, file.path(outDir, paste0(ct$name, "_samplesheet.csv")), row.names = FALSE)
print(sheet[, c("SampleID", "Condition", "Replicate")])

# ------------------------------------------------------------------------------
# DiffBind
# ------------------------------------------------------------------------------
db   <- dba(sampleSheet = sheet)
db   <- dba.count(db, summits = SUMMITS[[ct$mark]], bParallel = FALSE)
db   <- dba.normalize(db)
db   <- dba.contrast(db, design = "~Condition",
                     contrast = c("Condition", ct$site1, ct$site2))
db   <- dba.analyze(db, method = DBA_EDGER, bBlacklist = FALSE, bGreylist = FALSE)
print(dba.show(db, bContrasts = TRUE))

res   <- dba.report(db, method = DBA_EDGER, contrast = 1, th = 1)
resDF <- as.data.frame(res)
names(resDF)[names(resDF) == "seqnames"] <- "chr"     # getGeneDistances expects 'chr'
write.csv(resDF, outCSV, row.names = FALSE)
saveRDS(db, file.path(outDir, paste0(ct$name, "_analysis.RDS")))

pdf(file.path(outDir, paste0(ct$name, "_plots.pdf")), width = 7, height = 7)
try(dba.plotPCA(db, attributes = DBA_CONDITION, label = DBA_ID))
try(dba.plotPCA(db, contrast = 1, th = 1, attributes = DBA_CONDITION, label = DBA_ID))
try(dba.plotMA(db, contrast = 1, method = DBA_EDGER))
try(dba.plotVolcano(db, contrast = 1, method = DBA_EDGER))
dev.off()

writeLines(capture.output(sessionInfo()), file.path(outDir, paste0(ct$name, "_sessionInfo.txt")))
cat(sprintf("Done: %s  sites %d  FDR<=0.1: %d  FDR<=0.05: %d\n", outCSV, nrow(resDF),
            sum(resDF$FDR <= 0.1), sum(resDF$FDR <= 0.05)))
