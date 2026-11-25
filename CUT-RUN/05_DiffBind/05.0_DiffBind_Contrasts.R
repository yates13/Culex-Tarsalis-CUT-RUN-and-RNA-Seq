################################################################################
# STEP 5.0: DiffBind differential peaks (edgeR), one contrast per run
#
# Usage:   Rscript 05.0_DiffBind_Contrasts.R <contrast number 1-24>
#          Rscript 05.0_DiffBind_Contrasts.R list      # print the 24 contrasts
#          (on Alpine this is run by 05.0_DiffBind_array.sh, one array task each)
#
# Contrasts (12 per day = 3 marks x 4 pairs; built from DAYS x MARKS x PAIRS):
#   for day D1, D3; for mark AC, ME, SR:  ZH vs BF,  MP vs BF,  BF vs SF,  ZH vs MP
#   The first 3 pairs are the ones used in the getGeneDistances runs; ZH vs MP
#   was added 2026-11-25 (commit d744a85).
#   Fold = log2(group1) - log2(group2): positive Fold = higher in the FIRST site.
#   Task numbers follow the order of DAYS, MARKS and PAIRS. Adding or
#   reordering any of them renumbers the tasks: always check with "list".
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
#
# ENV:     diffbind_env (DiffBind 3.20.0, edgeR 4.8.2; see CUT-RUN/envs/)
#
# USER SETTINGS (lines tagged [EDIT]; list them: grep -n "\[EDIT\]" 05.0_DiffBind_Contrasts.R)
#   Paths
#     RUN              run folder; every input and output folder is under it
#     BAM_DIR, INPUT_DIR, POOLED_PEAK_DIR, REP_PEAK_DIR   input folders under RUN
#     OUT_BASE         output folder. For a NEW run use a NEW dated name; a
#                      contrast whose CSV already exists stops (nothing is
#                      overwritten), so reusing a folder only adds new contrasts.
#     NAS block        commented-out paths for running on the NAS (see its note)
#   What is compared
#     DAYS, MARKS      days and marks to run
#     PAIRS            site pairs; the FIRST site of each pair is group 1 (Fold > 0).
#                      Changing DAYS, MARKS or PAIRS changes the number of tasks:
#                      update --array in 05.0_DiffBind_array.sh, and add any new
#                      pair to the two comparison lists in 05.1_DiffBind_Summary.R.
#   Analysis choices. Changing these departs from the original sheets and the
#   Aedes paper (Ogg et al. 2026): decide with the PI first.
#     SUMMITS          FALSE = full peak widths (original scripts); a number
#                      re-centres every peak to +/- that many bp
#     MIN_REPS         fewest replicates a site may have before the contrast stops
#     PEAK_MODE        "pooled" = Step 3.1 pooled peaks (current); "per_replicate"
#                      = Step 3.3 per-replicate peaks (not run yet; this is what
#                      the Aedes paper's GitHub sample sheets use)
#     USE_CONTROLS     FALSE = no input controls (original sheets)
#     peakExt          broadPeak for ME, narrowPeak for AC and SR; must match how
#                      Step 3.1 called each mark
#     peakCaller       "macs" (original sheets)
#   FDR is NOT set here. The CSV keeps every site (th = 1) and the cutoff is
#   applied later: 05.1 (2nd argument), getGeneDistances, Gviz. The paper used
#   FDR <= 0.10. The FDR counts printed at the end of the log are a quick look only.
#   Genome: nothing to set here; DiffBind works only from the BAMs and peak files.
#
# FIXED (lines tagged [FIXED]; change only for a bug or wrong data):
#   dba.normalize(db) default (library size); method = DBA_EDGER;
#   bBlacklist/bGreylist = FALSE (no Culex blacklist); th = 1 in dba.report.
#
# KNOWN ISSUE (not fixed in this script): the 2nd dba.plotPCA call (contrast
#   PCA) has no method = DBA_EDGER, so that page is missing from every
#   _plots.pdf (log: "DESeq2 analysis has not been run for this contrast").
#   The CSVs are not affected. 05_PCA_all18_FDR01_03_2026-11-25.R redraws the
#   PCAs from the saved _analysis.RDS files.
#
# NOTE ON INTERPRETATION (2026-11-25 diagnostics): with pooled per-group peaks
#   the consensus is the union of the two groups' peak sets, and almost all
#   significant sites lie in peaks called in only one group, higher in that
#   group (diag_peak_origin_all18_2026-11-25.R). This pattern fits both real
#   differences and selection bias (Lun & Smyth 2014, NAR 42:e95); see the
#   diagnostics before interpreting large contrasts such as BF vs SF SR D3.
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
RUN        <- "/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates"   # [EDIT]
BAM_DIR    <- file.path(RUN, "02_Cxt_Alignment/02.3.1_Cxt_BamFiltered")   # [EDIT]
INPUT_DIR  <- file.path(RUN, "02_Cxt_Alignment/02.4_Cxt_BamMerged_with_replicate_count")   # [EDIT]
POOLED_PEAK_DIR <- file.path(RUN, "03_MACS2_Peaks_New_Peaks/03.1_MacsPeaks_with_replicate_count")   # [EDIT]
REP_PEAK_DIR    <- file.path(RUN, "03_MACS2_Peaks_New_Peaks/03.3_MacsPeaks_per_replicate")   # [EDIT]
OUT_BASE   <- file.path(RUN, "05_DiffBind_2026-11-24")   # [EDIT] new dated name for a new run

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
DAYS     <- c("D1", "D3")   # [EDIT]
MARKS    <- c("AC", "ME", "SR")   # [EDIT]
PAIRS    <- list(c("ZH", "BF"), c("MP", "BF"), c("BF", "SF"), c("ZH", "MP"))   # [EDIT] first site = group 1 (Fold > 0)
SUMMITS  <- list(AC = FALSE, ME = FALSE, SR = FALSE)   # [EDIT] full peak widths, as in the original DiffBind scripts
MIN_REPS <- 2          # [EDIT] stop if a site has fewer replicates than this
PEAK_MODE <- "pooled"  # [EDIT] "pooled" = Step 3.1 (default); "per_replicate" = Step 3.3 (not run)
USE_CONTROLS <- FALSE  # [EDIT] original sheets had no bamControl; TRUE adds the same-day pooled input

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
peakExt    <- if (ct$mark == "ME") "broadPeak" else "narrowPeak"   # [EDIT] must match Step 3.1 peak type per mark
peakCaller <- "macs"   # [EDIT] as in the original sample sheets

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
db   <- dba.normalize(db)   # [FIXED] library-size normalization (DiffBind default)
db   <- dba.contrast(db, design = "~Condition",
                     contrast = c("Condition", ct$site1, ct$site2))
db   <- dba.analyze(db, method = DBA_EDGER, bBlacklist = FALSE, bGreylist = FALSE)   # [FIXED] edgeR; no blacklist/greylist
print(dba.show(db, bContrasts = TRUE))

res   <- dba.report(db, method = DBA_EDGER, contrast = 1, th = 1)   # [FIXED] th = 1 keeps every site; FDR is filtered later
resDF <- as.data.frame(res)
names(resDF)[names(resDF) == "seqnames"] <- "chr"     # getGeneDistances expects 'chr'
write.csv(resDF, outCSV, row.names = FALSE)
saveRDS(db, file.path(outDir, paste0(ct$name, "_analysis.RDS")))

pdf(file.path(outDir, paste0(ct$name, "_plots.pdf")), width = 7, height = 7)
try(dba.plotPCA(db, attributes = DBA_CONDITION, label = DBA_ID))
try(dba.plotPCA(db, contrast = 1, th = 1, attributes = DBA_CONDITION, label = DBA_ID))   # KNOWN ISSUE: needs method = DBA_EDGER (see header)
try(dba.plotMA(db, contrast = 1, method = DBA_EDGER))
try(dba.plotVolcano(db, contrast = 1, method = DBA_EDGER))
dev.off()

writeLines(capture.output(sessionInfo()), file.path(outDir, paste0(ct$name, "_sessionInfo.txt")))
cat(sprintf("Done: %s  sites %d  FDR<=0.1: %d  FDR<=0.05: %d\n", outCSV, nrow(resDF),   # [EDIT] FDR counts for the log only (not the analysis)
            sum(resDF$FDR <= 0.1), sum(resDF$FDR <= 0.05)))
