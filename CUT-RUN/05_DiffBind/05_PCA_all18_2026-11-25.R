################################################################################
# FIGURES (2026-11-25): one PDF with the DiffBind summary (FDR 0.1) as page 1,
#   then one edgeR PCA page per contrast.
#   SUPERSEDED by 05_PCA_all18_FDR01_03_2026-11-25.R (same, plus the FDR 0.3
#   summary); kept because it made figures_2026-11-25/.
#
# Usage:   Rscript 05_PCA_all18_2026-11-25.R
# OUTPUT:  <DB>/figures_2026-11-25/DiffBind_summary_plus_PCA_2026-11-25.pdf
#
# Reads the saved *_analysis.RDS files only (no DiffBind rerun) and redraws the
#   contrast PCA with method = DBA_EDGER (the page missing from every
#   _plots.pdf, see the KNOWN ISSUE in 05.0_DiffBind_Contrasts.R).
# PCA uses th = 1 = every consensus site, so it does NOT depend on the FDR
#   cutoff ([FIXED]: using only significant sites would be circular).
# Every *_analysis.RDS in DB is included (now 24 contrasts, despite "all18" in
#   the name). Order: D1 then D3; within a day AC, ME, SR; then by name.
# A PCA that fails gets a page with the error instead of stopping the run.
# ENV:     diffbind_env (DiffBind, png). Run on a compute node.
#
# USER SETTINGS (lines tagged [EDIT]; list them: grep -n "\[EDIT\]" <this file>)
#   DB          DiffBind folder (05.0 OUT_BASE)
#   OUT         NEW output folder (the script stops if it exists)
#   PDF name    output file name
#   summary PNG folder: the PNG lines read from DB. For the 24-contrast
#               summaries, change file.path(DB, ...) to the 05.1 output folder
#               (e.g. summary_24contrasts_2026-11-25)
#   summary FDR the "0.1" in the PNG name and page title
################################################################################
suppressPackageStartupMessages(library(DiffBind))
DB  <- "/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/05_DiffBind_2026-11-24"   # [EDIT]
OUT <- file.path(DB, "figures_2026-11-25")   # [EDIT] new folder
if (dir.exists(OUT)) stop("Output folder exists, not overwriting: ", OUT)
dir.create(OUT)
f <- Sys.glob(file.path(DB, c("d1", "d3"), "*_analysis.RDS"))
nm <- sub("-OUT_analysis\\.RDS$", "", basename(f))
f <- f[order(sub(".*_(d[13])_.*", "\\1", nm), sub("^[A-Z]+_([A-Z]+)_.*", "\\1", nm), nm)]
pdf(file.path(OUT, "DiffBind_summary_plus_PCA_2026-11-25.pdf"), width = 7, height = 7)   # [EDIT] PDF name
img <- png::readPNG(file.path(DB, "DiffBind_summary_FDR0.1.png"))   # [EDIT] summary PNG folder and FDR
grid::grid.newpage(); grid::grid.raster(img)
grid::grid.text("DiffBind summary, FDR <= 0.1", y = grid::unit(0.985, "npc"), gp = grid::gpar(fontsize = 11, fontface = "bold"))   # [EDIT] page title FDR
for (x in f) {
  lab <- sub("-OUT_analysis\\.RDS$", "", basename(x))
  ok <- try(dba.plotPCA(readRDS(x), contrast = 1, method = DBA_EDGER, th = 1,   # [FIXED] edgeR, th = 1 (all sites)
                        attributes = DBA_CONDITION, label = DBA_ID), silent = TRUE)
  if (inherits(ok, "try-error")) { plot.new(); text(0.5, 0.5, paste(lab, "\nPCA failed:", ok)) }
  grid::grid.text(paste(lab, "(edgeR, all sites)"), y = grid::unit(0.985, "npc"), gp = grid::gpar(fontsize = 11, fontface = "bold"))
}
dev.off()
