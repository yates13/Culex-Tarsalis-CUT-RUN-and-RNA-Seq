## One PDF: page 1 = DiffBind summary FDR 0.1, page 2 = FDR 0.3 (existing PNGs); then edgeR PCA for all 18 contrasts
## Reads saved RDS only (no DiffBind rerun). Writes to a new folder; refuses to overwrite.
suppressPackageStartupMessages(library(DiffBind))
DB  <- "/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/05_DiffBind_2026-11-24"
OUT <- file.path(DB, "figures_FDR0.1_0.3_2026-11-25")
if (dir.exists(OUT)) stop("Output folder exists, not overwriting: ", OUT)
for (fdr in c("0.1", "0.3")) if (!file.exists(file.path(DB, paste0("DiffBind_summary_FDR", fdr, ".png")))) stop("Missing summary PNG for FDR ", fdr)
dir.create(OUT)
f <- Sys.glob(file.path(DB, c("d1", "d3"), "*_analysis.RDS"))
nm <- sub("-OUT_analysis\\.RDS$", "", basename(f))
f <- f[order(sub(".*_(d[13])_.*", "\\1", nm), sub("^[A-Z]+_([A-Z]+)_.*", "\\1", nm), nm)]
pdf(file.path(OUT, "DiffBind_summary_FDR0.1_0.3_plus_PCA_2026-11-25.pdf"), width = 7, height = 7)
for (fdr in c("0.1", "0.3")) {
  grid::grid.newpage(); grid::grid.raster(png::readPNG(file.path(DB, paste0("DiffBind_summary_FDR", fdr, ".png"))))
  grid::grid.text(paste0("DiffBind summary, FDR <= ", fdr), y = grid::unit(0.985, "npc"), gp = grid::gpar(fontsize = 11, fontface = "bold"))
}
for (x in f) {
  lab <- sub("-OUT_analysis\\.RDS$", "", basename(x))
  ok <- try(dba.plotPCA(readRDS(x), contrast = 1, method = DBA_EDGER, th = 1,
                        attributes = DBA_CONDITION, label = DBA_ID), silent = TRUE)
  if (inherits(ok, "try-error")) { plot.new(); text(0.5, 0.5, paste(lab, "\nPCA failed:", ok)) }
  grid::grid.text(paste(lab, "(edgeR, all sites)"), y = grid::unit(0.985, "npc"), gp = grid::gpar(fontsize = 11, fontface = "bold"))
}
dev.off()
