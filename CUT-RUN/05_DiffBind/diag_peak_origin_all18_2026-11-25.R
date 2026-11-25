## Peak-origin check for all 18 contrasts (read-only; outputs to a new folder)
suppressPackageStartupMessages({ library(DiffBind); library(GenomicRanges) })
DB  <- "/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates/05_DiffBind_2026-11-24"
OUT <- file.path(DB, "diagnostics_peak_origin_all18_2026-11-25")
if (dir.exists(OUT)) stop("Output folder exists, not overwriting: ", OUT)
dir.create(OUT)
readPk <- function(f) { x <- read.table(f, sep = "\t", stringsAsFactors = FALSE)[, 1:3]; GRanges(x$V1, IRanges(x$V2 + 1, x$V3)) }
res <- do.call(rbind, lapply(Sys.glob(file.path(DB, c("d1", "d3"), "*_analysis.RDS")), function(f) {
  nm <- sub("_analysis\\.RDS$", "", basename(f))
  s1 <- sub("_.*", "", nm); s2 <- sub(".*_vs_([A-Z]+)-OUT$", "\\1", nm)
  sh <- read.csv(sub("_analysis\\.RDS$", "_samplesheet.csv", f), stringsAsFactors = FALSE)
  p1 <- readPk(unique(sh$Peaks[sh$Condition == s1])); p2 <- readPk(unique(sh$Peaks[sh$Condition == s2]))
  r  <- dba.report(readRDS(f), method = DBA_EDGER, contrast = 1, th = 1, DataType = DBA_DATA_FRAME)
  gr <- GRanges(r$Chr, IRanges(r$Start, r$End))
  in1 <- suppressWarnings(overlapsAny(gr, p1)); in2 <- suppressWarnings(overlapsAny(gr, p2))
  sig <- r$FDR <= 0.1; only1 <- in1 & !in2; only2 <- in2 & !in1
  data.frame(contrast = nm, sites = nrow(r),
    ov_1in2 = round(mean(suppressWarnings(overlapsAny(p1, p2))), 3), ov_2in1 = round(mean(suppressWarnings(overlapsAny(p2, p1))), 3),
    sig = sum(sig), sig_only1 = sum(sig & only1), sig_both = sum(sig & in1 & in2), sig_only2 = sum(sig & only2),
    sig_dir_matches_origin = round(sum(sig & ((only1 & r$Fold > 0) | (only2 & r$Fold < 0))) / max(1, sum(sig & (only1 | only2))), 3))
}))
write.table(res, file.path(OUT, "peak_origin_all18.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
