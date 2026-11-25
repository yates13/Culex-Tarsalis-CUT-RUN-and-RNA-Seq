## Diagnostics for BF vs SF SR D3 (read-only; outputs to a new folder)
## Run: Rscript diag_BF_SR_d3_vs_SF_2026-11-25.R   (diffbind_env, compute node)
suppressPackageStartupMessages({ library(DiffBind); library(Rsamtools); library(GenomicRanges) })
RUN <- "/scratch/alpine/c832500103@colostate.edu/Cxt_Cut_Run_Pipeline/Round1_Round2_With_Duplicates"
DB  <- file.path(RUN, "05_DiffBind_2026-11-24")
OUT <- file.path(DB, "diagnostics_BF_SR_d3_vs_SF_2026-11-25")
if (dir.exists(OUT)) stop("Output folder exists, not overwriting: ", OUT)
dir.create(OUT)
sink(file.path(OUT, "report.txt"), split = TRUE)
isMt <- function(x) grepl("^Mt_", as.character(x))

## 1. Mito share of mapped reads, all filtered BAMs (idxstats counts reads, not pairs)
cat("== 1. Mt fraction of mapped reads\n")
bams <- Sys.glob(file.path(RUN, "02_Cxt_Alignment/02.3.1_Cxt_BamFiltered/*-Cxt.bam"))
mt <- do.call(rbind, lapply(bams, function(b) {
  s <- idxstatsBam(b); tot <- sum(as.numeric(s$mapped))
  data.frame(sample = sub("-Cxt\\.bam$", "", basename(b)), mapped = tot,
             Mt_frac = sum(as.numeric(s$mapped[isMt(s$seqnames)])) / tot,
             nonTig_nonMt_frac = sum(as.numeric(s$mapped[!isMt(s$seqnames) & !grepl("^tig", s$seqnames)])) / tot)
}))
write.table(mt, file.path(OUT, "Mt_fraction_all_bams.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("seqnames that are not tig*:", paste(grep("^tig", idxstatsBam(bams[1])$seqnames, invert = TRUE, value = TRUE), collapse = ", "), "\n")
mt$group <- sub("-R[0-9]+$", "", mt$sample)
print(aggregate(cbind(Mt_frac, nonTig_nonMt_frac) ~ group, data = mt, FUN = function(v) round(mean(v), 3)), row.names = FALSE)
cat("\nSR-D3 per sample:\n"); print(mt[grepl("-SR-D3-", mt$sample), 1:4], row.names = FALSE)

## 2. Pooled SR-D3 peak sets: size, Mt peaks, width, strength, overlap with BF
cat("\n== 2. Pooled SR-D3 peaks\n")
P <- file.path(RUN, "03_MACS2_Peaks_New_Peaks/03.1_MacsPeaks_with_replicate_count")
readNP <- function(f) { x <- read.table(f, sep = "\t", stringsAsFactors = FALSE)
  GRanges(x$V1, IRanges(x$V2 + 1, x$V3), fe = x$V7, q = x$V9) }
sites <- c("BF", "SF", "MP", "ZH")
pk <- lapply(setNames(sites, sites), function(s)
  readNP(file.path(P, paste0(s, "-SR-D3-rep123"), paste0(s, "-SR-D3-rep123_peaks.narrowPeak"))))
for (n in names(pk)) { g <- pk[[n]]
  cat(sprintf("%s: %d peaks | Mt peaks %d | width q10/50/90 %s | -log10q q10/50/90 %s | FE median %.2f\n",
      n, length(g), sum(isMt(seqnames(g))), paste(quantile(width(g), c(.1, .5, .9)), collapse = "/"),
      paste(round(quantile(g$q, c(.1, .5, .9)), 1), collapse = "/"), median(g$fe))) }
for (n in c("SF", "MP", "ZH"))
  cat(sprintf("BF peaks overlapping %s: %.3f | %s peaks overlapping BF: %.3f\n", n,
      mean(overlapsAny(pk$BF, pk[[n]])), n, mean(overlapsAny(pk[[n]], pk$BF))))

## 3. Saved DiffBind object: normalization, Mt sites, DB sites by peak origin
cat("\n== 3. DiffBind object\n")
db <- readRDS(file.path(DB, "d3/BF_SR_d3_vs_SF-OUT_analysis.RDS"))
print(dba.normalize(db, bRetrieve = TRUE))
rep <- dba.report(db, method = DBA_EDGER, contrast = 1, th = 1, bCounts = TRUE,
                  bNormalized = FALSE, DataType = DBA_DATA_FRAME)
cnt <- rep[, (which(names(rep) == "FDR") + 1):ncol(rep)]
cat("consensus sites:", nrow(rep), "| on Mt:", sum(isMt(rep$Chr)), "| Mt sites FDR<=0.1:", sum(isMt(rep$Chr) & rep$FDR <= 0.1), "\n")
cat("fraction of consensus-site reads on Mt, per sample:\n"); print(round(colSums(cnt[isMt(rep$Chr), , drop = FALSE]) / colSums(cnt), 3))
gr  <- GRanges(rep$Chr, IRanges(rep$Start, rep$End))
inB <- overlapsAny(gr, pk$BF); inS <- overlapsAny(gr, pk$SF)
origin <- ifelse(inB & inS, "both", ifelse(inB, "BF_only", ifelse(inS, "SF_only", "neither")))
cat("\nsites by peak origin x FDR<=0.1 x direction (Fold>0 = higher in BF):\n")
print(ftable(table(origin = origin, sig = rep$FDR <= 0.1, dir = ifelse(rep$Fold > 0, "up_BF", "up_SF"))))

## 4. PCA / correlation with edgeR (preview of the missing-PCA fix; this contrast only)
pdf(file.path(OUT, "PCA_corr_edgeR.pdf"), width = 7, height = 7)
try(dba.plotPCA(db, attributes = DBA_CONDITION, label = DBA_ID))
try(dba.plotPCA(db, contrast = 1, method = DBA_EDGER, th = 1, attributes = DBA_CONDITION, label = DBA_ID))
try(plot(db))
dev.off()
sink(); writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
