################################################################################
# STEP 5.1: Summarize the DiffBind contrasts (table + figure)
#
# Usage:   Rscript 05.1_DiffBind_Summary.R <DiffBind output folder> [FDR cutoff, default 0.1] [output folder, default = input folder]
#          e.g. D=.../05_DiffBind_2026-11-24; S=$D/summary_24contrasts_2026-11-25
#               mkdir $S && Rscript 05.1_DiffBind_Summary.R $D 0.1 $S && Rscript 05.1_DiffBind_Summary.R $D 0.3 $S
#
# Reads:   <folder>/d1/*-OUT.csv and <folder>/d3/*-OUT.csv (Step 5.0 output, th = 1)
#          File names must look like <SITE1>_<MARK>_<day>_vs_<SITE2>-OUT.csv
# Writes (into the output folder, which must already exist; stops if the .tsv
#   for this FDR is already there, so nothing is overwritten):
#          DiffBind_summary_FDR<cut>.tsv   one row per contrast: sites tested, significant,
#                                          higher in site 1, higher in site 2, % significant
#          DiffBind_summary_FDR<cut>.pdf / .png   diverging bars per contrast,
#                                          rows = day, columns = mark
# Fold > 0 = higher in the first site (DiffBind: group1 - group2).
# Significant = FDR <= cutoff. The paper used FDR <= 0.10 (Ogg et al. 2026, Fig. 3B).
# ENV:     diffbind_env (ggplot2, scales)
#
# USER SETTINGS
#   On the command line (no editing needed):
#     arg 1  input folder = OUT_BASE of 05.0_DiffBind_Contrasts.R
#     arg 2  FDR cutoff (0.1 = paper; 0.3 was used only as an exploratory view)
#     arg 3  output folder; use a new dated folder (default = input folder)
#   Code lines tagged [EDIT] (list them: grep -n "\[EDIT\]" 05.1_DiffBind_Summary.R):
#     FDR_CUT default      used only when arg 2 is not given
#     comparison lists     TWO lists (table row order, figure order). Every pair
#                          in 05.0 PAIRS must appear in both, written exactly
#                          "SITE1 vs SITE2"; a pair missing from the figure list
#                          becomes NA and its bars are lost.
#     mark order/labels    facet columns and their titles (must match 05.0 MARKS)
#     y-axis breaks        pseudo-log ticks; the top tick is 10,000 (BF vs SF SR D3
#                          already goes past it, 23,962); add +/-100000 to show more
#     colours, titles      figure look only
#     figure size, dpi     ggsave width/height (inches) and PNG resolution
#   Genome: nothing to set here.
################################################################################

suppressPackageStartupMessages(library(ggplot2))

args    <- commandArgs(trailingOnly = TRUE)
DIR     <- if (length(args) >= 1) args[1] else "."
FDR_CUT <- if (length(args) >= 2) as.numeric(args[2]) else 0.1   # [EDIT] default FDR cutoff if no 2nd argument
OUT_DIR <- if (length(args) >= 3) args[3] else DIR

files <- Sys.glob(file.path(DIR, c("d1", "d3"), "*-OUT.csv"))
if (length(files) == 0) stop("No *-OUT.csv files under ", DIR, "/d1 or /d3")

rows <- lapply(files, function(f) {
  nm <- sub("-OUT\\.csv$", "", basename(f))
  p  <- regmatches(nm, regexec("^([A-Z]+)_([A-Z]+)_(d[0-9])_vs_([A-Z]+)$", nm))[[1]]
  if (length(p) != 5) stop("Unexpected file name: ", basename(f))
  r   <- read.csv(f)
  sig <- !is.na(r$FDR) & r$FDR <= FDR_CUT
  data.frame(contrast = nm, site1 = p[2], mark = p[3], day = toupper(p[4]), site2 = p[5],
             sites = nrow(r), sig = sum(sig),
             higher_site1 = sum(sig & r$Fold > 0), higher_site2 = sum(sig & r$Fold < 0),
             pct_sig = round(100 * sum(sig) / nrow(r), 2), stringsAsFactors = FALSE)
})
s <- do.call(rbind, rows)
s$comparison <- paste(s$site1, "vs", s$site2)
s <- s[order(s$day, s$mark, match(s$comparison, c("ZH vs BF", "MP vs BF", "BF vs SF", "ZH vs MP"))), ]   # [EDIT] table row order of the pairs

out_tsv <- file.path(OUT_DIR, sprintf("DiffBind_summary_FDR%s.tsv", FDR_CUT))
if (file.exists(out_tsv)) stop("Output already exists, nothing overwritten: ", out_tsv)
write.table(s, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
print(s[, c("contrast", "sites", "sig", "higher_site1", "higher_site2", "pct_sig")], row.names = FALSE)

long <- rbind(
  data.frame(s[, c("comparison", "mark", "day")], direction = "Higher in first site",  n =  s$higher_site1),
  data.frame(s[, c("comparison", "mark", "day")], direction = "Higher in second site", n = -s$higher_site2))
long$comparison <- factor(long$comparison, levels = c("ZH vs BF", "MP vs BF", "BF vs SF", "ZH vs MP"))   # [EDIT] figure order of the pairs
long$mark <- factor(long$mark, levels = c("AC", "ME", "SR"),   # [EDIT] mark order
                    labels = c("H3K27ac (AC)", "H3K9me3 (ME)", "SREBP (SR)"))   # [EDIT] mark titles
long$lab  <- ifelse(long$n == 0, "", format(abs(long$n), big.mark = ",", trim = TRUE))
long$vj   <- ifelse(long$n > 0, -0.3, 1.3)
zeros <- unique(long[ave(abs(long$n), long$comparison, long$mark, long$day, FUN = sum) == 0, c("comparison", "mark", "day")])   # contrasts tested with 0 sites

p <- ggplot(long, aes(x = comparison, y = n, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(aes(label = lab, vjust = vj), size = 2.6) +
  geom_text(data = zeros, aes(x = comparison, y = 0, label = "0"), inherit.aes = FALSE, vjust = -0.4, size = 2.6, colour = "grey40") +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(-10000, -1000, -100, -10, 0, 10, 100, 1000, 10000), expand = expansion(mult = 0.12),   # [EDIT] y-axis ticks
                     labels = function(b) format(abs(b), big.mark = ",", scientific = FALSE, trim = TRUE)) +
  facet_grid(day ~ mark) +
  scale_fill_manual(values = c("Higher in first site" = "#3B6FB6", "Higher in second site" = "#A8A8A8")) +   # [EDIT] colours
  labs(x = NULL, fill = NULL,
       y = sprintf("Differential sites, FDR <= %s (pseudo-log scale)", FDR_CUT),
       title = "DiffBind (edgeR) differential sites per contrast",   # [EDIT] title
       subtitle = "Up = higher in the first site of the pair; down = higher in the second") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUT_DIR, sprintf("DiffBind_summary_FDR%s.pdf", FDR_CUT)), p, width = 9, height = 6)   # [EDIT] size
ggsave(file.path(OUT_DIR, sprintf("DiffBind_summary_FDR%s.png", FDR_CUT)), p, width = 9, height = 6, dpi = 200)   # [EDIT] size, dpi
cat("Wrote:", out_tsv, "and the .pdf/.png next to it\n")
