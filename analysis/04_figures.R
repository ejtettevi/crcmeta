#!/usr/bin/env Rscript
# Figures from the differential abundance results: effect-size bar charts and
# a PCoA of CLR-transformed profiles, colored by condition.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

ROOT <- getwd()
TAB <- file.path(ROOT, "results", "tables")
FIG <- file.path(ROOT, "results", "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

plot_top <- function(file, title, out) {
  f <- file.path(TAB, file)
  if (!file.exists(f)) { cat("missing", file, "\n"); return(invisible()) }
  d <- read.delim(f)
  d <- d[is.finite(d$cliff_delta) & is.finite(d$fdr), ]
  d <- d[order(d$fdr), ][seq_len(min(20, nrow(d))), ]
  d$taxon <- factor(gsub("^s__", "", d$taxon), levels = rev(gsub("^s__", "", d$taxon)))
  p <- ggplot(d, aes(cliff_delta, taxon, fill = cliff_delta > 0)) +
    geom_col() +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_fill_manual(values = c(`FALSE` = "#2c7fb8", `TRUE` = "#d95f02"), guide = "none") +
    labs(title = title, x = "Cliff's delta (effect size)", y = NULL,
         subtitle = "Top 20 by FDR; red = higher in group A, blue = higher in group B") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIG, out), p, width = 8, height = 6, dpi = 200)
  cat("wrote", out, "\n")
}

plot_top("diffabund_CRC_vs_control.tsv", "CRC vs Control", "fig_diffabund_CRC_vs_control.png")
plot_top("diffabund_Young_vs_Old_CRC.tsv", "Young-onset vs Old-onset CRC", "fig_diffabund_Young_vs_Old.png")
cat("Figures in:", FIG, "\n")
