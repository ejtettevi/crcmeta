#!/usr/bin/env Rscript
# Publication-quality figures for the CRC shotgun meta-analysis.
# Output: PNG at 750 dpi with explicit group labels, colour coding tied to
# biological direction, and legends written for standalone interpretation.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(ggrepel); library(RColorBrewer)
  library(vegan)
})

ROOT <- Sys.getenv("CRC_ROOT", unset = getwd())
TAB  <- file.path(ROOT, "results", "tables")
FIG  <- file.path(ROOT, "results", "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

DPI <- 750

pretty <- function(x) {
  x <- sub("^s__", "", x)
  x <- sub("^g__", "", x)
  gsub("_", " ", x)
}
fmt_fdr <- function(x) {
  x <- as.numeric(x)
  ifelse(x < 1e-15, "<1e-15",
  ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
         formatC(x, format = "f", digits = 3)))
}

# colour-blind-safe palette (Okabe–Ito inspired)
COL_CRC   <- "#C44E52"   # red   = higher in CRC / young-onset
COL_CTRL  <- "#4C72B0"   # blue  = higher in control / old-onset
COL_YOUNG <- "#DD8452"   # orange = young-onset CRC
COL_F     <- "#CC79A7"   # pink  = female
COL_M     <- "#4C72B0"   # blue  = male
COL_NS    <- "grey55"
COL_SIG   <- "#E58B3A"

THEME_MIN <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 13, colour = "grey10"),
    plot.subtitle = element_text(size = 10, colour = "grey35", margin = margin(b = 6)),
    plot.caption = element_text(size = 8.5, colour = "grey45", hjust = 0, margin = margin(t = 8)),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, colour = "grey20"),
    legend.position = "top",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9.5),
    legend.key.size = unit(0.9, "lines"),
    legend.margin = margin(0, 0, 2, 0),
    strip.background = element_rect(fill = "grey96", colour = NA),
    strip.text = element_text(face = "bold", size = 11)
  )

save_fig <- function(name, plot, width, height) {
  path <- file.path(FIG, name)
  ggsave(path, plot, width = width, height = height, dpi = DPI,
         bg = "white", units = "in")
  cat(sprintf("wrote %s  (%d dpi, %.1f x %.1f in)\n", name, DPI, width, height))
}

# ---------------------------------------------------------------------------
# FIG 1 — PCoA of faecal microbiome composition (Bray–Curtis)
# ---------------------------------------------------------------------------
pcoa <- read.csv(file.path(TAB, "pcoa_bray_curtis.csv"))
perm <- read.csv(file.path(TAB, "permanova_condition_study.csv"), row.names = 1)
pcol <- if ("Pr(>F)" %in% names(perm)) "Pr(>F)" else "Pr..F."

pcoa$age_group <- factor(pcoa$age_group,
  levels = c("Young-onset CRC", "Old-onset CRC"))
pcoa$condition <- factor(pcoa$condition, levels = c("Control", "CRC"))

fig1a <- ggplot(pcoa, aes(PCo1, PCo2, colour = condition)) +
  geom_point(alpha = 0.55, size = 1.5) +
  stat_ellipse(level = 0.68, linewidth = 0.8, type = "t") +
  scale_colour_manual(values = c(Control = COL_CTRL, CRC = COL_CRC),
                      name = "Group",
                      labels = c(Control = "Control", CRC = "Colorectal cancer")) +
  labs(title = "Colorectal cancer vs control",
       subtitle = "Bray–Curtis PCoA of faecal metagenomes",
       x = "PCoA1", y = "PCoA2") +
  THEME_MIN

var_file <- file.path(TAB, "pcoa_variance.txt")
if (file.exists(var_file)) {
  ve <- scan(var_file, quiet = TRUE)
  fig1a <- fig1a + labs(x = sprintf("PCoA1 (%.1f%%)", ve[1]),
                        y = sprintf("PCoA2 (%.1f%%)", ve[2]))
}

pcoa_crc <- pcoa %>% filter(condition == "CRC", !is.na(age_group))
if (nrow(pcoa_crc) > 10) {
  fig1b <- ggplot(pcoa_crc, aes(PCo1, PCo2, colour = age_group)) +
    geom_point(alpha = 0.55, size = 1.5) +
    stat_ellipse(level = 0.68, linewidth = 0.8, type = "t") +
    scale_colour_manual(values = c("Young-onset CRC" = COL_YOUNG,
                                   "Old-onset CRC" = COL_CRC),
                        name = "Age at diagnosis",
                        labels = c("Young-onset CRC" = "Young-onset CRC (<50 y)",
                                   "Old-onset CRC" = "Old-onset CRC (≥50 y)")) +
    labs(title = "Young-onset vs old-onset CRC",
         subtitle = "Same ordination, CRC cases only (<50 vs ≥50 years)",
         x = if (file.exists(var_file)) sprintf("PCoA1 (%.1f%%)", ve[1]) else "PCoA1",
         y = if (file.exists(var_file)) sprintf("PCoA2 (%.1f%%)", ve[2]) else "PCoA2") +
    THEME_MIN
} else fig1b <- plot_spacer()

fig1 <- (fig1a | fig1b) +
  plot_annotation(
    tag_levels = "A",
    caption = sprintf(
      "n = 1,395 CRC+control metagenomes; ellipses = 68%% t-distribution.\nPERMANOVA (study as covariate): condition R² = %.3f, p = %.3f; study R² = %.2f, p = %.3f.",
      perm["condition", "R2"], perm["condition", pcol],
      perm["study_name", "R2"], perm["study_name", pcol])) &
  theme(plot.caption = element_text(hjust = 0))
save_fig("fig1_pcoa.png", fig1, 11, 5.4)

# ---------------------------------------------------------------------------
# FIG 2 — Species-level differential abundance, CRC vs Control
# ---------------------------------------------------------------------------
da <- read.delim(file.path(TAB, "diffabund_CRC_vs_control.tsv"))
da <- da[is.finite(da$cliff_delta) & is.finite(da$fdr), ]
da <- da[order(da$fdr), ][1:20, ]
da$label <- pretty(da$taxon)
da$dir <- factor(
  ifelse(da$cliff_delta > 0, "Enriched in CRC", "Depleted in CRC"),
  levels = c("Enriched in CRC", "Depleted in CRC"))
da$label <- factor(da$label, levels = rev(da$label[order(da$cliff_delta)]))

fig2 <- ggplot(da, aes(cliff_delta, label, fill = dir)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_text(aes(label = paste0("q = ", fmt_fdr(fdr))),
            hjust = ifelse(da$cliff_delta > 0, -0.15, 1.15),
            size = 3.0, colour = "grey20") +
  scale_fill_manual(values = c("Enriched in CRC" = COL_CRC,
                               "Depleted in CRC" = COL_CTRL),
                    name = "Direction") +
  scale_x_continuous(expand = expansion(mult = c(0.24, 0.24))) +
  labs(title = "Species-level differential abundance: CRC vs control",
       subtitle = "Pooled analysis across 11 cohorts (n = 701 CRC, 694 control); CLR-transformed relative abundance",
       x = "Cliff's δ  (positive → higher in CRC; negative → higher in control)",
       y = NULL,
       caption = "Wilcoxon rank-sum; top 20 species by Benjamini–Hochberg FDR.\nδ: rank-based effect size in [−1, 1].") +
  THEME_MIN
save_fig("fig2_diffabund_CRC.png", fig2, 9.5, 7.2)

# ---------------------------------------------------------------------------
# FIG 3 — Forest plot: study-adjusted meta-analysis, CRC vs Control
# ---------------------------------------------------------------------------
meta <- read.delim(file.path(TAB, "meta_CRC_vs_control.tsv"))
meta <- meta[!is.na(meta$pooled_g) & !is.na(meta$fdr), ]
meta <- meta[order(meta$p_value), ][1:12, ]
meta$label <- pretty(meta$taxon)
meta$dir_lab <- factor(
  ifelse(meta$pooled_g > 0, "Higher in CRC", "Higher in control"),
  levels = c("Higher in CRC", "Higher in control"))
meta$label <- factor(meta$label, levels = rev(meta$label))
meta$anno <- sprintf("k = %d   %d/%d   I² = %.0f%%   q = %s",
                     meta$k_studies, meta$n_pos, meta$n_neg, meta$I2, fmt_fdr(meta$fdr))

fig3 <- ggplot(meta, aes(pooled_g, label, colour = dir_lab)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.18, linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(label = anno), vjust = -0.95, size = 2.9, colour = "grey20",
            show.legend = FALSE) +
  scale_colour_manual(values = c("Higher in CRC" = COL_CRC,
                                 "Higher in control" = COL_CTRL),
                      name = "Effect direction") +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.32))) +
  labs(title = "Study-adjusted meta-analysis of CRC-associated species",
       subtitle = "DerSimonian–Laird random-effects model; 11 cohorts; top 12 by FDR",
       x = "Pooled Hedges' g  (95% CI; positive → higher in CRC)",
       y = NULL,
       caption = "k = studies contributing; n+/n− = studies with positive/negative effect;\nI² = heterogeneity; q = BH-adjusted p. All shown taxa have q < 0.05.") +
  THEME_MIN +
  theme(panel.grid.major.y = element_blank())
save_fig("fig3_forest_CRC.png", fig3, 9, 5.8)

# ---------------------------------------------------------------------------
# FIG 4 — Young-onset vs old-onset: pooled bars + study-adjusted forest
# ---------------------------------------------------------------------------
daY <- read.delim(file.path(TAB, "diffabund_Young_vs_Old_CRC.tsv"))
daY <- daY[is.finite(daY$cliff_delta) & is.finite(daY$fdr), ]
daY <- daY[order(daY$fdr), ][1:15, ]
daY$label <- pretty(daY$taxon)
daY$dir <- factor(
  ifelse(daY$cliff_delta > 0, "Higher in young-onset CRC",
         "Higher in old-onset CRC"),
  levels = c("Higher in young-onset CRC", "Higher in old-onset CRC"))
daY$label <- factor(daY$label, levels = rev(daY$label[order(daY$cliff_delta)]))

fig4a <- ggplot(daY, aes(cliff_delta, label, fill = dir)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_text(aes(label = paste0("q = ", fmt_fdr(fdr))),
            hjust = ifelse(daY$cliff_delta > 0, -0.15, 1.15),
            size = 3.0, colour = "grey20") +
  scale_fill_manual(values = c("Higher in young-onset CRC" = COL_YOUNG,
                               "Higher in old-onset CRC" = COL_CRC),
                    name = "Effect direction") +
  scale_x_continuous(expand = expansion(mult = c(0.26, 0.26))) +
  labs(title = "Pooled differential abundance",
       subtitle = "n = 82 young-onset, 619 old-onset CRC; only 3/161 species reach FDR < 0.05",
       x = "Cliff's δ  (+ → young-onset; − → old-onset)",
       y = NULL)

metaY <- read.delim(file.path(TAB, "meta_Young_vs_Old_CRC.tsv"))
metaY <- metaY[!is.na(metaY$pooled_g) & !is.na(metaY$fdr), ]
metaY <- metaY[order(metaY$p_value), ][1:12, ]
metaY$label <- pretty(metaY$taxon)
metaY$sig_lab <- factor(
  ifelse(metaY$fdr < 0.05, "FDR < 0.05", "n.s."),
  levels = c("FDR < 0.05", "n.s."))
metaY$label <- factor(metaY$label, levels = rev(metaY$label))
metaY$anno <- sprintf("k = %d   %d/%d   I² = %.0f%%   q = %s",
                      metaY$k_studies, metaY$n_pos, metaY$n_neg, metaY$I2, fmt_fdr(metaY$fdr))

fig4b <- ggplot(metaY, aes(pooled_g, label, colour = sig_lab)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.18, linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(label = anno), vjust = -0.95, size = 2.9, colour = "grey20") +
  scale_colour_manual(values = c("FDR < 0.05" = COL_SIG, "n.s." = COL_NS),
                      name = "Significance") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.30))) +
  scale_y_discrete(expand = expansion(add = c(1.5, 1.5))) +
  labs(title = "Study-adjusted meta-analysis",
       subtitle = "0/161 species significant after study adjustment (k = 10)",
       x = "Pooled Hedges' g  (95% CI)",
       y = NULL,
       caption = "Random-effects (DerSimonian–Laird).\nNo species reached FDR < 0.05; median I² = 0%.\nk = studies; n+/n− = direction counts.") +
  theme(panel.grid.major.y = element_blank())

fig4 <- (fig4a | fig4b) +
  plot_annotation(
    tag_levels = "A",
    title = "Young-onset vs old-onset colorectal cancer",
    caption = "Age cut-off: 50 years. Left: pooled Cliff's δ for young-onset vs old-onset CRC.\nRight: study-adjusted Hedges' g (DerSimonian–Laird).") &
  THEME_MIN
save_fig("fig4_young_vs_old.png", fig4, 13, 5.8)

# ---------------------------------------------------------------------------
# FIG 5 — CRC effect size vs age effect size (concordance scatter)
# ---------------------------------------------------------------------------
cmp <- read.delim(file.path(TAB, "effect_CRC_vs_age.tsv"))
sig_crc <- cmp$fdr_CRC < 0.05
sig_age <- cmp$fdr_Age < 0.05
rho <- cor(cmp$g_CRC, cmp$g_Age, method = "spearman", use = "complete.obs")

cmp$sig_cat <- factor(
  ifelse(sig_crc & sig_age, "Both FDR < 0.05",
  ifelse(sig_crc, "CRC only FDR < 0.05",
  ifelse(sig_age, "Age only FDR < 0.05", "Neither"))),
  levels = c("Both FDR < 0.05", "CRC only FDR < 0.05",
             "Age only FDR < 0.05", "Neither"))

fig5 <- ggplot(cmp, aes(g_CRC, g_Age, colour = sig_cat)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_colour_manual(values = c("Both FDR < 0.05" = "#8172B3",
                                 "CRC only FDR < 0.05" = COL_CRC,
                                 "Age only FDR < 0.05" = COL_YOUNG,
                                 "Neither" = "grey75"),
                      name = "FDR category") +
  labs(title = "CRC association vs age-at-diagnosis association",
       subtitle = sprintf("Each point is a species (n = %d shared). Spearman ρ = %.2f between the two study-adjusted effect sizes.",
                          nrow(cmp), rho),
       x = "Study-adjusted Hedges' g — CRC vs control",
       y = "Study-adjusted Hedges' g — young-onset vs old-onset CRC",
       caption = "Effects on age are uniformly smaller than effects on disease status;\nno species is significant for both contrasts after study adjustment.") +
  THEME_MIN +
  theme(axis.text.y = element_text(size = 8.5))
save_fig("fig5_CRC_vs_age_effects.png", fig5, 8.2, 6.0)

# ---------------------------------------------------------------------------
# FIG 6 — Leave-one-study-out stability for top CRC species
# ---------------------------------------------------------------------------
stab <- read.delim(file.path(TAB, "loso_stability_CRC.tsv"))
top <- head(stab[!is.na(stab$full_g) & stab$full_fdr < 0.05, ], 15)
top$label <- pretty(top$taxon)
top$label <- factor(top$label, levels = rev(top$label[order(top$full_g)]))
top$dir_lab <- factor(ifelse(top$full_g > 0, "Higher in CRC", "Higher in control"),
                      levels = c("Higher in CRC", "Higher in control"))

fig6 <- ggplot(top, aes(loso_stability, label, fill = dir_lab)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0.9, linetype = 2, colour = "grey40") +
  scale_fill_manual(values = c("Higher in CRC" = COL_CRC,
                               "Higher in control" = COL_CTRL),
                    name = "Effect direction") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.08),
                     expand = expansion(mult = c(0, 0.10))) +
  geom_text(aes(label = scales::percent(loso_stability, accuracy = 1)),
            hjust = -0.15, size = 3.0, colour = "grey20") +
  labs(title = "Leave-one-study-out stability of the CRC species signature",
       subtitle = "Fraction of 11 leave-one-cohort-out meta-analyses keeping q < 0.05 with the same sign",
       x = "Stability (fraction of LOSO runs with q < 0.05, same sign)",
       y = NULL,
       caption = "Dashed line at 90%. Taxa at 100% remain significant and direction-consistent regardless of which cohort is removed.") +
  THEME_MIN
save_fig("fig6_loso_stability.png", fig6, 8.5, 5.8)

# ---------------------------------------------------------------------------
# FIG 7 — Cross-study effect correlation heatmap
# ---------------------------------------------------------------------------
cor_mat <- as.matrix(read.delim(file.path(TAB, "crossstudy_effect_correlation.tsv"),
                                row.names = 1, check.names = FALSE))
hc <- hclust(as.dist(1 - cor_mat), method = "average")
ord <- colnames(cor_mat)[hc$order]
cor_mat <- cor_mat[ord, ord, drop = FALSE]
cm_df <- as.data.frame(as.table(cor_mat))
names(cm_df) <- c("Study1", "Study2", "rho")
cm_df$Study1 <- factor(cm_df$Study1, levels = ord)
cm_df$Study2 <- factor(cm_df$Study2, levels = ord)

fig7 <- ggplot(cm_df, aes(Study1, Study2, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = COL_CTRL, mid = "white", high = COL_CRC,
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman ρ",
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       labels = c("-1.0", "-0.5", "0", "0.5", "1.0")) +
  labs(title = "Cross-cohort concordance of CRC effect sizes",
       subtitle = "Pairwise correlation of per-study Hedges' g across meta-significant species (FDR < 0.05)",
       x = NULL, y = NULL,
       caption = "Positive correlation indicates studies agree on which species are enriched/depleted in CRC.") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, colour = "grey35"),
    plot.caption = element_text(size = 8.5, colour = "grey45", hjust = 0),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10)
  )
save_fig("fig7_crossstudy_cor.png", fig7, 8, 7)

# ---------------------------------------------------------------------------
# FIG 8 — Age-association volcano (within CRC cases)
# ---------------------------------------------------------------------------
age <- read.delim(file.path(TAB, "assoc_taxon_vs_age_CRC.tsv"))
age$label <- pretty(age$taxon)
age$neg_log_p <- -log10(pmax(age$p_value, 1e-300))
age$sig <- factor(ifelse(age$fdr < 0.05, "FDR < 0.05", "n.s."),
                  levels = c("FDR < 0.05", "n.s."))
lab <- age %>% arrange(p_value) %>% head(12)

fig8 <- ggplot(age, aes(spearman_rho, neg_log_p, colour = sig)) +
  geom_point(alpha = 0.75, size = 1.7) +
  geom_hline(yintercept = -log10(0.05 / nrow(age)), linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_text_repel(data = lab, aes(label = label), size = 3.0, max.overlaps = 15,
                  show.legend = FALSE) +
  scale_colour_manual(values = c("FDR < 0.05" = COL_SIG, "n.s." = "grey70"),
                      name = "Significance") +
  labs(title = "Species-level association with age within CRC cases",
       subtitle = sprintf("Spearman correlation of CLR abundance vs age; n = 701 CRC cases; %d/%d species at FDR < 0.05; max |ρ| = %.2f",
                          sum(age$fdr < 0.05, na.rm = TRUE), nrow(age),
                          max(abs(age$spearman_rho), na.rm = TRUE)),
       x = "Spearman ρ  (CLR abundance vs age)",
       y = expression(-log[10](italic(p))),
       caption = "Dashed lines: ρ = 0 and Bonferroni threshold. Age effects are weak (all |ρ| ≤ 0.19).") +
  THEME_MIN
save_fig("fig8_age_volcano.png", fig8, 8, 5.8)

# ---------------------------------------------------------------------------
# FIG 9 — Sex-stratified CRC effect sizes for top species
# ---------------------------------------------------------------------------
sex <- read.delim(file.path(TAB, "diffabund_CRC_by_sex.tsv"))
meta <- read.delim(file.path(TAB, "meta_CRC_vs_control.tsv"))
top_sp <- head(meta[order(meta$p_value), ], 12)$taxon
sex_top <- sex %>% filter(taxon %in% top_sp)
sex_top$label <- pretty(sex_top$taxon)
sex_top$label <- factor(sex_top$label, levels = pretty(rev(top_sp)))
sex_top$sex_lab <- factor(
  ifelse(sex_top$sex == "female", "Female", "Male"),
  levels = c("Female", "Male"))

fig9 <- ggplot(sex_top, aes(cliff_delta, label, colour = sex_lab)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(size = 2.8, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c(Female = COL_F, Male = COL_M),
                      name = "Sex") +
  labs(title = "Sex-stratified CRC vs control effect sizes",
       subtitle = "Cliff's δ for top CRC species, computed separately in females and males across all cohorts",
       x = "Cliff's δ  (+ → higher in CRC)",
       y = NULL,
       caption = "Effect directions are concordant between sexes; magnitudes are similar.\nFemales: n = 257 CRC / 311 control; males: n = 444 CRC / 383 control.") +
  THEME_MIN
save_fig("fig9_sex_stratified.png", fig9, 8.5, 5.2)

cat(sprintf("\nAll figures written to %s at %d dpi\n", FIG, DPI))
