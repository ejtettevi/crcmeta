#!/usr/bin/env Rscript
# Study-adjusted analysis: per-cohort effect sizes + random-effects meta-analysis.
#   Combines the 11 CRC/adenoma studies WITHOUT naively pooling samples, so that
#   cross-study batch effects are modelled as between-study heterogeneity rather
#   than contaminating the within-sample signal.
#
# For each comparison, per study & per species we compute Hedges' g (SMD) and its
# variance, then combine via DerSimonian-Laird random effects (tau^2 by DL, I^2).
# A Fisher's-combined-p and direction-consistency count are added as robustness.
#
# Comparisons:
#   (A) CRC vs Control
#   (B) Adenoma vs Control
#   (C) Young (<50) vs Old (>=50) CRC   [CRC cases only]

suppressPackageStartupMessages({ library(dplyr) })

ROOT <- Sys.getenv("CRC_ROOT", unset = getwd())
PROF <- file.path(ROOT, "data", "profiles")
TAB  <- file.path(ROOT, "results", "tables")
FIG  <- file.path(ROOT, "results", "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

species_key <- function(rn) sub("^.*\\|", "", rn)

# ---- load SE objects (species-level, per study) ----
rds <- list.files(PROF, pattern = "\\.relative_abundance\\.rds$", full.names = TRUE)
if (!length(rds)) stop("No .rds profiles. Run analysis/02_fetch_profiles.R first.")

mats <- list(); cds <- list()
for (f in rds) {
  se <- readRDS(f); st <- sub("\\..*$", "", basename(f))
  m  <- as.matrix(SummarizedExperiment::assay(se))
  agg <- rowsum(m, group = species_key(rownames(m)))
  cd  <- as.data.frame(SummarizedExperiment::colData(se)); cd$study_name <- st
  colnames(agg) <- rownames(cd)
  mats[[st]] <- agg; cds[[st]] <- cd
}

taxa <- Reduce(intersect, lapply(mats, rownames))
cat("Shared species-level taxa:", length(taxa), "across", length(mats), "studies\n")

# ---- CLR per study (compositional transform is study-internal) ----
clr_mat <- function(m) {
  x <- m; x[x <= 0] <- 1e-6
  sweep(log(x), 2, colMeans(log(x)), "-")
}
clr_list <- lapply(mats, function(m) clr_mat(m[taxa, , drop = FALSE]))

# ---- per-study Hedges' g for a two-group contrast ----
hedges_g <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  n1 <- length(a); n2 <- length(b)
  if (n1 < 2 || n2 < 2) return(c(g = NA, v = NA, n1 = n1, n2 = n2))
  s1 <- sd(a); s2 <- sd(b)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(c(g = NA, v = NA, n1 = n1, n2 = n2))
  d  <- (mean(a) - mean(b)) / sp
  J  <- 1 - 3 / (4 * (n1 + n2) - 9)
  g  <- J * d
  v  <- (n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2 - 3.94))
  c(g = g, v = v, n1 = n1, n2 = n2)
}

# ---- DerSimonian-Laird random-effects combination ----
dl_meta <- function(g, v) {
  ok <- is.finite(g) & is.finite(v) & v > 0
  g <- g[ok]; v <- v[ok]; k <- length(g)
  if (k == 0) return(list(k = 0, est = NA, se = NA, lo = NA, hi = NA, p = NA, tau2 = NA, I2 = NA))
  w  <- 1 / v
  fe <- sum(w * g) / sum(w)
  Q  <- sum(w * (g - fe)^2)
  C  <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - (k - 1)) / C)
  I2 <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
  ws <- 1 / (v + tau2)
  est <- sum(ws * g) / sum(ws); se <- sqrt(1 / sum(ws))
  z <- est / se; p <- 2 * pnorm(-abs(z))
  list(k = k, est = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se,
       p = p, tau2 = tau2, I2 = I2)
}

# ---- build group labels ----
std_cond <- function(cond) ifelse(grepl("CRC", cond, ignore.case = TRUE), "CRC",
                           ifelse(grepl("control", cond, ignore.case = TRUE), "Control",
                           ifelse(grepl("adenoma", cond, ignore.case = TRUE), "Adenoma", NA)))
cond_all <- lapply(cds, function(d) std_cond(d$study_condition))
age_all  <- lapply(cds, function(d) suppressWarnings(as.numeric(d$age)))

# ---- meta-analyse a contrast across studies ----
meta_contrast <- function(A, B, crc_only = FALSE, label, outfile) {
  G <- rep(NA_real_, length(taxa)); rowsums_p <- rep(1, length(taxa)); pos <- rep(0, length(taxa)); neg <- rep(0, length(taxa))
  per_study <- list(); k_study <- 0
  gmat <- matrix(NA_real_, length(taxa), length(mats), dimnames = list(taxa, names(mats)))
  vmat <- gmat
  for (st in names(mats)) {
    cond <- cond_all[[st]]; age <- age_all[[st]]
    if (crc_only) {
      keep <- which(!is.na(cond) & cond == "CRC" & !is.na(age))
      grp  <- ifelse(age[keep] < 50, "Young", "Old")
      a_idx <- keep[grp == A]; b_idx <- keep[grp == B]
    } else {
      a_idx <- which(cond == A); b_idx <- which(cond == B)
    }
    if (length(a_idx) < 2 || length(b_idx) < 2) next
    k_study <- k_study + 1
    X <- clr_list[[st]]
    gs <- vapply(taxa, function(t) {
      r <- hedges_g(X[t, a_idx], X[t, b_idx]); c(g = r[["g"]], v = r[["v"]])
    }, numeric(2))
    gmat[, st] <- gs["g", ]; vmat[, st] <- gs["v", ]
  }
  stopifnot(k_study > 0)
  res <- do.call(rbind, lapply(seq_along(taxa), function(i) {
    m <- dl_meta(gmat[i, ], vmat[i, ])
    g <- gmat[i, ]; v <- vmat[i, ]; ok <- is.finite(g)
    pv <- 2 * pnorm(-abs(g / sqrt(v)))[ok]
    fisher_p <- if (length(pv)) 1 - pchisq(-2 * sum(log(pv)), 2 * length(pv)) else NA
    data.frame(taxon = taxa[i], k_studies = m$k, pooled_g = m$est, se = m$se,
               ci_lb = m$lo, ci_ub = m$hi, z = m$est / m$se, p_value = m$p,
               fisher_p = fisher_p, tau2 = m$tau2, I2 = m$I2,
               n_pos = sum(g > 0, na.rm = TRUE), n_neg = sum(g < 0, na.rm = TRUE))
  }))
  res$fdr <- p.adjust(res$p_value, "BH"); res$comparison <- label
  res <- res[order(res$p_value), ]
  write.table(res, file.path(TAB, outfile), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%-22s studies=%d  q<0.05 = %d / %d  (median I2=%.0f%%)\n",
              label, k_study, sum(res$fdr < 0.05, na.rm = TRUE), nrow(res), median(res$I2, na.rm = TRUE)))
  list(res = res, gmat = gmat, vmat = vmat, studies = colnames(gmat)[colSums(!is.na(gmat)) > 0])
}

cat("\n=== Study-adjusted random-effects meta-analysis ===\n")
mC <- meta_contrast("CRC", "Control", FALSE, "CRC_vs_Control", "meta_CRC_vs_control.tsv")
mA <- meta_contrast("Adenoma", "Control", FALSE, "Adenoma_vs_Control", "meta_Adenoma_vs_control.tsv")
mY <- meta_contrast("Young", "Old", TRUE, "Young_vs_Old_CRC", "meta_Young_vs_Old_CRC.tsv")

# ---- forest plot for top taxa (CRC vs Control) ----
forest <- function(m, label, file, ntop = 12) {
  r <- head(m$res[!is.na(m$res$pooled_g), ], ntop)
  if (!nrow(r)) return(invisible())
  hpx <- max(500, 90 + 34 * nrow(r))
  png(file.path(FIG, file), width = 1500, height = hpx, units = "px", res = 150)
  op <- par(mar = c(4, 16, 3, 4)); on.exit({par(op); dev.off()}, add = TRUE)
  y <- rev(seq_len(nrow(r)))
  plot(r$pooled_g, y, xlim = range(c(r$ci_lb, r$ci_ub), na.rm = TRUE),
       pch = 19, col = ifelse(r$fdr < 0.05, "firebrick", "grey50"),
       yaxt = "n", ylab = "", xlab = "Hedges' g (random effects)",
       main = paste0(label, "  (k studies; red = FDR<0.05)"))
  segments(r$ci_lb, y, r$ci_ub, y, col = ifelse(r$fdr < 0.05, "firebrick", "grey50"))
  abline(v = 0, lty = 2, col = "grey40")
  axis(2, at = y, labels = sub("^s__", "", r$taxon), las = 1, cex.axis = 0.8)
  invisible()
}
forest(mC, "CRC vs Control", "fig_forest_CRC_vs_control.png")
forest(mY, "Young vs Old CRC", "fig_forest_Young_vs_Old.png")

cat("\nTables:", TAB, "\nFigures:", FIG, "\n")
