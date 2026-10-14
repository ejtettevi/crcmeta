#!/usr/bin/env Rscript
# Extended rigour analyses on the 11-cohort CRC shotgun resource:
#   1. Leave-one-study-out (LOSO) sensitivity for study-adjusted meta-analysis
#   2. Cross-cohort effect concordance (pairwise correlation of Hedges' g)
#   3. Overlap between CRC signature and age-associated taxa
#   4. Ordination (PCoA) of CLR profiles by condition, sex, and study
#   5. Sex-stratified CRC vs control (pooled CLR Wilcoxon + Cliff's delta)
#   6. Association of top CRC taxa with tumour location / stage where available
#   7. Prevalence of top taxa by group
# Outputs: results/tables/*.tsv (extended), results/figures/*.png

suppressPackageStartupMessages({
  library(dplyr); library(metafor); library(vegan)
  library(ggplot2); library(patchwork); library(ggrepel)
  library(RColorBrewer); library(tidyr)
})

ROOT <- Sys.getenv("CRC_ROOT", unset = getwd())
PROF <- file.path(ROOT, "data", "profiles")
TAB  <- file.path(ROOT, "results", "tables")
FIG  <- file.path(ROOT, "results", "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

set.seed(1)
species_key <- function(rn) sub("^.*\\|", "", rn)
pretty <- function(x) {
  x <- sub("^s__", "", x)
  x <- sub("^g__", "", x)
  gsub("_", " ", x)
}

# ---- load profiles ----
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
  cat(sprintf("loaded %-20s %4d x %4d\n", st, nrow(agg), ncol(agg)))
}
taxa <- Reduce(intersect, lapply(mats, rownames))
cat("Shared species:", length(taxa), " studies:", length(mats), "\n")

clr_mat <- function(m) {
  x <- m; x[x <= 0] <- 1e-6
  sweep(log(x), 2, colMeans(log(x)), "-")
}
clr_list <- lapply(mats, function(m) clr_mat(m[taxa, , drop = FALSE]))

std_cond <- function(cond) {
  ifelse(grepl("CRC", cond, ignore.case = TRUE), "CRC",
  ifelse(grepl("control", cond, ignore.case = TRUE), "Control",
  ifelse(grepl("adenoma", cond, ignore.case = TRUE), "Adenoma", NA)))
}
cond_all <- lapply(cds, function(d) std_cond(d$study_condition))
age_all  <- lapply(cds, function(d) suppressWarnings(as.numeric(as.character(d$age))))
sex_all  <- lapply(cds, function(d) {
  g <- as.character(d$gender)
  ifelse(g %in% c("female", "F", "Female"), "female",
  ifelse(g %in% c("male", "M", "Male"), "male", NA))
})

hedges_g <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  n1 <- length(a); n2 <- length(b)
  if (n1 < 2 || n2 < 2) return(c(g = NA, v = NA))
  s1 <- sd(a); s2 <- sd(b)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(c(g = NA, v = NA))
  d <- (mean(a) - mean(b)) / sp
  J <- 1 - 3 / (4 * (n1 + n2) - 9)
  g <- J * d
  v <- (n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2 - 3.94))
  c(g = g, v = v)
}

dl_meta <- function(g, v) {
  ok <- is.finite(g) & is.finite(v) & v > 0
  g <- g[ok]; v <- v[ok]; k <- length(g)
  if (k < 2) return(list(k = k, est = NA, se = NA, lo = NA, hi = NA, p = NA, tau2 = NA, I2 = NA))
  w <- 1 / v; fe <- sum(w * g) / sum(w)
  Q <- sum(w * (g - fe)^2); C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - (k - 1)) / max(C, 1e-12))
  I2 <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
  ws <- 1 / (v + tau2)
  est <- sum(ws * g) / sum(ws); se <- sqrt(1 / sum(ws))
  z <- est / se; p <- 2 * pnorm(-abs(z))
  list(k = k, est = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se,
       p = p, tau2 = tau2, I2 = I2)
}

# per-study Hedges' g matrices for CRC vs Control and Young vs Old
build_gmat <- function(A, B, crc_only = FALSE) {
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
    X <- clr_list[[st]]
    gs <- vapply(taxa, function(t) hedges_g(X[t, a_idx], X[t, b_idx]), numeric(2))
    gmat[, st] <- gs["g", ]; vmat[, st] <- gs["v", ]
  }
  list(g = gmat, v = vmat)
}

cat("\n=== Building per-study effect matrices ===\n")
mC <- build_gmat("CRC", "Control")
mY <- build_gmat("Young", "Old", crc_only = TRUE)

# ============================================================
# 1. Leave-one-study-out (LOSO) for CRC vs Control
# ============================================================
cat("\n=== 1. Leave-one-study-out sensitivity ===\n")
meta_full <- read.delim(file.path(TAB, "meta_CRC_vs_control.tsv"))
top_sp <- head(meta_full[order(meta_full$p_value), ], 20)$taxon

studies <- colnames(mC$g)
losos <- list()
for (omit in studies) {
  keep_st <- setdiff(studies, omit)
  g <- mC$g[, keep_st, drop = FALSE]; v <- mC$v[, keep_st, drop = FALSE]
  est <- se <- p <- rep(NA_real_, length(taxa))
  for (i in seq_along(taxa)) {
    m <- dl_meta(g[i, ], v[i, ])
    est[i] <- m$est; se[i] <- m$se; p[i] <- m$p
  }
  fdr <- p.adjust(p, "BH")
  losos[[omit]] <- data.frame(taxon = taxa, omitted = omit,
                              pooled_g = est, p_value = p, fdr = fdr)
}
loso_df <- bind_rows(losos)
write.table(loso_df, file.path(TAB, "loso_CRC_vs_control.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# stability: fraction of LOSO runs in which taxon stays FDR<0.05 with same sign
full_g <- setNames(meta_full$pooled_g, meta_full$taxon)
full_fdr <- setNames(meta_full$fdr, meta_full$taxon)
stability <- sapply(taxa, function(t) {
  sub <- loso_df[loso_df$taxon == t, ]
  if (!nrow(sub) || is.na(full_g[t])) return(NA_real_)
  sig <- sub$fdr < 0.05 & sign(sub$pooled_g) == sign(full_g[t])
  mean(sig, na.rm = TRUE)
})
stab_df <- data.frame(taxon = taxa,
                      full_g = as.numeric(full_g[taxa]),
                      full_fdr = as.numeric(full_fdr[taxa]),
                      loso_stability = stability)
stab_df <- stab_df[order(-stab_df$loso_stability, stab_df$full_fdr), ]
write.table(stab_df, file.path(TAB, "loso_stability_CRC.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Top-20 mean LOSO stability: %.0f%%\n",
            100 * mean(stab_df$loso_stability[match(top_sp, stab_df$taxon)], na.rm = TRUE)))

# ============================================================
# 2. Cross-cohort effect concordance
# ============================================================
cat("\n=== 2. Cross-cohort effect concordance ===\n")
# correlation of per-study Hedges' g among studies, for taxa significant in meta
sig_sp <- meta_full$taxon[meta_full$fdr < 0.05]
g_sig <- mC$g[sig_sp, , drop = FALSE]
# only studies with any data
st_use <- studies[colSums(!is.na(g_sig)) >= 5]
cor_mat <- cor(g_sig[, st_use, drop = FALSE], use = "pairwise.complete.obs",
               method = "spearman")
write.table(round(cor_mat, 3), file.path(TAB, "crossstudy_effect_correlation.tsv"),
            sep = "\t", quote = FALSE)
cat(sprintf("Median pairwise study correlation (sig taxa): %.2f\n",
            median(cor_mat[lower.tri(cor_mat)], na.rm = TRUE)))

# ============================================================
# 3. Overlap: CRC signature vs age-associated taxa
# ============================================================
cat("\n=== 3. CRC signature vs age association overlap ===\n")
age_res <- read.delim(file.path(TAB, "assoc_taxon_vs_age_CRC.tsv"))
crc_sig <- meta_full$taxon[meta_full$fdr < 0.05]
age_sig <- age_res$taxon[age_res$fdr < 0.05]
both <- intersect(crc_sig, age_sig)
only_crc <- setdiff(crc_sig, age_sig)
only_age <- setdiff(age_sig, crc_sig)
overlap_tab <- data.frame(
  category = c("CRC-significant", "Age-significant", "Both", "CRC only", "Age only"),
  n = c(length(crc_sig), length(age_sig), length(both), length(only_crc), length(only_age))
)
write.table(overlap_tab, file.path(TAB, "overlap_CRC_age.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(taxon = both),
            file.path(TAB, "overlap_taxa_CRC_age.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("CRC sig=%d, Age sig=%d, overlap=%d\n",
            length(crc_sig), length(age_sig), length(both)))

# ============================================================
# 4. Ordination (PCoA) of CLR profiles
# ============================================================
cat("\n=== 4. PCoA ordination ===\n")
# pooled CLR for samples with CRC/Control only
X_pool <- do.call(cbind, lapply(mats, function(m) m[taxa, , drop = FALSE]))
md <- do.call(rbind, lapply(cds, function(d)
  d[, intersect(c("study_name", "study_condition", "age", "gender", "country"),
                names(d)), drop = FALSE]))
rownames(md) <- colnames(X_pool)
md$condition <- std_cond(md$study_condition)
md$age_num <- suppressWarnings(as.numeric(as.character(md$age)))
md$sex <- ifelse(md$gender %in% c("female", "F", "Female"), "female",
          ifelse(md$gender %in% c("male", "M", "Male"), "male", NA))
md$age_group <- ifelse(!is.na(md$age_num) & md$condition == "CRC",
                ifelse(md$age_num < 50, "Young-onset CRC",
                ifelse(md$age_num >= 50, "Old-onset CRC", NA)), NA)

keep <- md$condition %in% c("CRC", "Control")
Xk <- X_pool[, keep]; mdk <- md[keep, ]
Xc <- clr_mat(Xk)
# Bray-Curtis on proportions for PCoA (more standard for microbiome) + CLR Euclidean
prop <- sweep(Xk, 2, colSums(Xk), "/")
bc <- vegdist(t(prop), method = "bray")
pcoa_bc <- cmdscale(bc, k = 2, eig = TRUE)
pcoa_df <- data.frame(
  PCo1 = pcoa_bc$points[, 1], PCo2 = pcoa_bc$points[, 2],
  condition = mdk$condition, study = mdk$study_name,
  sex = mdk$sex, age_group = mdk$age_group, country = mdk$country
)
var_exp <- 100 * pcoa_bc$eig[1:2] / sum(pcoa_bc$eig[pcoa_bc$eig > 0])
write.csv(pcoa_df, file.path(TAB, "pcoa_bray_curtis.csv"), row.names = FALSE)
writeLines(as.character(round(var_exp, 2)), file.path(TAB, "pcoa_variance.txt"))
cat(sprintf("PCoA variance: PCo1=%.1f%% PCo2=%.1f%%\n", var_exp[1], var_exp[2]))

# PERMANOVA: condition (controlling for study) and study (controlling for condition)
set.seed(1)
adonis_cond <- adonis2(bc ~ study_name + condition, data = mdk, permutations = 999, by = "margin")
write.csv(as.data.frame(adonis_cond), file.path(TAB, "permanova_condition_study.csv"))
print(adonis_cond)

# ============================================================
# 5. Sex-stratified CRC vs Control
# ============================================================
cat("\n=== 5. Sex-stratified analysis ===\n")
cliffs <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  r <- rank(c(a, b)); n1 <- length(a); n2 <- length(b)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2) * 2 - 1
}
sex_res <- list()
Xfull <- clr_mat(X_pool)
for (sx in c("female", "male")) {
  sel <- which(md$condition %in% c("CRC", "Control") & !is.na(md$sex) & md$sex == sx)
  if (length(sel) < 20) next
  g <- md$condition[sel]
  Xs <- Xfull[, sel, drop = FALSE]
  p <- apply(Xs, 1, function(v) tryCatch(wilcox.test(v ~ g)$p.value, error = function(e) NA_real_))
  d <- apply(Xs, 1, function(v) cliffs(v[g == "CRC"], v[g == "Control"]))
  sex_res[[sx]] <- data.frame(taxon = rownames(Xs), sex = sx,
                              n_CRC = sum(g == "CRC"), n_Control = sum(g == "Control"),
                              cliff_delta = d, p_value = p,
                              fdr = p.adjust(p, "BH"))
  cat(sprintf("%s: nCRC=%d nCtrl=%d  q<0.05=%d/%d\n", sx, sum(g == "CRC"),
              sum(g == "Control"), sum(sex_res[[sx]]$fdr < 0.05, na.rm = TRUE), nrow(Xs)))
}
sex_df <- bind_rows(sex_res)
sex_df <- sex_df[order(sex_df$sex, sex_df$p_value), ]
write.table(sex_df, file.path(TAB, "diffabund_CRC_by_sex.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================
# 6. Clinical correlates of top CRC taxa (Yachida only: location, stage)
# ============================================================
cat("\n=== 6. Clinical correlates (YachidaS_2019) ===\n")
y_st <- "YachidaS_2019"
cdY <- cds[[y_st]]; Xy <- clr_list[[y_st]]
condY <- std_cond(cdY$study_condition)
crcY <- which(condY == "CRC")
clin_out <- list()
# disease_location
if ("disease_location" %in% names(cdY)) {
  loc <- as.character(cdY$disease_location)
  tab <- table(loc[crcY], useNA = "ifany")
  cat("disease_location in CRC:\n"); print(tab)
  # proximal vs distal if codes allow
  loc_crc <- loc[crcY]
  if (any(grepl("proximal|right|colon", loc_crc, ignore.case = TRUE), na.rm = TRUE) ||
      any(grepl("distal|left|rectum", loc_crc, ignore.case = TRUE), na.rm = TRUE)) {
    grp <- ifelse(grepl("proximal|right", loc_crc, ignore.case = TRUE), "proximal",
           ifelse(grepl("distal|left|rectum", loc_crc, ignore.case = TRUE), "distal", NA))
    for (t in head(sig_sp, 15)) {
      a <- Xy[t, crcY][grp == "proximal"]; b <- Xy[t, crcY][grp == "distal"]
      if (length(a) >= 5 && length(b) >= 5) {
        pv <- tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA)
        clin_out[[length(clin_out) + 1]] <- data.frame(
          taxon = t, contrast = "proximal_vs_distal",
          n1 = length(a), n2 = length(b), cliff_delta = cliffs(a, b), p_value = pv)
      }
    }
  }
}
# AJCC stage
if ("ajcc" %in% names(cdY)) {
  stg <- as.character(cdY$ajcc)
  cat("AJCC in CRC:\n"); print(table(stg[crcY], useNA = "ifany"))
  # early (I/II) vs late (III/IV)
  early <- grepl("^I$|^II$|stage I|stage II", stg[crcY], ignore.case = TRUE)
  late  <- grepl("^III$|^IV$|stage III|stage IV", stg[crcY], ignore.case = TRUE)
  cat(sprintf("early(I/II)=%d late(III/IV)=%d\n", sum(early, na.rm=TRUE), sum(late, na.rm=TRUE)))
  if (sum(early, na.rm = TRUE) >= 5 && sum(late, na.rm = TRUE) >= 5) {
    for (t in head(sig_sp, 15)) {
      a <- Xy[t, crcY][early]; b <- Xy[t, crcY][late]
      pv <- tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA)
      clin_out[[length(clin_out) + 1]] <- data.frame(
        taxon = t, contrast = "early_vs_late_stage",
        n1 = sum(early, na.rm=TRUE), n2 = sum(late, na.rm=TRUE),
        cliff_delta = cliffs(a, b), p_value = pv)
    }
  }
}
if (length(clin_out)) {
  clin_df <- bind_rows(clin_out)
  clin_df$fdr <- p.adjust(clin_df$p_value, "BH")
  write.table(clin_df, file.path(TAB, "clinical_correlates_Yachida.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Clinical tests: %d, FDR<0.05: %d\n", nrow(clin_df), sum(clin_df$fdr < 0.05, na.rm = TRUE)))
}

# ============================================================
# 7. Prevalence of top taxa by group
# ============================================================
cat("\n=== 7. Prevalence of top CRC taxa ===\n")
prev_out <- list()
for (t in head(sig_sp, 20)) {
  for (cond_lab in c("CRC", "Control")) {
    vals <- c()
    for (st in names(mats)) {
      idx <- which(cond_all[[st]] == cond_lab)
      if (!length(idx)) next
      vals <- c(vals, mats[[st]][t, idx] > 0)
    }
    prev_out[[length(prev_out) + 1]] <- data.frame(
      taxon = t, condition = cond_lab,
      prevalence = mean(vals), n = length(vals))
  }
}
prev_df <- bind_rows(prev_out)
write.table(prev_df, file.path(TAB, "prevalence_top_taxa.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================
# 8. Effect-size comparison: CRC vs Control vs Young vs Old
# ============================================================
cat("\n=== 8. CRC effect vs age effect ===\n")
meta_young <- read.delim(file.path(TAB, "meta_Young_vs_Old_CRC.tsv"))
cmp <- merge(
  meta_full[, c("taxon", "pooled_g", "fdr")],
  meta_young[, c("taxon", "pooled_g", "fdr")],
  by = "taxon", suffixes = c("_CRC", "_Age")
)
names(cmp) <- c("taxon", "g_CRC", "fdr_CRC", "g_Age", "fdr_Age")
cmp <- cmp[order(cmp$fdr_CRC), ]
write.table(cmp, file.path(TAB, "effect_CRC_vs_age.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cor_eff <- cor(cmp$g_CRC, cmp$g_Age, method = "spearman", use = "complete.obs")
cat(sprintf("Spearman cor(CRC g, age g) across %d taxa: %.3f\n", nrow(cmp), cor_eff))

# ============================================================
# Summary counts for manuscript
# ============================================================
cat("\n=== SUMMARY COUNTS ===\n")
cat(sprintf("Shared species: %d\n", length(taxa)))
cat(sprintf("CRC meta FDR<0.05: %d\n", sum(meta_full$fdr < 0.05, na.rm = TRUE)))
cat(sprintf("Young-Old meta FDR<0.05: %d\n", sum(meta_young$fdr < 0.05, na.rm = TRUE)))
cat(sprintf("LOSO stability (top20 mean): %.0f%%\n",
            100 * mean(stab_df$loso_stability[match(top_sp, stab_df$taxon)], na.rm = TRUE)))
cat(sprintf("Cross-study median rho: %.2f\n", median(cor_mat[lower.tri(cor_mat)], na.rm = TRUE)))
cat(sprintf("Overlap CRC∩Age: %d\n", length(both)))
cat(sprintf("Sex-stratified FDR<0.05: F=%d M=%d\n",
            sum(sex_res$female$fdr < 0.05, na.rm = TRUE),
            sum(sex_res$male$fdr < 0.05, na.rm = TRUE)))

cat("\nDone. Extended tables in", TAB, "\n")
