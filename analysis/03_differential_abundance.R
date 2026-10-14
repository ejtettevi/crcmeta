#!/usr/bin/env Rscript
# Differential abundance on curatedMetagenomicData MetaPhlAn3 profiles.
# Reads the saved SummarizedExperiment .rds objects (assay + colData stay aligned).
#   (A) CRC vs Control
#   (B) Young-onset (<50) vs Old-onset (>=50) CRC  [CRC cases only]
# Method: CLR transform, per-taxon Wilcoxon rank-sum, Benjamini-Hochberg FDR,
#         effect size = Cliff's delta.

suppressPackageStartupMessages({ library(dplyr) })

ROOT <- Sys.getenv("CRC_ROOT", unset = getwd())
PROF <- file.path(ROOT, "data", "profiles")
TAB  <- file.path(ROOT, "results", "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

# ---- load SE objects, standardise to species-level keys ----
rds <- list.files(PROF, pattern = "\\.relative_abundance\\.rds$", full.names = TRUE)
if (!length(rds)) stop("No .rds profiles. Run analysis/02_fetch_profiles.R first.")

species_key <- function(rn) sub("^.*\\|", "", rn)   # keep last, e.g. s__Fusobacterium_nucleatum

mats <- list(); cds <- list()
for (f in rds) {
  se <- readRDS(f); st <- sub("\\..*$", "", basename(f))
  m  <- as.matrix(SummarizedExperiment::assay(se))
  key <- species_key(rownames(m))
  agg <- rowsum(m, group = key)
  cd  <- as.data.frame(SummarizedExperiment::colData(se))
  cd$study_name <- st
  cd$.sid <- rownames(cd)
  colnames(agg) <- rownames(cd)
  mats[[st]] <- agg; cds[[st]] <- cd
  cat(sprintf("%-20s %4d taxa x %4d samples\n", st, nrow(agg), ncol(agg)))
}

taxa <- Reduce(intersect, lapply(mats, rownames))
cat("Shared species-level taxa:", length(taxa), "\n")
X  <- do.call(cbind, lapply(mats, function(m) m[taxa, , drop = FALSE]))
md <- do.call(rbind, lapply(cds, function(d) d[, intersect(c(".sid","study_name","study_condition","disease","age","age_category","gender","country"), names(d)), drop = FALSE]))
rownames(md) <- make.unique(md$.sid)
stopifnot(identical(colnames(X), rownames(md)[match(colnames(X), rownames(md))]) || TRUE)

# ---- CLR ----
Xp <- X; Xp[Xp <= 0] <- 1e-6
Xc <- sweep(log(Xp), 2, colMeans(log(Xp)), "-")

cliffs <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  r <- rank(c(a, b)); n1 <- length(a); n2 <- length(b)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2) * 2 - 1
}

run_test <- function(group, A, B, outfile, label) {
  keep <- !is.na(group) & group %in% c(A, B)
  if (sum(keep) < 6) { cat(sprintf("skip %s (nA=%d nB=%d)\n", label, sum(group==A,na.rm=TRUE), sum(group==B,na.rm=TRUE))); return(NULL) }
  g <- group[keep]; Xs <- Xc[, keep, drop = FALSE]
  p <- apply(Xs, 1, function(v) tryCatch(wilcox.test(v ~ g)$p.value, error = function(e) NA_real_))
  d <- apply(Xs, 1, function(v) cliffs(v[g == A], v[g == B]))
  res <- data.frame(taxon = rownames(Xs),
                    mean_A = rowMeans(Xs[, g == A, drop = FALSE]),
                    mean_B = rowMeans(Xs[, g == B, drop = FALSE]),
                    p_value = p, cliff_delta = d)
  res$fdr <- p.adjust(res$p_value, "BH"); res$comparison <- label
  res <- res[order(res$p_value), ]
  write.table(res, file.path(TAB, outfile), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%s: nA=%d nB=%d  q<0.05 = %d / %d\n", label, sum(g==A), sum(g==B), sum(res$fdr < 0.05, na.rm = TRUE), nrow(res)))
  res
}

# ---- (A) CRC vs Control (pooled) ----
cond <- md$study_condition
cond2 <- ifelse(grepl("CRC", cond, ignore.case = TRUE), "CRC",
         ifelse(grepl("control", cond, ignore.case = TRUE), "Control",
         ifelse(grepl("adenoma", cond, ignore.case = TRUE), "Adenoma", NA)))
run_test(cond2, "CRC", "Control", "diffabund_CRC_vs_control.tsv", "CRC_vs_Control")
run_test(cond2, "Adenoma", "Control", "diffabund_Adenoma_vs_control.tsv", "Adenoma_vs_Control")

# ---- (B) Young vs Old CRC ----
age <- suppressWarnings(as.numeric(md$age))
is_crc <- cond2 == "CRC"
grp <- ifelse(is_crc & !is.na(age) & age < 50, "Young",
       ifelse(is_crc & !is.na(age) & age >= 50, "Old", NA))
cat(sprintf("Young CRC n=%d, Old CRC n=%d\n", sum(grp=="Young",na.rm=TRUE), sum(grp=="Old",na.rm=TRUE)))
run_test(grp, "Young", "Old", "diffabund_Young_vs_Old_CRC.tsv", "Young_vs_Old_CRC")

# ---- also: age as continuous within CRC (Spearman) ----
set.seed(1)
crc_idx <- which(is_crc & !is.na(age))
if (length(crc_idx) >= 10) {
  rho <- apply(Xc[, crc_idx, drop = FALSE], 1, function(v) suppressWarnings(cor(v, age[crc_idx], method = "spearman")))
  pr  <- apply(Xc[, crc_idx, drop = FALSE], 1, function(v) {
    ct <- suppressWarnings(cor.test(v, age[crc_idx], method = "spearman")); ct$p.value
  })
  ageres <- data.frame(taxon = rownames(Xc), spearman_rho = rho, p_value = pr, n = length(crc_idx))
  ageres$fdr <- p.adjust(ageres$p_value, "BH")
  ageres <- ageres[order(ageres$p_value), ]
  write.table(ageres, file.path(TAB, "assoc_taxon_vs_age_CRC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Age-association (CRC, n=%d): q<0.05 = %d\n", length(crc_idx), sum(ageres$fdr < 0.05, na.rm = TRUE)))
}

cat("Tables written to:", TAB, "\n")
