#!/usr/bin/env Rscript
# Export curatedMetagenomicData MetaPhlAn3 relative abundance + curated metadata
# for studies that contain colorectal cancer cases. Fast, low-disk, reproducible.

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  for (p in c("curatedMetagenomicData", "SummarizedExperiment", "dplyr")) {
    if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
  }
})

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(dplyr)
})

ROOT <- Sys.getenv("CRC_ROOT", unset = getwd())
if (!dir.exists(file.path(ROOT, "data"))) ROOT <- getwd()
OUT <- file.path(ROOT, "data", "profiles")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

md <- as.data.frame(sampleMetadata)
cat("curatedMetagenomicData samples:", nrow(md), " studies:", length(unique(md$study_name)), "\n")

is_crc_row <- grepl("CRC|adenoma|colorectal|cancer", md$study_condition, ignore.case = TRUE) |
              grepl("CRC|adenoma|colorectal|cancer", md$disease, ignore.case = TRUE)
crc_studies <- sort(unique(md$study_name[is_crc_row]))
cat("CRC-related studies found:", length(crc_studies), "\n")
print(crc_studies)

write.table(md, file.path(OUT, "sampleMetadata_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

for (st in crc_studies) {
  df <- tryCatch(
    curatedMetagenomicData(paste0(st, ".relative_abundance"), dryrun = FALSE),
    error = function(e) { message("skip ", st, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(df)) next
  se <- df[[1]]
  mat <- as.matrix(assay(se))
  saveRDS(se, file.path(OUT, paste0(st, ".relative_abundance.rds")))
  write.table(mat, file.path(OUT, paste0(st, ".relative_abundance.tsv")), sep = "\t", quote = FALSE)
  cat("saved", st, ":", nrow(mat), "taxa x", ncol(mat), "samples\n")
}

cat("Profiles written to:", OUT, "\n")
