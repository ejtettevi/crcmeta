suppressPackageStartupMessages({library(SummarizedExperiment)})
files <- list.files("data/profiles", pattern = "\\.relative_abundance\\.rds$", full.names = TRUE)
tot <- data.frame()
for (f in files) {
  se <- readRDS(f); st <- sub("\\..*$", "", basename(f))
  cd <- as.data.frame(colData(se)); cond <- as.character(cd$study_condition)
  c2 <- rep(NA_character_, length(cond))
  c2[grepl("CRC", cond, ignore.case = TRUE)] <- "CRC"
  c2[grepl("control", cond, ignore.case = TRUE) & is.na(c2)] <- "Control"
  c2[grepl("adenoma", cond, ignore.case = TRUE) & is.na(c2)] <- "Adenoma"
  age <- suppressWarnings(as.numeric(as.character(cd$age)))
  is_crc <- !is.na(c2) & c2 == "CRC"
  has_age <- !is.na(age)
  y <- sum(is_crc & has_age & age < 50)
  o <- sum(is_crc & has_age & age >= 50)
  ct <- if ("country" %in% names(cd)) paste(unique(na.omit(as.character(cd$country))), collapse = ", ") else ""
  tot <- rbind(tot, data.frame(
    study = st, n_samples = ncol(se),
    n_CRC = sum(c2 == "CRC", na.rm = TRUE),
    n_Control = sum(c2 == "Control", na.rm = TRUE),
    n_Adenoma = sum(c2 == "Adenoma", na.rm = TRUE),
    n_youngCRC = y, n_oldCRC = o,
    n_CRC_age_missing = sum(is_crc & !has_age),
    country = ct, stringsAsFactors = FALSE))
}
tot <- tot[order(tot$study), ]
write.table(tot, "results/tables/study_characteristics.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(tot)
cat("\nTOTALS n=", sum(tot$n_samples), " CRC=", sum(tot$n_CRC), " Ctrl=", sum(tot$n_Control),
    " Aden=", sum(tot$n_Adenoma), " yCRC=", sum(tot$n_youngCRC), " oCRC=", sum(tot$n_oldCRC),
    " ageMiss=", sum(tot$n_CRC_age_missing), "\n", sep = "")
