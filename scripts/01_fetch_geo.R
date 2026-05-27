# scripts/01_fetch_geo.R

suppressPackageStartupMessages(library(GEOquery))

dir.create("data_raw", showWarnings = FALSE)
dir.create("data_intermediate", showWarnings = FALSE)

### GSE17542 (Bulk Microarray, MPTP) ###
bulk_out <- "data_intermediate/gse17542_eset.rds"
bulk_compat_out <- "data_intermediate/gse17542_eset_unzipped.rds"

if (!file.exists(bulk_out) || !file.exists(bulk_compat_out)) {
  gse17542 <- getGEO("GSE17542", GSEMatrix = TRUE)
  gse17542 <- gse17542[[1]]

  saveRDS(gse17542, bulk_out)
  saveRDS(gse17542, bulk_compat_out)

  cat("GSE17542 downloaded and saved to:\n")
  cat(" - ", bulk_out, "\n", sep = "")
  cat(" - ", bulk_compat_out, "\n", sep = "")
} else {
  cat("GSE17542 cached files already exist; skipping download.\n")
}


### GSE187012 (snRNA-seq, Paraquat/Maneb) ###
# GEO stores this as a supplementary TAR file
supp_dir <- file.path("data_raw", "GSE187012")
supp_tar <- file.path(supp_dir, "GSE187012_RAW.tar")
extract_dir <- file.path(supp_dir, "extracted")
expected_files <- file.path(
  extract_dir,
  c(
    "GSM5667021_CTL_barcodes.tsv.gz",
    "GSM5667021_CTL_features.tsv.gz",
    "GSM5667021_CTL_matrix.mtx.gz",
    "GSM5667022_EXP_barcodes.tsv.gz",
    "GSM5667022_EXP_features.tsv.gz",
    "GSM5667022_EXP_matrix.mtx.gz"
  )
)

if (!file.exists(supp_tar)) {
  getGEOSuppFiles("GSE187012", baseDir = "data_raw")
  cat("GSE187012 supplementary TAR downloaded.\n")
} else {
  cat("GSE187012 supplementary TAR already exists; skipping download.\n")
}

dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

if (!all(file.exists(expected_files))) {
  utils::untar(supp_tar, exdir = extract_dir)
}

missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  stop(
    "Expected extracted GSE187012 files are still missing after untar: ",
    paste(missing_files, collapse = ", ")
  )
}

cat("GSE187012 extracted files are ready in ", extract_dir, "\n", sep = "")
