#!/usr/bin/env Rscript

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required; restore renv.lock first")
}

root <- normalizePath(".")
inputs <- data.frame(
  name = c(
    "GSE196190 raw archive", "GSE229460 raw archive",
    "GSE287941 count matrix", "GSE4773 series matrix"
  ),
  relative_path = c(
    "data_raw/gse196190/GSE196190_RAW.tar",
    "data_raw/gse229460/GSE229460_RAW.tar",
    "data_raw/gse287941/GSE287941_RawCounts_AllSamples.txt.gz",
    "data_raw/gse4773/GSE4773_series_matrix.txt.gz"
  ),
  url = c(
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/",
      "GSE196190/suppl/GSE196190_RAW.tar"
    ),
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/",
      "GSE229460/suppl/GSE229460_RAW.tar"
    ),
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/",
      "GSE287941/suppl/GSE287941_RawCounts_AllSamples.txt.gz"
    ),
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE4nnn/",
      "GSE4773/matrix/GSE4773_series_matrix.txt.gz"
    )
  ),
  sha256 = c(
    "1a20222c1836188b8d0ee9617460372dcf538a5d0bd7007b306eeb7293184b1d",
    "d3c5e7a93c12c4aa0ca91389be94fd30f70ae4997e288ff0d55fdd5ab0988eb2",
    "54bc3c7e6a272e7f2c476c066ed1a686dfc95d8156fd7a7f6584af472b5c54b9",
    "5b944e95d019f4fad73dc8a6175050a17d64abc54ac6843c26989f2863b7a638"
  ),
  stringsAsFactors = FALSE
)

file_sha256 <- function(path) {
  digest::digest(path, algo = "sha256", serialize = FALSE, file = TRUE)
}

for (i in seq_len(nrow(inputs))) {
  target <- file.path(root, inputs$relative_path[[i]])
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(target)) {
    observed <- file_sha256(target)
    if (!identical(observed, inputs$sha256[[i]])) {
      stop(
        inputs$name[[i]], " exists but has the wrong SHA256. Move it aside and rerun: ",
        target
      )
    }
    message("Verified existing ", inputs$name[[i]])
    next
  }

  partial <- paste0(target, ".download")
  if (file.exists(partial)) unlink(partial)
  message("Downloading ", inputs$name[[i]], " ...")
  download.file(inputs$url[[i]], partial, mode = "wb", quiet = FALSE)
  observed <- file_sha256(partial)
  if (!identical(observed, inputs$sha256[[i]])) {
    unlink(partial)
    stop(inputs$name[[i]], " download failed SHA256 verification")
  }
  if (!file.rename(partial, target)) {
    unlink(partial)
    stop("Could not move verified download into place: ", target)
  }
}

extract_if_needed <- function(archive, destination, sentinel) {
  if (file.exists(file.path(destination, sentinel))) return(invisible(NULL))
  message("Extracting ", basename(archive), " ...")
  utils::untar(archive, exdir = destination)
  if (!file.exists(file.path(destination, sentinel))) {
    stop("Archive extraction did not create expected file: ", sentinel)
  }
}

extract_if_needed(
  file.path(root, inputs$relative_path[[1L]]),
  file.path(root, "data_raw/gse196190"),
  "GSM5862314_16C_MPP_0_t74_S6.genes.results.gz"
)
extract_if_needed(
  file.path(root, inputs$relative_path[[2L]]),
  file.path(root, "data_raw/gse229460"),
  "GSM7163946_01_Ctrl_1_S1.R1.readcounts.tsv.gz"
)

reactome <- file.path(root, "data_raw/pathways/ReactomePathways.gmt")
if (!file.exists(reactome)) stop("Tracked Reactome GMT is missing: ", reactome)
if (!identical(
  file_sha256(reactome),
  "89983d5c1f0af11c52edfeee7323eb425580ac6281d387a528562ab1787ce56b"
)) {
  stop("Tracked Reactome GMT failed SHA256 verification")
}

message("All held-out inputs downloaded, extracted, and checksum-verified.")
