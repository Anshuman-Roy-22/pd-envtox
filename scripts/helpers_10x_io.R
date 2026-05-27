# Shared helpers for reading 10x-style gzipped feature/barcode/matrix triplets.

read_10x_features <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)

  x <- read.delim(con, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  genes <- if (ncol(x) >= 2) x[[2]] else x[[1]]
  toupper(trimws(as.character(genes)))
}

read_10x_barcodes <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)

  x <- read.delim(con, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  as.character(x[[1]])
}

read_10x_matrix <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)

  counts <- Matrix::readMM(con)
  methods::as(counts, "CsparseMatrix")
}

read_10x_counts <- function(feature_path, barcode_path, matrix_path, sample_name = NULL) {
  stopifnot(file.exists(feature_path), file.exists(barcode_path), file.exists(matrix_path))

  genes <- make.unique(read_10x_features(feature_path))
  barcodes <- read_10x_barcodes(barcode_path)
  counts <- read_10x_matrix(matrix_path)

  if (nrow(counts) == length(barcodes) && ncol(counts) == length(genes)) {
    counts <- t(counts)
  }

  stopifnot(nrow(counts) == length(genes), ncol(counts) == length(barcodes))

  rownames(counts) <- genes
  stopifnot(!any(duplicated(rownames(counts))))

  if (is.null(sample_name)) {
    colnames(counts) <- barcodes
  } else {
    colnames(counts) <- paste0(sample_name, "_", barcodes)
  }

  counts
}

make_10x_seurat_object <- function(
  feature_path,
  barcode_path,
  matrix_path,
  sample_name,
  condition,
  project = "GSE187012"
) {
  counts <- read_10x_counts(
    feature_path = feature_path,
    barcode_path = barcode_path,
    matrix_path = matrix_path,
    sample_name = sample_name
  )

  obj <- SeuratObject::CreateSeuratObject(counts = counts, project = project)
  obj$sample <- sample_name
  obj$condition <- condition
  obj
}
