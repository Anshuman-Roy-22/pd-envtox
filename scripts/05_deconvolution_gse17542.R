# scripts/05_deconvolution_gse17542.R
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(nnls)
})

bulk_path <- "data_intermediate/gse17542_clean.rds"
sig_path  <- "data_intermediate/gse187012_signatures.rds"
out_path  <- "data_intermediate/gse17542_deconv.rds"
out_csv   <- "results/deconv_gse17542_proportions.csv"

stopifnot(file.exists(bulk_path), file.exists(sig_path))
dir.create("results", showWarnings = FALSE)
dir.create("data_intermediate", showWarnings = FALSE)

bulk <- readRDS(bulk_path)
sigobj <- readRDS(sig_path)

bulk_expr <- bulk$expr_gene   # genes x samples
meta <- bulk$meta

sig <- sigobj$signature_matrix  # genes x (DA, Microglia)

# ---- harmonize gene symbols ----
genes_common <- intersect(rownames(bulk_expr), rownames(sig))
if (length(genes_common) < 500) stop("Too few shared genes: ", length(genes_common))

Y <- bulk_expr[genes_common, , drop = FALSE]  # genes x samples
S <- sig[genes_common, , drop = FALSE]        # genes x celltypes

# ---- scale genes (z-score across genes) to reduce platform bias ----
zscale_rows <- function(M) {
  mu <- rowMeans(M)
  sd <- apply(M, 1, sd)
  sd[sd == 0] <- 1
  (M - mu) / sd
}

Yz <- zscale_rows(Y)
Sz <- zscale_rows(S)

# ---- NNLS per sample: Y ~= S * w, w>=0 ----
celltypes <- colnames(Sz)
samples <- colnames(Yz)

W <- matrix(NA_real_, nrow = length(samples), ncol = length(celltypes),
            dimnames = list(samples, celltypes))

for (j in seq_along(samples)) {
  y <- Yz[, j]
  fit <- nnls::nnls(as.matrix(Sz), as.numeric(y))
  w <- pmax(coef(fit), 0)
  if (sum(w) == 0) w <- rep(1, length(w))
  w <- w / sum(w)
  W[j, ] <- w
}

W_dt <- as.data.table(W, keep.rownames = "sample")
W_dt <- cbind(W_dt, as.data.table(meta[W_dt$sample, , drop = FALSE], keep.rownames = "sample_id"))

# ---- quick group summary (if group exists) ----
if ("group" %in% colnames(meta)) {
  grp <- meta[rownames(meta) %in% samples, "group"]
  W_dt[, group := grp[match(sample, rownames(meta))]]
}

fwrite(W_dt, out_csv)

saveRDS(list(
  proportions = W,
  meta = meta,
  genes_used = genes_common,
  signature_info = sigobj$chosen_clusters
), out_path)

cat("Saved:\n -", out_path, "\n -", out_csv, "\n")

if ("group" %in% names(W_dt)) {
  cat("\nMean proportions by group:\n")
  print(W_dt[, lapply(.SD, mean), by = group, .SDcols = celltypes])
}

