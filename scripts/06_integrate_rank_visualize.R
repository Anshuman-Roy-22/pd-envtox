# scripts/06_integrate_rank_visualize.R
suppressPackageStartupMessages({
  library(data.table)
})

# Inputs
deg_path <- "data_intermediate/gse17542_deg_tables.rds"
pq_path  <- "data_intermediate/gse187012_pseudobulk.rds"
mark_path <- "data_intermediate/gse187012_marker_lists.rds"  # optional but recommended

stopifnot(file.exists(deg_path), file.exists(pq_path))

deg <- readRDS(deg_path)
pq  <- readRDS(pq_path)

tab2d  <- as.data.table(deg$mptp2d_vs_ctl)
tab10d <- as.data.table(deg$mptp10d_vs_ctl)
pq_dt  <- as.data.table(pq)

# Standardize column names
setnames(tab2d,  old = c("logFC","adj.P.Val","perturb","perturb_norm","gene"),
               new = c("mptp2d_logFC","mptp2d_fdr","mptp2d_perturb","mptp2d_norm","gene"))
setnames(tab10d, old = c("logFC","adj.P.Val","perturb","perturb_norm","gene"),
               new = c("mptp10d_logFC","mptp10d_fdr","mptp10d_perturb","mptp10d_norm","gene"))

# PQ table already has: gene, log2FC, perturb, perturb_norm
setnames(pq_dt, old = c("log2FC","perturb","perturb_norm"),
               new = c("pq_log2FC","pq_perturb","pq_norm"))

# Keep essentials only
tab2d  <- tab2d[,  .(gene, mptp2d_logFC,  mptp2d_fdr,  mptp2d_perturb,  mptp2d_norm)]
tab10d <- tab10d[, .(gene, mptp10d_logFC, mptp10d_fdr, mptp10d_perturb, mptp10d_norm)]
pq_dt  <- pq_dt[,  .(gene, pq_log2FC, pq_perturb, pq_norm)]

# Merge across datasets
dt <- merge(tab10d, pq_dt, by = "gene", all = FALSE)
dt <- merge(dt, tab2d, by = "gene", all.x = TRUE)

# Reactivity score (primary = MPTP 10d)
dt[, reactivity := 0.66 * mptp10d_norm + 0.34 * pq_norm]

# Shared / replication logic
# "Shared perturbed" means above a minimal norm threshold in BOTH
# (thresholds can be adjusted later if needed)
dt[, shared_loose := (mptp10d_norm >= 0.20 & pq_norm >= 0.20)]

# Direction consistency (10d vs PQ)
dt[, same_direction_10d := sign(mptp10d_logFC) == sign(pq_log2FC)]

# 2d replication: does it move in same direction too?
dt[, same_direction_2d := !is.na(mptp2d_logFC) & (sign(mptp2d_logFC) == sign(pq_log2FC))]

# "High-confidence MPTP" if FDR passes a standard cutoff
dt[, mptp10d_sig_fdr05 := (mptp10d_fdr <= 0.05)]
dt[, mptp2d_sig_fdr05  := (!is.na(mptp2d_fdr) & mptp2d_fdr <= 0.05)]

# Optional: cell-type attribution via marker lists
dt[, inferred_celltype := "Unknown"]

if (file.exists(mark_path)) {
  markers <- readRDS(mark_path)
  da_set <- unique(toupper(markers$DA_top50))
  mg_set <- unique(toupper(markers$Microglia_top50))

  dt[gene %in% da_set, inferred_celltype := "DA"]
  dt[gene %in% mg_set, inferred_celltype := "Microglia"]
  dt[gene %in% intersect(da_set, mg_set), inferred_celltype := "Both"]
}

# Final ranking + export
setorder(dt, -reactivity)

out_main <- "results/robust_sporadic_candidates_ranked.csv"
fwrite(dt, out_main)

# Also export a compact top table
top <- dt[, .(
  gene,
  reactivity,
  shared_loose,
  inferred_celltype,
  pq_log2FC,
  mptp10d_logFC,
  mptp10d_fdr,
  same_direction_10d,
  mptp2d_logFC,
  mptp2d_fdr,
  same_direction_2d
)]
fwrite(top[1:200], "results/top200_candidates_for_writeup.csv")

cat("Saved:\n - ", out_main, "\n - results/top200_candidates_for_writeup.csv\n", sep="")

# Quick summary counts
cat("\nSummary:\n")
cat("Genes merged:", nrow(dt), "\n")
cat("Shared_loose:", sum(dt$shared_loose, na.rm=TRUE), "\n")
cat("MPTP10d FDR<=0.05:", sum(dt$mptp10d_sig_fdr05, na.rm=TRUE), "\n")
cat("Same direction (10d vs PQ):", sum(dt$same_direction_10d, na.rm=TRUE), "\n")
cat("Replicated direction in 2d:", sum(dt$same_direction_2d, na.rm=TRUE), "\n")

cat("\nDONE: 06_integrate_rank_visualize\n")
