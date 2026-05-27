# scripts/15_results_pack_original20.R
suppressPackageStartupMessages({
  library(data.table)
})

dir.create("results", showWarnings = FALSE)

rank_path <- "results/robust_sporadic_candidates_ranked.csv"
ct_path   <- "results/original20_celltype_specific_logFC.csv"
mptp10_path <- "results/gse17542_limma_mptp10d_vs_control.csv"

stopifnot(file.exists(rank_path), file.exists(ct_path), file.exists(mptp10_path))

rank <- fread(rank_path); rank[, gene := toupper(gene)]
ct <- fread(ct_path); ct[, gene := toupper(gene)]
mptp <- fread(mptp10_path); mptp[, gene := toupper(gene)]
mptp <- mptp[, .(gene, mptp10d_logFC = logFC, mptp10d_FDR = adj.P.Val, mptp10d_t = t)]

genes20 <- unique(ct$gene)

out <- rank[gene %in% genes20]
out <- merge(out, ct, by="gene", all.x=TRUE)
out <- merge(out, mptp, by="gene", all.x=TRUE)

# add a readable "call" column
out[, call := fifelse(as.logical(same_direction_10d) & as.logical(same_direction_2d), "Conserved + Replicated",
               fifelse(as.logical(same_direction_10d), "Conserved (10d)", "Not conserved"))]

setorder(out, -reactivity)
fwrite(out, "results/original20_results_pack.csv")
cat("Saved: results/original20_results_pack.csv\n")
