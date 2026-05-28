# scripts/06a_extract_original20.R
suppressPackageStartupMessages({
  library(data.table)
})

in_rank <- "results/robust_sporadic_candidates_ranked.csv"
stopifnot(file.exists(in_rank))

dir.create("results", showWarnings = FALSE)

# Original 20-gene candidate panel
panel20 <- c(
  "SNX3","VPS26B","PSENEN","APH1A","SNX2","APH1B","EIF4E","EIF4B","EIF4A2","NCSTN",
  "PABPC1","DCTN2","EIF4A1","VPS26A","VPS29","APP","CLIP1","DYNC1H1","ACTR1B","DCTN3"
)

dt <- fread(in_rank)
dt[, gene := toupper(gene)]
panel20 <- toupper(panel20)

sub <- dt[gene %in% panel20]

# Add rank within full table
dt[, global_rank := .I]
sub <- merge(sub, dt[, .(gene, global_rank)], by="gene", all.x=TRUE)

# Ensure panel order is preserved
sub[, panel_order := match(gene, panel20)]
setorder(sub, panel_order)

# Identify missing panel genes
missing <- setdiff(panel20, sub$gene)

# Compact summary table (Table 1)
keep_cols <- c(
  "gene","global_rank","reactivity","shared_loose","inferred_celltype",
  "pq_log2FC","pq_norm",
  "mptp10d_logFC","mptp10d_fdr","mptp10d_norm","same_direction_10d",
  "mptp2d_logFC","mptp2d_fdr","mptp2d_norm","same_direction_2d"
)

avail_cols <- intersect(keep_cols, names(sub))
table1 <- sub[, ..avail_cols]

fwrite(table1, "results/original20_table1.csv")

cat("Saved: results/original20_table1.csv\n")
cat("Found", nrow(table1), "of 20 genes.\n")
if (length(missing) > 0) {
  cat("Missing genes (not present in merged table):\n")
  print(missing)
}
