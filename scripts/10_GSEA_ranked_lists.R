# scripts/10_GSEA_ranked_lists.R
suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(msigdbr)
})

dir.create("results", showWarnings = FALSE)

# Load ranked list from MPTP10d limma table
in_path <- "results/gse17542_limma_mptp10d_vs_control.csv"
stopifnot(file.exists(in_path))

dt <- fread(in_path)

# Expect columns: gene, perturb (from 04 script). If missing, fall back.
if (!("gene" %in% names(dt))) stop("Input must have a 'gene' column.")
if (!("perturb" %in% names(dt))) {
  # fallback: use abs(logFC) * -log10(adj.P.Val)
  if (!all(c("logFC", "adj.P.Val") %in% names(dt))) {
    stop("Input must contain 'perturb' OR both 'logFC' and 'adj.P.Val'.")
  }
  dt[, perturb := abs(logFC) * (-log10(pmax(adj.P.Val, 1e-300)))]
}

# Clean gene IDs and keep plausible symbol-ish strings
dt[, gene := toupper(trimws(as.character(gene)))]
dt <- dt[gene != ""]
dt <- dt[grepl("^[A-Z0-9._-]+$", gene)]         # remove weird junk
dt <- dt[!grepl("^\\d+$", gene)]                # drop pure numeric
dt <- dt[!is.na(perturb) & is.finite(perturb)]

# If duplicates exist, keep the max perturbation per gene
dt <- dt[, .(perturb = max(perturb, na.rm = TRUE)), by = gene]

# Build geneList (named numeric vector)
geneList <- dt$perturb
names(geneList) <- dt$gene
geneList <- sort(geneList, decreasing = TRUE)

# Save the actual rank list used (for transparency)
rank_out <- data.table(gene = names(geneList), rank_score = as.numeric(geneList))
fwrite(rank_out, "results/gsea_ranklist_mptp10d_used.csv")

cat("Ranklist size (unique genes):", length(geneList), "\n")

# MSigDB gene sets for Mouse
# Reactome (C2:CP:REACTOME) + GO Biological Process (C5:BP)
ms_react <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
ms_bp    <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")

term2gene_react <- ms_react[, c("gs_name", "gene_symbol")]
term2gene_bp    <- ms_bp[, c("gs_name", "gene_symbol")]

# Ensure uppercase gene symbols
term2gene_react$gene_symbol <- toupper(term2gene_react$gene_symbol)
term2gene_bp$gene_symbol    <- toupper(term2gene_bp$gene_symbol)

# Run GSEA (symbol-based, no Entrez mapping)
# exploratory cutoffs; tighten later if needed
gsea_react <- tryCatch(
  GSEA(
    geneList = geneList,
    TERM2GENE = term2gene_react,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.25,
    verbose = FALSE
  ),
  error = function(e) NULL
)

gsea_bp <- tryCatch(
  GSEA(
    geneList = geneList,
    TERM2GENE = term2gene_bp,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.25,
    verbose = FALSE
  ),
  error = function(e) NULL
)

saveRDS(list(reactome = gsea_react, gobp = gsea_bp), "results/gsea_mptp10d_msigdb.rds")
cat("Saved: results/gsea_mptp10d_msigdb.rds\n")

# Also export result tables (if present)
export_res <- function(obj, out_csv) {
  if (is.null(obj)) {
    fwrite(data.table(), out_csv)
    return(invisible(FALSE))
  }
  res <- as.data.frame(obj@result)
  if (nrow(res) == 0) {
    fwrite(data.table(), out_csv)
    return(invisible(FALSE))
  }
  fwrite(as.data.table(res), out_csv)
  invisible(TRUE)
}

ok1 <- export_res(gsea_react, "results/gsea_mptp10d_reactome_table.csv")
ok2 <- export_res(gsea_bp,    "results/gsea_mptp10d_gobp_table.csv")

cat("Reactome enriched terms found? ", ok1, "\n", sep = "")
cat("GO:BP enriched terms found?     ", ok2, "\n", sep = "")
