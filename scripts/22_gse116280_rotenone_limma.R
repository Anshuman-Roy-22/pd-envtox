#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))

root <- normalizePath(".")
out <- file.path(root, "results/v2/gse116280_rotenone")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

lines <- readLines(gzfile(file.path(root, "data_raw/gse116280/GSE116280_series_matrix.txt.gz")), warn = FALSE)
begin <- grep("^!series_matrix_table_begin", lines) + 1L
end <- grep("^!series_matrix_table_end", lines) - 1L
tab <- read.delim(text = paste(lines[begin:end], collapse = "\n"), check.names = FALSE,
                  quote = "\"", stringsAsFactors = FALSE)
rownames(tab) <- tab[[1]]
x <- as.matrix(tab[-1]); storage.mode(x) <- "numeric"

ann_lines <- readLines(file.path(root, "data_raw/gse116280/GPL17077_table.txt"), warn = FALSE)
header <- grep("^ID\\t", ann_lines)[1]
ann <- read.delim(text = paste(ann_lines[header:(grep("^!platform_table_end", ann_lines)[1]-1L)], collapse = "\n"),
                  check.names = FALSE, quote = "", fill = TRUE, stringsAsFactors = FALSE)
ann <- ann[, c("ID", "GENE_SYMBOL")]
names(ann) <- c("probe", "gene")
ann <- ann[ann$probe %in% rownames(x) & !is.na(ann$gene) & nzchar(ann$gene) &
           !grepl("///|//|;|,", ann$gene), ]
ann$gene <- toupper(trimws(ann$gene))
ann <- ann[grepl("^[A-Z0-9_.-]+$", ann$gene), ]

# Series matrix is already Entrez-collapsed; retain one row per unambiguous gene
# by greatest all-sample mean, with probe ID deterministic tie break.
ann$AveExpr <- rowMeans(x[ann$probe, , drop = FALSE])
ann <- ann[order(ann$gene, -ann$AveExpr, ann$probe), ]
ann <- ann[!duplicated(ann$gene), ]
x <- x[ann$probe, , drop = FALSE]

group <- factor(rep(c("D8_C", "D8_50_12", "D8_50_24", "D8_100_24",
                      "D15_C", "D15_50_12", "D15_50_24", "D15_100_24"), each = 3))
design <- model.matrix(~0 + group); colnames(design) <- levels(group)
contr <- makeContrasts(
  D15_50nM_24h = D15_50_24 - D15_C,
  D15_50nM_12h = D15_50_12 - D15_C,
  D15_100nM_24h = D15_100_24 - D15_C,
  D8_50nM_24h = D8_50_24 - D8_C,
  D8_50nM_12h = D8_50_12 - D8_C,
  D8_100nM_24h = D8_100_24 - D8_C,
  levels = design)
fit <- eBayes(contrasts.fit(lmFit(x, design), contr))
results <- do.call(rbind, lapply(colnames(contr), function(z) {
  tt <- topTable(fit, coef = z, number = Inf, sort.by = "none", adjust.method = "BH")
  tt$probe <- rownames(tt); tt <- merge(tt, ann[, c("probe", "gene")], by = "probe")
  tt$contrast <- z; tt
}))
write.table(results, file.path(out, "GSE116280_rotenone_limma_all_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ann, file.path(out, "GSE116280_probe_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
