#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(limma))

root <- normalizePath(".")
out <- file.path(root, "results", "v2", "limma_confirmation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines) + 1L
  end <- grep("^!series_matrix_table_end", lines) - 1L
  tab <- read.delim(text = paste(lines[begin:end], collapse = "\n"), check.names = FALSE,
                    quote = "\"", stringsAsFactors = FALSE)
  rownames(tab) <- tab[[1]]
  values <- as.matrix(tab[-1])
  suppressWarnings(storage.mode(values) <- "numeric")
  values[complete.cases(values), , drop = FALSE]
}

read_annotation <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  header <- grep("^ID\\t", lines)[1]
  ann <- read.delim(text = paste(lines[header:length(lines)], collapse = "\n"),
                    check.names = FALSE, quote = "", fill = TRUE, comment.char = "",
                    stringsAsFactors = FALSE)
  ann <- ann[, c("ID", "Gene.symbol")]
  names(ann) <- c("probe", "gene")
  ann <- ann[!is.na(ann$gene) & nzchar(ann$gene) &
             !grepl("///|//|;", ann$gene), ]
  ann$gene <- toupper(trimws(ann$gene))
  ann[grepl("^[A-Z0-9_.-]+$", ann$gene), ]
}

bh_table <- function(fit, coef_name, probe_map, collapse = "highest_mean") {
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none", adjust.method = "BH")
  tt$probe <- rownames(tt)
  tt <- merge(tt, probe_map, by = "probe")
  if (collapse == "highest_mean") {
    tt <- tt[order(tt$gene, -tt$AveExpr, tt$probe), ]
    tt <- tt[!duplicated(tt$gene), ]
  }
  tt$contrast <- coef_name
  tt
}

# GSE46798: same one-probe-per-gene, all-sample-mean mapping frozen in Python.
g46798_path <- file.path(out, "GSE46798_limma_all_genes.tsv")
if (!file.exists(g46798_path)) {
x <- read_matrix(file.path(root, "data_raw/gse46798/GSE46798_series_matrix.txt.gz"))
if (quantile(x, .99, na.rm = TRUE) > 100) x <- log2(pmax(x, 0) + 1)
map46798 <- read.delim(file.path(root, "results/v2/gse46798/GSE46798_frozen_probe_map.tsv"),
                       stringsAsFactors = FALSE)
map46798 <- map46798[map46798$probe %in% rownames(x), ]
x <- x[map46798$probe, , drop = FALSE]
genotype <- factor(rep(c("Corrected", "A53T", "Corrected", "A53T"), each = 3),
                   levels = c("Corrected", "A53T"))
exposure <- factor(rep(c("Vehicle", "Vehicle", "PQMB", "PQMB"), each = 3),
                   levels = c("Vehicle", "PQMB"))
design <- model.matrix(~ genotype * exposure)
colnames(design) <- make.names(colnames(design))
fit <- eBayes(lmFit(x, design))
contr <- makeContrasts(
  exposure_corrected = exposurePQMB,
  exposure_A53T = exposurePQMB + genotypeA53T.exposurePQMB,
  genotype_by_exposure = genotypeA53T.exposurePQMB,
  A53T_vehicle = genotypeA53T,
  levels = design
)
fitc <- eBayes(contrasts.fit(fit, contr))
res46798 <- do.call(rbind, lapply(colnames(contr), function(z) bh_table(fitc, z, map46798)))
write.table(res46798, g46798_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# GSE17542: use exact sample-wise median gene collapse frozen in Python.
g17542_path <- file.path(out, "GSE17542_limma_all_genes.tsv")
if (!file.exists(g17542_path)) {
x <- read_matrix(file.path(root, "data_raw/gse17542/GSE17542_series_matrix.txt.gz"))
manifest <- read.delim(file.path(root, "results/v2/gse17542_sn_vta/GSE17542_gene_probe_manifest.tsv"),
                       stringsAsFactors = FALSE)
gene_expr <- t(vapply(seq_len(nrow(manifest)), function(i) {
  probes <- strsplit(manifest$eligible_probes[i], ";", fixed = TRUE)[[1]]
  apply(x[probes, , drop = FALSE], 2, median)
}, numeric(ncol(x))))
rownames(gene_expr) <- manifest$gene
group <- factor(rep(c("SN_C", "SN_2", "SN_10", "VTA_C", "VTA_2", "VTA_10"), each = 3),
                levels = c("SN_C", "SN_2", "SN_10", "VTA_C", "VTA_2", "VTA_10"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contr <- makeContrasts(
  SN_MPTP10_vs_control = SN_10 - SN_C,
  interaction10_SNminusVTA = (SN_10 - SN_C) - (VTA_10 - VTA_C),
  VTA_MPTP10_vs_control = VTA_10 - VTA_C,
  SN_MPTP2_vs_control = SN_2 - SN_C,
  interaction2_SNminusVTA = (SN_2 - SN_C) - (VTA_2 - VTA_C),
  VTA_MPTP2_vs_control = VTA_2 - VTA_C,
  baseline_SN_minus_VTA = SN_C - VTA_C,
  levels = design
)
fitc <- eBayes(contrasts.fit(lmFit(gene_expr, design), contr))
res17542 <- do.call(rbind, lapply(colnames(contr), function(z) {
  tt <- topTable(fitc, coef = z, number = Inf, sort.by = "none", adjust.method = "BH")
  tt$gene <- rownames(tt); tt$contrast <- z; tt
}))
write.table(res17542, g17542_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Human postmortem cohorts: reuse outcome-blind frozen probe choices.
human_results <- list()
for (accession in c("GSE20141", "GSE7621")) {
  cohort_path <- file.path(out, paste0(accession, "_limma_all_genes.tsv"))
  if (file.exists(cohort_path)) {
    human_results[[accession]] <- read.delim(cohort_path, stringsAsFactors = FALSE)
    next
  }
  x <- read_matrix(file.path(root, "data_raw/human_validation", paste0(accession, "_series_matrix.txt.gz")))
  if (quantile(x, .99, na.rm = TRUE) > 100) x <- log2(pmax(x, 0) + 1)
  frozen <- read.delim(file.path(root, "results/v2/human_validation", paste0(accession, "_gene_effects.tsv")),
                       stringsAsFactors = FALSE)
  map <- unique(frozen[, c("probe", "gene")])
  map <- map[map$probe %in% rownames(x), ]
  x <- x[map$probe, , drop = FALSE]
  if (accession == "GSE20141") {
    condition <- factor(c(rep("Control", 8), rep("PD", 10)), levels = c("Control", "PD"))
  } else {
    condition <- factor(c(rep("Control", 9), rep("PD", 16)), levels = c("Control", "PD"))
  }
  design <- model.matrix(~ condition)
  fit <- eBayes(lmFit(x, design))
  tt <- topTable(fit, coef = "conditionPD", number = Inf, sort.by = "none", adjust.method = "BH")
  tt$probe <- rownames(tt)
  tt <- merge(tt, map, by = "probe")
  tt$accession <- accession
  write.table(tt, cohort_path,
              sep = "\t", quote = FALSE, row.names = FALSE)
  human_results[[accession]] <- tt
}

a <- human_results[["GSE20141"]][, c("gene", "logFC", "P.Value", "adj.P.Val")]
b <- human_results[["GSE7621"]][, c("gene", "logFC", "P.Value", "adj.P.Val")]
names(a)[-1] <- paste0(names(a)[-1], "_20141")
names(b)[-1] <- paste0(names(b)[-1], "_7621")
meta <- merge(a, b, by = "gene")
meta$z_20141 <- sign(meta$logFC_20141) * qnorm(meta$P.Value_20141 / 2, lower.tail = FALSE)
meta$z_7621 <- sign(meta$logFC_7621) * qnorm(meta$P.Value_7621 / 2, lower.tail = FALSE)
meta$meta_z <- (meta$z_20141 + meta$z_7621) / sqrt(2)
meta$meta_p <- 2 * pnorm(abs(meta$meta_z), lower.tail = FALSE)
meta$meta_FDR <- p.adjust(meta$meta_p, method = "BH")
meta <- meta[order(meta$gene), ]
write.table(meta, file.path(out, "human_PD_limma_directional_meta.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
