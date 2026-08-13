#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

root <- normalizePath(".")
counts_path <- file.path(
  root,
  "data_raw/gse150005/GSE150005_Tong_et_al_DopaNeuron_PD_Tox_RAW_COUNT_MATRIX.txt.gz"
)
gmt_path <- file.path(root, "data_raw/pathways/ReactomePathways.gmt")
out <- file.path(root, "results/v2/gse150005_multitoxicant")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(counts_path)) stop("Missing raw count matrix: ", counts_path)
if (!file.exists(gmt_path)) stop("Missing Reactome GMT: ", gmt_path)

raw <- read.delim(
  gzfile(counts_path),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_metadata <- c(
  "Ensembl_Gene_ID", "Gene name", "HGNC symbol", "Gene description"
)
stopifnot(identical(colnames(raw)[seq_along(required_metadata)], required_metadata))

symbols <- trimws(raw[["HGNC symbol"]])
keep_symbol <- !is.na(symbols) & nzchar(symbols)
count_mat <- as.matrix(raw[keep_symbol, -(seq_along(required_metadata)), drop = FALSE])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- toupper(symbols[keep_symbol])

# HGNC symbols duplicated across Ensembl rows are combined before expression
# filtering so every pathway is evaluated at one row per gene.
count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)
genes_before_expression_filter <- nrow(count_mat)
keep_expression <- rowSums(count_mat >= 10L) >= 3L
count_mat_main <- count_mat[keep_expression, , drop = FALSE]

sample_names <- colnames(count_mat_main)
group_labels <- rep(NA_character_, length(sample_names))
group_labels[grepl("DMSO", sample_names)] <- "DMSO"
group_labels[grepl("MS275", sample_names)] <- "MS275"
group_labels[grepl("Paraquat", sample_names)] <- "Paraquat"
group_labels[grepl("_MPP_", sample_names)] <- "MPP"
group_labels[grepl("Ziram", sample_names)] <- "Ziram"
group_labels[grepl("Rotenone", sample_names)] <- "Rotenone"
group_labels[grepl("6-HD", sample_names)] <- "SixOHDA"
group_labels[grepl("Methyl_mercury", sample_names)] <- "Methylmercury"

if (anyNA(group_labels)) {
  stop("Unclassified samples: ", paste(sample_names[is.na(group_labels)], collapse = ", "))
}

group_levels <- c(
  "DMSO", "Rotenone", "MPP", "Paraquat", "SixOHDA", "Ziram",
  "Methylmercury", "MS275"
)
groups <- factor(group_labels, levels = group_levels)
stopifnot(all(table(groups) == 3L))

sample_manifest <- data.frame(
  sample = sample_names,
  group = as.character(groups),
  treatment_role = ifelse(
    groups == "DMSO",
    "reference",
    ifelse(groups == "MS275", "HDAC_inhibitor_comparator", "PD_environmental_toxicant")
  )
)
write.table(
  sample_manifest,
  file.path(out, "GSE150005_sample_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)
contrasts <- makeContrasts(
  Rotenone = Rotenone - DMSO,
  MPP = MPP - DMSO,
  Paraquat = Paraquat - DMSO,
  SixOHDA = SixOHDA - DMSO,
  Ziram = Ziram - DMSO,
  Methylmercury = Methylmercury - DMSO,
  MS275 = MS275 - DMSO,
  levels = design
)

gmt_fields <- strsplit(readLines(gmt_path, warn = FALSE), "\t", fixed = TRUE)
pathway_names <- vapply(gmt_fields, `[[`, character(1), 1L)
pathway_genes <- lapply(gmt_fields, function(x) unique(toupper(x[-c(1L, 2L)])))
names(pathway_genes) <- pathway_names

make_pathway_index <- function(expression_matrix) {
  idx <- lapply(pathway_genes, function(x) which(rownames(expression_matrix) %in% x))
  idx[lengths(idx) >= 10L & lengths(idx) <= 500L]
}

fit_dataset <- function(counts, filter_label) {
  dge <- DGEList(counts = counts, group = groups)
  dge <- calcNormFactors(dge, method = "TMM")
  voom_object <- voom(dge, design, plot = FALSE)
  fit <- eBayes(contrasts.fit(lmFit(voom_object, design), contrasts), robust = TRUE)
  list(
    counts = counts,
    dge = dge,
    voom = voom_object,
    fit = fit,
    pathway_index = make_pathway_index(voom_object$E),
    filter_label = filter_label
  )
}

main <- fit_dataset(count_mat_main, "counts_ge_10_in_ge_3_samples")

screen_rows <- list()
screen_tables <- list()
for (contrast_name in colnames(contrasts)) {
  camera_tab <- camera(
    main$voom,
    index = main$pathway_index,
    design = design,
    contrast = contrasts[, contrast_name],
    inter.gene.cor = NA,
    sort = FALSE
  )
  camera_tab$Pathway <- rownames(camera_tab)
  camera_tab$Contrast <- contrast_name
  camera_tab <- camera_tab[, c(
    "Contrast", "Pathway", "NGenes", "Correlation", "Direction", "PValue", "FDR"
  )]
  camera_tab <- camera_tab[order(camera_tab$PValue, camera_tab$Pathway), ]
  rownames(camera_tab) <- NULL
  screen_tables[[contrast_name]] <- camera_tab
  write.table(
    camera_tab,
    file.path(out, paste0("GSE150005_", contrast_name, "_Reactome_CAMERA.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  gene_tab <- topTable(
    main$fit,
    coef = contrast_name,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )
  gene_tab$Gene <- rownames(gene_tab)
  gene_tab$Contrast <- contrast_name
  gene_tab <- gene_tab[, c(
    "Contrast", "Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"
  )]
  write.table(
    gene_tab,
    gzfile(file.path(out, paste0("GSE150005_", contrast_name, "_limma_all_genes.tsv.gz"))),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  screen_rows[[contrast_name]] <- data.frame(
    Contrast = contrast_name,
    Treatment_role = ifelse(
      contrast_name == "MS275", "HDAC_inhibitor_comparator", "PD_environmental_toxicant"
    ),
    Genes_FDR_0.05 = sum(gene_tab$adj.P.Val <= 0.05, na.rm = TRUE),
    Pathways_FDR_0.05 = sum(camera_tab$FDR <= 0.05, na.rm = TRUE),
    Top_pathway = camera_tab$Pathway[[1]],
    Top_direction = camera_tab$Direction[[1]],
    Top_PValue = camera_tab$PValue[[1]],
    Top_FDR = camera_tab$FDR[[1]],
    stringsAsFactors = FALSE
  )
}

screen_summary <- do.call(rbind, screen_rows)
rownames(screen_summary) <- NULL
write.table(
  screen_summary,
  file.path(out, "GSE150005_contrast_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

mechanism_pathways <- c(
  "Microtubule-dependent trafficking of connexons from Golgi to the plasma membrane",
  "Intraflagellar transport",
  "Cargo trafficking to the periciliary membrane",
  "Carboxyterminal post-translational modifications of tubulin",
  "Hedgehog 'off' state",
  "Aggrephagy",
  "Ribosome-associated quality control",
  "Proteasome assembly",
  "ATF6 (ATF6-alpha) activates chaperone genes",
  "Response of EIF2AK1 (HRI) to heme deficiency"
)
mechanism_summary <- do.call(rbind, lapply(names(screen_tables), function(z) {
  tab <- screen_tables[[z]]
  tab[tab$Pathway %in% mechanism_pathways, , drop = FALSE]
}))
mechanism_summary$Mechanism_family <- ifelse(
  grepl("connexon|Intraflagellar|periciliary|tubulin|Hedgehog", mechanism_summary$Pathway),
  "cilium_Hedgehog_microtubule_trafficking",
  ifelse(
    mechanism_summary$Pathway == "Aggrephagy",
    "aggrephagy",
    ifelse(
      grepl("Ribosome-associated|Proteasome", mechanism_summary$Pathway),
      "proteostasis",
      "integrated_ER_stress_response"
    )
  )
)
mechanism_summary <- mechanism_summary[, c(
  "Mechanism_family", "Contrast", "Pathway", "NGenes", "Correlation",
  "Direction", "PValue", "FDR"
)]
write.table(
  mechanism_summary,
  file.path(out, "GSE150005_exploratory_mechanism_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

proteasome_main <- do.call(rbind, lapply(names(screen_tables), function(z) {
  tab <- screen_tables[[z]]
  ans <- tab[tab$Pathway == "Proteasome assembly", , drop = FALSE]
  ans$PValue_one_sided_down <- ifelse(
    ans$Direction == "Down", ans$PValue / 2, 1 - ans$PValue / 2
  )
  ans
}))
pd_toxicants <- proteasome_main$Contrast != "MS275"
proteasome_main$BH_six_PD_environmental_toxicants <- NA_real_
proteasome_main$BH_six_PD_environmental_toxicants[pd_toxicants] <- p.adjust(
  proteasome_main$PValue_one_sided_down[pd_toxicants], method = "BH"
)
write.table(
  proteasome_main,
  file.path(out, "GSE150005_proteasome_CAMERA_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Sensitivity analysis: replace the declared simple count filter with edgeR's
# design-aware filterByExpr and rerun only the proteasome CAMERA test.
dge_unfiltered <- DGEList(counts = count_mat, group = groups)
keep_filter_by_expr <- filterByExpr(dge_unfiltered, design = design)
sensitivity <- fit_dataset(
  count_mat[keep_filter_by_expr, , drop = FALSE],
  "edgeR_filterByExpr_design"
)
proteasome_sensitivity <- do.call(rbind, lapply(colnames(contrasts), function(z) {
  tab <- camera(
    sensitivity$voom,
    index = list("Proteasome assembly" = sensitivity$pathway_index[["Proteasome assembly"]]),
    design = design,
    contrast = contrasts[, z],
    inter.gene.cor = NA,
    sort = FALSE
  )
  tab$Contrast <- z
  tab$PValue_one_sided_down <- ifelse(
    tab$Direction == "Down", tab$PValue / 2, 1 - tab$PValue / 2
  )
  tab[, c(
    "Contrast", "NGenes", "Correlation", "Direction", "PValue",
    "PValue_one_sided_down"
  )]
}))
rownames(proteasome_sensitivity) <- NULL
write.table(
  proteasome_sensitivity,
  file.path(out, "GSE150005_proteasome_filterByExpr_sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

filter_manifest <- data.frame(
  quantity = c(
    "samples", "symbol_rows_before_expression_filter", "genes_main_filter",
    "genes_filterByExpr_sensitivity", "Reactome_pathways_main"
  ),
  value = c(
    ncol(count_mat), genes_before_expression_filter, nrow(main$counts),
    nrow(sensitivity$counts), length(main$pathway_index)
  )
)
write.table(
  filter_manifest,
  file.path(out, "GSE150005_analysis_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
