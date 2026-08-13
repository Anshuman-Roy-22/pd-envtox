#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(edgeR)
  library(limma)
  library(org.Hs.eg.db)
})

root <- normalizePath(".")
gmt_path <- file.path(root, "data_raw/pathways/ReactomePathways.gmt")
out <- file.path(root, "results/v2/heldout_external_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  line <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(line)) stop("Could not compute SHA256 for ", path)
  strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

expected_checksums <- c(
  Reactome_GMT = "89983d5c1f0af11c52edfeee7323eb425580ac6281d387a528562ab1787ce56b",
  GSE229460_RAW_tar = "d3c5e7a93c12c4aa0ca91389be94fd30f70ae4997e288ff0d55fdd5ab0988eb2",
  GSE287941_counts = "54bc3c7e6a272e7f2c476c066ed1a686dfc95d8156fd7a7f6584af472b5c54b9"
)
raw_229460 <- file.path(root, "data_raw/gse229460/GSE229460_RAW.tar")
counts_287941_path <- file.path(
  root,
  "data_raw/gse287941/GSE287941_RawCounts_AllSamples.txt.gz"
)
stopifnot(
  file.exists(gmt_path), file.exists(raw_229460), file.exists(counts_287941_path),
  identical(sha256(gmt_path), expected_checksums[["Reactome_GMT"]]),
  identical(sha256(raw_229460), expected_checksums[["GSE229460_RAW_tar"]]),
  identical(sha256(counts_287941_path), expected_checksums[["GSE287941_counts"]])
)

gmt_fields <- strsplit(readLines(gmt_path, warn = FALSE), "\t", fixed = TRUE)
pathway_names <- vapply(gmt_fields, `[[`, character(1), 1L)
pathway_genes <- lapply(gmt_fields, function(x) unique(toupper(x[-c(1L, 2L)])))
names(pathway_genes) <- pathway_names
primary_pathways <- c(
  "Microtubule-dependent trafficking of connexons from Golgi to the plasma membrane",
  "Intraflagellar transport",
  "Cargo trafficking to the periciliary membrane",
  "Carboxyterminal post-translational modifications of tubulin",
  "Hedgehog 'off' state"
)
stopifnot(all(primary_pathways %in% names(pathway_genes)))
primary_union_genes <- unique(unlist(pathway_genes[primary_pathways], use.names = FALSE))

map_ensembl_counts <- function(counts, ensembl_ids) {
  base_ids <- sub("\\..*$", "", ensembl_ids)
  counts_by_ensembl <- rowsum(counts, group = base_ids, reorder = TRUE)
  mapping <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db,
    keys = rownames(counts_by_ensembl),
    columns = "SYMBOL",
    keytype = "ENSEMBL"
  ))
  mapping <- mapping[!is.na(mapping$SYMBOL) & nzchar(trimws(mapping$SYMBOL)), ]
  mapping$SYMBOL <- toupper(trimws(mapping$SYMBOL))
  mapping <- unique(mapping)
  mapping <- mapping[order(mapping$ENSEMBL, mapping$SYMBOL), ]
  ambiguous <- names(which(table(mapping$ENSEMBL) > 1L))
  mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
  symbols <- mapping$SYMBOL[match(rownames(counts_by_ensembl), mapping$ENSEMBL)]
  keep <- !is.na(symbols) & nzchar(symbols)
  counts_by_symbol <- rowsum(
    counts_by_ensembl[keep, , drop = FALSE],
    group = symbols[keep],
    reorder = TRUE
  )
  list(
    counts = counts_by_symbol,
    input_rows = length(ensembl_ids),
    unique_ensembl = nrow(counts_by_ensembl),
    mapped_ensembl = sum(keep),
    ambiguous_ensembl = length(ambiguous)
  )
}

read_gse229460 <- function() {
  gsm <- c(
    "GSM7163946", "GSM7163947", "GSM7163948",
    "GSM7163949", "GSM7163950", "GSM7163951"
  )
  files <- c(
    "GSM7163946_01_Ctrl_1_S1.R1.readcounts.tsv.gz",
    "GSM7163947_02_Ctrl_2_S2.R1.readcounts.tsv.gz",
    "GSM7163948_03_Ctrl_3_S3.R1.readcounts.tsv.gz",
    "GSM7163949_09_MPP_1_S9.R1.readcounts.tsv.gz",
    "GSM7163950_06_MPP_2_S6.R1.readcounts.tsv.gz",
    "GSM7163951_07_MPP_3_S7.R1.readcounts.tsv.gz"
  )
  paths <- file.path(root, "data_raw/gse229460", files)
  if (any(!file.exists(paths))) stop("Extract GSE229460_RAW.tar before analysis")
  tabs <- lapply(paths, function(path) {
    read.delim(
      gzfile(path), header = FALSE, col.names = c("Ensembl", "count"),
      stringsAsFactors = FALSE
    )
  })
  ids <- tabs[[1L]]$Ensembl
  stopifnot(all(vapply(tabs, function(x) identical(x$Ensembl, ids), logical(1))))
  counts <- do.call(cbind, lapply(tabs, `[[`, "count"))
  colnames(counts) <- gsm
  storage.mode(counts) <- "double"
  mapped <- map_ensembl_counts(counts, ids)
  manifest <- data.frame(
    accession = "GSE229460",
    GSM = gsm,
    matrix_column = gsm,
    group = c(rep("control", 3L), rep("MPP", 3L)),
    treatment = c(rep("untreated", 3L), rep("10 micromolar MPP+ for 48 hours", 3L)),
    genotype = "wild type",
    model = "differentiated LUHMES",
    stringsAsFactors = FALSE
  )
  list(mapped = mapped, manifest = manifest)
}

read_gse287941 <- function() {
  raw <- read.delim(
    gzfile(counts_287941_path),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  stopifnot(identical(colnames(raw)[1:2], c("Ensembl Gene ID", "Gene Symbol")))
  columns <- c(
    "1-V300071532_L01_114", "2-V300071532_L01_32",
    "3-V300071532_L01_33", "4-V300071532_L01_34",
    "5-V300071532_L01_35", "6-V300071532_L01_36",
    "7-V300071532_L01_37", "8-V300071532_L01_38"
  )
  stopifnot(all(columns %in% colnames(raw)))
  gsm <- c(
    "GSM8755842", "GSM8755843", "GSM8755844", "GSM8755845",
    "GSM8755850", "GSM8755851", "GSM8755852", "GSM8755853"
  )
  counts <- as.matrix(raw[, columns, drop = FALSE])
  storage.mode(counts) <- "double"
  colnames(counts) <- gsm
  mapped <- map_ensembl_counts(counts, raw[["Ensembl Gene ID"]])
  manifest <- data.frame(
    accession = "GSE287941",
    GSM = gsm,
    matrix_column = columns,
    group = c(rep("control", 4L), rep("SixOHDA", 4L)),
    treatment = c(
      rep("untreated", 4L),
      rep("6-OHDA; dose and duration not supplied in GEO metadata", 4L)
    ),
    genotype = "wild type",
    model = "differentiated LUHMES",
    stringsAsFactors = FALSE
  )
  list(mapped = mapped, manifest = manifest)
}

one_sided_down <- function(direction, two_sided_p) {
  ifelse(direction == "Down", two_sided_p / 2, 1 - two_sided_p / 2)
}

analyze_dataset <- function(accession, input) {
  counts <- input$mapped$counts[, input$manifest$GSM, drop = FALSE]
  groups <- factor(input$manifest$group, levels = unique(input$manifest$group))
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  rownames(design) <- input$manifest$GSM
  contrast <- makeContrasts(contrasts = paste0(levels(groups)[2L], "-", levels(groups)[1L]),
                            levels = design)

  dge <- DGEList(counts = counts, group = groups)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  voom_object <- voom(dge, design, plot = FALSE)
  fit <- eBayes(
    contrasts.fit(lmFit(voom_object, design), contrast),
    robust = TRUE
  )

  index_for <- function(genes) which(rownames(voom_object$E) %in% genes)
  constituent_indices <- lapply(pathway_genes[primary_pathways], index_for)
  union_index <- index_for(primary_union_genes)
  union_camera <- camera(
    voom_object,
    index = list(cilium_Hedgehog_microtubule_trafficking = union_index),
    design = design,
    contrast = contrast[, 1L],
    inter.gene.cor = NA,
    sort = FALSE
  )
  constituent_camera <- camera(
    voom_object,
    index = constituent_indices,
    design = design,
    contrast = contrast[, 1L],
    inter.gene.cor = NA,
    sort = FALSE
  )
  constituent <- data.frame(
    accession = accession,
    pathway = primary_pathways,
    analyzed_genes = constituent_camera$NGenes,
    intergene_correlation = constituent_camera$Correlation,
    CAMERA_direction = constituent_camera$Direction,
    CAMERA_p_two_sided = constituent_camera$PValue,
    CAMERA_p_one_sided_down = one_sided_down(
      constituent_camera$Direction,
      constituent_camera$PValue
    ),
    mean_moderated_t = vapply(
      constituent_indices,
      function(idx) mean(fit$t[idx, 1L]),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  constituent$negative_mean_moderated_t <- constituent$mean_moderated_t < 0
  negative_constituents <- sum(constituent$negative_mean_moderated_t)
  p_down <- one_sided_down(union_camera$Direction[[1L]], union_camera$PValue[[1L]])
  primary_style_criteria <-
    union_camera$NGenes[[1L]] >= 10L &&
    union_camera$Direction[[1L]] == "Down" &&
    p_down < 0.05 &&
    negative_constituents >= 3L
  summary <- data.frame(
    accession = accession,
    contrast = paste0(levels(groups)[2L], "_minus_", levels(groups)[1L]),
    model = input$manifest$model[[1L]],
    n_control = sum(groups == levels(groups)[1L]),
    n_exposed = sum(groups == levels(groups)[2L]),
    genes_after_filterByExpr = nrow(voom_object$E),
    analyzed_mechanism_genes = union_camera$NGenes[[1L]],
    intergene_correlation = union_camera$Correlation[[1L]],
    CAMERA_direction = union_camera$Direction[[1L]],
    CAMERA_p_two_sided = union_camera$PValue[[1L]],
    CAMERA_p_one_sided_down = p_down,
    constituent_pathways_with_negative_mean_t = negative_constituents,
    passes_primary_style_criteria_as_sensitivity = primary_style_criteria,
    stringsAsFactors = FALSE
  )
  gene_results <- topTable(
    fit, coef = 1L, number = Inf, sort.by = "none", adjust.method = "BH"
  )
  gene_results$HGNC_symbol <- rownames(gene_results)
  gene_results$in_primary_union <- gene_results$HGNC_symbol %in% primary_union_genes
  gene_results <- gene_results[, c(
    "HGNC_symbol", "in_primary_union", "logFC", "AveExpr", "t",
    "P.Value", "adj.P.Val", "B"
  )]
  mapping <- data.frame(
    accession = accession,
    quantity = c(
      "input_rows", "unique_unversioned_Ensembl_rows", "mapped_Ensembl_rows",
      "Ensembl_rows_with_multiple_current_symbols", "unique_HGNC_symbols"
    ),
    value = c(
      input$mapped$input_rows, input$mapped$unique_ensembl,
      input$mapped$mapped_ensembl, input$mapped$ambiguous_ensembl,
      nrow(input$mapped$counts)
    ),
    stringsAsFactors = FALSE
  )
  qc <- data.frame(
    accession = accession,
    GSM = input$manifest$GSM,
    group = input$manifest$group,
    library_size_after_filter = dge$samples$lib.size,
    TMM_normalization_factor = dge$samples$norm.factors,
    stringsAsFactors = FALSE
  )
  list(
    summary = summary, constituent = constituent, genes = gene_results,
    mapping = mapping, qc = qc
  )
}

inputs <- list(
  GSE229460 = read_gse229460(),
  GSE287941 = read_gse287941()
)
analyses <- lapply(names(inputs), function(x) analyze_dataset(x, inputs[[x]]))
names(analyses) <- names(inputs)

write.table(
  do.call(rbind, lapply(inputs, `[[`, "manifest")),
  file.path(out, "external_sensitivity_sample_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  do.call(rbind, lapply(analyses, `[[`, "summary")),
  file.path(out, "external_sensitivity_primary_union_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  do.call(rbind, lapply(analyses, `[[`, "constituent")),
  file.path(out, "external_sensitivity_constituent_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  do.call(rbind, lapply(analyses, `[[`, "mapping")),
  file.path(out, "external_sensitivity_identifier_mapping.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  do.call(rbind, lapply(analyses, `[[`, "qc")),
  file.path(out, "external_sensitivity_sample_QC.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
for (accession in names(analyses)) {
  write.table(
    analyses[[accession]]$genes,
    gzfile(file.path(out, paste0(accession, "_limma_all_genes.tsv.gz"))),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

orgdb_path <- system.file("extdata", "org.Hs.eg.sqlite", package = "org.Hs.eg.db")
manifest <- data.frame(
  quantity = c(
    "analysis_role", "primary_result_can_be_rescued", "Reactome_GMT_SHA256",
    "GSE229460_RAW_tar_SHA256", "GSE287941_counts_SHA256",
    "org.Hs.eg.db_version", "org.Hs.eg.sqlite_SHA256", "limma_version",
    "edgeR_version"
  ),
  value = c(
    "declared_external_sensitivity", "FALSE", sha256(gmt_path),
    sha256(raw_229460), sha256(counts_287941_path),
    as.character(packageVersion("org.Hs.eg.db")), sha256(orgdb_path),
    as.character(packageVersion("limma")), as.character(packageVersion("edgeR"))
  ),
  stringsAsFactors = FALSE
)
write.table(
  manifest,
  file.path(out, "external_sensitivity_analysis_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
session_lines <- sub("[[:space:]]+$", "", capture.output(sessionInfo()))
writeLines(session_lines, file.path(out, "sessionInfo.txt"))

print(do.call(rbind, lapply(analyses, `[[`, "summary")))
