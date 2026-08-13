#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(edgeR)
  library(limma)
  library(org.Hs.eg.db)
})

root <- normalizePath(".")
raw_dir <- file.path(root, "data_raw/gse196190")
gmt_path <- file.path(root, "data_raw/pathways/ReactomePathways.gmt")
out <- file.path(root, "results/v2/gse196190_heldout")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

reactome_sha256_expected <-
  "89983d5c1f0af11c52edfeee7323eb425580ac6281d387a528562ab1787ce56b"
raw_tar_sha256_expected <-
  "1a20222c1836188b8d0ee9617460372dcf538a5d0bd7007b306eeb7293184b1d"

sha256 <- function(path) {
  line <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(line)) stop("Could not compute SHA256 for ", path)
  strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

raw_tar <- file.path(raw_dir, "GSE196190_RAW.tar")
if (!file.exists(raw_tar)) stop("Missing raw archive: ", raw_tar)
if (!file.exists(gmt_path)) stop("Missing Reactome GMT: ", gmt_path)
stopifnot(identical(sha256(raw_tar), raw_tar_sha256_expected))
stopifnot(identical(sha256(gmt_path), reactome_sha256_expected))

samples <- data.frame(
  GSM = c(
    "GSM5862314", "GSM5862315", "GSM5862316",
    "GSM5862317", "GSM5862318", "GSM5862319"
  ),
  file = c(
    "GSM5862314_16C_MPP_0_t74_S6.genes.results.gz",
    "GSM5862315_17C_MPP_0_t74_S24.genes.results.gz",
    "GSM5862316_18C_MPP_0_t74_S42.genes.results.gz",
    "GSM5862317_19C_MPP_100uM_t74_S7.genes.results.gz",
    "GSM5862318_20C_MPP_100uM_t74_S25.genes.results.gz",
    "GSM5862319_21C_MPP_100uM_t74_S43.genes.results.gz"
  ),
  group = c(rep("control_t74", 3L), rep("MPP_100uM_t74", 3L)),
  replicate_series = c("S6_S7", "S24_S25", "S42_S43", "S6_S7", "S24_S25", "S42_S43"),
  role = c(rep("reference", 3L), rep("primary_treatment", 3L)),
  GEO_endpoint_label = "t74",
  paper_endpoint_label = "72 hours",
  stringsAsFactors = FALSE
)
samples$path <- file.path(raw_dir, samples$file)
if (any(!file.exists(samples$path))) {
  stop("Missing locked sample file(s): ", paste(samples$file[!file.exists(samples$path)], collapse = ", "))
}

read_rsem_counts <- function(path) {
  tab <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  required <- c(
    "gene_id", "transcript_id(s)", "length", "effective_length",
    "expected_count", "TPM", "FPKM"
  )
  if (!identical(colnames(tab), required)) stop("Unexpected RSEM columns in ", path)
  if (anyDuplicated(tab$gene_id)) stop("Duplicated Ensembl rows in ", path)
  tab[, c("gene_id", "expected_count")]
}

rsem <- lapply(samples$path, read_rsem_counts)
reference_ids <- rsem[[1L]]$gene_id
if (!all(vapply(rsem, function(x) identical(x$gene_id, reference_ids), logical(1)))) {
  stop("RSEM files do not contain identical ordered gene identifiers")
}

counts_ensembl <- do.call(cbind, lapply(rsem, `[[`, "expected_count"))
rownames(counts_ensembl) <- reference_ids
colnames(counts_ensembl) <- samples$GSM
storage.mode(counts_ensembl) <- "double"
if (anyNA(counts_ensembl) || any(counts_ensembl < 0)) stop("Invalid expected counts")

# Resolve any one-to-many current OrgDb mappings deterministically by choosing
# the lexicographically first uppercase HGNC symbol for each Ensembl identifier.
mapping_all <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(reference_ids),
  columns = "SYMBOL",
  keytype = "ENSEMBL"
))
mapping_all <- mapping_all[!is.na(mapping_all$SYMBOL) & nzchar(trimws(mapping_all$SYMBOL)), ]
mapping_all$SYMBOL <- toupper(trimws(mapping_all$SYMBOL))
mapping_all <- unique(mapping_all)
mapping_all <- mapping_all[order(mapping_all$ENSEMBL, mapping_all$SYMBOL), ]
ambiguous_ensembl <- names(which(table(mapping_all$ENSEMBL) > 1L))
mapping_one <- mapping_all[!duplicated(mapping_all$ENSEMBL), ]
mapped_symbols <- mapping_one$SYMBOL[match(reference_ids, mapping_one$ENSEMBL)]
has_symbol <- !is.na(mapped_symbols) & nzchar(mapped_symbols)

counts_symbol <- rowsum(
  counts_ensembl[has_symbol, , drop = FALSE],
  group = mapped_symbols[has_symbol],
  reorder = TRUE
)

groups <- factor(samples$group, levels = c("control_t74", "MPP_100uM_t74"))
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)
rownames(design) <- samples$GSM
contrast <- makeContrasts(
  MPP_100uM_t74_minus_control_t74 = MPP_100uM_t74 - control_t74,
  levels = design
)

fit_counts <- function(counts, keep, filter_name) {
  dge <- DGEList(counts = counts[, samples$GSM, drop = FALSE], group = groups)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  voom_object <- voom(dge, design, plot = FALSE)
  fit <- eBayes(
    contrasts.fit(lmFit(voom_object, design), contrast),
    robust = TRUE
  )
  list(dge = dge, voom = voom_object, fit = fit, filter_name = filter_name)
}

dge_unfiltered <- DGEList(counts = counts_symbol, group = groups)
keep_primary <- filterByExpr(dge_unfiltered, design = design)
main <- fit_counts(counts_symbol, keep_primary, "edgeR_filterByExpr_locked_design")

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
secondary_pathway <- "Aggrephagy"
required_pathways <- c(primary_pathways, secondary_pathway)
if (!all(required_pathways %in% names(pathway_genes))) {
  stop("Frozen Reactome pathway name(s) absent from GMT")
}

index_for <- function(gene_set, expression_matrix = main$voom$E) {
  which(rownames(expression_matrix) %in% unique(toupper(gene_set)))
}

run_camera <- function(fitted, index_list) {
  camera(
    fitted$voom,
    index = index_list,
    design = design,
    contrast = contrast[, 1L],
    inter.gene.cor = NA,
    sort = FALSE
  )
}

one_sided_down <- function(direction, two_sided_p) {
  ifelse(direction == "Down", two_sided_p / 2, 1 - two_sided_p / 2)
}

primary_union_genes <- unique(unlist(pathway_genes[primary_pathways], use.names = FALSE))
primary_union_index <- index_for(primary_union_genes)
union_camera <- run_camera(
  main,
  list(cilium_Hedgehog_microtubule_trafficking = primary_union_index)
)
union_direction <- union_camera$Direction[[1L]]
union_p_two <- union_camera$PValue[[1L]]
union_p_down <- one_sided_down(union_direction, union_p_two)

constituent_indices <- lapply(pathway_genes[primary_pathways], index_for)
constituent_camera <- run_camera(main, constituent_indices)
constituent_results <- data.frame(
  pathway = primary_pathways,
  membership_genes = lengths(pathway_genes[primary_pathways]),
  mapped_before_expression_filter = vapply(
    pathway_genes[primary_pathways],
    function(x) sum(rownames(counts_symbol) %in% x),
    integer(1)
  ),
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
    function(idx) mean(main$fit$t[idx, 1L]),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)
constituent_results$CAMERA_BH_five_two_sided <- p.adjust(
  constituent_results$CAMERA_p_two_sided,
  method = "BH"
)
constituent_results$negative_mean_moderated_t <- constituent_results$mean_moderated_t < 0
negative_constituents <- sum(constituent_results$negative_mean_moderated_t)

primary_union_significant <-
  union_camera$NGenes[[1L]] >= 10L &&
  union_direction == "Down" &&
  union_p_down < 0.05
primary_confirmed <- primary_union_significant && negative_constituents >= 3L
primary_outcome <- if (primary_confirmed) {
  "CONFIRMED"
} else if (
  (union_direction == "Down" && union_p_down >= 0.05 && union_p_down < 0.10) ||
  (primary_union_significant && negative_constituents < 3L)
) {
  "PARTIAL_SUPPORT"
} else {
  "NOT_CONFIRMED"
}

primary_results <- data.frame(
  accession = "GSE196190",
  contrast = "MPP_100uM_t74_minus_control_t74",
  mechanism = "cilium_Hedgehog_microtubule_trafficking",
  membership_genes = length(primary_union_genes),
  mapped_before_expression_filter = sum(rownames(counts_symbol) %in% primary_union_genes),
  analyzed_genes = union_camera$NGenes[[1L]],
  intergene_correlation = union_camera$Correlation[[1L]],
  CAMERA_direction = union_direction,
  CAMERA_p_two_sided = union_p_two,
  CAMERA_p_one_sided_down = union_p_down,
  constituent_pathways_with_negative_mean_t = negative_constituents,
  criterion_at_least_10_genes = union_camera$NGenes[[1L]] >= 10L,
  criterion_direction_down = union_direction == "Down",
  criterion_one_sided_p_lt_0_05 = union_p_down < 0.05,
  criterion_at_least_3_negative_constituents = negative_constituents >= 3L,
  frozen_outcome = primary_outcome,
  stringsAsFactors = FALSE
)

# The secondary test is run only when the primary frozen gate passes.
if (primary_confirmed) {
  secondary_index <- index_for(pathway_genes[[secondary_pathway]])
  secondary_camera <- run_camera(main, list(Aggrephagy = secondary_index))
  secondary_p_down <- one_sided_down(
    secondary_camera$Direction[[1L]],
    secondary_camera$PValue[[1L]]
  )
  secondary_confirmed <-
    secondary_camera$NGenes[[1L]] >= 10L &&
    secondary_camera$Direction[[1L]] == "Down" &&
    secondary_p_down < 0.05
  secondary_results <- data.frame(
    gate_status = "OPEN_PRIMARY_CONFIRMED",
    mechanism = "Aggrephagy",
    analyzed_genes = secondary_camera$NGenes[[1L]],
    intergene_correlation = secondary_camera$Correlation[[1L]],
    CAMERA_direction = secondary_camera$Direction[[1L]],
    CAMERA_p_two_sided = secondary_camera$PValue[[1L]],
    CAMERA_p_one_sided_down = secondary_p_down,
    frozen_outcome = ifelse(secondary_confirmed, "CONFIRMED", "NOT_CONFIRMED"),
    stringsAsFactors = FALSE
  )
} else {
  secondary_results <- data.frame(
    gate_status = "CLOSED_PRIMARY_NOT_CONFIRMED",
    mechanism = "Aggrephagy",
    analyzed_genes = NA_integer_,
    intergene_correlation = NA_real_,
    CAMERA_direction = NA_character_,
    CAMERA_p_two_sided = NA_real_,
    CAMERA_p_one_sided_down = NA_real_,
    frozen_outcome = "NOT_TESTED_GATE_CLOSED",
    stringsAsFactors = FALSE
  )
}

# Frozen leave-one-pathway-out union sensitivity.
leave_one_out_indices <- lapply(primary_pathways, function(omitted) {
  retained <- setdiff(primary_pathways, omitted)
  index_for(unique(unlist(pathway_genes[retained], use.names = FALSE)))
})
names(leave_one_out_indices) <- primary_pathways
leave_one_out_camera <- run_camera(main, leave_one_out_indices)
leave_one_out_results <- data.frame(
  omitted_pathway = primary_pathways,
  analyzed_genes = leave_one_out_camera$NGenes,
  intergene_correlation = leave_one_out_camera$Correlation,
  CAMERA_direction = leave_one_out_camera$Direction,
  CAMERA_p_two_sided = leave_one_out_camera$PValue,
  CAMERA_p_one_sided_down = one_sided_down(
    leave_one_out_camera$Direction,
    leave_one_out_camera$PValue
  ),
  stringsAsFactors = FALSE
)

# Frozen alternative count-filter sensitivity.
keep_simple <- rowSums(counts_symbol >= 10) >= 3L
simple <- fit_counts(counts_symbol, keep_simple, "counts_ge_10_in_ge_3_samples")
simple_union_index <- index_for(primary_union_genes, simple$voom$E)
simple_camera <- run_camera(
  simple,
  list(cilium_Hedgehog_microtubule_trafficking = simple_union_index)
)
simple_filter_results <- data.frame(
  filter = simple$filter_name,
  retained_genes = nrow(simple$voom$E),
  analyzed_mechanism_genes = simple_camera$NGenes[[1L]],
  intergene_correlation = simple_camera$Correlation[[1L]],
  CAMERA_direction = simple_camera$Direction[[1L]],
  CAMERA_p_two_sided = simple_camera$PValue[[1L]],
  CAMERA_p_one_sided_down = one_sided_down(
    simple_camera$Direction[[1L]],
    simple_camera$PValue[[1L]]
  ),
  stringsAsFactors = FALSE
)

# Frozen alternative valid handling of pairing/batch. The submitter filenames
# group the samples into three consistent series (S6/S7, S24/S25, S42/S43),
# although GEO does not expose a formal pairing field. The locked unpaired model
# remains primary; this metadata-derived block model is sensitivity-only.
replicate_series <- factor(samples$replicate_series)
design_block <- model.matrix(~ replicate_series + groups)
rownames(design_block) <- samples$GSM
block_contrast <- rep(0, ncol(design_block))
names(block_contrast) <- colnames(design_block)
block_coefficient <- grep("^groupsMPP_100uM_t74$", colnames(design_block))
if (length(block_coefficient) != 1L) stop("Could not identify block-model treatment coefficient")
block_contrast[block_coefficient] <- 1
dge_block_unfiltered <- DGEList(counts = counts_symbol, group = groups)
keep_block <- filterByExpr(dge_block_unfiltered, design = design_block)
dge_block <- dge_block_unfiltered[keep_block, , keep.lib.sizes = FALSE]
dge_block <- calcNormFactors(dge_block, method = "TMM")
voom_block <- voom(dge_block, design_block, plot = FALSE)
fit_block <- eBayes(
  contrasts.fit(lmFit(voom_block, design_block), block_contrast),
  robust = TRUE
)
block_union_index <- which(rownames(voom_block$E) %in% primary_union_genes)
block_constituent_indices <- lapply(
  pathway_genes[primary_pathways],
  function(x) which(rownames(voom_block$E) %in% x)
)
block_camera <- camera(
  voom_block,
  index = list(cilium_Hedgehog_microtubule_trafficking = block_union_index),
  design = design_block,
  contrast = block_contrast,
  inter.gene.cor = NA,
  sort = FALSE
)
block_p_down <- one_sided_down(
  block_camera$Direction[[1L]],
  block_camera$PValue[[1L]]
)
block_negative_constituents <- sum(vapply(
  block_constituent_indices,
  function(idx) mean(fit_block$t[idx, 1L]) < 0,
  logical(1)
))
block_sensitivity_results <- data.frame(
  model = "replicate_series_block_inferred_from_submitter_filenames",
  retained_genes = nrow(voom_block$E),
  analyzed_mechanism_genes = block_camera$NGenes[[1L]],
  intergene_correlation = block_camera$Correlation[[1L]],
  CAMERA_direction = block_camera$Direction[[1L]],
  CAMERA_p_two_sided = block_camera$PValue[[1L]],
  CAMERA_p_one_sided_down = block_p_down,
  constituent_pathways_with_negative_mean_t = block_negative_constituents,
  sensitivity_only = TRUE,
  stringsAsFactors = FALSE
)

gene_results <- topTable(
  main$fit,
  coef = 1L,
  number = Inf,
  sort.by = "none",
  adjust.method = "BH"
)
gene_results$HGNC_symbol <- rownames(gene_results)
gene_results$in_primary_union <- gene_results$HGNC_symbol %in% primary_union_genes
gene_results <- gene_results[, c(
  "HGNC_symbol", "in_primary_union", "logFC", "AveExpr", "t",
  "P.Value", "adj.P.Val", "B"
)]

pca <- prcomp(t(main$voom$E), scale. = FALSE)
pca_scores <- data.frame(
  GSM = rownames(pca$x),
  group = samples$group[match(rownames(pca$x), samples$GSM)],
  PC1 = pca$x[, 1L],
  PC2 = pca$x[, 2L],
  stringsAsFactors = FALSE
)

sample_qc <- data.frame(
  GSM = samples$GSM,
  group = samples$group,
  expected_count_library_size_before_filter = colSums(counts_symbol),
  expected_count_library_size_after_filter = main$dge$samples$lib.size,
  TMM_normalization_factor = main$dge$samples$norm.factors,
  stringsAsFactors = FALSE
)

orgdb_path <- system.file("extdata", "org.Hs.eg.sqlite", package = "org.Hs.eg.db")
orgdb_metadata <- AnnotationDbi::metadata(org.Hs.eg.db)
orgdb_ensource_date <- orgdb_metadata$value[orgdb_metadata$name == "ENSOURCEDATE"]
mapping_stats <- data.frame(
  quantity = c(
    "input_Ensembl_rows", "mapped_Ensembl_rows", "unmapped_Ensembl_rows",
    "Ensembl_rows_with_multiple_symbols", "unique_HGNC_symbols",
    "genes_after_filterByExpr", "org.Hs.eg.db_version",
    "org.Hs.eg.db_Ensembl_source_date", "org.Hs.eg.sqlite_SHA256"
  ),
  value = c(
    nrow(counts_ensembl), sum(has_symbol), sum(!has_symbol),
    length(ambiguous_ensembl), nrow(counts_symbol), nrow(main$voom$E),
    as.character(packageVersion("org.Hs.eg.db")), orgdb_ensource_date,
    sha256(orgdb_path)
  ),
  stringsAsFactors = FALSE
)

analysis_manifest <- data.frame(
  quantity = c(
    "analysis_status", "primary_frozen_outcome", "samples",
    "control_samples", "MPP_samples", "primary_filter",
    "genes_before_filter", "genes_after_filter", "Reactome_GMT_SHA256",
    "GSE196190_RAW_tar_SHA256", "limma_version", "edgeR_version",
    "AnnotationDbi_version", "org.Hs.eg.db_version"
  ),
  value = c(
    "completed", primary_outcome, nrow(samples), sum(groups == "control_t74"),
    sum(groups == "MPP_100uM_t74"), main$filter_name, nrow(counts_symbol),
    nrow(main$voom$E), sha256(gmt_path), sha256(raw_tar),
    as.character(packageVersion("limma")), as.character(packageVersion("edgeR")),
    as.character(packageVersion("AnnotationDbi")),
    as.character(packageVersion("org.Hs.eg.db"))
  ),
  stringsAsFactors = FALSE
)

samples_out <- samples[, setdiff(colnames(samples), "path")]
write.table(samples_out, file.path(out, "GSE196190_sample_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_qc, file.path(out, "GSE196190_sample_QC.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pca_scores, file.path(out, "GSE196190_PCA_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mapping_stats, file.path(out, "GSE196190_identifier_mapping_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(primary_results, file.path(out, "GSE196190_primary_mechanism_result.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(constituent_results, file.path(out, "GSE196190_constituent_pathway_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(secondary_results, file.path(out, "GSE196190_gated_aggrephagy_result.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(leave_one_out_results, file.path(out, "GSE196190_leave_one_pathway_out.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(simple_filter_results, file.path(out, "GSE196190_simple_filter_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(block_sensitivity_results,
            file.path(out, "GSE196190_replicate_block_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(analysis_manifest, file.path(out, "GSE196190_analysis_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  gene_results,
  gzfile(file.path(out, "GSE196190_limma_all_genes.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
session_lines <- sub("[[:space:]]+$", "", capture.output(sessionInfo()))
writeLines(session_lines, file.path(out, "sessionInfo.txt"))

message("Primary frozen outcome: ", primary_outcome)
message(
  "Union CAMERA: ", union_direction,
  ", one-sided-down p=", signif(union_p_down, 5),
  ", negative constituent means=", negative_constituents, "/5"
)
