#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(limma)
  library(org.Hs.eg.db)
})

root <- normalizePath(".")
matrix_path <- file.path(root, "data_raw/gse4773/GSE4773_series_matrix.txt.gz")
gmt_path <- file.path(root, "data_raw/pathways/ReactomePathways.gmt")
out <- file.path(root, "results/v2/heldout_external_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  line <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(line)) stop("Could not compute SHA256 for ", path)
  strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
stopifnot(
  file.exists(matrix_path), file.exists(gmt_path),
  identical(
    sha256(matrix_path),
    "5b944e95d019f4fad73dc8a6175050a17d64abc54ac6843c26989f2863b7a638"
  ),
  identical(
    sha256(gmt_path),
    "89983d5c1f0af11c52edfeee7323eb425580ac6281d387a528562ab1787ce56b"
  )
)

lines <- readLines(gzfile(matrix_path), warn = FALSE)
table_begin <- match("!series_matrix_table_begin", lines)
table_end <- match("!series_matrix_table_end", lines)
if (is.na(table_begin) || is.na(table_end) || table_end <= table_begin + 1L) {
  stop("Could not locate series matrix table")
}
matrix_tab <- read.delim(
  textConnection(lines[(table_begin + 1L):(table_end - 1L)]),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(colnames(matrix_tab)[[1L]] == "ID_REF")
probe_ids <- matrix_tab[["ID_REF"]]
expression_raw <- as.matrix(matrix_tab[, -1L, drop = FALSE])
storage.mode(expression_raw) <- "double"
if (anyNA(expression_raw) || any(expression_raw < 0)) stop("Invalid MAS5 signal values")

# GEO documents these as MAS5/GCOS signal values. Log2(x + 1) is applied
# before limma. Ambiguous probe-to-gene mappings are excluded.
expression_log2 <- log2(expression_raw + 1)
rownames(expression_log2) <- probe_ids
probe_to_entrez <- suppressMessages(AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  columns = "ENTREZID",
  keytype = "PROBEID"
))
probe_to_entrez <- unique(probe_to_entrez[!is.na(probe_to_entrez$ENTREZID), ])
unambiguous_probes <- names(which(table(probe_to_entrez$PROBEID) == 1L))
probe_to_entrez <- probe_to_entrez[probe_to_entrez$PROBEID %in% unambiguous_probes, ]
entrez_to_symbol <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(probe_to_entrez$ENTREZID),
  columns = "SYMBOL",
  keytype = "ENTREZID"
))
entrez_to_symbol <- unique(
  entrez_to_symbol[!is.na(entrez_to_symbol$SYMBOL) &
                     nzchar(trimws(entrez_to_symbol$SYMBOL)), ]
)
entrez_to_symbol$SYMBOL <- toupper(trimws(entrez_to_symbol$SYMBOL))
unambiguous_entrez <- names(which(table(entrez_to_symbol$ENTREZID) == 1L))
entrez_to_symbol <- entrez_to_symbol[
  entrez_to_symbol$ENTREZID %in% unambiguous_entrez,
]

probe_map <- merge(probe_to_entrez, entrez_to_symbol, by = "ENTREZID")
probe_map <- probe_map[match(probe_ids, probe_map$PROBEID, nomatch = 0L), ]
probe_index <- match(probe_map$PROBEID, rownames(expression_log2))
probe_mean <- rowMeans(expression_log2[probe_index, , drop = FALSE])
selection <- data.frame(
  probe = probe_map$PROBEID,
  symbol = probe_map$SYMBOL,
  all_sample_mean_log2_signal = probe_mean,
  index = probe_index,
  stringsAsFactors = FALSE
)
selection <- selection[order(
  selection$symbol,
  -selection$all_sample_mean_log2_signal,
  selection$probe
), ]
selection <- selection[!duplicated(selection$symbol), ]
expression_gene <- expression_log2[selection$index, , drop = FALSE]
rownames(expression_gene) <- selection$symbol

samples <- data.frame(
  GSM = c(
    "GSM107850", "GSM107851", "GSM107852",
    "GSM107853", "GSM107854", "GSM107855",
    "GSM107856", "GSM107857", "GSM107858",
    "GSM107859", "GSM107860", "GSM107861",
    "GSM107862", "GSM107863", "GSM107864",
    "GSM107865", "GSM107866", "GSM107867"
  ),
  group = rep(c("control_1w", "rotenone_1w", "control_2w", "rotenone_2w",
                "control_4w", "rotenone_4w"), each = 3L),
  exposure = rep(c("vehicle", "5 nM rotenone"), each = 3L, times = 3L),
  duration_weeks = rep(c(1L, 1L, 2L, 2L, 4L, 4L), each = 3L),
  model = "SK-N-MC neuroblastoma line described by submitter as dopaminergic",
  confidence = "borderline sensitivity only",
  stringsAsFactors = FALSE
)
stopifnot(all(samples$GSM %in% colnames(expression_gene)))
expression_selected <- expression_gene[, samples$GSM, drop = FALSE]

group <- factor(samples$group, levels = c(
  "control_1w", "rotenone_1w", "control_2w", "rotenone_2w",
  "control_4w", "rotenone_4w"
))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- samples$GSM
contrasts <- makeContrasts(
  rotenone_1w_minus_control_1w = rotenone_1w - control_1w,
  rotenone_2w_minus_control_2w = rotenone_2w - control_2w,
  rotenone_4w_minus_control_4w = rotenone_4w - control_4w,
  levels = design
)
fit <- eBayes(
  contrasts.fit(lmFit(expression_selected, design), contrasts),
  robust = TRUE
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
primary_union_genes <- unique(unlist(pathway_genes[primary_pathways], use.names = FALSE))
index_for <- function(genes) which(rownames(expression_selected) %in% genes)
constituent_indices <- lapply(pathway_genes[primary_pathways], index_for)
union_index <- index_for(primary_union_genes)

one_sided_down <- function(direction, two_sided_p) {
  ifelse(direction == "Down", two_sided_p / 2, 1 - two_sided_p / 2)
}

summary_rows <- list()
constituent_rows <- list()
for (contrast_name in colnames(contrasts)) {
  union_camera <- camera(
    expression_selected,
    index = list(cilium_Hedgehog_microtubule_trafficking = union_index),
    design = design,
    contrast = contrasts[, contrast_name],
    inter.gene.cor = NA,
    sort = FALSE
  )
  constituent_camera <- camera(
    expression_selected,
    index = constituent_indices,
    design = design,
    contrast = contrasts[, contrast_name],
    inter.gene.cor = NA,
    sort = FALSE
  )
  mean_t <- vapply(
    constituent_indices,
    function(idx) mean(fit$t[idx, contrast_name]),
    numeric(1)
  )
  negative_constituents <- sum(mean_t < 0)
  p_down <- one_sided_down(union_camera$Direction[[1L]], union_camera$PValue[[1L]])
  summary_rows[[contrast_name]] <- data.frame(
    accession = "GSE4773",
    contrast = contrast_name,
    confidence = "borderline sensitivity only",
    analyzed_mechanism_genes = union_camera$NGenes[[1L]],
    intergene_correlation = union_camera$Correlation[[1L]],
    CAMERA_direction = union_camera$Direction[[1L]],
    CAMERA_p_two_sided = union_camera$PValue[[1L]],
    CAMERA_p_one_sided_down = p_down,
    constituent_pathways_with_negative_mean_t = negative_constituents,
    passes_primary_style_criteria_as_sensitivity =
      union_camera$NGenes[[1L]] >= 10L &&
      union_camera$Direction[[1L]] == "Down" && p_down < 0.05 &&
      negative_constituents >= 3L,
    stringsAsFactors = FALSE
  )
  constituent_rows[[contrast_name]] <- data.frame(
    accession = "GSE4773",
    contrast = contrast_name,
    pathway = primary_pathways,
    analyzed_genes = constituent_camera$NGenes,
    intergene_correlation = constituent_camera$Correlation,
    CAMERA_direction = constituent_camera$Direction,
    CAMERA_p_two_sided = constituent_camera$PValue,
    CAMERA_p_one_sided_down = one_sided_down(
      constituent_camera$Direction,
      constituent_camera$PValue
    ),
    mean_moderated_t = mean_t,
    negative_mean_moderated_t = mean_t < 0,
    stringsAsFactors = FALSE
  )
  gene_results <- topTable(
    fit,
    coef = contrast_name,
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
  write.table(
    gene_results,
    gzfile(file.path(out, paste0("GSE4773_", contrast_name, "_limma_all_genes.tsv.gz"))),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

write.table(
  samples,
  file.path(out, "GSE4773_sample_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  do.call(rbind, summary_rows),
  file.path(out, "GSE4773_primary_union_sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  do.call(rbind, constituent_rows),
  file.path(out, "GSE4773_constituent_sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

hgu_db <- system.file("extdata", "hgu133plus2.sqlite", package = "hgu133plus2.db")
org_db <- system.file("extdata", "org.Hs.eg.sqlite", package = "org.Hs.eg.db")
manifest <- data.frame(
  quantity = c(
    "analysis_role", "primary_result_can_be_rescued", "input_probes",
    "unambiguously_mapped_symbols_after_probe_selection",
    "signal_processing", "probe_selection", "series_matrix_SHA256",
    "Reactome_GMT_SHA256", "hgu133plus2.db_version",
    "hgu133plus2.sqlite_SHA256", "org.Hs.eg.db_version",
    "org.Hs.eg.sqlite_SHA256", "limma_version"
  ),
  value = c(
    "borderline_external_sensitivity", "FALSE", nrow(expression_raw),
    nrow(expression_selected), "deposited MAS5 or GCOS signal; log2(x + 1)",
    "highest all-sample mean per current HGNC symbol; lexical probe tie-break",
    sha256(matrix_path), sha256(gmt_path),
    as.character(packageVersion("hgu133plus2.db")), sha256(hgu_db),
    as.character(packageVersion("org.Hs.eg.db")), sha256(org_db),
    as.character(packageVersion("limma"))
  ),
  stringsAsFactors = FALSE
)
write.table(
  manifest,
  file.path(out, "GSE4773_analysis_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
session_lines <- sub("[[:space:]]+$", "", capture.output(sessionInfo()))
writeLines(session_lines, file.path(out, "GSE4773_sessionInfo.txt"))

print(do.call(rbind, summary_rows))
