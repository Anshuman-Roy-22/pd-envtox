# scripts/02_clean_gse17542_microarray.R

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
})

infile_candidates <- c(
  "data_intermediate/gse17542_eset_unzipped.rds",
  "data_intermediate/gse17542_eset.rds"
)
infile  <- infile_candidates[file.exists(infile_candidates)][1]
outfile <- "data_intermediate/gse17542_clean.rds"

if (is.na(infile)) {
  stop(
    "Missing required GSE17542 input. Expected one of: ",
    paste(infile_candidates, collapse = ", ")
  )
}

message("Using GSE17542 input: ", infile)
eset <- readRDS(infile)

# 1) Expression + metadata
expr_probe <- exprs(eset)                 # probes x samples
meta <- pData(eset)                       # samples x fields
fdat <- tryCatch(fData(eset), error = function(e) NULL)

# 2) Probe -> Gene Symbol mapping
#    Try fData first; else use GPL table
probe_ids <- rownames(expr_probe)

get_symbol_from_fdata <- function(fdat) {
  if (is.null(fdat) || nrow(fdat) == 0) return(NULL)
  cn <- colnames(fdat)

  # Prefer common symbol columns
  preferred <- c("Gene Symbol", "GENE_SYMBOL", "Symbol", "SYMBOL", "gene_symbol",
                 "Gene symbol", "GENESYMBOL", "GeneSymbol", "gene assignment")

  pick <- preferred[preferred %in% cn]
  if (length(pick) > 0) {
    sym <- as.character(fdat[[pick[1]]])
    return(sym)
  }

  # Fallback: any column containing "symbol"
  idx <- grep("symbol", cn, ignore.case = TRUE)
  if (length(idx) > 0) {
    sym <- as.character(fdat[[cn[idx[1]]]])
    return(sym)
  }

  return(NULL)
}

symbols <- get_symbol_from_fdata(fdat)

# If fData didn't have symbols, pull GPL annotation table and map
if (is.null(symbols) || all(is.na(symbols)) || sum(nzchar(symbols), na.rm = TRUE) == 0) {
  gpl_id <- annotation(eset)
  if (is.na(gpl_id) || gpl_id == "") {
    # sometimes stored in metadata
    gpl_id <- unique(meta$platform_id)
    gpl_id <- gpl_id[!is.na(gpl_id)][1]
  }
  stopifnot(!is.na(gpl_id), gpl_id != "")

  message("No gene symbols in fData; fetching GPL annotation: ", gpl_id)
  gpl <- getGEO(gpl_id)
  gpl_tab <- Table(gpl)
  gpl_dt <- as.data.table(gpl_tab)

  # Identify probe ID column (usually "ID")
  id_col <- if ("ID" %in% names(gpl_dt)) "ID" else names(gpl_dt)[1]
  setnames(gpl_dt, id_col, "PROBE_ID")

  # Identify a symbol-like column
  sym_col_candidates <- names(gpl_dt)[grep("symbol|gene.?symbol|genesymbol", names(gpl_dt), ignore.case = TRUE)]
  if (length(sym_col_candidates) == 0) {
    stop("Could not find a gene symbol column in GPL table. Inspect Table(getGEO('", gpl_id, "')).")
  }
  sym_col <- sym_col_candidates[1]
  setnames(gpl_dt, sym_col, "GENE_SYMBOL_RAW")

  map_dt <- gpl_dt[, .(PROBE_ID, GENE_SYMBOL_RAW)]
  map_dt <- unique(map_dt)

  # Clean symbol: take first token if multiple (e.g., "TH|..." or "TH /// ...")
  map_dt[, GENE_SYMBOL := toupper(trimws(sub(" ///.*$", "", sub("\\|.*$", "", as.character(GENE_SYMBOL_RAW)))))]

  # Match to expr rows
  m <- match(probe_ids, map_dt$PROBE_ID)
  symbols <- map_dt$GENE_SYMBOL[m]
}

# Standardize symbols
symbols <- toupper(trimws(as.character(symbols)))
symbols[symbols %in% c("", "NA", "N/A", "---")] <- NA

# 3) Collapse probes -> genes (median)
keep <- !is.na(symbols)
expr_probe <- expr_probe[keep, , drop = FALSE]
symbols_kept <- symbols[keep]

dt <- as.data.table(expr_probe, keep.rownames = "PROBE_ID")
dt[, GENE_SYMBOL := symbols_kept]

# median per gene per sample
gene_cols <- setdiff(names(dt), c("PROBE_ID", "GENE_SYMBOL"))
expr_gene_dt <- dt[, lapply(.SD, median, na.rm = TRUE), by = GENE_SYMBOL, .SDcols = gene_cols]

expr_gene <- as.matrix(expr_gene_dt[, ..gene_cols])
rownames(expr_gene) <- expr_gene_dt$GENE_SYMBOL
colnames(expr_gene) <- gene_cols

# 4) Filter to Substantia Nigra samples (robust match)
meta_text <- apply(meta, 1, function(r) paste(r, collapse = " | "))
is_sn <- grepl("substantia\\s*nigra|\\bSN\\b", meta_text, ignore.case = TRUE)

if (sum(is_sn) < 2) {
  message("SN filter found <2 samples; leaving all samples in place (inspect meta to refine).")
  is_sn <- rep(TRUE, nrow(meta))
}

meta_sn <- meta[is_sn, , drop = FALSE]
expr_gene_sn <- expr_gene[, rownames(meta_sn), drop = FALSE]

# 5) Label groups: control / mptp_2d / mptp_10d (robust parse)
txt <- apply(meta_sn, 1, function(r) tolower(paste(r, collapse = " | ")))

group <- rep(NA_character_, length(txt))

# Control first
group[grepl("control|saline|vehicle|untreated", txt)] <- "control"

# MPTP timepoints
group[grepl("mptp", txt) & grepl("2\\s*d|2\\s*day|day\\s*2|\\b2d\\b", txt)]  <- "mptp_2d"
group[grepl("mptp", txt) & grepl("10\\s*d|10\\s*day|day\\s*10|\\b10d\\b", txt)] <- "mptp_10d"

# Fallback: if MPTP present but no time parsed
group[is.na(group) & grepl("mptp", txt)] <- "mptp_unknown"

meta_sn$group <- factor(group, levels = c("control", "mptp_2d", "mptp_10d", "mptp_unknown"))

# Basic sanity
message("Group counts:")
print(table(meta_sn$group, useNA = "ifany"))

# 6) Save
dir.create("data_intermediate", showWarnings = FALSE)

saveRDS(
  list(
    expr_gene = expr_gene_sn,  # genes x samples
    meta = meta_sn,            # samples x fields incl group
    notes = list(
      probe_to_gene = "fData preferred, else GPL Table mapping; probes collapsed by median"
    )
  ),
  outfile
)

message("Saved: ", outfile)
