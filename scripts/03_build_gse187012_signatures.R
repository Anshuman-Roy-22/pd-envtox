# scripts/03_build_gse187012_signatures.R
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
})

in_dir <- "data_raw/GSE187012/extracted"
dir.create("data_intermediate", showWarnings = FALSE)
source("scripts/helpers_10x_io.R")

# ---- file paths (you listed these exact names) ----
ctl_feat <- file.path(in_dir, "GSM5667021_CTL_features.tsv.gz")
ctl_bc   <- file.path(in_dir, "GSM5667021_CTL_barcodes.tsv.gz")
ctl_mtx  <- file.path(in_dir, "GSM5667021_CTL_matrix.mtx.gz")

exp_feat <- file.path(in_dir, "GSM5667022_EXP_features.tsv.gz")
exp_bc   <- file.path(in_dir, "GSM5667022_EXP_barcodes.tsv.gz")
exp_mtx  <- file.path(in_dir, "GSM5667022_EXP_matrix.mtx.gz")

ctl <- make_10x_seurat_object(ctl_feat, ctl_bc, ctl_mtx, "CTL", "control")
exp <- make_10x_seurat_object(exp_feat, exp_bc, exp_mtx, "EXP", "exposed")

obj <- merge(ctl, y = exp)

cat("Merged cells:", ncol(obj), " genes:", nrow(obj), "\n")

# ---- standard clustering ----
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)

# ---- "smart" labeling without web searches: module scores ----
DA_MARKERS <- c("TH","SLC6A3","DDC","SLC18A2","NR4A2","FOXA2","LMX1A","PITX3")
MG_MARKERS <- c("P2RY12","TMEM119","CX3CR1","AIF1","C1QA","C1QB","C1QC","TYROBP","TREM2")

# keep markers present
DA_MARKERS <- intersect(DA_MARKERS, rownames(obj))
MG_MARKERS <- intersect(MG_MARKERS, rownames(obj))

if (length(DA_MARKERS) < 3) stop("Too few DA markers found in matrix (", length(DA_MARKERS), ").")
if (length(MG_MARKERS) < 3) stop("Too few Microglia markers found in matrix (", length(MG_MARKERS), ").")

obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"
obj <- AddModuleScore(obj, features = list(DA_MARKERS), name = "DA_SCORE", verbose = FALSE)
obj <- AddModuleScore(obj, features = list(MG_MARKERS), name = "MG_SCORE", verbose = FALSE)

# module score columns end with "1"
da_col <- "DA_SCORE1"
mg_col <- "MG_SCORE1"

# cluster-level average scores
cl <- as.character(Idents(obj))
meta <- obj@meta.data
meta$cluster <- cl

avg <- as.data.table(meta)[, .(
  DA = mean(get(da_col), na.rm = TRUE),
  MG = mean(get(mg_col), na.rm = TRUE),
  n = .N
), by = cluster]

# label rules: pick clusters with highest DA and highest MG
da_cluster <- avg[which.max(DA), cluster]
mg_cluster <- avg[which.max(MG), cluster]

cat("DA cluster:", da_cluster, " Microglia cluster:", mg_cluster, "\n")
print(avg[order(-DA)][1:min(8,.N)])
print(avg[order(-MG)][1:min(8,.N)])

# assign cell type label
meta$celltype <- "Other"
meta$celltype[meta$cluster == da_cluster] <- "DA"
meta$celltype[meta$cluster == mg_cluster] <- "Microglia"

# ---- export barcode->celltype (for your records + report) ----
annot <- data.table(
  barcode = rownames(meta),
  sample = meta$sample,
  condition = meta$condition,
  cluster = meta$cluster,
  celltype = meta$celltype
)
fwrite(annot, "data_intermediate/gse187012_cell_annotations.csv")

# ---- build signature matrix (average log-normalized expression) ----
expr <- SeuratObject::LayerData(obj[["RNA"]], layer = "data") # log-normalized
da_cells <- rownames(meta)[meta$celltype == "DA"]
mg_cells <- rownames(meta)[meta$celltype == "Microglia"]

stopifnot(length(da_cells) > 20, length(mg_cells) > 20)

sig <- cbind(
  DA = Matrix::rowMeans(expr[, da_cells, drop = FALSE]),
  Microglia = Matrix::rowMeans(expr[, mg_cells, drop = FALSE])
)
sig <- as.matrix(sig)

saveRDS(list(
  signature_matrix = sig,
  chosen_clusters = list(DA = da_cluster, Microglia = mg_cluster),
  cluster_score_table = avg
), "data_intermediate/gse187012_signatures.rds")

# ---- marker lists (for no-tool fallback & reporting) ----
Idents(obj) <- factor(meta$celltype)
da_mark <- FindMarkers(obj, ident.1 = "DA", ident.2 = "Other", only.pos = TRUE, logfc.threshold = 0.25)
mg_mark <- FindMarkers(obj, ident.1 = "Microglia", ident.2 = "Other", only.pos = TRUE, logfc.threshold = 0.25)

markers <- list(
  DA_top50 = rownames(head(da_mark[order(-da_mark$avg_log2FC), , drop=FALSE], 50)),
  Microglia_top50 = rownames(head(mg_mark[order(-mg_mark$avg_log2FC), , drop=FALSE], 50))
)
saveRDS(markers, "data_intermediate/gse187012_marker_lists.rds")

cat("Saved:\n",
    " - data_intermediate/gse187012_cell_annotations.csv\n",
    " - data_intermediate/gse187012_signatures.rds\n",
    " - data_intermediate/gse187012_marker_lists.rds\n")

