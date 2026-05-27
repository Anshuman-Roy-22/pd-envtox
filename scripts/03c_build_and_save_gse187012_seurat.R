# scripts/03c_build_and_save_gse187012_seurat.R
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

in_dir <- "data_raw/GSE187012/extracted"
dir.create("data_intermediate", showWarnings = FALSE)
source("scripts/helpers_10x_io.R")

ctl_feat <- file.path(in_dir, "GSM5667021_CTL_features.tsv.gz")
ctl_bc   <- file.path(in_dir, "GSM5667021_CTL_barcodes.tsv.gz")
ctl_mtx  <- file.path(in_dir, "GSM5667021_CTL_matrix.mtx.gz")

exp_feat <- file.path(in_dir, "GSM5667022_EXP_features.tsv.gz")
exp_bc   <- file.path(in_dir, "GSM5667022_EXP_barcodes.tsv.gz")
exp_mtx  <- file.path(in_dir, "GSM5667022_EXP_matrix.mtx.gz")

ctl <- make_10x_seurat_object(ctl_feat, ctl_bc, ctl_mtx, "CTL", "control")
exp <- make_10x_seurat_object(exp_feat, exp_bc, exp_mtx, "EXP", "exposed")
obj <- merge(ctl, y = exp)

# Standard workflow
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)
obj <- JoinLayers(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)

# Label DA and Microglia clusters via marker module scores
DA_MARKERS <- intersect(c("TH","SLC6A3","DDC","SLC18A2","NR4A2","FOXA2","LMX1A","PITX3"), rownames(obj))
MG_MARKERS <- intersect(c("P2RY12","TMEM119","CX3CR1","AIF1","C1QA","C1QB","C1QC","TYROBP","TREM2"), rownames(obj))

stopifnot(length(DA_MARKERS) >= 3, length(MG_MARKERS) >= 3)

DefaultAssay(obj) <- "RNA"
obj <- AddModuleScore(obj, features = list(DA_MARKERS), name = "DA_SCORE", assay = "RNA", layer = "data", verbose = FALSE)
obj <- AddModuleScore(obj, features = list(MG_MARKERS), name = "MG_SCORE", assay = "RNA", layer = "data", verbose = FALSE)

meta <- obj@meta.data
meta$cluster <- as.character(Idents(obj))

avg <- aggregate(meta[, c("DA_SCORE1","MG_SCORE1")], by = list(cluster = meta$cluster), FUN = mean)
da_cluster <- as.character(avg$cluster[which.max(avg$DA_SCORE1)])
mg_cluster <- as.character(avg$cluster[which.max(avg$MG_SCORE1)])

meta$celltype <- "Other"
meta$celltype[meta$cluster == da_cluster] <- "DA"
meta$celltype[meta$cluster == mg_cluster] <- "Microglia"
obj@meta.data <- meta

saveRDS(obj, "data_intermediate/gse187012_seurat_object.rds")
cat("Saved: data_intermediate/gse187012_seurat_object.rds\n")
cat("DA cluster:", da_cluster, "Microglia cluster:", mg_cluster, "\n")
print(table(obj@meta.data$celltype))
