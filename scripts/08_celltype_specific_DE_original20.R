# scripts/08_celltype_specific_DE_original20.R
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

# Load Seurat object (from step 03)
obj <- readRDS("data_intermediate/gse187012_seurat_object.rds")
meta <- obj@meta.data

# Original 20 genes
genes20 <- c(
  "SNX3","VPS26B","PSENEN","APH1A","SNX2","APH1B","EIF4E","EIF4B","EIF4A2","NCSTN",
  "PABPC1","DCTN2","EIF4A1","VPS26A","VPS29","APP","CLIP1","DYNC1H1","ACTR1B","DCTN3"
)

genes20 <- intersect(genes20, rownames(obj))
stopifnot(length(genes20) > 10)

DefaultAssay(obj) <- "RNA"

# Helper: celltype-specific DE
celltype_DE <- function(celltype) {
  cells <- rownames(obj@meta.data)[obj@meta.data$celltype == celltype]
  if (length(cells) == 0) stop("No cells found for celltype: ", celltype)
  sub <- subset(obj, cells = cells)
  if (!all(c("control", "exposed") %in% unique(sub$condition))) {
    stop("Both control and exposed are required in celltype: ", celltype)
  }
  Idents(sub) <- "condition"  # CTL vs EXP
  de <- FindMarkers(
    sub,
    ident.1 = "exposed",
    ident.2 = "control",
    features = genes20,
    logfc.threshold = 0,
    min.pct = 0
  )
  de$gene <- rownames(de)
  de[, c("gene","avg_log2FC","pct.1","pct.2")]
}

da_de <- celltype_DE("DA")
mg_de <- celltype_DE("Microglia")

setnames(da_de, old=c("avg_log2FC","pct.1","pct.2"),
         new=c("DA_log2FC","DA_pct_EXP","DA_pct_CTL"))
setnames(mg_de, old=c("avg_log2FC","pct.1","pct.2"),
         new=c("MG_log2FC","MG_pct_EXP","MG_pct_CTL"))

dt <- merge(da_de, mg_de, by="gene", all=TRUE)

fwrite(dt, "results/original20_celltype_specific_logFC.csv")
cat("Saved: results/original20_celltype_specific_logFC.csv\n")
