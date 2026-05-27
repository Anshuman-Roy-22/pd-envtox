# scripts/11_fig4_celltype_heatmap_and_dotplot.R
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(Seurat)
})

dir.create("results", showWarnings = FALSE)

# Inputs
celltype_path <- "results/original20_celltype_specific_logFC.csv"
mptp10_path   <- "results/gse17542_limma_mptp10d_vs_control.csv"
seurat_path   <- "data_intermediate/gse187012_seurat_object.rds"

stopifnot(file.exists(celltype_path), file.exists(mptp10_path), file.exists(seurat_path))

ct <- fread(celltype_path)
ct[, gene := toupper(gene)]

mptp <- fread(mptp10_path)
mptp[, gene := toupper(gene)]
mptp <- mptp[, .(gene, mptp10d_logFC = logFC, mptp10d_fdr = adj.P.Val)]

dt <- merge(ct, mptp, by = "gene", all.x = TRUE)

# Order genes by |MG_log2FC| + |DA_log2FC| (cell-type salience) - stable ordering for visual
num <- function(x) suppressWarnings(as.numeric(x))
dt[, `:=`(
  DA_log2FC = num(DA_log2FC),
  MG_log2FC = num(MG_log2FC),
  mptp10d_logFC = num(mptp10d_logFC)
)]
dt[, order_score := abs(MG_log2FC) + abs(DA_log2FC)]
setorder(dt, -order_score)
dt[, gene_f := factor(gene, levels = rev(dt$gene))]

# ---- FIG 4A: Heatmap of cell-type logFC (DA, Microglia) + MPTP bulk logFC ----
hm <- melt(
  dt[, .(gene, gene_f, DA_log2FC, MG_log2FC, mptp10d_logFC)],
  id.vars = c("gene","gene_f"),
  variable.name = "source",
  value.name = "logFC"
)
hm[, source := factor(
  source,
  levels = c("DA_log2FC","MG_log2FC","mptp10d_logFC"),
  labels = c("Mb/Pq: DA (EXP-CTL)", "Mb/Pq: Microglia (EXP-CTL)", "MPTP 10d bulk (T-C)")
)]

# symmetric limits for color interpretability
lim <- max(abs(hm$logFC), na.rm = TRUE)
if (!is.finite(lim) || lim == 0) lim <- 1

p4a <- ggplot(hm, aes(x = source, y = gene_f, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  coord_cartesian(expand = FALSE) +
  labs(
    title = "Cell-type-specific perturbation of the original 20 candidates",
    x = NULL, y = NULL, fill = "logFC"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5,face="bold", size=16, margin=margin(b=6)),
    plot.subtitle = element_text(size=11, margin=margin(b=10)),
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(color="black"),
    panel.grid = element_blank(),
    plot.margin = margin(12, 16, 12, 16)
  ) +
  scale_fill_gradient2(
    low = "#BFDFFF",
    mid = "#F7FBFF",
    high = "#1E3A8A",
    midpoint = 0,
    limits = c(-lim, lim)
  )

ggsave("results/fig4A_original20_celltype_logFC_heatmap.png", p4a, width = 10, height = 6.5, dpi = 300)

# ---- FIG 4B: DotPlot across Mb/Pq cell types (DA/Microglia/Other) ----
obj <- readRDS(seurat_path)
DefaultAssay(obj) <- "RNA"

# ensure celltype exists
stopifnot("celltype" %in% colnames(obj@meta.data))
Idents(obj) <- factor(obj@meta.data$celltype, levels = c("DA","Microglia","Other"))

genes20 <- intersect(dt$gene, rownames(obj))
stopifnot(length(genes20) >= 10)

dp <- DotPlot(obj, features = genes20, assay = "RNA") +
  labs(
    title = "Expression context of the original 20 across cell types (Mb/Pq)",
    subtitle = NULL,
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 6)),
    plot.subtitle = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 16, 12, 16)
  )

ggsave("results/fig4B_original20_dotplot_celltypes.png", dp, width = 12, height = 5.5, dpi = 300)

cat("Saved:\n",
    " - results/fig4A_original20_celltype_logFC_heatmap.png\n",
    " - results/fig4B_original20_dotplot_celltypes.png\n", sep="")


