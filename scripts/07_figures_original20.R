# scripts/07_figures_original20.R
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

dir.create("results", showWarnings = FALSE)

tab_path <- "results/original20_table1.csv"
stopifnot(file.exists(tab_path))
t20 <- fread(tab_path)
t20[, gene := toupper(gene)]

# Safe numeric conversion
num <- function(x) suppressWarnings(as.numeric(x))
t20[, `:=`(
  reactivity = num(reactivity),
  pq_norm = num(pq_norm),
  mptp2d_norm = num(mptp2d_norm),
  mptp10d_norm = num(mptp10d_norm),
  pq_log2FC = num(pq_log2FC),
  mptp10d_logFC = num(mptp10d_logFC),
  mptp2d_logFC = num(mptp2d_logFC),
  mptp10d_fdr = num(mptp10d_fdr),
  mptp2d_fdr = num(mptp2d_fdr)
)]

# Order genes by reactivity (high -> low) for consistent visuals
t20 <- t20[order(-reactivity)]
t20[, gene_f := factor(gene, levels = rev(gene))]  # for coord_flip aesthetics

# Global theme (poster-like)
poster_theme <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 11, margin = margin(b = 10)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(12, 16, 12, 16)
    )
}

# ----------------------------
# FIG 1 - Reactivity barplot (clean)
# ----------------------------
p1 <- ggplot(t20, aes(x = gene_f, y = reactivity)) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    title = "Original 20 Candidate Panel: Integrated Reactivity",
    x = NULL,
    y = "Integrated Reactivity (0-1)"
  ) +
  poster_theme() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 6))
  )

ggsave("results/fig1_original20_reactivity_bar.png", p1, width = 9, height = 6, dpi = 300)

# ----------------------------
# FIG 2 - Scatter: MPTP10d_norm vs PQ_norm with non-overlapping labels
# ----------------------------
p2 <- ggplot(t20, aes(x = mptp10d_norm, y = pq_norm, label = gene)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  labs(
    title = "Cross-model Perturbation: MPTP 10d vs Mb/Pq",
    x = "MPTP 10d perturbation (0-1)",
    y = "Mb/Pq perturbation (0-1)"
  ) +
  poster_theme() +
  theme(panel.grid.major = element_line(linewidth = 0.3),
  plot.title = element_text(hjust = 0.5,face="bold", size=15, margin=margin(b=6)))

ggsave("results/fig2_original20_scatter_mptp10d_vs_pq.png", p2, width = 8, height = 6, dpi = 300)

# ----------------------------
# FIG 3 - Direction concordance tile (clean legend + spacing)
# ----------------------------
t20[, same_direction_10d := as.logical(same_direction_10d)]
t20[, same_direction_2d := as.logical(same_direction_2d)]

dir_dt <- melt(
  t20[, .(gene, same_direction_10d, same_direction_2d)],
  id.vars = "gene",
  variable.name = "comparison",
  value.name = "same_direction"
)
dir_dt[, gene_f := factor(gene, levels = levels(t20$gene_f))]
dir_dt[, comparison := fifelse(comparison == "same_direction_10d", "MPTP 10d vs Mb/Pq", "MPTP 2d vs Mb/Pq")]

p3 <- ggplot(dir_dt, aes(x = comparison, y = gene_f, fill = same_direction)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_x_discrete(expand = expansion(add = 0.2)) +
  scale_fill_manual(values = c("FALSE" = "#D6EAF8", "TRUE" = "#1E3A8A")) +
  labs(
    title = "Directional Concordance Across Models",
    x = NULL, y = NULL, fill = "Same Direction"
  ) +
  poster_theme() +
  theme(
    axis.text.x = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5,face="bold", size=15, margin=margin(b=6))
  )

ggsave("results/fig3_original20_direction_tile.png", p3, width = 8.2, height = 6.2, dpi = 300)

# ----------------------------
# FIG 4 - Cell-type attribution: gene-level strip (interpretable)
# ----------------------------
if ("inferred_celltype" %in% names(t20)) {
  t20[, inferred_celltype := fifelse(is.na(inferred_celltype) | inferred_celltype == "", "Unknown", inferred_celltype)]

  p4 <- ggplot(t20, aes(x = 1, y = gene_f, fill = inferred_celltype)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_x_continuous(breaks = NULL) +
    labs(
      title = "Cell-type Attribution (Mb/Pq single-cell reference)",
      x = NULL, y = NULL, fill = "Cell type"
    ) +
    poster_theme() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5,face="bold", size=15, margin=margin(b=6))
    )

  ggsave("results/fig4_original20_celltype_strip.png", p4, width = 7.8, height = 6.2, dpi = 300)
}

# ----------------------------
# FIG 5 - Heatmap: PQ_norm, MPTP2d_norm, MPTP10d_norm (proper layout)
# ----------------------------
hm <- t20[, .(gene, gene_f, pq_norm, mptp2d_norm, mptp10d_norm)]
hm_long <- melt(
  hm,
  id.vars = c("gene", "gene_f"),
  variable.name = "dataset",
  value.name = "norm_perturb"
)
hm_long[, dataset := factor(
  dataset,
  levels = c("pq_norm", "mptp2d_norm", "mptp10d_norm"),
  labels = c("Mb/Pq", "MPTP 2d", "MPTP 10d")
)]

p5 <- ggplot(hm_long, aes(x = dataset, y = gene_f, fill = norm_perturb)) +
  geom_tile(color = "white", linewidth = 0.5) +
  labs(
    title = "Perturbation Pattern Across Conditions (Normalized)",
    x = NULL, y = NULL, fill = "Perturb\n(0-1)"
  ) +
  poster_theme() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 6))
  )

ggsave("results/fig5_original20_norm_heatmap.png", p5, width = 7.8, height = 6.2, dpi = 300)

cat("Saved figures to results/:\n",
    " - fig1_original20_reactivity_bar.png\n",
    " - fig2_original20_scatter_mptp10d_vs_pq.png\n",
    " - fig3_original20_direction_tile.png\n",
    " - fig4_original20_celltype_strip.png (if inferred_celltype exists)\n",
    " - fig5_original20_norm_heatmap.png\n", sep="")


