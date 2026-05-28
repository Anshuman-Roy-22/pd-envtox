# scripts/12_compare_candidate20_vs_top20.R
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create("results", showWarnings = FALSE)

dt <- fread("results/robust_sporadic_candidates_ranked.csv")
dt[, gene := toupper(gene)]

candidate20 <- c(
  "SNX3","VPS26B","PSENEN","APH1A","SNX2","APH1B","EIF4E","EIF4B","EIF4A2","NCSTN",
  "PABPC1","DCTN2","EIF4A1","VPS26A","VPS29","APP","CLIP1","DYNC1H1","ACTR1B","DCTN3"
)
candidate20 <- intersect(candidate20, dt$gene)

# Define TOP20 by the final integrated score
top20 <- dt[order(-reactivity)][1:20]$gene

dt[, set := fifelse(gene %in% candidate20, "Candidate20",
             fifelse(gene %in% top20, "Top20_by_Reactivity", "Background"))]

# Metrics table
metrics <- dt[set != "Background", .(
  n = .N,
  mean_reactivity = mean(reactivity, na.rm=TRUE),
  median_reactivity = median(reactivity, na.rm=TRUE),
  mean_mptp10d_norm = mean(mptp10d_norm, na.rm=TRUE),
  mean_pq_norm = mean(pq_norm, na.rm=TRUE),
  frac_same_direction_10d = mean(as.logical(same_direction_10d), na.rm=TRUE),
  frac_same_direction_2d  = mean(as.logical(same_direction_2d), na.rm=TRUE)
), by = set]

fwrite(metrics, "results/compare_candidate20_vs_top20_metrics.csv")

# Visual: distributions
sub <- dt[set != "Background", .(gene, set, reactivity, mptp10d_norm, pq_norm)]
sub <- melt(sub, id.vars = c("gene","set"), variable.name = "metric", value.name = "value")
sub[, metric := factor(metric, levels=c("reactivity","mptp10d_norm","pq_norm"),
                       labels=c("Integrated Reactivity","MPTP 10d perturb (norm)","Mb/Pq perturb (norm)"))]

p <- ggplot(sub, aes(x = set, y = value)) +
  geom_boxplot(color="blue", fill="lightblue", outlier.alpha = 0.2) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7) +
  labs(
    title = "Candidate 20 vs Top 20 Perturbed (Integrated Score)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold", size=15, margin=margin(b=6)),
    plot.subtitle = element_text(size=10, margin=margin(b=10)),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face="bold"),
    plot.margin = margin(12, 16, 12, 16)
  ) +
  facet_wrap(~metric, scales = "free_y")

ggsave("results/fig6_candidate20_vs_top20_boxplots.png", p, width = 10, height = 5.8, dpi = 300)

cat("Saved:\n - results/compare_candidate20_vs_top20_metrics.csv\n - results/fig6_candidate20_vs_top20_boxplots.png\n")
