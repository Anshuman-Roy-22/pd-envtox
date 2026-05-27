# scripts/14_plot_gsea_results.R
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create("results", showWarnings = FALSE)

in_rds <- "results/gsea_mptp10d_msigdb.rds"
stopifnot(file.exists(in_rds))

g <- readRDS(in_rds)

plot_gsea <- function(gsea_obj, title, out_png, out_csv, n_terms = 15) {
  if (is.null(gsea_obj)) {
    fwrite(data.table(), out_csv)
    p_empty <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "no results object", size = 6, color = "#1E3A8A") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = title, x = NULL, y = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 8)),
        plot.title.position = "plot",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(12, 16, 12, 16)
      )
    ggsave(out_png, p_empty, width = 10, height = 6.2, dpi = 300)
    cat("Saved placeholder plot:", out_png, "\n")
    return(invisible(FALSE))
  }

  res <- as.data.frame(gsea_obj@result)
  if (nrow(res) == 0) {
    fwrite(data.table(), out_csv)
    p_empty <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "no enriched terms found", size = 6, color = "#1E3A8A") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = title, x = NULL, y = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 8)),
        plot.title.position = "plot",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(12, 16, 12, 16)
      )
    ggsave(out_png, p_empty, width = 10, height = 6.2, dpi = 300)
    cat("Saved placeholder plot:", out_png, "\n")
    return(invisible(FALSE))
  }

  res <- as.data.table(res)
  setorder(res, p.adjust)
  top <- res[1:min(n_terms, .N)]
  fwrite(top, out_csv)

  # Use NES for direction/strength; label with cleaned pathway names
  top[, term := gsub("^REACTOME_", "", Description)]
  top[, term := gsub("^GOBP_", "", term)]
  top[, term := gsub("_", " ", term)]
  top[, term := tolower(term)]
  top[, term := gsub("\\s+", " ", trimws(term))]
  top[, term := ifelse(
    nchar(term) > 0,
    paste0(toupper(substr(term, 1, 1)), substr(term, 2, nchar(term))),
    term
  )]
  top[, term := vapply(term, function(x) paste(strwrap(x, width = 42), collapse = "\n"), character(1))]
  top[, term := factor(term, levels = rev(term))]

  p <- ggplot(top, aes(x = term, y = NES)) +
    geom_hline(yintercept = 0, linewidth = 0.45, color = "gray75") +
    geom_col(width = 0.72, fill = "#D6EAF8", color = "#1E3A8A", linewidth = 0.55) +
    coord_flip() +
    labs(
      title = title,
      x = NULL,
      y = "normalized enrichment score (nes)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 8)),
      plot.title.position = "plot",
      plot.subtitle = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray88", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 16, 12, 16)
    )

  ggsave(out_png, p, width = 10, height = 6.2, dpi = 300)
  cat("Saved:", out_png, "\n")
  invisible(TRUE)
}

plot_gsea(g$reactome,
          "Fig 7A: GSEA (MPTP 10d ranked perturbation) - Reactome",
          "results/fig7A_gsea_mptp10d_reactome.png",
          "results/fig7A_gsea_mptp10d_reactome_table.csv")

plot_gsea(g$gobp,
          "Fig 7B: GSEA (MPTP 10d ranked perturbation) - GO: Biological Process",
          "results/fig7B_gsea_mptp10d_gobp.png",
          "results/fig7B_gsea_mptp10d_gobp_table.csv")


