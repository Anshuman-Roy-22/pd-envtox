# scripts/13b_plot_permutation_matched_original20.R
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create("results", showWarnings = FALSE)

# Shared styling to match the report aesthetic
fill_blue <- "#D6EAF8"
line_blue <- "#1E3A8A"
vline_blue <- "#0B3C8C"
clean_theme <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = margin(b = 6)),
    plot.subtitle = element_text(hjust = 0.5, size = 10.5, margin = margin(b = 10)),
    plot.title.position = "plot",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray88", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    plot.margin = margin(12, 16, 12, 16)
  )

# Load the permutation summary
sum_path <- "results/permutation_matched_original20.csv"
stopifnot(file.exists(sum_path))
s <- fread(sum_path)

# If the full null vector is available, plot it directly.
# Otherwise, reconstruct a summary plot from the reported moments.

# Try to load null vector if present:
null_path <- "results/permutation_matched_original20_null.csv"

if (file.exists(null_path)) {
  null <- fread(null_path)$null_mean
  obs <- s$observed_mean_reactivity[1]

  p <- ggplot(data.table(null = null), aes(x = null)) +
    geom_histogram(bins = 35, fill = fill_blue, color = line_blue, linewidth = 0.4) +
    geom_vline(xintercept = obs, linewidth = 1.1, color = vline_blue, linetype = "dashed") +
    labs(
      title = "Fig 5: Matched Permutation Test (Candidate20 Mean Reactivity)",
      subtitle = paste0("Observed mean = ", round(obs, 5),
                        " | Empirical p = ", s$empirical_p_value[1],
                        " | B = ", s$B[1]),
      x = "Null distribution: mean reactivity of matched random 20-gene panels",
      y = "Count"
    ) +
    clean_theme

  ggsave("results/fig5_permutation_matched_original20.png", p, width = 10, height = 6, dpi = 300)
  cat("Saved: results/fig5_permutation_matched_original20.png\n")

} else {
  # Fallback plot using only summary statistics
  obs <- s$observed_mean_reactivity[1]
  mu  <- s$null_mean[1]
  sd  <- s$null_sd[1]

  # Create a normal approximation curve to show where obs lies
  x <- seq(mu - 4*sd, mu + 4*sd, length.out = 400)
  y <- dnorm(x, mean = mu, sd = sd)
  dt <- data.table(x = x, y = y)

  p <- ggplot(dt, aes(x = x, y = y)) +
    geom_area(fill = fill_blue, alpha = 0.65) +
    geom_line(linewidth = 1, color = line_blue) +
    geom_vline(xintercept = obs, linewidth = 1.1, color = vline_blue, linetype = "dashed") +
    labs(
      title = "Matched Permutation Test (Normal Approximation)",
      subtitle = paste0("Observed mean = ", round(obs, 5),
                        " | Null mean = ", round(mu, 5),
                        " | Null SD = ", round(sd, 5),
                        " | Empirical p = ", s$empirical_p_value[1]),
      x = "Mean reactivity (matched random panels)",
      y = "Density"
    ) +
    clean_theme

  ggsave("results/fig5_permutation_matched_original20.png", p, width = 10, height = 6, dpi = 300)
  cat("Saved (approx): results/fig5_permutation_matched_original20.png\n")
  cat("Save the null draws in the permutation script to generate the full histogram.\n")
}
