# scripts/13_permutation_matched_original20.R
suppressPackageStartupMessages({
  library(data.table)
})

dt <- fread("results/robust_sporadic_candidates_ranked.csv")
dt[, gene := toupper(gene)]
required_cols <- c("gene", "pq_norm", "mptp10d_norm", "reactivity")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

candidate20 <- c(
  "SNX3","VPS26B","PSENEN","APH1A","SNX2","APH1B","EIF4E","EIF4B","EIF4A2","NCSTN",
  "PABPC1","DCTN2","EIF4A1","VPS26A","VPS29","APP","CLIP1","DYNC1H1","ACTR1B","DCTN3"
)
candidate20 <- intersect(candidate20, dt$gene)

obs <- dt[gene %in% candidate20]
obs_mean <- mean(obs$reactivity, na.rm=TRUE)

# Bin genes by pq_norm (proxy for detectability / scale)
dt <- dt[!is.na(pq_norm) & !is.na(mptp10d_norm) & !is.na(reactivity)]
q_breaks <- unique(quantile(dt[["pq_norm"]], probs = seq(0, 1, 0.1), na.rm = TRUE))
if (length(q_breaks) < 2) {
  dt[, bin := factor("all")]
} else {
  dt[, bin := cut(pq_norm, breaks = q_breaks, include.lowest = TRUE)]
}

# For each candidate gene, sample a random gene from the same bin (without replacement per draw)
set.seed(42)
B <- 2000
null <- numeric(B)

cand_bins <- dt[gene %in% candidate20, .(gene, bin)]
stopifnot(nrow(cand_bins) == length(candidate20))

for (i in seq_len(B)) {
  chosen <- character(0)
  for (b in cand_bins$bin) {
    pool <- dt[bin == b & !(gene %in% chosen), gene]
    if (length(pool) == 0) pool <- dt[bin == b, gene]
    chosen <- c(chosen, sample(pool, 1))
  }
  null[i] <- mean(dt[gene %in% chosen]$reactivity, na.rm=TRUE)
}

p_emp <- mean(null >= obs_mean)
out <- data.table(
  observed_mean_reactivity = obs_mean,
  null_mean = mean(null),
  null_sd = sd(null),
  empirical_p_value = p_emp,
  B = B
)
fwrite(out, "results/permutation_matched_original20.csv")

cat("Matched permutation test saved: results/permutation_matched_original20.csv\n")
cat("Empirical p-value:", p_emp, "\n")
