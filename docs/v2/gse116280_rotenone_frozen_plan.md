# Frozen GSE116280 rotenone replication plan

Frozen before any candidate or pathway outcome was calculated.

GSE116280 contains 24 human LUHMES dopaminergic-neuron profiles, three replicates
per condition, spanning differentiation days 8 and 15 and three rotenone exposure
regimens. The published series matrix is already background-corrected, quantile-
normalized, batch-corrected, filtered, and summarized to Entrez genes.

## Contrasts

Primary external replication:

- Day 15, 50 nM rotenone for 24 h versus day-15 DMSO.

Secondary sensitivity:

- Day 15, 50 nM for 12 h versus day-15 DMSO.
- Day 15, 100 nM for 24 h versus day-15 DMSO.
- Corresponding three day-8 contrasts.

The day-15 50 nM/24 h contrast is primary because it tests a sustained, lower-
dose mitochondrial toxicant response in the more differentiated culture.

## Inference

- Fit a two-factor cell-means limma model to all 24 samples.
- Use the GPL17077 annotation supplied by GEO and exclude ambiguous mappings.
- Test the fixed 21 genes without redefining the panel.
- Run the frozen official Reactome preranked analysis on the primary moderated t
  statistic using 1,000 permutations.
- Proteasome assembly is the single prespecified replication pathway. A same-
  direction nominal p < 0.05 supports external replication; dataset-wide FDR is
  reported but is not required for this one-pathway confirmatory test.
