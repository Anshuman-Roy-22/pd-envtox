# Frozen GSE17542 SN-versus-VTA analysis plan

Frozen before candidate or module outcomes were calculated.

## Design and evidence hierarchy

GSE17542 contains three biological replicates in each of six groups: substantia
nigra (SN) and ventral tegmental area (VTA) dopamine neurons under control,
2-day MPTP, and 10-day MPTP conditions. There is one common control group per
region, not duration-specific controls.

Primary tests:

1. 10-day MPTP effect in SN.
2. 10-day region-by-treatment interaction: SN response minus VTA response.

Secondary tests:

3. 10-day MPTP effect in VTA.
4. Corresponding 2-day SN, VTA, and interaction contrasts.

Context test: baseline SN minus VTA.

## Fixed panel and modules

The same outcome-independent 21-gene panel and four modules frozen for GSE46798
are used without modification.

## Processing and inference

- Use the GEO series matrix, documented as log2 transformed and RMA normalized.
- Use the August 2016 GPL1261 mouse annotation.
- Uppercase symbols to map one-to-one human/mouse symbols; exclude probes mapped
  to multiple symbols.
- For genes with multiple eligible probes, collapse sample-wise by the median.
- Fit a six-group cell-means linear model and the seven declared contrasts.
- Apply genome-wide Benjamini-Hochberg FDR separately within each contrast.
- Test panel/modules with 10,000 fixed-seed random gene sets matched on eligible
  probe-count category and residual-variance decile. Residual variance is a
  precision nuisance estimate and neither effect direction nor significance.
- Use mean absolute t as the primary responsiveness metric and mean signed t as
  a secondary coordination metric.
- Confirm final claims with limma empirical-Bayes moderation in R.
