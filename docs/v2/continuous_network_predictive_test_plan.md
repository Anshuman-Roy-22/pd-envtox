# Frozen continuous network-predictive test

This is a transparently labeled secondary analysis conceived after fixed-panel
validation was null. It is frozen before network proximity is merged with any
transcriptomic outcome table.

## Question

Across all non-seed network-ranked genes eligible under the network-v2 GTEx
substantia nigra filter, does closer seed proximity predict a larger absolute
transcriptomic response?

## Prespecified outcomes

1. GSE46798: pesticide effect in genetically corrected human DA neurons.
2. GSE17542: 10-day MPTP effect in SN DA neurons.
3. GSE17542: 10-day SN-minus-VTA treatment interaction.
4. Human postmortem: absolute directional Stouffer meta-z across GSE20141 and
   GSE7621.

## Model

- Exclude the 14 seed genes.
- Merge by uppercase gene symbol.
- Exposure datasets use absolute ordinary-model t; postmortem uses absolute
  meta-z.
- Primary descriptive test: Spearman correlation between closeness
  (`-proximity`) and absolute response.
- Primary adjusted test: partial Spearman correlation after separately
  residualizing ranked closeness and ranked absolute response against ranked
  log1p STRING degree and ranked log1p GTEx SN TPM.
- Two-sided p-values; Benjamini-Hochberg correction across the four adjusted
  tests.
- Directional response is not tested because the network model predicts
  prioritization/proximity, not upregulation or downregulation.

This analysis evaluates predictive validity of the continuous model. It cannot
retroactively validate the fixed 21-gene panel.
