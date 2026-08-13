# AAN V2 scientific figure package

These figures summarize the revised project without hiding failed validations.
Every plotted value is read from a tracked V2 result table by
`scripts/28_aan_v2_figures.py`. Release-value assertions stop generation if a
primary result changes unexpectedly.

## Figures

### Figure 1. Evidence design

`fig1_evidence_design` separates the original network-panel hypothesis,
retrospective reanalysis, frozen external tests, human disease-tissue testing,
and independent replication or sensitivity analyses. Solid nodes indicate
frozen tests, dashed nodes indicate retrospective or nominated evidence, and
dotted nodes indicate sensitivity evidence.

### Figure 2. Proteasome cross-model convergence

`fig2_proteasome_cross_model` shows model-specific Reactome normalized
enrichment scores for corrected human DA neurons exposed to paraquat/maneb,
mouse substantia nigra DA neurons after MPTP, and human postmortem substantia
nigra PD meta-analysis. The figure reports the signed cross-model meta-FDR and
the eight shared leading-edge genes.

This is retrospective convergence. It is not labeled independent toxicant
replication.

### Figure 3. Frozen rotenone replication

`fig3_rotenone_replication` shows all six declared GSE116280 LUHMES contrasts.
The day-15, 50 nM, 24-hour contrast is the frozen primary. Day-8 nominal results
are labeled sensitivities and do not relabel the primary outcome.

### Figure 4. Human substantia nigra validation

`fig4_human_sn_validation` shows donor-level pathway scores for the frozen
GSE178265 primary and GSE243639 replication cohorts. Its forest panel reports
age-, sex-, and postmortem-interval-adjusted PD coefficients with two-sided 95%
confidence intervals and prespecified one-sided-down p values.

The frozen overall outcome is `PRIMARY_ONLY`.

### Figure 5. Validation outcome matrix

`fig5_validation_outcome_matrix` distinguishes nomination, retrospective signal,
mixed evidence, frozen confirmation, failed confirmation, and analyses that
were not performed. The matrix prevents a positive retrospective result or
sensitivity result from being presented as an independently replicated
mechanism.

## Formats and derived tables

Each figure is provided as:

- 300-DPI PNG for submission systems and slides.
- Vector PDF for publication-quality assembly.

The `fig*.tsv` files contain the exact derived values used in each figure.
`source_manifest.tsv` records the SHA256 digest of every source table.
`figure_manifest.tsv` records output checksums, sizes, and PNG dimensions.

## Reproduction

From the repository root:

```powershell
python scripts/28_aan_v2_figures.py
```

The generator overwrites only the named files inside
`results/v2/aan_figures/`.
