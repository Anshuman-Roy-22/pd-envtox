# Human PD validation checkpoint

## Frozen design

- Cohorts: GSE20141 (8 control, 10 PD laser-captured SNpc samples) and GSE7621
  (9 control, 16 PD bulk substantia nigra samples).
- Platform: GPL570 for both cohorts.
- Panel: 21 genes frozen by the primary network-v2 analysis before human outcomes.
- Preprocessing: complete-case series-matrix probes; log2(x + 1) when the 99th
  percentile exceeds 100.
- Mapping: August 2016 GEO GPL570 annotation; probes mapping to more than one
  symbol excluded; one probe per gene selected by highest mean expression across
  all samples, with probe ID as a deterministic tie-break.
- Cohort model in this recovery run: equal-variance two-group linear model.
- Meta-analysis: equal-weight directional Stouffer z.
- Multiple testing: Benjamini-Hochberg across all evaluable genes.
- Matched null: 10,000 fixed-seed draws, matched jointly on expression quintile
  in each cohort; panel genes excluded from the background.

## Result

- 20/21 genes were evaluable in both cohorts. HSPA8 was excluded because every
  GPL570 probe in the frozen GEO annotation maps jointly to HSPA8 and snoRNAs.
- 14/20 genes had matching PD-control directions.
- No panel gene passed FDR <= 0.05 in either cohort or the directional meta-analysis.
- Direction agreement was not greater than the matched background (p = 0.1974).
- Mean absolute meta-z was not greater than matched genes (p = 0.4579).
- The panel's mean meta-z was negative and borderline against the two-sided
  matched null (mean = -0.7508, p = 0.0522). This is exploratory, not confirmed.

## Important recovery limitation

The earlier frozen implementation used limma moderated statistics. The transient
workspace containing it was lost before results were delivered, and this runtime
does not contain R. This checkpoint reconstructs all outcome-blind preprocessing,
mapping, matching, and meta-analysis decisions, but uses ordinary linear-model
t statistics. Final competition claims must be confirmed with the frozen limma
implementation in an R environment before release.
