# Fixed functional-rescue reanalysis specification

Date: 2026-09-15. Source state: `444ace0`. Commit this specification and the
design audit before calculating compound effects or training a phenotype score.
This is a prespecified **retrospective reanalysis of published experiments**.
The paper's major rescue claims and selected compounds are already known.
It is neither a prospectively unseen experiment nor independent PD replication.

## Question and population

Does a small multivariate control-derived phenotype score predict a repeatable
joint alpha-synuclein/neuronal-integrity response better than alpha-synuclein
intensity alone, in the separate follow-up imaging experiment?

The primary population is the 31 directly catalogue-matched compounds in
`matched_compound_design.tsv`, at the nearest follow-up dose within 5% of the
primary screen dose. Alectinib is tested at 2.5 versus 2.51 micromolar; most
others at 5 versus 5.13 micromolar. No compound is substituted based on effects.
Dasatinib and Prostratin lack an exact catalogue match and are excluded from
this prediction comparison. Prostratin remains a named positive assay control.
Follow-up selection was performed by the original authors, so this population
is enriched and cannot estimate whole-library or new-chemical performance.

The raw author DRC notebook explicitly removes Alectinib without an explanation
in that code. Our analysis retains it as deposited, uses its actual dose,
and reports an exclusion sensitivity. Existing missing wells are retained as
missing. We do not invent the 30 rows absent from an ideal four-full-plate
follow-up layout. Existing author classifier distances are not used.

## Experimental units

Plate/well keys establish 4,800 primary-screen observations and 930 follow-up
observations. Technical wells are not independent donors. The 20 screen plates
form five library blocks with four plate suffixes each; naming does not prove
that the suffixes identify biological differentiations. Each eligible compound
has four wells on distinct plates at its matched dose. Report plate consistency
and a finite-panel prediction benchmark, not donor-level inference.

Biochemical exports/workbooks do not establish biological replicate membership.
No order-based grouping of adjacent rows, pooled-cell P values, or invented
standard errors is permitted. This limits mechanistic inference but does not
prevent the specified imaging benchmark and descriptive assay summaries.

## Model and calibration

Use the four fixed features in `analysis_config.json`: individual-cell SNCA
intensity, MAP2 and TH neurite length per nucleus, and TH mean intensity.
The two separately evaluated viability features are living-nuclei count and
living-nuclei fraction. All are present in the input headers.

For each plate and feature, divide each well by the median mutant-vehicle
control. At least eight finite mutant and corrected control wells per plate
and a positive mutant median are required. Reject negative/nonfinite feature
values as invalid, preserve explicit zero values, and use a ratio floor of
1e-6 only for logarithms. A zero count remains a failed retention test.
No outlier trimming or toxicity-based removal of compounds is applied.

Train solely on mutant-vehicle and corrected-vehicle wells in the separate
144-well Fig1 control dataset. Positive Prostratin controls are not training
data. Use log2 feature ratios. Let m0 and m1 be class means and S the pooled
within-class covariance of these log ratios. Use
Sreg = 0.5*S + 0.5*diag(S) + 1e-8*I, w = inverse(Sreg)*(m1-m0).
The multivariate score is (x-m0)'w / ((m1-m0)'w), with mean training mutant
score zero and corrected score one. This is an interpretable regularized
linear discriminant; it is not a cellular simulator.

The marker-only score uses the corresponding one-dimensional SNCA projection.
The training SNCA corrected-minus-mutant mean difference must be negative
and finite; otherwise this benchmark fails its calibration gate. Report
leave-one-Fig1-plate-out control AUC as assay separation QC, with normalization
using known vehicle controls on each held-out plate. This is not clinical
diagnostic accuracy. Do not optimize shrinkage or select new features.

Apply the frozen weights to all screen and follow-up wells after the same
plate-specific vehicle normalization. Taking a median across a compound's
eligible wells gives its predictor score. At least three usable wells are
required; exclusions and their reasons must be output.

## Follow-up response definition

A well passes the operational joint-response rule if all conditions hold:

1. SNCA marker-only score is at least 0.5, representing at least half the
   control-calibrated shift in log intensity. This is not equivalent to
   50% absolute protein lowering.
2. Each of the three structural/identity feature ratios is at least 0.8
   relative to mutant vehicle on its own plate.
3. Both viability ratios are at least 0.8 relative to mutant vehicle.

A compound passes if at least 75% of its valid wells pass, with at least
three valid wells. At n=4 this means at least three; at n=3 all three.
The 0.8 boundary is an operational 20% deterioration tolerance, not a validated
biological safety margin or proof of noninferiority. Improvement of structural
features is not required by this rule, so call it marker improvement with
preserved measured features, not established functional restoration.
Threshold sensitivities of 0.7 and 0.9 must be reported without replacing 0.8.

## Primary benchmark endpoint and decision

Use **average precision (AP)** in the 31-compound follow-up panel, predicting
the above binary response from primary-screen multivariate versus marker-only
scores. Primary endpoint is AP(multivariate) minus AP(marker-only).
Calculate non-interpolated AP with ties handled as score groups, not arbitrary
row order. Also report prevalence, ROC AUC and Spearman repeatability of each
score between matched screen/follow-up doses. These are finite-panel summaries,
without population P values or biological confidence intervals.

The benchmark requires at least 20 evaluable compounds and at least five
positive and five negative follow-up outcomes. If this gate fails, label
`BENCHMARK_NOT_EVALUABLE`; descriptive statistics may still be reported.
Otherwise label `WITHIN_STUDY_PREDICTIVE_GAIN` only if the primary AP gain
is >=0.10 and >0 in at least 80% of leave-one-compound-out recalculations.
All other evaluable outcomes are `NO_DEMONSTRATED_PREDICTIVE_GAIN`.
These are preselected engineering criteria, not statistical significance.
Even a pass cannot establish a new molecular mechanism or external replication.

Mandatory sensitivities: retention 0.7/0.9; covariance shrinkage 0.1/0.9;
exclude Alectinib; collapse the two Bazedoxifene salt entries by taking median
predictor scores and requiring both follow-up labels to pass. No tuning uses
these outcomes. Show positive-control performance and compound-level results.

## Biochemical secondary descriptions

Re-extract ROS and JC-1 workbooks from sheets, preserving sheet, original
row, column and asterisk flags. Do not use the processed exports as if they
preserved experimental units. Exclude flagged values in the main descriptive
summary and include parseable flagged values in a parallel sensitivity.
Malformed nonnumeric entries (including TIC10 H16 in ROS) are missing and
recorded; do not silently repair them. The workbook contains TIC10 even though
the tidy export omits it. Retain its descriptive result but do not invent a
matched imaging compound. Empty sheets are not assays.

Summarize medians within each source sheet/condition, normalized to its own
mutant control median. Control vectors repeated across sheets must not be
pooled as independent measurements. Report source-vector checksums and any
duplication, plus count-of-values, explicitly not n biological replicates.
Normalize Greek/Latin micro-molar unit variants when parsing column labels;
never default an unrecognized treated concentration to zero.

At 1 micromolar, compare biochemical ratios descriptively with the nearest
follow-up imaging dose within 5%. Exact stripped compound-name matches are
the main cross-assay subset. Report the plausible source spelling/salt aliases
(SLK2001/SKL2001, Tyrphostatin 9/Tyrphostin 9, Chloroquine/Chloroquine
Phosphate, Sertraline/Sertraline HCl) separately and label their join as inferred.
Do not claim manufacturer confirmation alone verifies a mislabeled source well.

Both ROS and JC-1 directions are reported against mutant and corrected
controls. A drop in JC-1 is not automatically improvement. Do not infer an
uncoupling mechanism or mediation from cross-assay correlation. No biochemical
P values, biological confidence intervals, or independent-rescue labels.

## Outputs, interpretation and stopping

Produce control/model audit, well and compound scores, the complete 31-compound
benchmark, dose-response summaries, all fixed sensitivities, biochemical
extraction audit and summaries, a compact figure set, and a report with exact
commands, source checksums, software versions and deterministic verification.
Preserve all prior null results. Stop after this bounded analysis and assess
whether it adds mechanistic information; do not switch thresholds or compounds
to obtain a favorable result.
