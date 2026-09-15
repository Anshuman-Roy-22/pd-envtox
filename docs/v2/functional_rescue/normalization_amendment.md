# Source-scale amendment, before model fitting or compound ranking

Date: 2026-09-15. Original plan commit: `61e73fd`.

The initial run stopped at follow-up input calibration before fitting a
model. Every follow-up plate had zero mutant-control wells nonnegative in
all six selected fields. Several nominal intensity, neurite-length and
living-fraction columns contain negative normalized values. Living-nuclei
counts retain an uncentered count scale. Thus ratios/logarithms of all
columns, and interpreting 0.8 as 80% of physical neurite length, are invalid.
`normalization_failure_audit.json` records the control-only evidence.

The first attempt normalized initial-screen features in memory before the
follow-up gate failed. No fitted model, compound score, ranking, statistical
result or target-effect table was inspected or exported. This amendment is
motivated by source units, not an unfavorable drug result. The original plan
and configuration remain preserved. No thresholds from completed earlier
studies are changed.

The author's primary-screen processing notebook supplies per-feature positive
affine normalization functions (`z_std`, `z_robust`, mean, median and min-max).
Its immutable URL/hash are recorded in `additional_method_source.json`.
It does not establish exact raw-to-follow-up transform parameters. We therefore
do not reconstruct absolute physical units from these normalized columns.

## Replacement normalization

Use `analysis_config_v2.json`. For every selected feature on each plate:

`z = (value - median_mutant_control) / pooled_within_genotype_control_SD`.

The pooled SD is the square root of the sum of within-mutant and within-
corrected squared residuals divided by n_mutant+n_corrected-2. It excludes
the between-genotype difference. It must be positive and finite. Keep the
minimum of eight finite controls per genotype. Accept negative centered
feature values; require living-cell counts to be nonnegative. Preserve zero
cell counts as possible depletion, not an outlier to discard.

This transformation is invariant to a common positive affine transform of
each feature within each plate. Its applicability assumes the deposited
features retain their original ordering, consistent with the author code.
It does not correct unknown nonlinear transforms, changed segmentation,
missing observations or assay-selection bias. Report these limits.

Fit the same four-feature regularized linear discriminant and marker-only
baseline to these standardized values in the separate Fig1 control data.
Keep feature membership, shrinkage, matched catalogue IDs/doses, validation
folds, minimum sample rules, AP endpoint and model-success criteria unchanged.
The score now describes a shift in control-standardized feature coordinates;
0.5 is half the training control separation, not a log fold-change or 50%
absolute protein reduction.

## Replacement preservation rule

At the primary setting a well must satisfy all of:

- Marker-only score >=0.5.
- Each of MAP2 neurite length, TH neurite length and TH intensity has z >=-0.5.
- Living-nuclei fraction has z >=-0.5.
- Living-nuclei count is >=0.8 times its plate's mutant-control median.

Keep the requirement that at least 75% of valid wells pass with >=3 wells.
This is an operational tolerance of half a pooled control SD for the
normalized features and 20% for counts. It does not establish biological
noninferiority, safety or restoration. Report every separate component.

The previously planned looser/stricter sensitivities become:

| Label | Minimum living-count ratio | Minimum normalized integrity/fraction z |
|---|---:|---:|
| Looser (0.7) | 0.7 | -0.75 |
| Primary (0.8) | 0.8 | -0.50 |
| Stricter (0.9) | 0.9 | -0.25 |

The numeric labels retained in output columns refer to the **living-count
ratio only**. Structural fields are reported as standardized changes, never
as physical ratios. No other primary or secondary policy is changed.

Commit this amendment before rerunning the analysis. Report the original
calibration failure and this source-scale correction alongside all results.
