# Functional-rescue reanalysis: completed, benchmark not evaluable

Date: 2026-09-15. **This extension does not establish a stronger molecular
mechanism or an independently validated rescue predictor.** The data do
contain repeatable marker responses, but the transferred score and combined
response definition fail important control checks. The primary label is
`BENCHMARK_NOT_EVALUABLE`, not a successful prediction result or proof that
the compounds cannot rescue neurons.

## What was done

We reanalyzed published drug-intervention measurements in human
SNCA-triplication dopaminergic neurons and genetically corrected controls.
The original experiments and published drug findings belong to
Gorgogietas et al., *Scientific Reports* (2025),
https://doi.org/10.1038/s41598-025-14735-0.
The source is an immutable author-repository snapshot, commit
`1ed064f3d7881e7f520ae6afefcc8aed75ffc6ba`.

The independent contribution here is a specified computational comparison
and its audit, using an alpha-synuclein-only score versus a four-feature
regularized linear discriminant trained on a separate control experiment.
No drug labels enter model fitting. Published follow-up selection remains
a source of bias: the 31-compound matched panel is not a random test set.
Technical wells and plate names are not treated as independent donors.

| Component | Analyzed observations |
|---|---:|
| Separate control experiment | 144 wells; 48 mutant and 48 corrected wells used for fitting |
| Initial compound screen | 4,800 wells; 1,020 compounds |
| Follow-up imaging | 930 wells; 33 compounds at six doses, 198 compound-dose combinations |
| Exact catalogue and matched-dose comparison | 31 compounds, four wells each per stage |
| ROS and JC-1 extraction | 15,201 parseable values across 20 compound sheets per assay |

The default matched dose is 5 versus 5.13 micromolar. Alectinib is matched
at 2.5 versus 2.51 micromolar and retained despite an unexplained exclusion
in the author notebook. Its exclusion is a fixed sensitivity. Dasatinib
and Prostratin have no exact catalogue match in the screen and do not enter
the 31-compound predictive comparison. Prostratin is retained as a control.

## Prespecification and the source-scale correction

Commit `61e73fd` froze the design before compound scoring. The first input
check then showed that nominal follow-up lengths, intensities and fractions
include normalized negative values. Physical ratios and logarithms would
therefore be invalid. That run stopped before fitting or ranking compounds.
It had normalized screen values in memory; we do not call it completely
value-blind. The control-only failure audit is preserved.

Commit `57cdfae` recorded a necessary amendment before fitting: center each
feature by the plate's mutant-control median and divide by its pooled
within-genotype control SD. This is invariant to a common positive affine
feature transformation. Raw living-nuclei counts retain count ratios.
Unknown nonlinear transforms, segmentation differences and source-version
changes cannot be repaired by this step. The original specification remains
preserved beside the amendment and active `analysis_config_v2.json`.

The fixed response requires, in at least 75% of at least three valid wells:
marker score >=0.5, each of three neuronal features and living fraction
>=-0.5 pooled control SD, and living count >=80% of mutant vehicle.
These are operational cutoffs, not validated safety, noninferiority or
functional-restoration criteria. Score 0.5 is half the training control
separation; it does not mean a 50% absolute protein reduction.

## Primary result

| Fixed endpoint | Result |
|---|---:|
| Eligible matched compounds | 31 |
| Follow-up joint positives | 0 |
| Required positives for AP comparison | At least 5 |
| AP difference, either AUC, leave-one-compound-out AP gain | Undefined |
| Decision | `BENCHMARK_NOT_EVALUABLE` |
| Marker rank repeatability, Spearman rho | 0.838306 |
| Multivariate rank repeatability, Spearman rho | 0.673387 |

The seven planned primary/sensitivity configurations all have zero joint
positives. These include the looser and stricter preservation rules,
shrinkage 0.1/0.9, excluding Alectinib, and collapsing Bazedoxifene salts.
Across all follow-up doses, none of the 198 compound-dose entries meets
the joint rule at any of the three preservation settings. This is a rule
output, not evidence of universal biological nonresponse.

No invalid well was removed for missing/nonfinite selected features.
There are no inferred biological P values or confidence intervals.
Undefined AP metrics are explicitly missing, including leave-one-out gain;
they are never reported as measured zero performance.

## Why zero positives cannot establish absence of rescue

The model separates genotype controls in the training experiment well:
leave-one-plate-out AUC is 0.9306-1.0000 for the marker and 0.9375-1.0000
for the multivariate model. That is assay QC, not clinical accuracy.
The calibration does not transfer reliably to the follow-up experiment:

| Corrected-vehicle controls | Mean marker score | Mean multivariate score | Jointly passing wells |
|---|---:|---:|---:|
| Training experiment | 1.0000 | 1.0000 | 4/48 |
| Initial screen | 0.8923 | 0.9476 | 14/240 |
| Follow-up | 0.2038 | 0.0937 | 0/48 |

The preservation rule itself is restrictive and directionally problematic.
It assumes that decreases from mutant levels are undesirable for every
neuronal feature. Yet corrected controls have lower values for some of
these features. In follow-up, corrected TH-neurite length averages
-3.10 pooled control SD relative to the mutant median. It cannot satisfy
the rule's -0.5 SD floor. Prostratin also has zero jointly passing wells
in follow-up (0/44), and only 2/44 pass the transferred marker threshold.

Thus the reanalysis exposes a limitation in **our transferred calibration
and endpoint definition**. It does not refute the source authors' rescue
experiments. We retain the locked outcome and do not retroactively choose
thresholds, features or a new correction to manufacture positives.
Any subsequent recalibration on these data would be exploratory and would
need a new untouched experiment to support a validation claim.

## The useful descriptive observation

Three compounds meet the fixed marker-only response threshold in all four
matched follow-up wells. Each also has fewer living nuclei than its own
plate's mutant-vehicle controls:

| Compound, 5.13 micromolar | Median marker score | Median living-count ratio | Count decrease |
|---|---:|---:|---:|
| SNS-032 | 1.0633 | 0.7309 | 26.9% |
| Flavopiridol | 1.0513 | 0.6168 | 38.3% |
| Halofuginone | 0.8273 | 0.4135 | 58.6% |

All four wells for each compound are below the 80% count floor. This
supports caution in interpreting marker lowering as beneficial rescue.
The marker is measured per individual cell, so reduced counts do not
automatically explain the lower intensity mathematically. The two changes
co-occur. We have not established whether death, detachment, segmentation,
changes in cell state, or another process causes the count difference.
These are descriptive measurements from the selected source experiments,
not newly established drug toxicities or a proteasome-dependent mechanism.

## Biochemical results and limits

Original workbooks, rather than melted exports alone, were parsed with
sheet, row, column and asterisk flags retained. There are 3,054 flagged
parseable values; main summaries exclude them and the parallel sensitivity
includes them. No 5-SD trimming or unrecorded corrections were applied.
The malformed ROS TIC10 H16 value `1..099436` remains missing. TIC10 is
preserved in the workbook summaries despite being absent from the tidy
export and unmatched to imaging. Repeated control vectors are audited by
checksum and never pooled as additional independent experiments.

There are 15 exact compound-name cross-assay joins and four explicitly
inferred spelling/salt joins. Comparisons use biochemical 1 micromolar
versus imaging 1.03 micromolar. Examples of descriptive assay changes:

| Source compound label | ROS / mutant median | JC-1 / mutant median | Mapping |
|---|---:|---:|---|
| FCCP | 0.6302 | 0.3883 | Exact name |
| Tyrphostatin 9 | 0.6543 | 0.3602 | Inferred spelling match to Tyrphostin 9 |

Lower ROS and JC-1 values do not on their own establish neuroprotection.
Biological replicate membership is not recoverable from these workbooks;
row counts are numbers of values, not independent cultures. No biochemical
significance tests, mediation estimates or causal chains are claimed.
The original authors already report mitochondrial-uncoupling findings;
their mechanism is not an original discovery of this project.

The two archived proteasome-activity files have separate MG132 and drug
columns, not drug-plus-MG132 combination arms. They cannot test whether
blocking proteasome function abolishes drug rescue. They were acquired
for feasibility but are not promoted to a dependency experiment here.

## Implication for AAN and the next scientific decision

This is a completed and reproducible analysis attempt with useful assay
diagnostics. It does **not** deliver the mechanistic breakthrough sought
for AAN. It does not change the earlier human `PRIMARY_ONLY` result, the
held-out toxicant nulls, or the nonsignificant neuronal survival comparison.
There is no defensible estimate of a maximum winning probability from
these results. The most rigorous presentation would retain this as a
qualified exploratory appendix, rather than elevate it to the central
mechanistic claim.

A further mechanism-directed extension must first supply the missing
experiment, not just a larger model. For the proteasome hypothesis the
minimum useful design is an intervention crossed with disease context,
with identified biological replicates, and an orthogonal neuronal endpoint.
A stronger dependency test additionally crosses the rescue intervention
with pathway blockade or genetic perturbation in the same system. The
key estimand is whether the intervention's benefit changes when that
pathway is perturbed, while accounting for the perturbation's own harm.

No such complete public matched dataset has been established in the
current audit. Raw experiment identifiers, matching normalization records
and combined-treatment measurements are concrete missing inputs; their
existence should be checked before another analysis is promised. Obtaining
them could enable a defensible experimental reanalysis within the available
month. Without them, additional fitting of the present files cannot resolve
proteasome dependence. This package stops at the fixed analysis and preserves
the outcome; it makes no new-data discovery or causal-validation promise.

## Figures and result inventory

- `fig1_benchmark_and_calibration.png/.pdf`: transferred control calibration,
  all 31 marker-score pairs, the three marker responders' living counts,
  and the failed joint-positive gate. Dots are measured wells or plate
  means, as labeled; no biological error bars are invented.
- `fig2_all_matched_response_components.png/.pdf`: all 31 compounds and
  each individual component's pass fraction. Sorting is descriptive.
- `screen_compounds.tsv`, `followup_doses.tsv`, `matched_benchmark.tsv`:
  complete outcomes, including impaired counts and failed criteria.
- `control_*`, `normalization_audit.tsv`, `matched_component_pass_rates.tsv`:
  calibration and component diagnostics.
- `sensitivity_summary.tsv`, `leave_one_compound_out.tsv`: all fixed checks.
- `biochemical_*`, `cross_assay_descriptive.tsv`: source-preserving secondary
  extraction and summaries, including duplication and parse issues.
- `analysis_summary.json`, `environment.json`, `verification.json`,
  `reproduction_verification.json`, `SHA256SUMS.txt`: conclusion and provenance.

## Reproduction

From the repository root, using the observed Python 3.12.14 Linux environment:

```bash
python -m pip install -r requirements-functional-rescue.txt
python scripts/run_functional_rescue.py --fetch
```

The runner verifies/acquires the sixteen pinned public source files,
rebuilds design joins, runs the fixed analysis and figures, then checks
release hashes and numerical edge cases. It never executes author notebooks.
Dependency versions are pinned to the observed installed closure; this is
not a promise of bitwise agreement on an untested operating system.
The additional normalization notebook is recorded by immutable URL/hash
for method provenance but is not a numerical execution dependency.

Results were checked by an isolated reproduction of this extension, recorded
in `reproduction_verification.json`. Exact source/result checksums, undefined
metrics, tied-score metrics, affine-invariant normalization, raw count-ratio
reconstruction and gzip integrity are covered by `verify_functional_rescue.py`.
This verification does not close the separate repository-wide master V2
runner and R clean-clone release gates. Those remain outstanding.
