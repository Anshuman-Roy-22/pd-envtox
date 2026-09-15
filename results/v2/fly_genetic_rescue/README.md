# Genetic interactions and intervention response: NOT_SUPPORTED

Date: 2026-09-15. **This extension does not establish a stronger molecular mechanism or a validated treatment-selection model.** The fixed predictor performed worse than the average-response baseline. It is an auditable supporting analysis, not a new headline mechanistic result for the AAN submission.

## Experimental foundation and original contribution

We used the public Source data and Supplementary Dataset 1 workbooks from Kaempf et al., *Nature Communications* 17, 3761 (2026), DOI https://doi.org/10.1038/s41467-026-70303-8. The authors performed the genetic perturbations and drug interventions. Their molecular subgroup and rescue findings are prior work.

The original publication and associated source data are credited under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). This project extracted the interaction matrix, aggregated treatment measurements, calculated response contrasts, and generated new prediction tables and a figure; these transformations are distinct from the authors' original analyses.

Our distinct question was whether an author-measured genetic-interaction profile predicts relative drug response when a genotype's drug outcomes are withheld from training. The predictor uses the published 24-locus interaction matrix; 18 genotypes have both Q10 and R55 arms with their respective vehicles. The interaction screen measures electrophysiology; the drug outcome measures TH-positive innervation area. These are different endpoints and, for several loci, different allele dosages. Exact genotypes are preserved in the mapping and result tables.

The source was selected after published findings were visible. Commit `5d6968a49890cf724566c762aae2a9299f878efe` fixed the numeric analysis before treatment effects were extracted or calculated. This is retrospective within-study cross-validation, not an independent confirmatory dataset. There was no search over models after seeing the result.

## Primary result

| Predictor | Genotype-withheld mean absolute error |
|---|---:|
| Average response of other genotypes | 0.449791 |
| Three neighbors based on vehicle severity | 0.491310 |
| Three neighbors based on genetic interactions | 0.500592 |

The genetic-interaction model had **11.29% higher error** than the average-response baseline and **1.89% higher error** than the severity baseline. Spearman correlation was 0.2479. The prespecified improvement statistic was -0.112944; the 10,000-permutation randomization p value was **0.486651**. The classification is **`NOT_SUPPORTED`**.

The frozen success rule required at least 10% improvement over both baselines, p at most 0.05, and positive correlation. The first two conditions failed. Correlation alone cannot support the model. This result pertains to the specified predictor and endpoint; it does not prove that every interaction-based model would fail or contradict the original authors' drug experiments.

## Methods in brief

For each genotype and drug, take the log2 ratio of mean treated to mean vehicle measurements, then subtract the corresponding treated-to-vehicle log2 ratio in controls. The prediction target is the adjusted Q10 response minus the adjusted R55 response. Both drug responses are reported separately so that relative preference cannot hide deterioration under both treatments.

For each withheld genotype, choose three training genotypes with the closest interaction profiles. Distances exclude both comparison loci as anchor coordinates; 22 anchors remain for every comparison. The outcome from the withheld genotype never enters its prediction. The pretreatment-severity comparator uses the two vehicle measurements relative to their respective controls. No published subgroup or best-drug label is a predictor.

The analysis unit is one genotype, equally weighted. Source observations do not retain vial/batch identifiers; consequently, we make no animal-row significance claim and provide no animal-bootstrap uncertainty interval. The randomization p value is a global label-exchangeability benchmark. It is not a conditional causal test, and it assumes more genotype exchangeability than can be independently established in this study.

## Prespecified sensitivity results

| Analysis | Genotypes | Error change versus average-response baseline |
|---|---:|---:|
| Primary, three neighbors | 18 | 11.29% worse |
| One neighbor | 18 | 4.07% worse |
| Five neighbors | 18 | 10.52% worse |
| Median arm values | 18 | 5.60% worse |
| Exclude four younger-age genotypes | 14 | 25.26% worse |

None improved on the average-response baseline. Sensitivities cannot replace the primary result. After the age restriction, the prediction correlation was also negative (-0.2885).

For the four genotypes with both climbing and innervation measurements, relative treatment-preference directions agreed in **4/4**: Q10 for Pink1 and nutcracker; R55 for DJ-1a/b and iPLA2-VIA. This was a prespecified descriptive check, shares a source study, uses a selected subset, and has no independent replication p value. It does not rescue the prediction claim or establish a new mechanism. Climbing table column names explicitly identify that endpoint.

## What this means for the project

Actual perturbation and intervention experiments offer a better mechanistic foundation than attaching an unvalidated structural model to transcriptomic enrichment. However, this reanalysis did not establish the proposed bridge between interaction profiles and differential rescue. The original proteasome results remain `PRIMARY_ONLY`; the prior functional-rescue benchmark remains `BENCHMARK_NOT_EVALUABLE`.

This addition supports transparency and tests a new question, but it should not displace the main study's coherent evidence narrative. It does not justify claiming that the project now has a demonstrated causal mechanism or a reliably higher probability of winning AAN. No computational method can supply an unmeasured biological intervention outcome.

A further mechanism-focused extension should be qualified against a specific causal chain before implementation:

1. A defined candidate perturbation changes a measured molecular intermediate.
2. A disease-relevant neuronal endpoint changes under that perturbation.
3. Reversal of the perturbation or an appropriate dependency intervention changes the effect as predicted.
4. The actual source includes control arms, independently identified biological replicates, and an independent perturbation or model for verification.

Reanalysis of such existing experiments is feasible without the student running a wet lab. The two workbooks used here do not supply this candidate-specific design for the original proteasome lead. A more flexible model or another favorable subset would not resolve that missing experiment. No additional dataset search was performed after this benchmark result.

## Reproduce this extension

The tested runtime was CPython 3.12.14. Direct dependency pins and the installed transitive dependency closure are under `docs/v2/fly_genetic_rescue/`. These are version pins, not hashes of platform-specific wheels. Linux byte identity was checked; identical rendering on every operating system is not asserted.

In an environment with those packages installed, from the repository root:

```bash
python scripts/run_fly_genetic_rescue.py --fetch
```

For a separate Windows environment, these commands do not modify the Git branch:

```powershell
py -3.12 -m venv ..\pd-fly-venv
..\pd-fly-venv\Scripts\python.exe -m pip install -r docs\v2\fly_genetic_rescue\environment_freeze.txt
..\pd-fly-venv\Scripts\python.exe scripts\run_fly_genetic_rescue.py --fetch
```

Verification alone:

```bash
python scripts/verify_fly_genetic_rescue.py
```

The input cache is ignored by Git. Downloads must match the two recorded SHA256 hashes. This standalone extension does not close the preexisting whole-project master-runner/R clean-clone gate.

## Outputs and figure caption

- `arm_aggregates.tsv`: all source arm observation counts, means and medians for both endpoints.
- `interaction_matrix.tsv`, `interaction_audit.json`: fixed author scores, matrix checks and supplementary summary-label coverage.
- `genotype_effects.tsv`: both responses, controls, vehicle severity and exact allele mapping.
- `heldout_predictions.tsv`, `prediction_neighbors.tsv`: every prediction, error and training neighbor.
- `sensitivity_summary.tsv`, `climbing_descriptive_comparison.tsv`: all fixed secondary analyses.
- `permutation_statistics.tsv`, `validation_summary.json`: the complete randomization reference and primary result.
- `genetic_interaction_prediction.png` and `.pdf`: 300-DPI raster and vector figure.
- `mathematical_checks.json`, `environment.json`, `source_inventory.json`, `reproduction_check.json`, `release_checksums.json`: verification records.

**Figure.** A: pooled, control-adjusted intervention responses for all 18 genotypes. No animal-level confidence intervals are shown because batch/vial information is unavailable. B: genotype-withheld prediction error; lower is better. C: the genetic-interaction model's improvement over the average-response baseline under each fixed analysis; all values are negative. All values are derived from tracked result tables. Original experiments: Kaempf et al. (2026); prediction benchmark: this project.
