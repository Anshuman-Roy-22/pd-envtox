# Neuronal survival follow-up: complete result

**Primary result: `NO_CONDITIONAL_RANK_SHIFT`, two-sided empirical P = 0.14558.**

The eight existing proteasome leading-edge genes moved toward greater relative
dependence in neurons without antioxidant support, but the shift was not
unusual enough under the fixed matched-gene reference model to pass the
prespecified threshold. This adds a transparent analysis of real experimental
measurements. It does not provide the independently validated mechanism or
successful AI benchmark sought for a major project improvement.

The analysis plan was committed at `6d493a3` before target-specific survival
phenotypes were inspected. The source is Supplementary Table 1 of Tian et al.
(2021), DOI 10.1038/s41593-021-00862-0. This is a reanalysis of their experiments,
not an experiment performed for this project.

## Results

| Fixed comparison | Target genes | Two-sided reference P | Interpretation |
|---|---:|---:|---|
| Primary: shared eight, CRISPRi, 20 baseline bins | 8 | 0.14558 | Primary threshold not met |
| Sensitivity: shared eight, CRISPRi, 10 bins | 8 | 0.13380 | Also nonsignificant |
| Sensitivity: shared eight, CRISPRi, 40 bins | 8 | 0.13694 | Also nonsignificant |
| Secondary: shared eight, CRISPRa | 8 | 0.58823 | Exploratory and nonsignificant |
| Secondary: full evaluable pathway, CRISPRi | 49 | Not evaluable | One matching stratum had 18 eligible controls, below the frozen minimum of 20 |

The primary statistic, mean change in percentile rank from +AO to -AO, was
-8.74 percentile points. Its exact matched-background expectation was +1.06
points, for an excess shift of **-9.80 percentile points**. Five of eight genes
had a negative raw rank change; seven had a negative change relative to their
stratum-specific reference expectation. Removing any one gene retained a
negative excess mean, ranging from -13.59 to -6.39 percentile points. These
directional diagnostics do not turn the primary null into confirmation.

For clarity, the 100,000 simulated reference sets had a mean of +1.09 percentile
points and a central 95% reference interval from -11.98 to +16.20 points. That
interval describes matched random sets. It is not a biological confidence
interval for the candidate set's effect.

The full pathway secondary comparison was retained as `NOT_EVALUABLE` under
the frozen matching rule. Its cutoff was not loosened after seeing the primary
result. No alternative gene set, contrast, one-sided threshold or model was
substituted to obtain a positive result.

## What the experiment and analysis can establish

The laboratory experiment perturbed genes in iPSC-derived neurons. The source
readout compares day-10 neurons with an earlier iPSC reference, under either
standard antioxidant support or antioxidant withdrawal. It can reflect effects
on differentiation as well as survival. Strong CRISPR knockdown or activation
does not necessarily reproduce the modest expression changes observed in PD.

The analysis matches random gene sets on baseline neuronal dependence and on
one versus multiple source rows per gene. It asks about **relative dependence
within these screens**. It does not compare neurons against another cell type,
estimate an absolute survival interaction, or establish selective dopaminergic
vulnerability. The summary workbook does not provide replicate-level errors
needed for those inferences.

The reference test assumes that targets and background genes are exchangeable
within the chosen strata. Correlated proteasome components, shared screen
errors, expression differences and variable knockdown efficiency remain
limitations. Its P value is a conditional gene-set statistic, not the P value
of a replicated biological interaction. A nonsignificant result also does not
demonstrate equivalence or absence of biological effects of individual genes.

Existing human and toxicant validation classifications are unchanged. No AI
model has been trained or validated by this analysis.

## Reproduction

From the repository root, using Python 3.12.14:

```bash
python -m pip install -r requirements-neuronal-survival.txt
python scripts/neuronal_survival_reanalysis.py --fetch
```

The command downloads the public workbook only if absent, requires the exact
SHA256, and generates the result tables and deterministic reference files.
Default input: `data_raw/neuronal_perturbation/tian2021_s1.xlsx`.
Default output: `results/v2/neuronal_survival`.

For an already downloaded input:

```bash
python scripts/neuronal_survival_reanalysis.py --input /path/to/tian2021_s1.xlsx
```

The package file pins the analysis packages and their required pure-Python
dependencies. Exact Python and package versions used are recorded in
`provenance.json`. Results were reproduced in a separate output directory and
then in a fresh virtual environment installed from these pins, using a newly
downloaded checksum-verified input. This is not a claim of full V2 clean-clone
or cross-platform verification.

Generated files include:

- Gene-level means, promoter row counts, percentiles and rank changes for both
  CRISPR modes, with raw input-derived values retained.
- Complete primary/secondary/sensitivity comparison summary, including the
  unevaluable comparison.
- Target membership and matching-stratum inventory for each evaluable test.
- Four files containing all 100,000 reference draws each.
- Leave-one-gene-out effect diagnostics, with no additional confirmation tests.
- Input, plan, script and environment provenance.

## Integrity audit

An initial compressed sensitivity file was incomplete when checked against the
independent rerun. The numerical summaries agreed. Reference serialization was
changed to atomic writes, with decompression verification before replacement
and a final 100,000-row check. The final outputs, not that incomplete file, are
the release artifacts. Numerical results and frozen decision rules were not
changed by this repair.

`verification.json` records the final reproducibility and independent
calculation checks. `SHA256SUMS.txt` covers the final results and supporting
analysis files; run it from the repository root with `sha256sum -c` on systems
that provide that utility. The repository-wide R clean-clone gate remains a
separate unfinished task.

## Source

Tian et al. (2021), *Genome-wide CRISPRi/a screens in human neurons link
lysosomal failure to ferroptosis*:
https://www.nature.com/articles/s41593-021-00862-0

The exact publisher download URL and checksum are recorded in `provenance.json`
and `docs/v2/neuronal_survival_frozen_plan.md`.
