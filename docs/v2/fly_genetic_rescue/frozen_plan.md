# Genetic interaction profiles and genotype-specific drug response

## Status and question

Frozen on 2026-09-15 before this project's extraction of numeric genetic-interaction profiles or calculation of treatment effects. This is a **retrospective, source-selected prediction benchmark**. Published findings, including treatment preferences and the authors' subgroup analysis, were visible during selection. A timestamped plan does not turn this into an independent confirmatory experiment.

Question: can a measured genetic-interaction profile predict the relative response to two interventions for a Parkinsonism genotype whose drug outcomes are withheld from training?

The experimental foundation is Kaempf et al., Nature Communications 17, 3761 (2026), DOI [10.1038/s41467-026-70303-8](https://doi.org/10.1038/s41467-026-70303-8). The public workbook supplies a 24-locus interaction matrix and both Q10 and R55 treatment arms for 18 mutant genotypes. The interventions and published pathways belong to that study. Our proposed contribution is an explicit prediction benchmark with genotype exclusion and simple competing explanations.

This is a new functional question. It neither replaces the original proteasome hypothesis nor changes any prior validation label. Positive prediction would be within-study experimental triangulation, not demonstration of a new molecular causal chain, human efficacy, or independent wet-lab replication.

## Inputs and eligibility

- Use only the two checksum-locked workbooks in `input_manifest.json`.
- Primary predictor: author-processed `Fig3d_GI_cluster` in the Source data workbook. Preserve all 24 anchor loci and the published continuous scores. This is not a reanalysis of raw pairwise electrophysiology.
- Primary outcomes: `Fig4b_Q10` and `Fig4d_R55`. Use all 18 genotypes having both drug and vehicle columns, with the exact locus/genotype mapping in `genotype_map.tsv`.
- Different allele dosages between interaction screening and treatment experiments remain visible in that mapping. Transfer is across locus-level profiles, not identical experimental genotypes.
- Blank cells are missing observations; nonblank cells must be finite numeric values. Do not pair observations by workbook row. Do not delete numerical outliers. Require at least three source observations per mutant arm and control arm, and positive arm means. No substitution if an eligibility requirement fails.
- Predictor distances require at least 18 jointly finite anchors after excluding the two loci being compared. Do not fill missing interactions with zero. Diagonal entries are always excluded. Report matrix asymmetry rather than silently symmetrizing it.
- A missing or nonpositive aggregate outcome, missing mapping, insufficient arm coverage, or incomplete prediction makes the primary `NOT_EVALUABLE`. Preserve the failure and stop; do not choose another outcome or dataset.

## Experimental units and limitations

Drug source tables do not preserve vial or experiment identifiers. Their rows cannot establish the number of independent treatment assignments. The prediction unit is the genotype: **18 units**, with equal weight. We will not conduct animal-row t tests, animal-level permutation tests, or animal bootstraps. Report source observation counts without representing them as independent biological replicates.

The outcome is normalized TH-positive innervation area. It is not a count of living neurons, a direct proteasome activity measurement, or the interaction-screen electrophysiology endpoint. Shared control normalization and related genotypes limit independence. Cross-validation separates drug outcomes within this study; it does not create a second cohort.

## Frozen outcome definition

For genotype g and drug d, let M(g,d) denote the arithmetic mean of all available source observations. Let v(d) be its own vehicle and WT be the corresponding control columns. Define:

`D(g,d) = log2[M(g,d) / M(g,v(d))] - log2[M(WT,d) / M(WT,v(d))]`.

The primary target is `S(g) = D(g,Q10) - D(g,R55)`.

Positive S indicates relatively greater Q10 response; negative S indicates relatively greater R55 response. Report both D values and the unadjusted within-genotype treatment ratios. Preference does not imply benefit if both interventions decrease the readout. Ratios use existing normalized source measurements, not reconstructed raw images. There is no pseudocount and no outcome-based gene exclusion.

## Prediction and baselines

1. Order the 18 loci lexicographically using `genotype_map.tsv`. For each held-out locus i, withhold its complete S outcome. Its measured interaction profile remains available. This is prediction for a genotype with an existing interaction profile, not zero-shot prediction for a new gene.
2. For each candidate training locus j, compute root-mean-square distance between their interaction profiles over jointly finite anchors, excluding anchors i and j. All scores have the same units, so no fitted rescaling is used. Break distance ties by lexicographic locus name.
3. Predict S(i) as the equally weighted mean S of its **three** closest training loci. No hyperparameter search, published subgroup labels, published best-drug calls, sequence model, or external pathway annotations enter the predictor.
4. Baseline A predicts the mean S of the other 17 loci.
5. Baseline B uses the same three-neighbor algorithm with a two-coordinate vehicle-severity profile: `log2[M(g,EtOH)/M(WT,EtOH)]` and `log2[M(g,H2O)/M(WT,H2O)]`. Euclidean distance; no fitted rescaling. These pretreatment measurements contain no active-drug arm values.

For each model, calculate mean absolute error (MAE) across the 18 predictions. The primary statistic is:

`T = min(1 - MAE_GI / MAE_mean, 1 - MAE_GI / MAE_severity)`.

Thus T measures the worse of the improvements over the two baselines. A zero baseline MAE makes improvement unevaluable. Also report Spearman correlation, all predictions, neighbors, and per-genotype absolute errors. Do not use correlation alone to declare success.

## Randomization benchmark and decision

Use NumPy `default_rng(20260915)` and 10,000 random permutations of S across the 18 genotype labels. Recompute all three sets of predictions and T for each permutation; predictor distances stay fixed. Compute `p = (1 + number[T_perm >= T_observed]) / 10001`.

This is an exploratory global label-exchangeability benchmark. It is not a conditional causal test of genetic interactions after adjustment for severity, and related loci, age differences and shared controls weaken its exchangeability assumption. There is one primary statistic; no pathway selection or multiplicity of primary outcomes.

Call the result `WITHIN_STUDY_PREDICTION_SUPPORTED` only if all of the following hold:

- T is at least 0.10 (at least 10% lower MAE than each baseline; a prespecified practical benchmark, not a clinical threshold).
- Randomization p is at most 0.05.
- The correlation between held-out predictions and S is positive.

Otherwise label an evaluable test `NOT_SUPPORTED`. These labels cannot be upgraded to independent replication or a proven mechanism.

## Sensitivities fixed before calculation

- One-neighbor and five-neighbor predictors, with matching neighbor counts in the severity baseline.
- Replace arm means by medians in both outcomes and the severity baseline, keeping three neighbors. If a required median is nonpositive, record that sensitivity as unevaluable.
- Exclude the four genotypes assayed at younger ages because of reduced survival: Pink1, nutcracker, iPLA2-VIA, gba. Recalculate the three-neighbor comparison on the remaining 14 genotypes. Excluded loci can remain predictor anchors, but their outcomes cannot train the model.
- Supplementary Dataset 1 contains author Bayesian interaction summaries. Audit coverage and identity only; do not fit or select a second primary interaction model from it.
- The four available climbing-assay genotypes in `Suppl Fig12b_Q10_SING` and `Suppl Fig12c_R55_SING` provide a descriptive comparison of the direction of treatment preference. Their source points represent groups; no replication p value is allowed. Do not use this subset to train or tune anything.

Sensitivities receive error and correlation summaries; they cannot rescue the primary label. No further dataset or model search follows the observed result within this extension.

## Reproduction and reporting commitments

Publish a source/arm audit, exact mapping, arm aggregates, both drug effects, genotype predictions and neighbors, all sensitivity results, permutation summary, one scientific figure, dependency versions and checksums. Verify numeric output stability in a second run and use a fresh checksum-checked input acquisition. Include mathematical checks that the held-out outcome cannot influence its own prediction and that genotype-specific common multiplicative scale cancels from the response contrast.

The previous full V2/R clean-clone release gate remains separate and unfinished. This extension cannot declare the complete submission technically finished.
