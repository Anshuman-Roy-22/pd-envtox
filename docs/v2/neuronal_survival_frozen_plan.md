# Fixed neuronal survival reanalysis plan

Date: 2026-09-14. This plan is committed before inspecting the target-specific
survival phenotypes. It is a secondary analysis of published experiments, not
a new wet-lab experiment or external preregistration.

## Question and scope

Do the existing eight shared proteasome leading-edge genes exhibit an unusual
change in **relative neuronal survival dependence** when antioxidant support is
removed, compared with other genes of similar baseline neuronal dependence?

The fixed eight are POMP, PSMA4, PSMB4, PSMC5, PSMD1, PSMD2, PSMD4, and PSMG1.
They were selected in the preceding PD analyses, not from these CRISPR effects.
The full 52-member Reactome set remains exactly as recorded in
`perturbation_feasibility/candidate_manifest.json`.

The metadata audit supports this bounded CPU analysis. A direct neuronal
transcriptome-prediction benchmark for these genes is not currently established.
The 2021 CROP-seq libraries contain none of the eight; the older study's full
guide mapping is not yet resolved. No State or GEARS result is implied.

## Fixed experimental source

Tian et al., Nature Neuroscience (2021), DOI 10.1038/s41593-021-00862-0,
Supplementary Table 1, publisher file `41593_2021_862_MOESM3_ESM.xlsx`.

URL:
https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41593-021-00862-0/MediaObjects/41593_2021_862_MOESM3_ESM.xlsx

SHA256:
`82785da6b6fd4f9732a1e6bd701ab7fc2aa4d12f0d33fddc676054e37841a464`

Use only these four sheets:

- `Survival_plusAO_CRISPRi`
- `Survival_noAO_CRISPRi`
- `Survival_plusAO_CRISPRa`
- `Survival_noAO_CRISPRa`

The paper's methods specify day-10 neurons grown with (+AO) or without (-AO)
antioxidants, compared with a day-minus-3 iPSC reference. Thus the readout can
include differentiation and survival effects. It is not a pure acute challenge
to otherwise identical mature neurons. Both conditions share the same nominal
endpoint day. The published Phenotype column is centered and scaled against
negative-control quasi-genes separately for each screen.

Only workbook structure, gene/TSS labels, library coverage and methods have
been examined before this plan. The workbook was downloaded for this metadata
audit, but numerical target phenotypes have not been inspected. The metadata
snapshot gives 49/52 pathway genes for CRISPRi and 51/52 for CRISPRa; all eight
lead genes are present in both modes. No aliases are silently substituted.

## Preprocessing

1. Require the locked workbook checksum. Read only `Gene`, `TSS`, and
   `Phenotype` for analysis. Do not use published P values to select promoters,
   targets or contrasts.
2. Trim surrounding whitespace in gene symbols. Exclude empty symbols and
   names beginning, case-insensitively, with `negative`, `non-target`,
   `nontarget`, or `control`. Retain other author-provided symbols as-is.
3. Require finite Phenotype values. Aggregate multiple rows for a gene by their
   unweighted mean. Record the number of rows and unique TSS labels. This is
   promoter aggregation, not biological replication.
4. Within each CRISPR mode, use the intersection of genes with finite means in
   both AO conditions. Require all eight primary genes; otherwise report
   `NOT_EVALUABLE`, with no replacement target set.
5. Compute each gene's percentile in the gene-level Phenotype distribution
   within each condition using `(average_rank - 0.5) / N`. Lower values indicate
   more deleterious perturbations relative to other genes in that screen.
6. Define `delta_rank = percentile_noAO - percentile_plusAO`. Negative values
   indicate relatively greater dependence without antioxidants. This is a
   relative rank change, not an absolute survival difference or a fitted
   treatment-by-gene interaction.

## Primary comparison and reference model

The primary statistic is the mean delta_rank across the eight genes for
CRISPRi. Compare it with 100,000 randomly sampled eight-gene sets matched to the
primary set's composition across:

- 20 fixed baseline percentile bins: `min(19, floor(20 * percentile_plusAO))`;
- one versus multiple source rows in the +AO sheet.

The control pool excludes all 52 proteasome pathway genes. Within each stratum,
sample without replacement within a draw, preserving the number of target genes
in that stratum. Draws are independent. Require at least max(20, 5 * number of
targets in the stratum) eligible control genes per occupied stratum; otherwise
report `NOT_EVALUABLE`, without adapting the binning to obtain significance.

Use NumPy PCG64 with seed 20260914. Calculate both empirical tail probabilities
with the plus-one correction `(1 + count)/(100000 + 1)`. The primary P value is
`min(1, 2 * min(lower_tail, upper_tail))`. The primary threshold is two-sided
P < 0.05. Report the observed statistic, reference mean, observed-minus-reference
effect, and central 95% reference interval. The interval is **not** a confidence
interval from biological replicates.

This is a conditional gene-set reference test. Its validity depends on
exchangeability of target and background genes within the matching strata.
Correlated proteasome genes, unmeasured expression/knockdown differences, and
shared screen errors can violate that assumption. It is not a replicate-level
experimental interaction test. Matching baseline neuronal dependence addresses
one form of generic essentiality; it does not establish neuron specificity
relative to other cell types.

## Prespecified secondary and sensitivity summaries

Use the same method for the following comparisons, with separate recorded RNG
seeds, and label their P values exploratory. They cannot replace the primary:

- full evaluable pathway, CRISPRi, 20 bins, seed 20260915;
- shared eight, CRISPRa, 20 bins, seed 20260916;
- shared eight, CRISPRi, 10 bins, seed 20260917;
- shared eight, CRISPRi, 40 bins, seed 20260918.

For the fixed eight, also report gene-level Phenotype values and percentiles in
all four conditions and leave-one-gene-out primary effect summaries against
the corresponding reference expectation. These are robustness diagnostics,
not eight new confirmatory tests. No target-specific P values will be newly
computed from single summary measurements.

## Interpretation and stopping rules

- Primary P < 0.05 supports an unusual conditional rank change under the stated
  gene-set reference model. A negative excess indicates relatively greater
  vulnerability; a positive excess indicates relative attenuation.
- A null does not establish equivalence. An opposite direction must be reported.
- Do not call this independent PD replication, neuronal selectivity, therapeutic
  rescue, or validation of an AI-generated transcriptome.
- Do not treat CRISPRi and CRISPRa as independent human cohorts. Strong knockdown
  and overexpression also need not mimic modest disease-associated expression
  changes or reverse one another.
- The summary tables do not expose the replicate-level uncertainty needed for
  a biological interaction test. No fabricated donor count, standard error,
  power estimate or experimental confidence interval will be supplied.
- No replacement pathways, new target genes or further outcome-driven searches
  follow this test. Existing frozen validation labels remain unchanged.

## Reproducibility deliverables

Provide a checksum-verified acquisition command, analysis script, dependency
versions, gene-level input-derived table, matching inventory, primary and
secondary summaries, saved deterministic reference draws, and an explicit
limitations report. Commit this plan separately before running the effects
analysis. Complete the analysis regardless of its direction or significance.
