# Neuronal perturbation extension: initial feasibility audit

**Update, 2026-09-15:** The user supplied the original 2019 Table S2. Its
183 guides name 89 genes plus controls, with zero of the fixed 52 pathway
genes and zero of the shared eight. In each neuronal mapping file, 55 guide
sequences map to the key (26 non-control genes); two additional 20-nt
sequences remain unidentified. Thus no fixed-target benchmark is established.
Do not infer the unidentified genes or claim their absence has been proven.
The source workbook, reproducible audit and complete mapping counts are now
recorded. See `followup_status.md` and
[`../functional_rescue/README.md`](../functional_rescue/README.md) for the
separately authorized redesign. Earlier records below remain historical.

Metadata retrieved: 2026-09-13. Reproduction verified: 2026-09-14. Existing scientific source: commit
`4b2c005bcc93ce28c3ba4c52fa71785fd7d05e1e` on `aan-v2-overhaul`.

**Update, 2026-09-14:** The 2021 survival screens cover all eight lead genes.
A separately committed protocol and completed CPU reanalysis are in
[`neuronal_survival_frozen_plan.md`](../neuronal_survival_frozen_plan.md) and
[`results/v2/neuronal_survival/README.md`](../../../results/v2/neuronal_survival/README.md).
The primary conditional rank-change test is nonsignificant (P = 0.14558).
The AI transcriptome benchmark remains unestablished. The record below
preserves the initial audit; follow-up details are in `followup_status.md`.

**Initial decision: the computational extension is in metadata feasibility. No new
biological result, trained model, final dataset selection, or frozen statistical
analysis plan is claimed.** The completed V2 results retain their original
interpretation, including the null independent human replication.

## Scientific objective

Determine whether the proteasome genes implicated by the existing PD analyses
have context-dependent functional effects in public neuronal perturbation
experiments, and, if suitable data exist, test how accurately a perturbation
model predicts their measured transcriptional consequences.

The 52-gene Reactome Proteasome assembly set and its existing eight shared
leading-edge genes are copied from committed V2 inputs into
`candidate_manifest.json`. The eight genes are POMP, PSMA4, PSMB4, PSMC5,
PSMD1, PSMD2, PSMD4, and PSMG1. This selection derives from earlier observed
results; it is not a de novo prospective discovery. No targets will be substituted
because their new experimental outcomes look more favorable.

Reanalyzing public laboratory measurements supplies experimental evidence in
their original cell models. It does not constitute an experiment performed by
the student, establish a causal mechanism in human PD, or show that an AI model
can replace a wet lab. Predicting transcriptomes and validating neuronal
survival are different endpoints and must be evaluated separately.

## What has been checked

| Item | Observed status | Consequence |
|---|---|---|
| Available local compute | Cgroup quota of 8 CPU cores, 14 GiB memory limit, about 31 GB free disk; no exposed NVIDIA device or nvidia-smi | Metadata, streamed preprocessing, pseudobulk statistics and modest baselines are feasible; GPU model training is not currently available here |
| GSE124703 | Author-deposited design includes iPSCs and day-7 neurons, with two 10X lanes for each CROP-seq sample | Potential within-study context comparison; lanes must not be counted as independent biological replicates |
| GSE124703 target mapping | Two neuronal sgRNA mapping files downloaded; they contain guide sequences and unprocessed labels, not a ready gene-symbol assignment | Requires the authors' guide-sequence-to-gene key and assignment rules before target coverage can be established |
| GSE152988 | Author-deposited design specifies 184 CRISPRi and 100 CRISPRa targets; focused CROP-seq targets were selected using preceding screen hits and neurodegeneration associations | This is a selected target panel, not an unbiased genome-wide transcriptome screen |
| GSE152988 CRISPRi assignments | Only PSMF1 from the fixed pathway is represented; none of the shared eight genes | Does not support a direct transcriptomic benchmark of the shared eight |
| GSE152988 CRISPRa assignments | No fixed pathway or shared-eight targets represented | Cannot directly validate activation predictions for these genes |
| GEARS documentation | Maintainers explicitly exclude cross-cell-type transfer and warn about sparse perturbation coverage | Do not train on cancer cells and describe neuronal predictions as supported GEARS use |
| State documentation | Provides genetic perturbation training and explicit context/perturbation splits | Candidate only; usable data, task fit, checkpoint provenance, compute and baselines still need checking |

GSE152988 coverage was computed from the author's six mapping files without
opening expression matrices or screen effect estimates. The files contain
99,433 CRISPRi and 32,129 CRISPRa cell assignments before any expression QC;
481 CRISPRi assignments are labeled PSMF1. These are assignment counts, not
independent experiments or expression-qualified sample sizes. Coverage refers
to targets actually represented in these mappings, not proof of absence from
every original guide library. Exact symbol matching was used after removing
the authors' `_iN` or `_aN` guide suffixes; no alias substitutions were made.

This rules out using GSE152988 as a direct benchmark of perturbing the shared
eight. It does not assess the separate genome-wide survival screen's coverage.
Switching the primary target to PSMF1 solely because it is available would
change the scientific question and is not done here.

## Next decisions before outcome access

1. Resolve GSE124703 guide identities, usable target coverage, perturbation
   contexts, control groups and biological replication. Reconcile original
   GEO metadata with harmonized data; some catalog descriptions mix the two
   Tian study accessions.
2. Audit the original genome-wide neuronal survival-screen library and its
   matched comparison contexts. Distinguish neuronal sensitivity from generic
   essentiality; a significant neuronal effect alone is not neuron specificity.
3. Require actual coverage of a defensible fixed candidate set. A model cannot
   compensate for missing experimental validation targets. If coverage is
   inadequate, report that and decide whether the survival-screen analysis
   alone provides a worthwhile extension.
4. If the prediction task is feasible, specify one primary endpoint and contrast,
   control selection, exclusions, multiplicity handling, model-selection rules,
   train/validation/test split, and failure criteria in a separate committed
   protocol before reading new target outcomes.
5. Audit model pretraining and benchmark exposure. GSE124703 appears among the
   GEARS paper's datasets; reusing it is not a novel independent model benchmark
   by default. Holding data out of our fine-tuning does not establish absence
   from a pretrained checkpoint's training set. Unknown overlap stays unknown.

The feasibility search should be short and based on design and target metadata,
not iterative inspection of candidate p values. Training a foundation model
from scratch is outside the proposed one-month scope.

## Evaluation requirements if a model proceeds

- Include no-change, average-training-perturbation and a suitable simple linear
  baseline. Select model settings using training/validation data only.
- Hold out entire perturbations or experimental contexts as appropriate to the
  scientific claim. A random split of cells from the same perturbation is not
  an unseen-perturbation test.
- Evaluate prediction errors in perturbation effects relative to controls.
  Correlation of absolute expression across all genes can be dominated by
  baseline expression and is insufficient.
- Estimate uncertainty at the available biological experimental unit. Cells,
  guides and sequencing lanes must not be silently treated as independent donors.
- Define neuron selectivity using a direct context comparison or interaction,
  with comparable controls and scales. Different significance labels across
  contexts do not demonstrate a context difference.
- Preserve null results and model failures. This extension cannot retroactively
  change the frozen human or toxicant validation labels.

## Division of work and time

The assistant can implement acquisition and checksums, metadata QC, preprocessing,
baselines, statistical analysis, model configuration, GPU execution notebooks,
evaluation, figures and repository documentation. Actual local execution is
limited by the resources above. No external GPU account, paid compute, model
checkpoint access, or completed training is assumed.

The student's contribution must include understanding and choosing the scientific
question, checking the assumptions and interpretation, and personally owning the
submission and its account of contributions. If a GPU becomes necessary, a
concrete notebook and resource estimate can be prepared before requesting access.

Planning budget, not a completion guarantee: 2–3 days for feasibility and protocol;
about one week for experimental-data analysis and baselines; at most one further
week for a suitable model and held-out evaluation; reserve the remaining time
for audits, figures and the submission. Do not spend the final week debugging a
new model at the expense of a defensible report.

## Reproducing the completed coverage audit

Python standard library only. Download the authors' `maps.zip` from the URL in
`GSE152988_target_coverage.json`; the script requires its exact SHA256:

```text
516d0b462e340bce4a08f0a5796b52c46d00f2af15cd0742dc1f857c22c20c45
```

From the repository root:

```bash
python scripts/audit_neuronal_perturbation_metadata.py /path/to/maps.zip
```

Expected coverage is CRISPRa 0/52 and 0/8, CRISPRi 1/52 and 0/8. The output
records checksums for each mapping file. Series metadata snapshots retain only
design, accession, publication and supplementary-file fields. The coverage
script does not open expression values or laboratory phenotypes.

## Sources

- Original GSE124703: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124703
- Original GSE152988: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152988
- Author's GSE152988 CROP-seq mapping download: https://kampmannlab.ucsf.edu/crop-seq
- Tian 2019 paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC6813890/
- Tian 2021 paper: https://www.nature.com/articles/s41593-021-00862-0
- GEARS maintained usage limitations: https://github.com/snap-stanford/GEARS
- GEARS original model paper and dataset list: https://www.nature.com/articles/s41587-023-01905-6
- State maintained implementation: https://github.com/ArcInstitute/state

Public methods, study descriptions and the existing PD project results have
been consulted. The new target-specific expression, survival, and model-test
outcomes have not been analyzed in this audit. This is a local feasibility
record, not a registration with an external preregistration service.
