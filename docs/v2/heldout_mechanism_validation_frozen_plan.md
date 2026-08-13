# Frozen held-out mechanism-validation plan

Frozen on 2026-08-13 before selecting an accession, downloading an expression
matrix, or inspecting outcomes from any new dataset.

## Objective

Test whether the strongest mechanism nominated by the completed exploratory
analyses, downregulation of cilium/Hedgehog/microtubule trafficking machinery,
generalizes to an independent PD-relevant neuronal toxicant dataset. Aggrephagy
is a gated secondary hypothesis.

Existing datasets used anywhere in discovery, auditing, or model development are
ineligible for the held-out test: GSE17542, GSE187012, GSE46798, GSE116280,
GSE150005, GSE20141, and GSE7621.

## Dataset eligibility

An eligible dataset must satisfy all of the following using metadata alone:

- publicly downloadable transcriptomic data
- human neuronal culture with dopaminergic differentiation or a clearly stated
  midbrain/dopaminergic neuronal identity
- exposure to rotenone, MPP+, paraquat, maneb, ziram, or
  6-hydroxydopamine
- a matched untreated or vehicle control
- at least three biological replicates in both the selected exposure and control
  groups
- gene-level raw counts or a documented normalized expression matrix with usable
  gene annotation
- samples independent of every discovery dataset listed above

The following are excluded before outcome inspection:

- non-neuronal or glial-only cultures
- genetic perturbation without a toxicant contrast
- fewer than three biological replicates per group
- pooled samples for which biological replication cannot be reconstructed
- single-cell data without enough independent biological samples for pseudobulk
  inference
- exposure and control groups confounded with donor, batch, sequencing run, or
  differentiation protocol
- the same biological samples deposited under a second accession

## Dataset search and selection rule

Record every accession screened and the metadata reason for inclusion or
exclusion. Do not inspect differential-expression tables, pathway results, or the
expression matrix before the primary accession and contrast are locked.

If several datasets are eligible, choose the primary dataset using this ordered
metadata-only priority:

1. human iPSC-derived midbrain dopaminergic neurons or mature LUHMES neurons
2. raw RNA-seq counts over microarray or preprocessed-only data
3. rotenone or MPP+ exposure over the other eligible toxicants
4. larger number of independent biological replicates in the primary contrast
5. sustained lower-dose exposure over an explicitly lethal or near-lethal dose
6. greater stated neuronal maturation at exposure
7. lexicographically smaller GEO accession as the final deterministic tie-break

All additional eligible datasets are analyzed as declared external sensitivity
datasets. A dataset may be replaced only for a documented metadata or technical
eligibility failure identified without looking at the mechanism outcomes.

## Frozen pathway definitions

All gene membership is fixed by the tracked
`data_raw/pathways/ReactomePathways.gmt` file with SHA256
`89983d5c1f0af11c52edfeee7323eb425580ac6281d387a528562ab1787ce56b`.
Gene symbols are converted to uppercase and intersected with the measured genes.

### Primary mechanism

The primary `cilium_Hedgehog_microtubule_trafficking` set is the unique union of
the following five Reactome pathways:

1. `Microtubule-dependent trafficking of connexons from Golgi to the plasma membrane`
2. `Intraflagellar transport`
3. `Cargo trafficking to the periciliary membrane`
4. `Carboxyterminal post-translational modifications of tubulin`
5. `Hedgehog 'off' state`

The declared direction is downregulated in toxicant-exposed neurons relative to
matched controls.

### Gated secondary mechanism

The secondary set is the exact Reactome pathway `Aggrephagy`, with a declared
downregulated direction.

Proteasome assembly, ribosome-quality control, individual genes, and every other
Reactome pathway are descriptive controls or exploratory results. They cannot
replace either frozen mechanism after outcome inspection.

## Dataset-specific design lock

After the accession is selected, sample metadata may be inspected to freeze the
sample IDs, primary contrast, and covariates before expression values are read.
The design must reflect donor, batch, pairing, or repeated measures when present.
Any sample exclusion requires a metadata-based reason recorded before outcome
analysis. No expression-based outlier removal is allowed in the primary result.

## Preprocessing

For RNA-seq counts:

1. map to one uppercase HGNC symbol per row and sum duplicated symbols
2. use edgeR `filterByExpr` with the locked design
3. apply TMM normalization
4. fit voom/limma with the locked design and contrast
5. use robust empirical Bayes moderation

For an eligible microarray dataset:

1. use the deposited documented normalized matrix
2. exclude ambiguous gene mappings
3. select one probe per gene by the greatest all-sample mean expression, with
   probe identifier as a deterministic tie-break
4. fit limma with the locked design and contrast

## Primary and secondary inference

Use limma CAMERA with pathway-specific intergene correlation
(`inter.gene.cor = NA`). The hypothesis direction is fixed as down, so the
one-sided-down p-value is half the two-sided p-value when CAMERA reports `Down`
and `1 - p/2` when it reports `Up`.

The primary mechanism is confirmed only if all of the following hold in the
locked primary dataset:

- the union set contains at least 10 measured genes
- CAMERA direction is `Down`
- one-sided-down p < 0.05
- at least three of the five constituent pathways have negative mean moderated
  t statistics

The secondary aggrephagy hypothesis is tested only after the primary mechanism
passes. It is confirmed if it contains at least 10 measured genes, CAMERA reports
`Down`, and its one-sided-down p < 0.05. This hierarchical gate keeps the total
confirmatory family-wise type-I error at 0.05.

Constituent-pathway CAMERA results are reported with Benjamini-Hochberg correction
across the five pathways, but they do not override the union-set primary result.
No individual-gene FDR hit is required for pathway confirmation.

## Sensitivity analyses

The following are declared secondary and cannot rescue a failed primary result:

- two-sided CAMERA p-values
- constituent-pathway FDR results
- leave-one-pathway-out versions of the union
- a simple count filter requiring at least 10 counts in at least three samples
  instead of `filterByExpr`
- alternative valid handling of donor or batch supported by metadata
- all additional eligible datasets found by the frozen search rule

## Outcome labels

- `CONFIRMED`: every primary criterion passes.
- `PARTIAL_SUPPORT`: the primary direction is down with one-sided p from 0.05 to
  less than 0.10, or the union passes while fewer than three constituent pathways
  have negative mean moderated t statistics.
- `NOT_CONFIRMED`: every other evaluable outcome, including an opposite direction.
- `NOT_EVALUABLE`: the locked dataset fails a technical eligibility condition
  identified without reference to the pathway result.

The accession search log, metadata lock, code, complete results, and failures will
be committed regardless of outcome.

Confirmation means that the frozen gene set shows a reproducible coordinated
transcriptional shift. It does not by itself establish altered Hedgehog signaling
activity, ciliary function, protein flux, or neuronal survival; those functional
claims would require experiments outside these transcriptomic datasets.
