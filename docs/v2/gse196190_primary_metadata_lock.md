# GSE196190 primary held-out metadata lock

Locked: 2026-08-13 EDT, before opening any expression-value file.

Parent protocol: `docs/v2/heldout_mechanism_validation_frozen_plan.md`.

## Metadata-only selection

GSE196190 is the primary held-out dataset under the frozen priority rules. It is public human RNA-seq from a named iPSC line differentiated for 54 days into dopaminergic neurons, includes MPP+ and matched vehicle groups with three biological samples per dose-time cell, and supplies gene-level RSEM output plus SRA raw reads. It outranks the other eligible candidates because it is an iPSC-derived dopaminergic model, uses RNA-seq, and tests MPP+.

The candidate search was run through NCBI GEO DataSets with this metadata query:

`(rotenone OR MPP OR paraquat OR maneb OR ziram OR 6-hydroxydopamine) AND Homo sapiens[Organism] AND gse[Entry Type] AND (neuron OR neuronal OR dopaminergic OR LUHMES OR midbrain)`

All 47 returned accessions and their eligibility decisions are recorded in `docs/v2/heldout_dataset_search_log.tsv`. No target-gene, target-pathway, differential-expression, or expression-matrix values were inspected during selection.

## Locked primary contrast

Treatment: 100 micromolar MPP+ at GEO time label `t74`:

- GSM5862317
- GSM5862318
- GSM5862319

Matched control at GEO time label `t74`:

- GSM5862314
- GSM5862315
- GSM5862316

The primary design is the two-group contrast `MPP_100uM_t74 - control_t74`, with three biological samples per group and no metadata-supported covariate. The submitter's paper describes the sustained endpoint as 72 hours, whereas the GEO sample titles and filenames label it `t74`; the six exact GSM accessions above govern the analysis. The lowest MPP+ dose at the sustained endpoint was selected before expression inspection, following the frozen preference for a sustained lower-dose exposure over a lethal or near-lethal exposure.

The remaining GSE196190 dose-time cells are excluded from the primary contrast. They may be reported only as clearly labeled within-study sensitivity analyses after the primary result is fixed.

## Locked processing and decision rule

The per-sample RSEM gene-level `expected_count` field will be used. Gene identifiers will be mapped to current HGNC symbols with the mapping source, retrieval date, and checksum recorded. Duplicate mapped symbols will be summed before filtering. The RNA-seq workflow remains the frozen edgeR `filterByExpr` plus TMM normalization and limma-voom gene model.

The primary pathway set, direction, CAMERA settings, minimum union size, and success labels are unchanged from the parent frozen protocol. Aggrephagy remains gated and will be tested only if the frozen primary cilium/Hedgehog/microtubule-trafficking union criterion passes.

## Declared external sensitivity datasets

- GSE229460: mature LUHMES, subtoxic MPP+, three control and three MPP+ samples, RNA-seq.
- GSE287941: differentiated wild-type LUHMES, untreated and 6-OHDA groups, four samples per group, RNA-seq.
- GSE4773: chronic rotenone in SK-N-MC cells, sensitivity-only and explicitly lower confidence because of tumor-cell identity and composite pooling.

These sensitivity datasets cannot replace or rescue the locked primary outcome.

## Sources consulted for metadata

- NCBI GEO accession GSE196190 and its sample records.
- GSE196190 supplementary file listing, without opening expression-value files.
- PMID 35725899 / DOI 10.1038/s41380-022-01663-y for culture maturity and the paper's 72-hour endpoint description.
