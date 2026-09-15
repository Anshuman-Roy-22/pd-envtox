# Why this next approach was selected

The previous functional-rescue benchmark was `BENCHMARK_NOT_EVALUABLE`: the calibration rule did not transfer. That is preserved. The next question should use identifiable experimental interventions and a neuronal-function connection, rather than add a structural prediction without a specific experimental test.

## Source audit, before numerical analysis

| Candidate | Verified public input | What it permits | Important constraint |
|---|---|---|---|
| Zhou et al., Cell Reports 2020, GSE154112 | 7,569 ordered guide pairs, 28 genes plus dummy controls, untreated/MPP+/rotenone columns | Comparison of double-gene and single-gene perturbations | Two biological replicates; CPM rather than raw counts; tumor-derived SK-N-MC cells; no targets from the original 52-gene set |
| Kaempf et al., Nature Communications 2026 | 24-locus interaction matrix; both drugs and corresponding vehicles for 18 mutant genotypes | Predict intervention response from an experimentally measured genetic-interaction profile | One fly study; author-processed interaction scores; missing vial/batch identifiers and different allele dosages |

Select the second candidate for the fixed benchmark in `frozen_plan.md`. Do not run both and select the favorable result. The first remains an audited alternative, not a new result. The published HSP90B1–HDAC2 combination in the first study and the subgroup/drug findings in the second are prior work, not discoveries by this project.

## Source provenance

- Zhou et al.: https://doi.org/10.1016/j.celrep.2020.108020 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154112
- ATCC identifies SK-N-MC HTB-10 as a neuroepithelioma cell line: https://www.atcc.org/products/htb-10
- Kaempf et al.: https://doi.org/10.1038/s41467-026-70303-8 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC13106710/
- Author code snapshot: https://github.com/verstrekenlab/drosophila-Parkinsonism-subgroups/tree/82f71c1adc356d1819f38e32f97dc0d42afce091

The author repository's `run.py` expects `raw/20250828_Raw_data.csv` and `raw/240627_list_gene_pairs.csv`, which were not present in the inspected tree. Its `models` file contains a local filesystem path. Consequently, we cannot claim to reproduce raw pairwise-electrophysiology fitting from the verified public files. The fixed published interaction matrix is the predictor.

The treatment workbook has numeric observations and blank cells without vial, animal or experiment identifiers. Workbook row position is not a pairing key. Source sheet labels use `Suppl Fig12` for the climbing data although the final article's corresponding supplementary numbering differs; scripts bind to the actual sheet names.

At the freeze point, the two workbooks had been downloaded, checksum-verified and inspected for sheets, labels, coverage and observation counts. Numeric effect estimates had not been calculated. Published results and the author's categorical treatment summary were already visible. This is a retrospective benchmark with outcome-withheld cross-validation, not a prospectively unseen dataset.

## What this could add

Successful prediction would connect measured genetic interactions to differential intervention response beyond baseline severity. That is stronger functional evidence than pathway enrichment alone. It would still leave the specific molecular mediator and relevance to human PD unresolved. If the fixed test fails, changing models until it passes would not solve that evidence gap.
