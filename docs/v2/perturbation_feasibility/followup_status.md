# Feasibility follow-up, 2026-09-14

## Completed survival branch

The original 2021 publisher workbook and guide-library supplement were obtained.
Only gene identities and workbook structure were examined for the initial
coverage check. `Tian2021_screen_library_coverage.json` records coverage and
checksums. The CRISPRi survival screens include 49 of the fixed 52 pathway
genes; the CRISPRa screens include 51. All eight shared leading-edge genes
occur in both modes. The focused CROP-seq libraries independently verify the
previous finding: zero of the shared eight occur in either transcriptomic
perturbation library.

The survival analysis was frozen at commit `6d493a3` and completed without
changing targets or thresholds. Its primary P value is 0.14558. Both binning
sensitivities and the activation secondary are also nonsignificant. The full
pathway secondary is unevaluable under its fixed control-count requirement.
See the complete report in `results/v2/neuronal_survival/README.md`.

## Older transcriptomic dataset: remaining mapping dependency

Original GEO accession: GSE124703, Tian et al. (2019). The two downloaded
neuronal sgRNA-enrichment files contain cell barcodes, guide sequences,
read counts and UMI counts. Sequence-to-gene identities are not supplied as
gene symbols in those files. Unprocessed sequence labels are also present.

The harmonizer's public processing notebook explicitly reads an additional
file called `mmc3.csv`, with columns `protospacer sequence` and
`sgRNA_short name`. It maps the original guide sequences using that file.
This verifies that the gene key is a real missing dependency, not something
that should be guessed from gene expression or inferred from the candidates.

Original study Table S2 contains short and long sgRNA names, protospacer
sequences, and subsequent phenotype columns. Its PMC attachment is named
`NIHMS1535621-supplement-8.xlsx`. Access attempts to the attachment and public
publisher/mirror endpoints returned challenge pages or HTTP errors here.
Two attempts at the harmonized neuronal H5AD download also did not yield a
usable file at the time of the initial follow-up. These are access limitations,
not evidence that the data do not exist.

Only row/column perturbation labels were extracted from the harmonizer's public
pairwise-distance tables; distance values were not inspected. Their retained
labels do not contain a fixed pathway target, including after splitting
multi-target labels. These are filtered analysis tables and **do not establish
the complete original library's absence of those genes**. They cannot substitute
for the missing guide key. Sources and source-file hashes are recorded in
`Tian2019_harmonizer_target_labels.json`.

The exact small file needed to finish the original-library check is the authors'
Table S2 workbook, or a faithful CSV export preserving its guide-name and
protospacer columns. The phenotype columns are not needed for that check.
Once available, the key can be joined to the already downloaded neuronal guide
files to establish represented candidate targets. No gene identities will be
invented and no model will be presented as validated on missing targets.

## Model decision

- State/GEARS training has not started. The binding issue is experimental task
  and target coverage; obtaining a GPU alone does not solve it.
- GEARS' maintained documentation excludes cross-cell-type transfer, and its
  paper already used GSE124703. Reusing that dataset would require an explicit
  distinction between a model benchmark and new disease-specific evidence.
- The survival tables contain actual laboratory measurements, but they are not
  transcriptomes and cannot directly score a model's predicted RNA vector.
- The completed secondary survival analysis supplies a bounded experimental-data
  check. It does not establish the stronger result sought for the competition.

## Sources

- Original 2019 paper and Table S2 description:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6813890/
- Harmonizer processing code, code cells only inspected:
  https://github.com/sanderlab/scPerturb/blob/master/dataset_processing/notebooks/data_processing_JS.ipynb
- Harmonized data record:
  https://zenodo.org/records/13350497
- 2021 original paper:
  https://www.nature.com/articles/s41593-021-00862-0
- 2021 guide-library supplement:
  https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41593-021-00862-0/MediaObjects/41593_2021_862_MOESM6_ESM.xlsx

Processing notebook SHA256 (retrieved 2026-09-14):
`09837dd5ce1cc98648a9153ca337598e940262a51b475be92141e5f8ec61ae16`.
