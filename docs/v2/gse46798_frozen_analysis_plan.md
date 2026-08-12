# Frozen GSE46798 analysis plan

Frozen before expression outcomes were inspected beyond confirming matrix scale.

## Design

GSE46798 contains 12 GPL10558 profiles from human iPSC-derived dopaminergic
neurons, with three replicates in each cell of a balanced design:

- genetically corrected, vehicle
- SNCA A53T, vehicle
- genetically corrected, paraquat/maneb
- SNCA A53T, paraquat/maneb

The linear model is `expression ~ genotype + exposure + genotype:exposure`.

## Prespecified contrasts

1. Primary: paraquat/maneb effect in genetically corrected neurons.
2. Secondary: paraquat/maneb effect in A53T neurons.
3. Secondary: genotype-by-exposure interaction.
4. Context: A53T effect under vehicle.

Positive log fold change means higher expression in the first named condition.

## Frozen candidate panel and modules

The boundary-inclusive 21-gene network-v2 panel was fixed without reference to
GSE46798 outcomes.

- Retromer/endosomal: SNX1, SNX2, SNX3, VPS26A, VPS26B, VPS29
- Mitochondrial dynamics/transport: MFN1, MFN2, RHOT1, TOMM20
- Endocytosis/vesicle trafficking: CLTC, CLTCL1, EPS15, ITSN1, SH3GL2
- Ubiquitin/proteostasis: UBB, UBC, UBA52, HSPA8

SNCAIP and WLS remain in the full panel but are not assigned post hoc to a module.

## Processing and statistics

- Apply log2(x + 1) because the matrix is on a linear intensity scale.
- Use the August 2016 GEO GPL10558 annotation.
- Exclude probes mapped to multiple gene symbols.
- Choose one probe per gene by greatest mean expression across all 12 samples,
  breaking ties by probe ID. This is outcome-blind.
- Fit ordinary least-squares factorial models in the recovery runtime. Benjamini-
  Hochberg FDR is calculated genome-wide for each contrast.
- The final competition pipeline must confirm results with limma empirical-Bayes
  moderation.
- Panel and module tests use 10,000 fixed-seed expression-matched gene-set draws.
  Primary module metric is mean absolute contrast t; signed mean t is secondary.
- No module or candidate will be selected based on significance in this dataset.
