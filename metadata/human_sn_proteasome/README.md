# GSE178265 frozen DA-cell metadata

`GSE178265_DA_cell_metadata.tsv.gz` contains the 22,048 author-annotated
dopaminergic nuclei used by the frozen GSE178265 analysis. It was created and
checksum-locked on 2026-08-13 before any expression value for the target
Proteasome assembly pathway was requested.

Columns are the exact SCP1768 cell identifier, author DA subtype, donor, disease
status, sex, age, postmortem interval, and per-cell total UMI count. The row
order is the public `UMAP: Human DA Neurons` cluster order and is checked against
SCP1768 whenever inputs are fetched.

Sources:

- Cell identifiers, subtype, donor, status, sex, age, and PMI: public SCP1768
  cluster and annotation API responses retrieved 2026-08-13.
- `total_umi`: `nUMI` from
  `Analyses/IntermediateObjects/individualobject_seuratmetadataframe.tar.gz`
  in the authors' official `tkamath1/Kamathetal2022` repository at commit
  `cae72e19290efe1ea8756e2032507905bdea52ba`.
- Official metadata archive SHA256:
  `97b7cdd379507e279ba7168dc32ea6d4e78271a820bcf4348dd36bcf2a5fed55`.
- Extracted `individualobject_seuratmetadataframe.qs` SHA256:
  `c4a7ece6a277574c93109c28f20a60a59573e4589fd9058cc590e4a56041d16f`.
- This tracked gzip file SHA256:
  `6732fb4323cc36bad0a23bb53391d3242a2d68a1ebef5dd3d36857a72c2ab916`.

The authors' `Pos` libraries correspond to `Nurr` in the intermediate object.
The following additional naming harmonizations were required and verified by
exact barcode matching: `5610Neg` to `5610DAPI`, `6173DAPIA` to `6173DAPI`, and
`2569PosA` to `2569Nurr`. All 22,048 cells matched uniquely. This file contains
metadata and library-size values only, with no target-gene expression.
