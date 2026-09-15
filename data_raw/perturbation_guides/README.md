# Original public guide key

`NIHMS1535621-supplement-8.xlsx` is the user-supplied copy of Table S2 from
Tian et al., *Neuron* (2019), DOI: 10.1016/j.neuron.2019.07.014.
Original paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC6813890/.

The workbook is preserved byte-for-byte; SHA256:
`149ec75dda0adaae6238c00d3ac5ecd7a308a4d910e4f23359247eb0372e1247`.
It was supplied on 2026-09-15 after public attachment downloads failed here.
The guide-coverage audit uses only its first three identity columns.
The remaining survival-phenotype columns were not analyzed.

Reproduce after downloading the two GEO mapping files at the URLs in the
coverage JSON:

```bash
python scripts/audit_tian2019_guide_key.py /path/to/mapping_files
```

This uses `openpyxl==3.1.5`, already pinned in
`requirements-neuronal-survival.txt`. Derived identities and coverage are in
`docs/v2/perturbation_feasibility/`.
