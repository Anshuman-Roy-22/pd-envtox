# Reproducibility Metadata

Store frozen metadata here so a release can be validated without re-running every step.

## Recommended Files

- `sessionInfo.txt`: captured from the release environment
- `raw_checksums.csv`: SHA256 hashes for committed raw inputs
- `results_checksums.csv`: SHA256 hashes for committed release outputs

## Capture Session Info

Run from R at the repository root:

```r
dir.create("metadata", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "metadata/sessionInfo.txt")
```

## Capture Checksums

Run from PowerShell at the repository root:

```powershell
New-Item -ItemType Directory -Force metadata
Get-ChildItem data_raw -Recurse -File |
  Get-FileHash -Algorithm SHA256 |
  Export-Csv metadata\raw_checksums.csv -NoTypeInformation

Get-ChildItem results -Recurse -File |
  Get-FileHash -Algorithm SHA256 |
  Export-Csv metadata\results_checksums.csv -NoTypeInformation
```

Regenerate these files for each tagged release.

