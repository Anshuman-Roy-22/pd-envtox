#!/usr/bin/env python3
"""Verify the extension snapshot and recompute its reported prediction errors."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/v2/fly_genetic_rescue'


def main():
    manifest = json.loads((OUT / 'release_checksums.json').read_text())
    for relative, expected in manifest['sha256'].items():
        path = Path(relative)
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('Unsafe release path')
        actual = hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f'Changed release file: {relative}')
    summary = json.loads((OUT / 'validation_summary.json').read_text())
    pred = pd.read_csv(OUT / 'heldout_predictions.tsv', sep='\t')
    y = pred.relative_response_Q10_minus_R55.to_numpy()
    errors = {}
    for model in ['GI', 'mean', 'severity']:
        errors[model] = np.mean(np.abs(pred[f'prediction_{model}'].to_numpy() - y))
        np.testing.assert_allclose(errors[model], summary[f'MAE_{model}'], rtol=0, atol=1e-10)
    assert summary['classification'] == 'NOT_SUPPORTED'
    assert len(pred) == 18 and pred.locus.nunique() == 18
    assert errors['GI'] > errors['mean']
    null = pd.read_csv(OUT / 'permutation_statistics.tsv', sep='\t')
    count = int((null['T'] >= summary['T']).sum())
    np.testing.assert_allclose((1 + count) / 10001, summary['randomization_p'], rtol=0, atol=1e-12)
    neighbors = pd.read_csv(OUT / 'prediction_neighbors.tsv', sep='\t')
    assert (neighbors.held_out_locus != neighbors.training_locus).all()
    assert neighbors.groupby(['held_out_locus', 'model']).size().eq(3).all()
    labels = pred.set_index('locus').relative_response_Q10_minus_R55
    for (heldout, model), group in neighbors.groupby(['held_out_locus', 'model']):
        expected = labels.loc[group.training_locus].mean()
        actual = pred.loc[pred.locus == heldout, f'prediction_{model}'].iloc[0]
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-10)
    print(f"PASS: {len(manifest['sha256'])} release checksums, genotype exclusion, neighbor predictions and primary statistics.")


if __name__ == '__main__':
    main()
