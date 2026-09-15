#!/usr/bin/env python3
"""Frozen, retrospective genotype-withheld prediction of intervention response."""
import argparse
import importlib.metadata
import json
import math
import platform
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy.stats import spearmanr

from acquire_fly_genetic_rescue import DOCS, ROOT, acquire

FREEZE = '5d6968a49890cf724566c762aae2a9299f878efe'
SEED = 20260915
N_PERM = 10000
DRUGS = [('Q10', 'EtOH'), ('R55', 'H2O')]


class NotEvaluable(ValueError):
    pass


def clean_json(value):
    if isinstance(value, dict):
        return {k: clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(v) for v in value]
    if isinstance(value, np.generic):
        return clean_json(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path, value):
    path.write_text(json.dumps(clean_json(value), indent=2, ensure_ascii=False,
                               allow_nan=False) + '\n')


def table(frame, path):
    frame.to_csv(path, sep='\t', index=False, float_format='%.12g',
                 na_rep='NA', lineterminator='\n')


def arm_aggregates(wb):
    records = []
    sheetmap = {'TH': {'Q10': 'Fig4b_Q10', 'R55': 'Fig4d_R55'},
                'SING': {'Q10': 'Suppl Fig12b_Q10_SING', 'R55': 'Suppl Fig12c_R55_SING'}}
    for endpoint, sheets in sheetmap.items():
        for drug, title in sheets.items():
            rows = list(wb[title].values)
            for j, header in enumerate(rows[0]):
                if header is None:
                    continue
                genotype, arm = header.rsplit('_', 1)
                raw = [row[j] for row in rows[1:] if row[j] is not None]
                if any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in raw):
                    raise NotEvaluable(f'Non-numeric source value: {title}/{header}')
                values = np.asarray(raw, dtype=float)
                if len(values) < 3 or not np.isfinite(values).all():
                    raise NotEvaluable(f'Insufficient/invalid observations: {title}/{header}')
                records.append(dict(endpoint=endpoint, sheet=title, drug=drug,
                                    genotype=genotype, arm=arm, source_observations=len(values),
                                    mean=float(values.mean()), median=float(np.median(values))))
    return pd.DataFrame(records)


def effects(arms, mapping, endpoint='TH', aggregate='mean'):
    subset = arms[arms.endpoint == endpoint]
    lookup = subset.set_index(['drug', 'genotype', 'arm'])[aggregate]
    if lookup.index.has_duplicates:
        raise NotEvaluable('Duplicate source arm')
    rows = []
    for item in mapping.to_dict('records'):
        genotype = item['drug_genotype']
        if endpoint == 'SING' and genotype not in set(subset.genotype):
            continue
        record = dict(item)
        for drug, vehicle in DRUGS:
            keys = [(drug, genotype, drug), (drug, genotype, vehicle),
                    (drug, 'control', drug), (drug, 'control', vehicle)]
            if any(key not in lookup for key in keys):
                raise NotEvaluable(f'Missing arm for {genotype}, {drug}, {endpoint}')
            treated, untreated, wt_t, wt_v = [float(lookup.loc[key]) for key in keys]
            if min(treated, untreated, wt_t, wt_v) <= 0:
                raise NotEvaluable(f'Nonpositive {aggregate}: {genotype}, {drug}, {endpoint}')
            record[f'{drug}_raw_log2_ratio'] = np.log2(treated / untreated)
            record[f'{drug}_control_log2_ratio'] = np.log2(wt_t / wt_v)
            record[f'{drug}_adjusted_log2_response'] = (np.log2(treated / untreated)
                                                       - np.log2(wt_t / wt_v))
            record[f'{drug}_vehicle_log2_severity'] = np.log2(untreated / wt_v)
        record['relative_response_Q10_minus_R55'] = (
            record['Q10_adjusted_log2_response'] - record['R55_adjusted_log2_response'])
        record['both_adjusted_responses_positive'] = (
            record['Q10_adjusted_log2_response'] > 0 and record['R55_adjusted_log2_response'] > 0)
        record['both_adjusted_responses_negative'] = (
            record['Q10_adjusted_log2_response'] < 0 and record['R55_adjusted_log2_response'] < 0)
        rows.append(record)
    return pd.DataFrame(rows)


def interaction_matrix(wb, mapping):
    rows = list(wb['Fig3d_GI_cluster'].values)
    labels = list(rows[0][1:])
    index = [row[0] for row in rows[1:]]
    if len(set(labels)) != 24 or len(set(index)) != 24 or set(labels) != set(index):
        raise NotEvaluable('Interaction matrix labels are not a unique 24-locus set')
    matrix = pd.DataFrame([row[1:] for row in rows[1:]], columns=labels, index=index)
    try:
        matrix = matrix.astype(float).reindex(index=labels, columns=labels)
    except (ValueError, TypeError) as exc:
        raise NotEvaluable('Interaction matrix is not numeric') from exc
    if not set(mapping.interaction_genotype).issubset(matrix.index):
        raise NotEvaluable('Missing genotype mapping in interaction matrix')
    a = matrix.to_numpy()
    comparable = np.isfinite(a) & np.isfinite(a.T)
    audit = dict(loci=24, nonfinite_entries=int((~np.isfinite(a)).sum()),
                 max_absolute_asymmetry=float(np.max(np.abs(a - a.T)[comparable])),
                 diagonal_used=False, source='author-processed Fig3d_GI_cluster')
    return matrix, audit


def gi_distances(matrix, data):
    names = list(data.interaction_genotype)
    n = len(names)
    distances = np.full((n, n), np.inf)
    anchors = np.zeros((n, n), dtype=int)
    for i, name_i in enumerate(names):
        for j, name_j in enumerate(names):
            if i == j:
                continue
            common = [a for a in matrix.columns if a not in {name_i, name_j}
                      and np.isfinite(matrix.loc[name_i, a]) and np.isfinite(matrix.loc[name_j, a])]
            if len(common) < 18:
                raise NotEvaluable(f'Insufficient interaction anchors: {name_i}, {name_j}')
            difference = matrix.loc[name_i, common].to_numpy() - matrix.loc[name_j, common].to_numpy()
            distances[i, j] = np.sqrt(np.mean(difference ** 2))
            anchors[i, j] = len(common)
    return distances, anchors


def neighbor_weights(distances, k):
    n = distances.shape[0]
    if n <= k or not np.isinf(np.diag(distances)).all():
        raise NotEvaluable('Invalid neighbor design')
    weights = np.zeros((n, n))
    # Data are lexicographically ordered. Stable sorting implements the frozen tie rule.
    nearest = np.argsort(distances, axis=1, kind='stable')[:, :k]
    for i, cols in enumerate(nearest):
        if i in cols or not np.isfinite(distances[i, cols]).all():
            raise NotEvaluable('Invalid or incomplete held-out prediction')
        weights[i, cols] = 1 / k
    return weights


def build_weights(matrix, data, k):
    gi, anchors = gi_distances(matrix, data)
    severity = data[['Q10_vehicle_log2_severity', 'R55_vehicle_log2_severity']].to_numpy()
    delta = severity[:, None, :] - severity[None, :, :]
    sd = np.sqrt(np.sum(delta ** 2, axis=2))
    np.fill_diagonal(sd, np.inf)
    n = len(data)
    mean = (np.ones((n, n)) - np.eye(n)) / (n - 1)
    return {'GI': neighbor_weights(gi, k), 'mean': mean,
            'severity': neighbor_weights(sd, k)}, gi, sd, anchors


def evaluate(y, weights):
    predictions = {key: weight @ y for key, weight in weights.items()}
    mae = {key: float(np.mean(np.abs(pred - y))) for key, pred in predictions.items()}
    if mae['mean'] <= 0 or mae['severity'] <= 0:
        raise NotEvaluable('Zero baseline MAE')
    improvement_mean = 1 - mae['GI'] / mae['mean']
    improvement_severity = 1 - mae['GI'] / mae['severity']
    rho = float(spearmanr(y, predictions['GI']).statistic)
    summary = dict(n_genotypes=len(y), MAE_GI=mae['GI'], MAE_mean=mae['mean'],
                   MAE_severity=mae['severity'], improvement_over_mean=improvement_mean,
                   improvement_over_severity=improvement_severity,
                   T=min(improvement_mean, improvement_severity), spearman_rho=rho)
    return summary, predictions


def permute(y, weights, observed):
    rng = np.random.default_rng(SEED)
    shuffled = np.vstack([rng.permutation(y) for _ in range(N_PERM)])
    errors = {key: np.mean(np.abs(shuffled @ w.T - shuffled), axis=1)
              for key, w in weights.items()}
    if (errors['mean'] <= 0).any() or (errors['severity'] <= 0).any():
        raise NotEvaluable('Degenerate permuted baseline')
    statistic = np.minimum(1 - errors['GI'] / errors['mean'],
                           1 - errors['GI'] / errors['severity'])
    exceed = int(np.sum(statistic >= observed))
    return (1 + exceed) / (N_PERM + 1), pd.DataFrame({
        'permutation': np.arange(1, N_PERM + 1), 'T': statistic})


def mathematical_checks(y, weights):
    for name, weight in weights.items():
        assert np.all(np.diag(weight) == 0), name
        np.testing.assert_allclose(weight.sum(axis=1), 1, atol=1e-14)
        prediction = weight @ y
        for i in range(len(y)):
            changed = y.copy()
            changed[i] += 10000
            np.testing.assert_allclose((weight @ changed)[i], prediction[i], rtol=0, atol=1e-12)
            expected = sum(weight[i, j] * y[j] for j in range(len(y)) if j != i)
            np.testing.assert_allclose(prediction[i], expected, rtol=0, atol=1e-12)
    # Known treatment/control ratios and invariance to common genotype scale.
    def response(t, v, wt_t, wt_v):
        return np.log2(t / v) - np.log2(wt_t / wt_v)
    np.testing.assert_allclose(response(4, 2, 2, 2), 1)
    np.testing.assert_allclose(response(4 * 7.9, 2 * 7.9, 2, 2), 1)
    return {'held_out_outcome_invariance': 'PASS', 'direct_loop_prediction_agreement': 'PASS',
            'known_ratio_and_common_scale_cancellation': 'PASS'}


def run(cache, out, fetch):
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / 'source_inventory.json', acquire(cache, fetch))
    mapping = pd.read_csv(DOCS / 'genotype_map.tsv', sep='\t').sort_values('locus').reset_index(drop=True)
    if len(mapping) != 18 or mapping.locus.duplicated().any():
        raise NotEvaluable('Expected 18 uniquely mapped loci')
    wb = openpyxl.load_workbook(cache / 'source_data.xlsx', read_only=True, data_only=True)
    arms = arm_aggregates(wb)
    matrix, audit = interaction_matrix(wb, mapping)
    wb.close()
    data = effects(arms, mapping)
    table(arms, out / 'arm_aggregates.tsv')
    table(data, out / 'genotype_effects.tsv')
    table(matrix.rename_axis('interaction_genotype').reset_index(), out / 'interaction_matrix.tsv')
    sw = openpyxl.load_workbook(cache / 'supplementary_dataset_1.xlsx', read_only=True, data_only=True)
    summary_rows = list(sw['bayesian_statistics'].values)
    audit['supplementary_summary_rows'] = len(summary_rows) - 1
    audit['supplementary_single_labels'] = sorted({str(r[c]) for r in summary_rows[1:] for c in [2, 3]})
    audit['supplementary_summary_role'] = 'label/coverage audit only; no second fitted model'
    sw.close()
    write_json(out / 'interaction_audit.json', audit)
    weights, gi, sd, anchors = build_weights(matrix, data, 3)
    y = data.relative_response_Q10_minus_R55.to_numpy()
    summary, predictions = evaluate(y, weights)
    write_json(out / 'mathematical_checks.json', mathematical_checks(y, weights))
    p, null = permute(y, weights, summary['T'])
    table(null, out / 'permutation_statistics.tsv')
    summary.update(randomization_p=p, n_permutations=N_PERM, seed=SEED, freeze_commit=FREEZE,
                   evidence_scope='retrospective within-study genotype-withheld prediction')
    summary['classification'] = ('WITHIN_STUDY_PREDICTION_SUPPORTED'
                                 if summary['T'] >= .10 and p <= .05 and summary['spearman_rho'] > 0
                                 else 'NOT_SUPPORTED')
    output = data[['locus', 'interaction_genotype', 'drug_genotype',
                   'relative_response_Q10_minus_R55']].copy()
    for name, pred in predictions.items():
        output[f'prediction_{name}'] = pred
        output[f'absolute_error_{name}'] = np.abs(pred - y)
    table(output, out / 'heldout_predictions.tsv')
    neighbors = []
    for i, locus in enumerate(data.locus):
        for name, distance in [('GI', gi), ('severity', sd)]:
            for rank, j in enumerate(np.argsort(distance[i], kind='stable')[:3], 1):
                neighbors.append(dict(held_out_locus=locus, model=name, rank=rank,
                                      training_locus=data.locus.iloc[j], distance=distance[i, j],
                                      interaction_anchors=int(anchors[i, j]) if name == 'GI' else None))
    table(pd.DataFrame(neighbors), out / 'prediction_neighbors.tsv')
    sensitivities = []
    for name, k, aggregate, age_filter in [('primary', 3, 'mean', False),
                                         ('one_neighbor', 1, 'mean', False),
                                         ('five_neighbors', 5, 'mean', False),
                                         ('median_arms', 3, 'median', False),
                                         ('exclude_younger_age', 3, 'mean', True)]:
        try:
            d = effects(arms, mapping, aggregate=aggregate)
            if age_filter:
                d = d[d.exclude_age_sensitivity == 0].reset_index(drop=True)
            w, _, _, _ = build_weights(matrix, d, k)
            result, _ = evaluate(d.relative_response_Q10_minus_R55.to_numpy(), w)
            sensitivities.append(dict(analysis=name, status='EVALUABLE', **result))
        except NotEvaluable as exc:
            sensitivities.append(dict(analysis=name, status='NOT_EVALUABLE', reason=str(exc)))
    table(pd.DataFrame(sensitivities), out / 'sensitivity_summary.tsv')
    sing = effects(arms, mapping, endpoint='SING')
    sing = sing[['locus', 'Q10_adjusted_log2_response', 'R55_adjusted_log2_response',
                 'relative_response_Q10_minus_R55']].rename(columns={
                     'Q10_adjusted_log2_response': 'SING_Q10_adjusted_log2_response',
                     'R55_adjusted_log2_response': 'SING_R55_adjusted_log2_response',
                     'relative_response_Q10_minus_R55': 'SING_relative_response'})
    compare = data[['locus', 'relative_response_Q10_minus_R55']].merge(sing, on='locus')
    compare['preference_direction_agrees'] = (np.sign(compare.relative_response_Q10_minus_R55)
                                             == np.sign(compare.SING_relative_response))
    table(compare, out / 'climbing_descriptive_comparison.tsv')
    summary['climbing_genotypes'] = len(compare)
    summary['climbing_direction_agreements'] = int(compare.preference_direction_agrees.sum())
    write_json(out / 'validation_summary.json', summary)
    write_json(out / 'environment.json', dict(python=platform.python_version(), packages={
        p: importlib.metadata.version(p) for p in ['numpy', 'pandas', 'scipy', 'openpyxl', 'matplotlib']}))
    print(json.dumps(clean_json(summary), indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache', type=Path, default=ROOT / 'data_raw/fly_genetic_rescue')
    p.add_argument('--output', type=Path, default=ROOT / 'results/v2/fly_genetic_rescue')
    p.add_argument('--fetch', action='store_true')
    args = p.parse_args()
    try:
        run(args.cache, args.output, args.fetch)
    except NotEvaluable as exc:
        args.output.mkdir(parents=True, exist_ok=True)
        write_json(args.output / 'validation_summary.json', dict(
            classification='NOT_EVALUABLE', reason=str(exc), freeze_commit=FREEZE))
        raise SystemExit(f'NOT_EVALUABLE: {exc}')


if __name__ == '__main__':
    main()
