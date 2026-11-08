"""Profile six-seed 1D held-out predictions without loading structural predictions."""

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
REVIEW = ROOT / 'analysis/reviewer_response_20261102'
PACKAGE = REVIEW / 'six_seed_nested_autodl_20261106'
OLD = REVIEW / 'nested_1d_received_20261106/nested_1d_results/runs'
NEW = REVIEW / 'six_seed_received_20261106/six_seed_new_results/runs'
SPLIT = PACKAGE / 'data/nested_folds.csv'
GEOMETRY = ROOT / 'analysis/structure_expert_20261107/training_geometry.csv'
EVIDENCE = ROOT / 'analysis/lowcys_background_20261027/all_samples_with_symmetric_family.json'
SEEDS = (42, 123, 2025, 314, 2718, 2026)
KEYS = ['sample_id', 'group_id', 'fold', 'sequence', 'label']
FEATURES = ['length', 'cys', 'mean_plddt', 'nonlocal_contact_density_8A', 'radius_gyration_ca',
            'net_charge_proxy', 'mean_hydropathy', 'fragment', 'conflict', 'uniform_flag', 'family_known']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def one_dimensional_scores():
    spec = importlib.util.spec_from_file_location('validated_comparison', PACKAGE / 'compare.py')
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    manifest = pd.read_csv(SPLIT)
    expected = {
        'manifest_sha256': sha(SPLIT),
        'config_sha256': sha(PACKAGE / 'config/model.json'),
        'source_sha256': sha(PACKAGE / 'src/ghxtox/lora_mamba_1d.py'),
        'prott5_weight_sha256': '7f51ba885541c7dc569d46b796af57cc7a2ba7945107dced4f19d1b5ec091157',
        'prott5_aux_sha256': {
            'config.json': 'c93cc8c31433fb062cad1825a9c16cde962a2015c8a3b10cfdc5f76a85128db3',
            'spiece.model': '74da7b4afcde53faa570114b530c726135bdfcdb813dec3abfb27f9d44db7324',
        },
    }
    scores, calibration, hashes = [], [], {}
    for fold in range(5):
        fold_rows = manifest[manifest.outer_fold.eq(fold)]
        role_groups = {role: set(rows.group_id) for role, rows in fold_rows.groupby('role')}
        assert sum(map(len, role_groups.values())) == len(set.union(*role_groups.values()))
        predictions = {}
        for role in ('calibration', 'test'):
            rows = fold_rows[fold_rows.role.eq(role)].copy()
            assert rows.sample_id.is_unique
            reference = rows.set_index('sample_id')
            columns = []
            for seed in SEEDS:
                base = OLD if seed in SEEDS[:3] else NEW
                run = prior.load_run(base, seed, fold, role, expected)
                assert set(run) == set(reference.index)
                for sid, pred in run.items():
                    row = reference.loc[sid]
                    assert (pred['sequence'], int(pred['label']), pred['group_id'], int(pred['fold']), int(pred['seed'])) == (
                        row.sequence, row.label, row.group_id, fold, seed)
                columns.append([float(run[sid]['probability']) for sid in rows.sample_id])
                path = base / f'seed{seed}_fold{fold}/{role}_predictions.csv'
                hashes[str(path.relative_to(ROOT))] = sha(path)
            raw = np.array(columns).T
            assert np.isfinite(raw).all() and ((raw >= 0) & (raw <= 1)).all()
            predictions[role] = (rows, raw)
        cal_rows, cal_raw = predictions['calibration']
        held, raw = predictions['test']
        slope, intercept = prior.fit(cal_raw.mean(axis=1), cal_rows.label.to_numpy())
        held = held.rename(columns={'outer_fold': 'fold'})[KEYS].copy()
        held['p_1d_raw_mean'] = raw.mean(axis=1)
        held['p_1d'] = prior.apply(held.p_1d_raw_mean.to_numpy(), (slope, intercept))
        held['seed_std_raw'] = raw.std(axis=1, ddof=1)
        raw_neutral = 1 / (1 + np.exp(intercept / slope))
        for j, seed in enumerate(SEEDS):
            held[f'p_raw_seed_{seed}'] = raw[:, j]
        seed_votes = (raw >= raw_neutral) == held.label.to_numpy()[:, None]
        held['seeds_support_label_side'] = seed_votes.sum(axis=1)
        scores.append(held)
        calibration.append({'fold': fold, 'calibration_n': len(cal_rows), 'held_n': len(held),
                            'slope': slope, 'intercept': intercept, 'raw_ensemble_neutral': raw_neutral})
    frame = pd.concat(scores, ignore_index=True)
    assert len(frame) == 6387 and frame.sample_id.is_unique and frame.group_id.nunique() == 5258
    assert frame.groupby('group_id').fold.nunique().eq(1).all()
    frame['p_true'] = np.where(frame.label.eq(1), frame.p_1d, 1 - frame.p_1d)
    frame['bce'] = -np.log(frame.p_true.clip(lower=1e-12))
    frame['distance_from_neutral'] = abs(frame.p_1d - .5)
    frame['band'] = np.select(
        [frame.p_true.le(.1), frame.p_true.lt(.4), frame.p_true.le(.6), frame.p_true.lt(.9)],
        ['strong_reverse', 'weaker_reverse', 'near_neutral', 'moderate_support'], default='high_support')
    frame['high_seed_disagreement'] = frame.seed_std_raw.ge(.15)
    return frame, calibration, hashes


def attach_context(frame):
    geometry = pd.read_csv(GEOMETRY, usecols=KEYS + ['length', 'cys', 'mean_plddt', 'radius_gyration_ca',
                                                    'nonlocal_contact_density_8A', 'cys_min_ca_distance'])
    frame = frame.merge(geometry, on=KEYS, validate='one_to_one')
    evidence = json.loads(EVIDENCE.read_text(encoding='utf-8-sig'))
    annotations = []
    for item in evidence:
        family = (item.get('current_family_sidecar') or {}).get('full_families') or item.get('curated_family') or []
        annotations.append({**{key: item[key] for key in KEYS},
                            **{key: item[key] for key in ['fragment', 'conflict', 'uniform_flag', 'unresolved_identity',
                                                          'original_database', 'original_endpoint', 'evidence_tier',
                                                          'strict_original_label_verified']},
                            'name': item.get('name'), 'family': ' | '.join(family), 'family_known': int(bool(family)),
                            'net_charge_proxy': item['features'].get('net_charge_proxy'),
                            'mean_hydropathy': item['features'].get('mean_hydropathy')})
    frame = frame.merge(pd.DataFrame(annotations), on=KEYS, validate='one_to_one')
    assert len(frame) == 6387 and frame.length.eq(frame.sequence.str.len()).all()
    assert frame.cys.eq(frame.sequence.str.count('C')).all()
    frame['cys_band'] = np.select([frame.cys.eq(0), frame.cys.ge(7)], ['0', '7+'], default='1–6')
    return frame


def cohorts(frame):
    return {**{name: frame.band.eq(name) for name in sorted(frame.band.unique())},
            'high_seed_disagreement': frame.high_seed_disagreement}


def summaries(frame):
    rows = []
    populations = [('all', frame)] + [(f'label_{label}', part) for label, part in frame.groupby('label')]
    populations += [(f'label_{label}_cys_{cys}', part) for (label, cys), part in frame.groupby(['label', 'cys_band'])]
    for name, part in populations:
        for band, mask in {'all': np.ones(len(part), dtype=bool), **cohorts(part)}.items():
            selected = part[mask]
            rows.append({'population': name, 'type': band, 'denominator': len(part), 'n': len(selected),
                         'groups': selected.group_id.nunique(), 'fraction': len(selected) / len(part),
                         'median_bce': selected.bce.median(), 'six_seeds_reverse_side': int(selected.seeds_support_label_side.eq(0).sum()),
                         **{f'median_{key}': selected[key].median() for key in FEATURES[:7]},
                         **{f'count_{key}': int(selected[key].sum()) for key in FEATURES[7:]},
                         'known_original_database': int(selected.original_database.ne('unknown').sum()),
                         'strict_label_verified': int(selected.strict_original_label_verified.sum())})
    return pd.DataFrame(rows)


def matched_controls(frame):
    controls = frame[frame.band.eq('high_support') & ~frame.high_seed_disagreement]
    pools = dict(tuple(controls.groupby(['fold', 'label', 'cys'])))
    comparisons = []
    for name, mask in cohorts(frame).items():
        if name in {'high_support', 'moderate_support'}:
            continue
        for row in frame[mask].itertuples(index=False):
            pool = pools.get((row.fold, row.label, row.cys))
            if pool is None:
                continue
            distance = abs(pool.length - row.length)
            candidates = pool[pool.group_id.ne(row.group_id) & distance.le(3)]
            if candidates.empty:
                continue
            nearest = candidates[distance.loc[candidates.index].eq(distance.loc[candidates.index].min())]
            comparisons.append({'sample_id': row.sample_id, 'type': name, 'fold': row.fold, 'label': row.label,
                                'cys_band': row.cys_band, 'n_controls': len(nearest),
                                'control_ids': '|'.join(nearest.sample_id),
                                **{f'delta_{key}': getattr(row, key) - nearest[key].median() for key in FEATURES}})
    paired = pd.DataFrame(comparisons)
    rows = []
    for population, part in [('all', frame), *[(f'label_{y}', p) for y, p in frame.groupby('label')],
                             *[(f'label_{y}_cys_{c}', p) for (y, c), p in frame.groupby(['label', 'cys_band'])]]:
        for name, mask in cohorts(part).items():
            if name in {'high_support', 'moderate_support'}:
                continue
            cases = part[mask]
            matched = paired[paired.type.eq(name) & paired.sample_id.isin(cases.sample_id)]
            rows.append({'population': population, 'type': name, 'cases': len(cases), 'matched': len(matched),
                         **{f'median_delta_{key}': matched[f'delta_{key}'].median() for key in FEATURES},
                         **{f'fold_medians_{key}': '|'.join(f'{fold}:{value:.4g}' for fold, value in
                            matched.groupby('fold')[f'delta_{key}'].median().items()) for key in FEATURES[2:5]}})
    return paired, pd.DataFrame(rows)


def main():
    scores, calibration, hashes = one_dimensional_scores()
    frame = attach_context(scores)
    summary = summaries(frame)
    pairs, matched = matched_controls(frame)
    families = []
    for name, mask in cohorts(frame).items():
        for (label, family), part in frame[mask & frame.family_known.eq(1)].groupby(['label', 'family']):
            families.append({'type': name, 'label': label, 'family_annotation': family, 'n': len(part),
                             'groups': part.group_id.nunique()})
    for name, table in [('sample_profiles', frame), ('population_summary', summary), ('matched_cases', pairs),
                        ('matched_summary', matched), ('annotated_families', pd.DataFrame(families))]:
        table.to_csv(OUT / f'{name}.csv', index=False, encoding='utf-8-sig')
    sensitivity = {
        'strong_reverse_at': {str(cut): int(frame.p_true.le(cut).sum()) for cut in (.05, .1, .2)},
        'near_neutral_radius': {str(cut): int(frame.distance_from_neutral.le(cut).sum()) for cut in (.05, .1, .15)},
    }
    assert int(summary[summary.population.eq('all') & ~summary.type.isin(['all','high_seed_disagreement'])].n.sum()) == len(frame)
    (OUT / 'provenance.json').write_text(json.dumps({'split_sha256': sha(SPLIT), 'geometry_sha256': sha(GEOMETRY),
        'evidence_sha256': sha(EVIDENCE), 'seed_files_sha256': hashes, 'calibration': calibration,
        'sensitivity': sensitivity, 'n': len(frame), 'group_count': int(frame.group_id.nunique()),
        'warning': 'descriptive held-out profiles only; not training weights'}, indent=2), encoding='utf-8')
    print(summary[summary.population.isin(['all','label_1_cys_0','label_0_cys_7+'])][
        ['population','type','n','denominator','groups','median_bce','six_seeds_reverse_side']].to_string(index=False))


if __name__ == '__main__':
    main()
