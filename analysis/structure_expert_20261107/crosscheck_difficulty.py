"""Compare saved 3D branch predictions on independently defined 1D cohorts."""

from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
PROFILE = OUT.parent / 'six_seed_1d_weakness_20261107/sample_profiles.csv'
PRED = OUT / 'esm_if_probe_predictions.csv'
KEYS = ['sample_id', 'fold', 'group_id', 'sequence', 'label', 'cys']


def summarize(part, cohort, fold, cutoff):
    old = (part.p_old3d >= cutoff).eq(part.label)
    new = (part.p_esmif >= cutoff).eq(part.label)
    result = dict(cohort=cohort, fold=fold, cutoff=cutoff, n=len(part),
                  groups=part.group_id.nunique(), toxic=int(part.label.sum()),
                  egnn_correct=int(old.sum()), esmif_correct=int(new.sum()),
                  both_correct=int((old & new).sum()), egnn_only=int((old & ~new).sum()),
                  esmif_only=int((~old & new).sum()), both_wrong=int((~old & ~new).sum()))
    assert sum(result[k] for k in ('both_correct', 'egnn_only', 'esmif_only', 'both_wrong')) == len(part)
    for name, right in [('egnn', old), ('esmif', new)]:
        result[f'{name}_false_positive'] = int((~right & part.label.eq(0)).sum())
        result[f'{name}_false_negative'] = int((~right & part.label.eq(1)).sum())
        result[f'{name}_mean_bce'] = part[f'{name}_bce'].mean()
    return result


def main():
    profile = pd.read_csv(PROFILE, usecols=KEYS + ['length', 'p_1d', 'band', 'high_seed_disagreement'])
    prediction = pd.read_csv(PRED, usecols=KEYS + ['p_old3d', 'p_esmif'])
    data = profile.merge(prediction, on=KEYS, validate='one_to_one')
    assert len(data) == len(profile) == len(prediction) == 6387
    assert data.sample_id.is_unique
    for name, column in [('egnn', 'p_old3d'), ('esmif', 'p_esmif')]:
        assert data[column].between(0, 1).all()
        data[f'{name}_p_true'] = np.where(data.label.eq(1), data[column], 1-data[column])
        data[f'{name}_bce'] = -np.log(data[f'{name}_p_true'].clip(lower=1e-12))
    hard = data[data.label.eq(1) & data.cys.eq(0) & data.band.eq('strong_reverse')]
    negative = data[data.label.eq(0) & data.cys.eq(0)]
    controls = {}
    matched_cases = {}
    for name, pool in [('matched_negative', negative), ('matched_negative_1d_le_0.1', negative[negative.p_1d.le(.1)])]:
        ids = set()
        matched_cases[name] = 0
        for case in hard.itertuples():
            candidates = pool[pool.fold.eq(case.fold) & pool.group_id.ne(case.group_id)
                              & (pool.length-case.length).abs().le(3)]
            ids.update(candidates.sample_id)
            matched_cases[name] += int(not candidates.empty)
        controls[name] = pool[pool.sample_id.isin(ids)]
    cohorts = {'cys0_toxic_strong_reverse': hard, 'weaker_reverse': data[data.band.eq('weaker_reverse')],
               'near_neutral': data[data.band.eq('near_neutral')],
               'high_seed_disagreement': data[data.high_seed_disagreement], **controls}
    assert [len(cohorts[k]) for k in list(cohorts)[:4]] == [67, 180, 166, 962]
    summaries, cases = [], []
    for name, part in cohorts.items():
        for cutoff in (.4, .5, .6):
            summaries.append(summarize(part, name, 'all', cutoff))
            for fold in range(5):
                summaries.append(summarize(part[part.fold.eq(fold)], name, fold, cutoff))
        listed = part.copy()
        listed['cohort'] = name
        old = listed.p_old3d.ge(.5).eq(listed.label)
        new = listed.p_esmif.ge(.5).eq(listed.label)
        listed['outcome'] = np.select([old & new, old, new], ['both_correct', 'egnn_only', 'esmif_only'],
                                      default='both_wrong')
        cases.append(listed)
    report = pd.DataFrame(summaries)
    report.to_csv(OUT / 'difficulty_crosscheck_summary.csv', index=False, encoding='utf-8-sig')
    pd.concat(cases).to_csv(OUT / 'difficulty_crosscheck_cases.csv', index=False, encoding='utf-8-sig')
    meta = {'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (PROFILE, PRED)},
            'aligned_rows': len(data), 'matched_positive_cases': matched_cases, 'main_cutoff': .5,
            'match_rule': 'Same fold, Cys=0, length difference <=3, different sequence group. Unique negative controls.',
            'sensitivity_cutoffs': [.4, .6], 'test1_test2_used': False}
    (OUT / 'difficulty_crosscheck_provenance.json').write_text(json.dumps(meta, indent=2), encoding='utf-8')
    print(report[report.cutoff.eq(.5)].to_string(index=False))


if __name__ == '__main__':
    main()
