"""Join the 67 target sequences to historical evidence and a primary assay panel."""
from pathlib import Path
import json
import re
import xml.etree.ElementTree as ET

import pandas as pd

OUT = Path(__file__).resolve().parent
ANALYSIS = OUT.parent


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def text(element):
    return ' '.join(''.join(element.itertext()).split())


def main():
    data = pd.read_csv(ANALYSIS / 'six_seed_1d_weakness_20261107/sample_profiles.csv')
    target = data[data.label.eq(1) & data.cys.eq(0) & data.band.eq('strong_reverse')].copy()
    assert len(target) == 67 and target.sample_id.is_unique
    semantic = read_json(ANALYSIS / 'membrane_class_feature_diagnosis_20261026/semantic_and_feature_rows.json')
    by_id = {x['sample_id']: x for x in semantic}
    target['prior_membrane_semantics'] = [
        by_id[x.sample_id]['semantic_state'] if x.sample_id in by_id else ''
        for x in target.itertuples()
    ]
    for row in target.itertuples():
        if row.sample_id in by_id:
            entry = by_id[row.sample_id]
            assert (entry['sequence'], int(entry['label'])) == (row.sequence, row.label)
    xml = ET.parse(OUT / 'sources/PMC10611374.xml')
    panel = []
    for row in xml.findall(".//table-wrap[@id='toxins-15-00591-t001']//tr"):
        cells = [text(cell) for cell in row if cell.tag in ('td', 'th')]
        match = next((i for i, cell in enumerate(cells) if re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWY]{10,30}(?:-NH2)?', cell)), None)
        if match is None:
            continue
        chemical = cells[match]
        seq = chemical.removesuffix('-NH2')
        tail = cells[match+1:]
        assert len(tail) == 9, cells
        panel.append({'peptide': cells[match-1], 'sequence': seq, 'form_in_table': chemical,
                      'human_rbc_hemolysis_320_ug_ml': tail[3], 'human_rbc_ec50_uM': tail[4],
                      'rat_rbc_hemolysis_320_ug_ml': tail[6], 'rat_rbc_ec50_uM': tail[7]})
    assert len(panel) == 55
    panel = pd.DataFrame(panel)
    assert panel.form_in_table.is_unique
    joined = panel.merge(data[['sample_id','sequence','label','band','p_1d']], on='sequence', how='left')
    hits = target.merge(panel[['sequence','peptide','form_in_table','human_rbc_hemolysis_320_ug_ml','human_rbc_ec50_uM']],
                        on='sequence', how='left', validate='one_to_many')
    assert hits.sample_id.nunique() == 67
    cols = ['sample_id','sequence','fold','group_id','p_1d','bce','name','family','fragment','conflict',
            'uniform_flag','prior_membrane_semantics','peptide','form_in_table',
            'human_rbc_hemolysis_320_ug_ml','human_rbc_ec50_uM']
    hits[cols].to_csv(OUT / 'target_evidence.csv', index=False, encoding='utf-8-sig')
    joined.to_csv(OUT / 'mastoparan_assay_panel.csv', index=False, encoding='utf-8-sig')
    for name in ('test1', 'test2'):
        fasta = (ANALYSIS.parent / 'dataset' / f'{name}.fasta').read_text()
        sequences = {''.join(record.splitlines()[1:]).strip() for record in fasta.split('>') if record.strip()}
        joined[f'{name}_exact_sequence'] = joined.sequence.isin(sequences)
    joined.to_csv(OUT / 'mastoparan_panel_sequence_overlap.csv', index=False, encoding='utf-8-sig')
    print('Target n, prior membrane annotation, fragment, flag, exact assay panel hits:',
          len(target), target.prior_membrane_semantics.ne('').sum(), target.fragment.sum(),
          target.uniform_flag.sum(), hits.loc[hits.peptide.notna(),'sample_id'].nunique())
    print(hits[hits.peptide.notna()][cols].to_string(index=False))
    print('Panel coverage in training', joined.sample_id.nunique())
    print(joined[joined.peptide.isin(['EMP-EM1', 'EMP-EM2', 'Eumenitin-R'])].to_string(index=False))


if __name__ == '__main__':
    main()
