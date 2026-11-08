import json
from pathlib import Path
from difflib import SequenceMatcher
from collections import Counter
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
audit=json.loads((ROOT/'analysis/cys0_common_error_assay_20261108/audit.json').read_text(encoding='utf-8'))
anns=json.loads((ROOT/'analysis/family_mechanism_comparison_20261026/all_training_annotation_index.json').read_text(encoding='utf-8'))
profiles=pd.read_csv(ROOT/'analysis/six_seed_1d_weakness_20261107/sample_profiles.csv').set_index('sample_id')
split=pd.read_csv(ROOT/'analysis/reviewer_response_20261102/six_seed_nested_autodl_20261106/data/nested_folds.csv')
archive=pd.read_csv(ROOT/'analysis/experimental_dataset_inventory_20261031/curated/experimental_records_all.csv',low_memory=False).fillna('')
index=pd.read_csv(ROOT/'analysis/error_evidence_linkage_20261025/all_training_evidence_index.csv').fillna('').set_index('sample_id')
targets=[r for r in audit['records'] if r['category']=='damage']
assert len(targets)==10 and sum(r['all_models_call_nontoxic'] for r in targets)==8
byseq={seq:g.to_dict('records') for seq,g in archive.groupby('sequence')}
fields=['source','source_id','n_terminal','c_terminal','target','endpoint','measurement_raw','concentration_raw','temperature','exposure_time','pmid','doi','table_figure','source_locator','extraction_method']
def evidence(sid,sequence):
    a=anns.get(sid,{})
    old=index.loc[sid] if sid in index.index else {}
    return {'name':a.get('name') or (old.get('name','') if len(old) else ''),
            'family':a.get('curated_family',[]),
            'prior_note':old.get('manual_note','') if len(old) else '',
            'prior_source_urls':old.get('source_urls','') if len(old) else '',
            'archived_assays':[{k:v for k,v in rec.items() if k in fields and v!=''} for rec in byseq.get(sequence,[])]}
def similarity(a,b):
    return round(SequenceMatcher(None,a,b,autojunk=False).ratio(),5)
def pack(rec,score):
    return {'sample_id':rec.sample_id,'sequence':rec.sequence,'group_id':rec.group_id,
            'original_label':int(rec.label),'similarity':score,**evidence(rec.sample_id,rec.sequence)}
results=[]
for t in targets:
    fold=t['fold']
    train=split[(split.outer_fold==fold)&(split.role=='train')]
    assert not (train.group_id==t['group_id']).any()
    candidates=[]
    for rec in train.itertuples(index=False):
        candidates.append((similarity(t['sequence'],rec.sequence),rec))
    candidates.sort(key=lambda x:(-x[0],abs(len(x[1].sequence)-len(t['sequence'])),x[1].sample_id))
    nearest={}
    for label in (0,1):
        chosen=[]; groups=set()
        for score,rec in candidates:
            if int(rec.label)==label and rec.group_id not in groups:
                chosen.append(pack(rec,score)); groups.add(rec.group_id)
                if len(chosen)==5: break
        nearest[str(label)]=chosen
    tf=set(anns.get(t['sample_id'],{}).get('curated_family',[]))
    family=[pack(rec,score) for score,rec in candidates if tf.intersection(anns.get(rec.sample_id,{}).get('curated_family',[]))]
    p=profiles
    controls=p[(p.fold==fold)&(p.label==1)&(p.cys==0)&(p.band=='high_support')&
               (~p.high_seed_disagreement.astype(bool))&(p.group_id!=t['group_id'])&
               ((p.length-len(t['sequence'])).abs()<=3)].copy()
    if len(controls):
        mind=(controls.length-len(t['sequence'])).abs().min()
        controls=controls[(controls.length-len(t['sequence'])).abs()==mind]
    control_rows=[]
    for sid,rec in controls.iterrows():
        control_rows.append({'sample_id':sid,'group_id':rec.group_id,'sequence':rec.sequence,
                             'p_1d':float(rec.p_1d),'similarity':similarity(t['sequence'],rec.sequence),
                             **evidence(sid,rec.sequence)})
    results.append({'target':t,'nearest_by_label':nearest,'same_annotated_family':family,'controls':control_rows})
assert len(results)==10
summary={'targets':len(results),'common_errors':sum(x['target']['all_models_call_nontoxic'] for x in results),
         'targets_with_controls':sum(bool(x['controls']) for x in results),
         'control_assignments':sum(len(x['controls']) for x in results),
         'unique_controls':len({c['sample_id'] for x in results for c in x['controls']}),
         'targets_with_archived_nearest':sum(any(n['archived_assays'] for ns in x['nearest_by_label'].values() for n in ns) for x in results)}
(OUT/'results.json').write_text(json.dumps({'summary':summary,'cases':results},ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(summary,ensure_ascii=False))
for x in results:
 t=x['target']
 print(t['sample_id'],t['prior_name'][:35], 'fold',t['fold'],'fam',[(z['sample_id'],z['original_label'],z['similarity']) for z in x['same_annotated_family']], 'top+',[(z['sample_id'],z['similarity'],len(z['archived_assays'])) for z in x['nearest_by_label']['1'][:2]],'top-',[(z['sample_id'],z['similarity'],len(z['archived_assays'])) for z in x['nearest_by_label']['0'][:2]],'controls',len(x['controls']))
