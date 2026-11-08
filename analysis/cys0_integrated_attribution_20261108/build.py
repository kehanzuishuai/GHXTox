import csv, json
from collections import Counter, defaultdict
from pathlib import Path

root = Path(__file__).resolve().parents[2]
out = Path(__file__).resolve().parent
audit = json.loads((root/'analysis/cys0_common_error_assay_20261108/audit.json').read_text(encoding='utf8'))
annotations = json.loads((root/'analysis/family_mechanism_comparison_20261026/all_training_annotation_index.json').read_text(encoding='utf8'))
manifest = root/'analysis/reviewer_response_20261102/six_seed_nested_autodl_20261106/data/nested_folds.csv'
with manifest.open(encoding='utf-8-sig', newline='') as f:
    folds = list(csv.DictReader(f))
train_by_fold = defaultdict(list)
for row in folds:
    if row['role'] == 'train':
        train_by_fold[int(row['outer_fold'])].append(row)
rows=[]
for item in audit['records']:
    ann=annotations.get(item['sample_id'], {})
    family=set(ann.get('curated_family', []))
    train=train_by_fold[item['fold']]
    hits=[x for x in train if family.intersection(annotations.get(x['sample_id'], {}).get('curated_family', []))] if family else []
    groups={x['group_id'] for x in hits}
    row={'sample_id':item['sample_id'], 'sequence':item['sequence'], 'outer_fold':item['fold'],
         'p_1d':item['p_prott5_six'], 'common_error':item['all_models_call_nontoxic'],
         'evidence_category':item['category'], 'note':item['note'], 'sources':item['references'],
         'original_label_trace':item['original_label_trace'],
         'curated_family':sorted(family), 'train_family_sample_count':len(hits) if family else None,
         'train_family_group_count':len(groups) if family else None,
         'train_family_positive_count':sum(x['label']=='1' for x in hits) if family else None,
         'family_coverage_scope':'外折 role=train，沿用旧 UniProt curated_family 注释，未核终点',
         'prior_assay_target':item['sample_id'] in {'peptide_102','peptide_251','peptide_347','peptide_424','peptide_44','peptide_463','peptide_59','peptide_97'}}
    rows.append(row)
assert len(rows)==67 and len({x['sample_id'] for x in rows})==67
assert sum(x['common_error'] for x in rows)==62
known=[x for x in rows if x['curated_family']]
summary={'n':len(rows),'common_errors':sum(x['common_error'] for x in rows),
         'categories':dict(Counter(x['evidence_category'] for x in rows)),
         'family_annotated':len(known),'family_train_has_match':sum(x['train_family_sample_count']>0 for x in known),
         'family_train_no_match':sum(x['train_family_sample_count']==0 for x in known),
         'family_unknown':len(rows)-len(known)}
(out/'attribution.json').write_text(json.dumps({'summary':summary,'rows':rows},ensure_ascii=False,indent=2),encoding='utf8')
with (out/'逐样本归因表.md').open('w',encoding='utf8') as f:
    f.write('# 六种子 Cys=0 强反向样本归因表\n\n复用已核查原文与旧家族注释。家族数量仅统计该样本外折 role=train 的同注释样本，未核对毒性终点。未知表示缺少家族注释。\n\n')
    f.write('|样本|折|六种子分数|共同漏判|证据类型|外折训练同家族 序列/组/正类|已核查的对象与证据|\n|---|---:|---:|---|---|---|---|\n')
    for x in rows:
        cov=f"{x['train_family_sample_count']}/{x['train_family_group_count']}/{x['train_family_positive_count']}" if x['curated_family'] else '未知'
        refs=' '.join(f'[来源{i+1}]({u})' for i,u in enumerate(x['sources']))
        note=x['note'].replace('|','／').replace('\n',' ')
        f.write(f"|{x['sample_id']}|{x['outer_fold']}|{x['p_1d']:.4f}|{'是' if x['common_error'] else '否'}|{x['evidence_category']}|{cov}|{note} {refs}|\n")
print(json.dumps(summary,ensure_ascii=False,indent=2))
