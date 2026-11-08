"""Does reusing a different outer-fold predictor expose target-fold groups?"""
import json
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
def main():
 m=pd.read_csv(R/'analysis/reviewer_response_20261102/six_seed_nested_autodl_20261106/data/nested_folds.csv')
 records=[]
 for f in range(5):
  held=set(m[(m.outer_fold==f)&(m.role=='test')].group_id)
  for g in range(5):
   trained=set(m[(m.outer_fold==g)&(m.role=='train')].group_id)
   records.append({'target_fold':f,'source_fold':g,'target_test_groups_seen_in_source_train':len(held&trained)})
 assert all(r['target_test_groups_seen_in_source_train']==0 for r in records if r['target_fold']==r['source_fold'])
 assert all(r['target_test_groups_seen_in_source_train']>0 for r in records if r['target_fold']!=r['source_fold'])
 result={'rows':records,'interpretation':'Different outer-fold model train groups overlap target-fold held groups. Naive reuse of global OOF difficulty for training another outer fold is not strictly outer-safe; group counts do not prove actual prediction impact.'}
 (O/'difficulty_lineage_audit.json').write_text(json.dumps(result,indent=2),encoding='utf8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
