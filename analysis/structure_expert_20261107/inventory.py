"""Descriptive, overlapping discovery populations; never used to select test labels."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

R=Path(__file__).resolve().parents[2]; O=Path(__file__).resolve().parent
G=O/'training_geometry.csv'; S=R/'analysis/six_seed_gate_probe_20261106/crossfit_outer_predictions.csv'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 d=pd.read_csv(G); s=pd.read_csv(S)[['sample_id','sequence','label','fold','group_id','seed_std_raw']]
 assert not d.sample_id.duplicated().any() and not s.sample_id.duplicated().any()
 d=d.merge(s,on=['sample_id','sequence','label','fold','group_id'],how='inner',validate='one_to_one')
 assert len(d)==6387 and np.isfinite(d.seed_std_raw).all()
 d['seq_correct']=d.seq_correct.astype(bool); d['structure_correct']=d.structure_correct.astype(bool)
 d['one_d_error']=~d.seq_correct; d['both_error']=~d.seq_correct & ~d.structure_correct
 d['uncertain_seq']=d.p_seq6.between(.3,.7,inclusive='both')
 d['six_seed_disagree']=d.seed_std_raw>=.15
 d['geometry_dense']=d.nonlocal_contact_density_8A>=.15
 d['plddt_confident']=d.mean_plddt>=.7
 d['cys_zero']=d.cys.eq(0); d['cys_sevenplus']=d.cys.ge(7)
 d['dense_quality']=d.geometry_dense & d.plddt_confident
 populations={'all':d.index==d.index,'1d_uncertain_0.3_0.7':d.uncertain_seq,
  'six_seed_std_at_least_0.15':d.six_seed_disagree,'nonlocal_density_at_least_0.15':d.geometry_dense,
  'dense_and_quality_at_least_0.7':d.dense_quality,'cys_zero':d.cys_zero,'cys_sevenplus':d.cys_sevenplus,
  'seq_wrong_evaluation_only':d.one_d_error,'both_wrong_evaluation_only':d.both_error}
 out=[]
 for name,mask in populations.items():
  for label in ('all',0,1):
   p=d[mask] if label=='all' else d[mask & d.label.eq(label)]
   if p.empty:continue
   out.append({'population':name,'label':label,'n':len(p),'group_count':p.group_id.nunique(),
    'seq_errors':int(p.one_d_error.sum()),'both_errors':int(p.both_error.sum()),
    'one_d_error_rate':float(p.one_d_error.mean()),'mean_length':float(p.length.mean()),
    'median_plddt':float(p.mean_plddt.median()),'median_nonlocal_density':float(p.nonlocal_contact_density_8A.median()),
    'fold_counts':' '.join(f'{f}:{len(v)}' for f,v in p.groupby('fold'))})
 result=pd.DataFrame(out);result.to_csv(O/'population_summary.csv',index=False,encoding='utf-8-sig')
 d.to_csv(O/'candidate_inventory.csv',index=False,encoding='utf-8-sig')
 (O/'inventory_manifest.json').write_text(json.dumps({'n':len(d),'geometry_sha256':digest(G),'seed_uncertainty_sha256':digest(S),
  'note':'Overlapping post-hoc descriptive groups; OOF seed_std_raw from 6 seeds, other prediction columns from another fold calibration; labels for description only; not proof of geometry usefulness.'},indent=2),encoding='utf8')
 print(result[result.label.eq('all')].to_string(index=False))
if __name__=='__main__':main()
