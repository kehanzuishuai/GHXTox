"""Exploratory fixed logistic probe of frozen ESM-IF1; nested held-out groups only."""
import csv,hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import matthews_corrcoef,average_precision_score,roc_auc_score

R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
sys.path.insert(0,str(R/'analysis/six_seed_gate_probe_20261106'))
from run_probe import collect,prior,prob,manifest,TH
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 fp=O/'esm_if_train_mean.npz';meta=json.loads((O/'esm_if_features_manifest.json').read_text());assert sha(fp)==meta['features_sha256']
 z=np.load(fp,allow_pickle=False); ids=z['sample_ids'];seqs=z['sequences'];features=z['features'];assert features.shape==(6387,512) and np.isfinite(features).all() and len(set(ids))==6387
 by={str(v):i for i,v in enumerate(ids)};out=[];folds=[]
 for fold in range(5):
  roles={role:[r for r in manifest if int(r['outer_fold'])==fold and r['role']==role] for role in ('train','validation','calibration','test')}
  groups={k:set(r['group_id'] for r in rows) for k,rows in roles.items()}
  for a in groups:
   for b in groups:
    if a!=b:assert groups[a].isdisjoint(groups[b])
  def xy(rows):
   for r in rows:assert r['sample_id'] in by and seqs[by[r['sample_id']]]==r['sequence']
   return features[[by[r['sample_id']] for r in rows]],np.array([int(r['label']) for r in rows])
  x,y=xy(roles['train']);model=make_pipeline(StandardScaler(),LogisticRegression(C=.1,max_iter=1000,random_state=20261107))
  model.fit(x,y);cal=collect(fold,'calibration'); test=collect(fold,'test');
  for role,d in [('calibration',cal),('test',test)]:assert [r['sample_id'] for r in roles[role]]==[r['sample_id'] for r in d['rows']]
  xc,yc=xy(roles['calibration']);xt,yt=xy(roles['test']);assert np.array_equal(yc,cal['y']) and np.array_equal(yt,test['y'])
  rawcal=model.predict_proba(xc)[:,1];rawtest=model.predict_proba(xt)[:,1]
  fit,_=next(StratifiedGroupKFold(n_splits=2,shuffle=True,random_state=20261106+fold).split(np.zeros(len(yc)),yc,cal['groups']))
  c=prior.fit(rawcal[fit],yc[fit]);p=prob(rawtest,c)
  for i,r in enumerate(roles['test']):
   out.append({'fold':fold,'sample_id':r['sample_id'],'group_id':r['group_id'],'sequence':r['sequence'],'label':int(yt[i]),'cys':r['sequence'].count('C'),
    'p_seq6':float(test['q'][i]),'p_old3d':float(test['p3'][i]),'p_esmif':float(p[i]),'p_fused_old':float(.72*test['q'][i]+.28*test['p3'][i]),'p_fused_esmif':float(.72*test['q'][i]+.28*p[i])})
  folds.append({'fold':fold,'train_n':len(x),'calibration_fit_n':len(fit),'test_n':len(xt),'platt':c})
 assert len(out)==6387 and len({r['sample_id'] for r in out})==6387
 with (O/'esm_if_probe_predictions.csv').open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
 y=np.array([r['label'] for r in out]);p0=np.array([r['p_seq6'] for r in out]);base=(p0>=TH)==y
 results={}
 for name in ['p_seq6','p_old3d','p_esmif','p_fused_old','p_fused_esmif']:
  p=np.array([r[name] for r in out]);correct=(p>=TH)==y
  results[name]={'mcc':float(matthews_corrcoef(y,p>=TH)),'auroc':float(roc_auc_score(y,p)),'auprc':float(average_precision_score(y,p)),'brier':float(np.mean((p-y)**2)),'errors':int((~correct).sum()),'corrects_vs_seq6':int((correct&~base).sum()),'harms_vs_seq6':int((~correct&base).sum())}
 group=np.array([r['group_id'] for r in out]);unique=np.unique(group);chunks=[np.flatnonzero(group==g) for g in unique];rng=np.random.default_rng(20261107)
 predictions={name:np.array([r[name] for r in out])>=TH for name in results}
 delta={name:[] for name in ['p_fused_old','p_fused_esmif']}
 for _ in range(1000):
  ix=np.concatenate([chunks[k] for k in rng.integers(len(chunks),size=len(chunks))]);base_mcc=matthews_corrcoef(y[ix],predictions['p_seq6'][ix])
  for name in delta:delta[name].append(float(matthews_corrcoef(y[ix],predictions[name][ix])-base_mcc))
 comparisons={name:{'delta_mcc':results[name]['mcc']-results['p_seq6']['mcc'],'group_bootstrap_ci95':np.quantile(vals,[.025,.975]).tolist()} for name,vals in delta.items()}
 for name in results:
  pred=predictions[name];results[name]['fold_mcc']=[float(matthews_corrcoef(y[np.array([r['fold']==f for r in out])],pred[np.array([r['fold']==f for r in out])])) for f in range(5)]
  c=np.array([r['cys'] for r in out]);results[name]['cys0_toxic_errors']=int(((c==0)&(y==1)&~pred).sum());results[name]['cys7plus_nontoxic_errors']=int(((c>=7)&(y==0)&pred).sum())
 res={'protocol':'PROTOCOL.md','description':'post hoc initial fixed frozen geometry probe, not optimized or confirmatory','feature_sha256':sha(fp),'folds':folds,'metrics':results,'comparisons_vs_seq6':comparisons,'test1_test2_used':False,'limitations':'sequence group OOF, no new independent test; six seed and old 3d upstream limitations; ESM-IF1 unsupervised pretraining overlap unknown; no geometry-free matched control; frozen linear probe only; fixed-prediction bootstrap excludes retraining uncertainty'}
 (O/'esm_if_probe_summary.json').write_text(json.dumps(res,indent=2),encoding='utf8');print(json.dumps(res,indent=2))
if __name__=='__main__':main()
