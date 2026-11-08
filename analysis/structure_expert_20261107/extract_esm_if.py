"""Frozen, centered-coordinate ESM-IF1 mean representations for training only."""
import hashlib,json,time
from pathlib import Path
import numpy as np
import torch,esm,esm.inverse_folding

O=Path(__file__).resolve().parent
MODEL=O/'esm_if1_gvp4_t16_142M_UR50.pt'
CACHE=Path('D:/dachuang2026/GHXTox-ziliao/train_chemical_sites_final_esm2_multisource_1_30_full2795.pt')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 assert not (O/'esm_if_train_mean.npz').exists()
 assert sha(CACHE)=='f8a6ddf7dd81198142ad8c30a0f7926de1967618b01720bd90cb483a45dfdb5b'
 assert sha(MODEL)=='be4ba36edec22a9bfaa4946ff6b2815f1f19d8a3d7e0eada8b796d5a0eae9fd4'
 device=torch.device('cuda' if torch.cuda.is_available() else 'cpu')
 model,alphabet=esm.pretrained.load_model_and_alphabet_local(str(MODEL)); model.eval().to(device)
 records=torch.load(CACHE,map_location='cpu',weights_only=False)['records']; assert len(records)==6387
 converter=esm.inverse_folding.util.CoordBatchConverter(alphabet); features=[];ids=[];seqs=[];t=time.perf_counter()
 with torch.inference_mode():
  for i,r in enumerate(records):
   x=r['backbone_coords'][:,:3,:].float();mask=r['backbone_mask'][:,:3];seq=r['sequence']
   assert tuple(x.shape)==(len(seq),3,3) and bool(mask.bool().all()) and torch.isfinite(x).all()
   x=x-x.mean(dim=(0,1)); coords,confidence,_,_,padding=converter([(x.numpy(),None,None)],device=device)
   emb=model.encoder(coords,padding,confidence,return_all_hiddens=False)['encoder_out'][0][1:-1,0]
   assert tuple(emb.shape)==(len(seq),512) and torch.isfinite(emb).all()
   features.append(emb.mean(0).cpu().numpy());ids.append(r['sample_id']);seqs.append(seq)
   if (i+1)%500==0:print('encoded',i+1,'elapsed',round(time.perf_counter()-t,1),flush=True)
 assert len(set(ids))==6387
 np.savez_compressed(O/'esm_if_train_mean.npz',features=np.stack(features),sample_ids=np.asarray(ids),sequences=np.asarray(seqs))
 manifest={'model_sha256':sha(MODEL),'cache_sha256':sha(CACHE),'features_sha256':sha(O/'esm_if_train_mean.npz'),'n':len(ids),'width':512,'seconds':time.perf_counter()-t,'device':str(device),'centering':'mean of all N CA C atoms per peptide','confidence_input':'1.0 default','training_only':True,'toxicity_labels_used':False,'limitation':'predicted geometry; preliminary rotation/translation checks not exact'}
 (O/'esm_if_features_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8');print(json.dumps(manifest,indent=2))
if __name__=='__main__':main()
