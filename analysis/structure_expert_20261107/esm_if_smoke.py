"""Eight training records: ESM-IF1 encoder technical check, not toxicity training."""
import argparse,hashlib,json,time
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
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=8);args=ap.parse_args()
 assert MODEL.stat().st_size==1700450121
 device=torch.device('cuda' if torch.cuda.is_available() else 'cpu')
 model,alphabet=esm.pretrained.load_model_and_alphabet_local(str(MODEL));model.eval().to(device)
 all_records=torch.load(CACHE,map_location='cpu',weights_only=False)['records']
 assert 1<=args.n<=len(all_records)
 indices=np.linspace(0,len(all_records)-1,args.n,dtype=int)
 records=[all_records[int(i)] for i in indices]
 rows=[];t=time.perf_counter()
 with torch.inference_mode():
  for r in records:
   x=r['backbone_coords'][:,:3,:];mask=r['backbone_mask'][:,:3]
   assert tuple(x.shape)==(len(r['sequence']),3,3) and bool(mask.bool().all()) and torch.isfinite(x).all()
   converter=esm.inverse_folding.util.CoordBatchConverter(alphabet)
   def encode(coords):
    coords,conf,_,_,padding=converter([(coords.numpy(),None,None)],device=device)
    return model.encoder(coords,padding,conf,return_all_hiddens=False)['encoder_out'][0][1:-1,0]
   emb=encode(x);repeat=encode(x);shift=encode(x+torch.tensor([7.,-3.,4.])); err=float((emb-shift).abs().max().item())
   repeat_err=float((emb-repeat).abs().max().item())
   xc=x-x.mean(dim=(0,1)); xs=x+torch.tensor([7.,-3.,4.]);xs=xs-xs.mean(dim=(0,1))
   center_embedding=encode(xc); shifted_embedding=encode(xs)
   centered_err=float((center_embedding-shifted_embedding).abs().max().item())
   rotation=torch.tensor([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
   rotated=xc @ rotation.T
   rotated_embedding=encode(rotated)
   rotate_err=float((center_embedding-rotated_embedding).abs().max().item())
   pooled_shift=float((center_embedding.mean(0)-shifted_embedding.mean(0)).abs().max().item())
   pooled_rotate=float((center_embedding.mean(0)-rotated_embedding.mean(0)).abs().max().item())
   assert emb.shape[0]==len(r['sequence']) and torch.isfinite(emb).all() and repeat_err<1e-5, (r['sample_id'],list(emb.shape),repeat_err)
   rows.append({'sample_id':r['sample_id'],'length':len(r['sequence']),'shape':list(emb.shape),'repeat_max_abs_difference':repeat_err,'translation_max_abs_difference':err,'centered_translation_max_abs_difference':centered_err,'centered_rotation_max_abs_difference':rotate_err,'pooled_centered_shift_max_abs_difference':pooled_shift,'pooled_centered_rotation_max_abs_difference':pooled_rotate})
 res={'model':'ESM-IF1','source':'https://dl.fbaipublicfiles.com/fair-esm/models/esm_if1_gvp4_t16_142M_UR50.pt','sha256':sha(MODEL),'device':str(device),'n':len(rows),'seconds':time.perf_counter()-t,'samples':rows,'note':'encoder sanity only, no labels, no evaluation or training'}
 (O/'esm_if_smoke.json').write_text(json.dumps(res,indent=2),encoding='utf8');print(json.dumps(res,indent=2))
if __name__=='__main__':main()
