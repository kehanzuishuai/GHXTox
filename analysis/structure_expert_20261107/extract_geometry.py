"""Training-only, descriptive geometric inventory from verified structure cache."""
import csv, hashlib, json
from pathlib import Path
import torch

R = Path(__file__).resolve().parents[2]
O = Path(__file__).resolve().parent
CACHE = Path('D:/dachuang2026/GHXTox-ziliao/train_chemical_sites_final_esm2_multisource_1_30_full2795.pt')
EXPECTED = 'f8a6ddf7dd81198142ad8c30a0f7926de1967618b01720bd90cb483a45dfdb5b'

def main():
    torch.set_num_threads(1)
    h = hashlib.sha256()
    with CACHE.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024), b''): h.update(b)
    assert h.hexdigest() == EXPECTED
    with (R/'analysis/3d_contribution_20261106/sample_analysis.csv').open(encoding='utf-8-sig',newline='') as f:
        reference = {r['sample_id']:r for r in csv.DictReader(f) if r['split']=='training_oof'}
    assert len(reference)==6387
    records = torch.load(CACHE,map_location='cpu',weights_only=False)['records']
    assert len(records)==6387
    output=[]; seen=set()
    for record in records:
        sid=record['sample_id']; seq=record['sequence']; ref=reference[sid]
        assert sid not in seen and seq==ref['sequence'] and int(record['label'])==int(ref['label'])
        seen.add(sid); x=record['coords'].float(); p=record['plddt'].float().flatten(); n=len(seq)
        assert tuple(x.shape)==(n,3) and len(p)==n and torch.isfinite(x).all() and torch.isfinite(p).all()
        dist=torch.cdist(x,x); idx=torch.arange(n); sep=(idx[:,None]-idx[None,:]).abs()
        eligible=torch.triu(sep>2,diagonal=1); denominator=int(eligible.sum())
        contacts=eligible & (dist<8.0)
        ci=torch.tensor([i for i,c in enumerate(seq) if c=='C'],dtype=torch.long)
        cysdist=dist[ci[:,None],ci[None,:]] if len(ci)>1 else None
        cysmin=float(cysdist[torch.triu(torch.ones_like(cysdist,dtype=torch.bool),diagonal=1)].min()) if cysdist is not None else None
        rg=float(((x-x.mean(dim=0)).square().sum(dim=1).mean()).sqrt())
        output.append({'sample_id':sid,'group_id':ref['group_id'],'fold':int(float(ref['fold'])),'sequence':seq,'label':int(ref['label']),
          'length':n,'cys':seq.count('C'),'mean_plddt':float(p.mean()),'radius_gyration_ca':rg,
          'nonlocal_contact_density_8A':int(contacts.sum())/denominator if denominator else 0.0,
          'nonlocal_contact_pairs_8A':int(contacts.sum()),'cys_min_ca_distance':cysmin,
          'p_seq6':float(ref['p_seq6']),'p_3d':float(ref['p_3d']),
          'seq_correct':ref['seq_correct'],'structure_correct':ref['structure_correct']})
    assert seen==set(reference)
    with (O/'training_geometry.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
    (O/'geometry_provenance.json').write_text(json.dumps({'cache':str(CACHE),'sha256':h.hexdigest(),'n':len(output),'test1_test2_used':False,'notes':'C-alpha predicted structures; descriptive proxies, not disulfide or experimental evidence'},indent=2),encoding='utf8')
    print('Extracted',len(output),'verified training records')

if __name__=='__main__': main()
