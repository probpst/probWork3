import pickle,numpy as np,sys,time,json; sys.path.insert(0,'/home/claude/hg')
import mb
from mb import FACES,fax
L=pickle.load(open('final_L.pkl','rb')); R=pickle.load(open('final_R.pkl','rb'))
blocks=L['blocks']+R['blocks']; inout=L['inout']+R['inout']
# ---- auto-tag inlet/outlet faces on every block (cells whose 4 nodes lie in an inlet/outlet disc)
ntag=0
for B in blocks:
    B.tags=[t for t in B.tags]  # keep core end-cap tags
    for f in FACES:
        Fc=B.face(f)
        for io in inout:
            if io['dom']!=B.dom: continue
            dist=(Fc-io['P'])@io['n']; rad=np.linalg.norm(Fc-io['P']-dist[...,None]*io['n'],axis=-1)
            on=(np.abs(dist)<1e-6)&(rad<io['r']*1.05+1e-6)
            cell=on[:-1,:-1]&on[1:,:-1]&on[:-1,1:]&on[1:,1:]
            if cell.any():
                for pr,qr in mb.rects(cell):
                    if not any(t['face']==f and t['pr']==pr and t['qr']==qr for t in B.tags):
                        B.tag(f,io['type'],io['label'],pr,qr); ntag+=1
print('auto-tagged faces',ntag)
t=time.time(); ifc,part=mb.connectivity(blocks); print('connectivity %.0fs'%(time.time()-t),flush=True)
bcs=mb.bc_surfaces(blocks,part)
# ---- checks
err=max(I['err'] for I in ifc); nsame=sum(I['order']=='same' for I in ifc); ntr=len(ifc)-nsame
ambig=[I for I in ifc if I['order']=='transposed' and (I['pr'][1]-I['pr'][0])==(I['qr'][1]-I['qr'][0])]
cover_ok=True; totcells=0
for b,B in enumerate(blocks):
    for f in FACES:
        n=part[(b,f)].size; cov=(part[(b,f)]>=0).sum()
        for s in bcs.get(b,[]):
            if s['f']==f: cov+=(s['pr'][1]-s['pr'][0]+1)*(s['qr'][1]-s['qr'][0]+1)
        if cov!=n: cover_ok=False; print('coverage mismatch',B.name,f,cov,n)
negs=sum(int((mb.hex_jac(B.X)<=0).sum()) for B in blocks)
vol=sum(float(mb.hex_vol(B.X).sum()) for B in blocks if B.dom=='L'),sum(float(mb.hex_vol(B.X).sum()) for B in blocks if B.dom=='R')
cells=sum(int(np.prod(np.array(B.dims)-1)) for B in blocks)
from collections import Counter
bcc=Counter((s['type'],s['label']) for S in bcs.values() for s in S)
print(f'blocks {len(blocks)} cells {cells} interfaces {len(ifc)} (same {nsame}, transposed {ntr}, ambiguous-transposed {len(ambig)}) max point err {err:.2e} coverage {cover_ok} neg {negs} vol L {vol[0]/1000:.1f} mL R {vol[1]/1000:.1f} mL')
print(sorted(set((k[0],k[1]) for k in bcc if k[0]!='wall')))
pickle.dump(dict(blocks=blocks,ifc=ifc,bcs=bcs,inout=inout),open('assembled.pkl','wb'))
