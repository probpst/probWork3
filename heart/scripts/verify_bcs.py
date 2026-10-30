"""Independent verification of heart_grid.grd + heart_grid.bcs (reads only the written files)."""
import numpy as np, struct, sys
def read_grd(fn):
    tok=open(fn).read().split(); nb=int(tok[0]); dims=[tuple(int(t) for t in tok[1+3*b:4+3*b]) for b in range(nb)]
    arr=np.array(tok[1+3*nb:],dtype=float); G=[]; q=0
    for (ni,nj,nk) in dims:
        n=ni*nj*nk; G.append(np.stack([arr[q+c*n:q+(c+1)*n].reshape(nk,nj,ni).transpose(2,1,0) for c in range(3)],-1)); q+=3*n
    assert q==len(arr); return dims,G
def read_bcs(fn,nb):
    L=[l.split('!')[0].split() for l in open(fn).read().replace('\r','').split('\n')]
    L=[[int(x) for x in l] for l in L if l]; k=0; BC=[]
    for b in range(nb):
        ns=L[k][0];k+=1; S=[]
        for s in range(ns):
            e=dict(type=L[k][0],dir=L[k+1][0],rng=[L[k+2],L[k+3],L[k+4]]); k+=6
            if e['type']==92: e.update(donor=L[k][0],ddir=L[k+1][0],drng=[L[k+2],L[k+3],L[k+4]]); k+=5
            S.append(e)
        BC.append(S)
    assert k==len(L); return BC
def face(X,rng,dirn):
    sl=[]
    for s,e in rng:
        st=1 if e>=s else -1; sl.append(slice(s-1,(e-1+st) if (e-1+st)>=0 else None,st))
    return np.squeeze(X[tuple(sl)],axis=abs(dirn)-1)
dims,G=read_grd(sys.argv[1]); BC=read_bcs(sys.argv[2],len(dims))
errs=[]; tr=0; bad=0; cov_ok=True; types={}
for b,S in enumerate(BC):
    ni,nj,nk=dims[b]; cnt={d:np.zeros([ (ni,nj,nk)[a]-1 for a in range(3) if a!=abs(d)-1],int) for d in (1,-1,2,-2,3,-3)}
    for e in S:
        types[e['type']]=types.get(e['type'],0)+1
        a=abs(e['dir'])-1; fixed=e['rng'][a][0]; expect=1 if e['dir']>0 else (ni,nj,nk)[a]
        assert e['rng'][a][0]==e['rng'][a][1]==expect, ('dir/range mismatch',b+1,e)
        r=[x for i,x in enumerate(e['rng']) if i!=a]
        cnt[e['dir']][r[0][0]-1:r[0][1]-1,r[1][0]-1:r[1][1]-1]+=1
        if e['type']==92:
            R=face(G[b],e['rng'],e['dir']); D=face(G[e['donor']-1],e['drng'],e['ddir'])
            if R.shape==D.shape: errs.append(np.abs(R-D).max())
            elif R.shape==D.transpose(1,0,2).shape: errs.append(np.abs(R-D.transpose(1,0,2)).max()); tr+=1
            else: bad+=1
    for d,c in cnt.items():
        if not (c==1).all(): cov_ok=False; print('coverage problem block',b+1,'face',d,'min',c.min(),'max',c.max())
print(f'blocks {len(dims)}  cells {sum((a-1)*(b-1)*(c-1) for a,b,c in dims)}  BC entries by type {types}')
print(f'type-92 entries {len(errs)+bad}: max point mismatch {max(errs):.3e} mm, transposed {tr}, shape errors {bad}; every face cell covered exactly once: {cov_ok}')
