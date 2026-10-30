#!/usr/bin/env python3
"""Per case: (1) resample the wall-normal (j) lines of every j-min wall block to a
COMMON point count (the case maximum) with the geometric law from each line's own
first-spacing fraction - identical lines on both sides of a block-to-block face
give identical nodes; (2) reorder blocks so that added cap/lid blocks precede the
refinement/background blocks (background last); rewrite grd, bcs, csv."""
import numpy as np, sys, re, os
from scipy.interpolate import CubicSpline
case, cdir, order = sys.argv[1], sys.argv[2], sys.argv[3]      # order: e.g. "1-72,74-76,73" (1-based, original numbering)
txt=open(f'{cdir}/{case}.grd').read().split(); nb=int(txt[0]); dims=np.array(txt[1:1+3*nb],int).reshape(nb,3); p=1+3*nb; B=[]
for b in range(nb):
    ni,nj,nk=dims[b]; n=ni*nj*nk; a=np.array(txt[p:p+3*n],float).reshape(3,nk,nj,ni).transpose(0,3,2,1); p+=3*n; B.append(a)
# ---- bcs: split into block sections ----
bcs=open(f'{cdir}/{case}.bcs').read()
head=bcs[:bcs.index('! ==================== block 1:')]
secs=re.split(r'(?m)^(?=! ==================== block \d+:)', bcs[bcs.index('! ==================== block 1:'):])
secs=[s for s in secs if s.strip()]
assert len(secs)==nb, (len(secs), nb)
def jmin_wall(sec):
    lines=[l for l in sec.split('\n') if l.strip() and not l.strip().startswith('!')]
    for n,l in enumerate(lines):
        v=l.split('!')[0].split()
        if len(v)==1 and v[0]=='2' and n+6<len(lines):        # direction line of a j-min face
            f=[lines[n+4+m].split('!')[0].split() for m in range(3)]
            if all(len(x)>=2 and x[0]=='1' and float(x[1])==0.0 for x in f): return True
    return False
def solve_r(g1,n):
    lo=np.full_like(g1,1.0+1e-12); hi=np.full_like(g1,6.0)
    for _ in range(300):
        r=0.5*(lo+hi); s=g1*(r**n-1)/(r-1); hi=np.where(s>1,r,hi); lo=np.where(s>1,lo,r)
    return 0.5*(lo+hi)
wall=[b for b in range(nb) if jmin_wall(secs[b]) and B[b].shape[2]>=10]
m=max(B[b].shape[2] for b in wall) if wall else 0
print(case, 'j-min wall blocks resampled to nj =', m, ':', [b+1 for b in wall])
for b in wall:
    a=B[b]; ni,nj,nk=a.shape[1:]
    d=np.sqrt(((a[:,:,1:,:]-a[:,:,:-1,:])**2).sum(0)); t=np.concatenate([np.zeros((ni,1,nk)),np.cumsum(d,axis=1)],axis=1)
    L=t[:,-1,:]; g1=np.clip(d[:,0,:]/L,1e-12,0.95/(m-1)); r=solve_r(g1,m-1)
    frac=np.zeros((ni,m,nk))
    for k in range(1,m): frac[:,k,:]=g1*(r**k-1)/(r-1)
    frac[:,-1,:]=1.0
    new=np.empty((3,ni,m,nk))
    for i in range(ni):
        for kk in range(nk):
            tt=np.maximum.accumulate(t[i,:,kk].copy()); tt[1:]=np.where(np.diff(tt)<=0,tt[:-1]+1e-14,tt[1:])
            tg=frac[i,:,kk]*L[i,kk]
            for c in range(3): new[c,i,:,kk]=CubicSpline(tt,a[c,i,:,kk])(tg)
            new[:,i,0,kk]=a[:,i,0,kk]; new[:,i,-1,kk]=a[:,i,-1,kk]
    B[b]=new
# ---- reorder ----
newidx=[]
for part in order.split(','):
    if '-' in part: lo,hi=part.split('-'); newidx+=list(range(int(lo)-1,int(hi)))
    else: newidx.append(int(part)-1)
assert sorted(newidx)==list(range(nb)), newidx
B2=[B[b] for b in newidx]; secs2=[secs[b] for b in newidx]
out=[]
for newb,sec in enumerate(secs2):
    sec=re.sub(r'^! ==================== block \d+:', '! ==================== block %d:'%(newb+1), sec, count=1)
    sec=re.sub(r'(! number of boundary surfaces of block )\d+', r'\g<1>%d'%(newb+1), sec)
    # update the "x x x nodes" text of resampled blocks
    ni,nj,nk=B2[newb].shape[1:]
    sec=re.sub(r'(\d+) x (\d+) x (\d+) nodes', lambda mm: '%d x %d x %d nodes'%(ni,nj,nk) if int(mm.group(1))==ni and int(mm.group(3))==nk else mm.group(0), sec, count=1)
    # fix j-ranges of k/i-face records (1..nj) inside this block
    lines=sec.split('\n')
    for i,l in enumerate(lines):
        if 'j-range (nodes)' in l:
            v=l.split('!')[0].split()
            if len(v)>=2 and v[0]=='1' and int(v[1])>1: lines[i]="    %4d %4d           ! j-range (nodes)"%(1,nj)
            elif len(v)>=2 and v[0]==v[1] and int(v[0])>1: lines[i]="    %4d %4d           ! j-range (nodes)"%(nj,nj)
    out.append('\n'.join(lines))
open(f'{cdir}/{case}.bcs','w').write(head+''.join(out) if head.endswith('\n') else head+'\n'+''.join(out))
with open(f'{cdir}/{case}.grd','w') as f:
    f.write("%8d\n"%nb)
    for a in B2: f.write("%8d%8d%8d\n"%a.shape[1:])
    for a in B2:
        for c in range(3):
            v=a[c].transpose(2,1,0).ravel(); n4=(len(v)//4)*4
            if n4: np.savetxt(f, v[:n4].reshape(-1,4), fmt='%22.15e')
            if n4<len(v): np.savetxt(f, v[n4:][None,:], fmt='%22.15e')
hdr="probCFD polyMesh Input"+"".join(",block-%d"%(b+1) for b in range(nb))
big=[b for b in range(nb) if (B2[b].shape[1]-1)*(B2[b].shape[2]-1)*(B2[b].shape[3]-1)>200000]
sp=[4 if b in big else 1 for b in range(nb)]
open(f'{cdir}/Input_polyMesh_probCFD.csv','w').write("\n".join([hdr,"i_split,"+",".join(str(v) for v in sp),
 "j_split,"+",".join("1" for _ in range(nb)),"k_split,"+",".join(("2" if b in big else "1") for b in range(nb)),
 "n_of_BLK,"+",".join(str(sp[b]*(2 if b in big else 1)) for b in range(nb))])+"\n")
print('  written: %d blocks, cells %.2f M, order %s'%(nb, sum((a.shape[1]-1)*(a.shape[2]-1)*(a.shape[3]-1) for a in B2)/1e6, order))
