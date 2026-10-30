#!/usr/bin/env python3
"""drop given blocks (1-based, current numbering) from <case>.grd / .bcs / csv"""
import numpy as np, sys, re
case, cdir, drop = sys.argv[1], sys.argv[2], [int(v) for v in sys.argv[3].split(',')]
txt=open(f'{cdir}/{case}.grd').read().split(); nb=int(txt[0]); dims=np.array(txt[1:1+3*nb],int).reshape(nb,3); p=1+3*nb; B=[]
for b in range(nb):
    ni,nj,nk=dims[b]; n=ni*nj*nk; B.append(np.array(txt[p:p+3*n],float).reshape(3,nk,nj,ni).transpose(0,3,2,1)); p+=3*n
bcs=open(f'{cdir}/{case}.bcs').read()
head=bcs[:bcs.index('! ==================== block 1:')]
secs=[s for s in re.split(r'(?m)^(?=! ==================== block \d+:)', bcs[bcs.index('! ==================== block 1:'):]) if s.strip()]
assert len(secs)==nb, (len(secs),nb)
keep=[b for b in range(nb) if b+1 not in drop]
B2=[B[b] for b in keep]; out=[]
for newb,b in enumerate(keep):
    sec=re.sub(r'^! ==================== block \d+:', '! ==================== block %d:'%(newb+1), secs[b], count=1)
    sec=re.sub(r'(! number of boundary surfaces of block )\d+', r'\g<1>%d'%(newb+1), sec)
    out.append(sec)
head=re.sub(r'^(\s*)%d(\s+! number of grid blocks)'%nb, r'\g<1>%d\g<2>'%len(keep), head, flags=re.M)
open(f'{cdir}/{case}.bcs','w').write(head+''.join(out))
with open(f'{cdir}/{case}.grd','w') as f:
    f.write("%8d\n"%len(B2))
    for a in B2: f.write("%8d%8d%8d\n"%a.shape[1:])
    for a in B2:
        for c in range(3):
            v=a[c].transpose(2,1,0).ravel(); n4=(len(v)//4)*4
            if n4: np.savetxt(f, v[:n4].reshape(-1,4), fmt='%22.15e')
            if n4<len(v): np.savetxt(f, v[n4:][None,:], fmt='%22.15e')
n=len(B2); big=[b for b in range(n) if (B2[b].shape[1]-1)*(B2[b].shape[2]-1)*(B2[b].shape[3]-1)>200000]; sp=[4 if b in big else 1 for b in range(n)]
open(f'{cdir}/Input_polyMesh_probCFD.csv','w').write("\n".join(["probCFD polyMesh Input"+"".join(",block-%d"%(b+1) for b in range(n)),"i_split,"+",".join(str(v) for v in sp),"j_split,"+",".join("1" for _ in range(n)),"k_split,"+",".join(("2" if b in big else "1") for b in range(n)),"n_of_BLK,"+",".join(str(sp[b]*(2 if b in big else 1)) for b in range(n))])+"\n")
print(case, 'dropped', drop, '->', n, 'blocks; grid blocks header:', head.strip().split('\n')[-1].strip()[:40])
