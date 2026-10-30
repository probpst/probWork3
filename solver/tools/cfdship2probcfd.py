#!/usr/bin/env python3
"""cfdship2probcfd.py - convert a CFDShip-Iowa example (binary multi-block
PLOT3D grid + CFDShip .bcs) into the probCFD input format:
    <case>.grd  ASCII multi-block PLOT3D (optionally coarsened by striding)
    <case>.bcs  probCFD boundary-condition file
Type map: 20 no-slip -> Dirichlet u=v=w=0 (wall, used by the overset assembly
to define the solid body); 10 inlet -> Dirichlet (1,0,0); 11/12 exit ->
p = 0 outlet; 13/28/40/61/9 farfield/zero-gradient/symmetry -> zero gradient
(61: normal velocity 0); 92 patched -> omitted (auto-detected b2b).
Usage: cfdship2probcfd.py grid.xyz|p3dxyzi case.bcs outcase [stride]
"""
import sys, numpy as np
def rec(f):
    n=np.frombuffer(f.read(4),'<i4')[0]; d=f.read(n); f.read(4); return d
def read_grid(fn):
    f=open(fn,'rb'); nb=np.frombuffer(rec(f),'<i4')[0]
    dims=np.frombuffer(rec(f),'<i4')
    if len(dims)==4*nb: dims=dims.reshape(nb,4)[:,:3]
    else: dims=dims.reshape(nb,3)
    blocks=[]
    for b in range(nb):
        ni,nj,nk=dims[b]; n=ni*nj*nk; d=rec(f)
        if len(d)==24*n or len(d)==32*n:   # doubles (3 or 4 vars)
            a=np.frombuffer(d[:24*n],'<f8').reshape(3,nk,nj,ni).transpose(0,3,2,1)
        else:                              # singles + int iblank
            a=np.frombuffer(d[:12*n],'<f4').reshape(3,nk,nj,ni).transpose(0,3,2,1).astype(np.float64)
        blocks.append(a)
    return blocks
def stride_idx(n,s):
    idx=list(range(0,n,s))
    if idx[-1]!=n-1: idx.append(n-1)
    return np.array(idx)
def read_bcs(fn):
    toks=[l.split('!')[0].split() for l in open(fn)]
    toks=[t for t in toks if t]
    p=0; blocks=[]
    while p<len(toks):
        nsurf=int(toks[p][0]); p+=1; surfs=[]
        for _ in range(nsurf):
            typ=int(toks[p][0]); d=int(toks[p+1][0])
            ir=list(map(int,toks[p+2][:2])); jr=list(map(int,toks[p+3][:2])); kr=list(map(int,toks[p+4][:2])); p+=6
            if typ==92: p+=5   # donor mesh, direction, 3 ranges
            surfs.append((typ,d,ir,jr,kr))
        blocks.append(surfs)
    return blocks
gridfile,bcsfile,out=sys.argv[1:4]; s=int(sys.argv[4]) if len(sys.argv)>4 else 1
G=read_grid(gridfile); B=read_bcs(bcsfile)
assert len(G)==len(B),(len(G),len(B))
newG=[]; maps=[]
for a in G:
    ii=stride_idx(a.shape[1],s); jj=stride_idx(a.shape[2],s); kk=stride_idx(a.shape[3],s)
    newG.append(a[:,ii][:,:,jj][:,:,:,kk]); maps.append((ii,jj,kk))
with open(out+'.grd','w') as f:
    f.write(f"{len(newG)}\n")
    for a in newG: f.write(f"{a.shape[1]} {a.shape[2]} {a.shape[3]}\n")
    for a in newG:
        for c in range(3):
            v=a[c].transpose(2,1,0).ravel()          # i fastest, then j, then k
            n4=(len(v)//4)*4
            if n4>0: np.savetxt(f, v[:n4].reshape(-1,4), fmt='%.16g')
            if n4<len(v): np.savetxt(f, v[n4:][None,:], fmt='%.16g')
def newidx(old,idxarr):
    # old 1-based node index -> new 1-based index of the nearest kept node
    j=np.searchsorted(idxarr, old-1)
    j=min(j,len(idxarr)-1)
    return int(j)+1
lines=[f"! probCFD bcs converted from {bcsfile} (stride {s}); type 92 (b2b) omitted", f"{len(newG)}"]
for b,(surfs,(ii,jj,kk)) in enumerate(zip(B,maps)):
    keep=[sf for sf in surfs if sf[0]!=92]
    lines.append(f"{len(keep)}     ! block {b+1}")
    for typ,d,ir,jr,kr in keep:
        r=[[newidx(ir[0],ii),newidx(ir[1],ii)],[newidx(jr[0],jj),newidx(jr[1],jj)],[newidx(kr[0],kk),newidx(kr[1],kk)]]
        lines.append(f"{d}     ! type {typ}")
        for a in range(3): lines.append(f"{min(r[a])} {max(r[a])}")
        ax=abs(d)-1
        if typ==20: fl=["1 0.0","1 0.0","1 0.0","0 0.0"]
        elif typ==10: fl=["1 1.0","1 0.0","1 0.0","0 0.0"]
        elif typ in (11,12): fl=["0 0.0","0 0.0","0 0.0","1 0.0"]
        elif typ==61:
            fl=["0 0.0","0 0.0","0 0.0","0 0.0"]; fl[ax]="1 0.0"
        else: fl=["0 0.0","0 0.0","0 0.0","0 0.0"]
        lines+=fl
open(out+'.bcs','w').write("\n".join(lines)+"\n")
print("blocks", [a.shape[1:] for a in newG])
