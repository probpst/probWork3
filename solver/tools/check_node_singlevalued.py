#!/usr/bin/env python3
"""Check that node values of a multi-block PLOT3D .q file are single-valued
on the block-boundary nodes shared by several blocks (junction lines, O-cuts).
Usage: python3 check_node_singlevalued.py Case.x Case_00000100.q"""
import sys, numpy as np
from scipy.io import FortranFile
def readx(fn):
    f=FortranFile(fn,'r'); nb=f.read_ints()[0]; dims=f.read_ints().reshape(nb,3); out=[]
    for b in range(nb):
        ni,nj,nk=dims[b]; out.append(f.read_reals('f4').reshape(3,nk,nj,ni))
    return out
def readq(fn):
    f=FortranFile(fn,'r'); nb=f.read_ints()[0]; dims=f.read_ints().reshape(nb,4); out=[]
    for b in range(nb):
        ni,nj,nk,nv=dims[b]; out.append(f.read_reals('f4').reshape(nv,nk,nj,ni))
    return out
X=readx(sys.argv[1]); Q=readq(sys.argv[2])
pts=[];vals=[]
for b in range(len(X)):
    x=X[b]; q=Q[b]; nk,nj,ni=x.shape[1:]
    m=np.zeros((nk,nj,ni),bool); m[0]=m[-1]=1; m[:,0]=m[:,-1]=1; m[:,:,0]=m[:,:,-1]=1
    pts.append(np.stack([x[c][m] for c in range(3)],1)); vals.append(np.stack([q[c][m] for c in range(q.shape[0])],1))
P=np.concatenate(pts); V=np.concatenate(vals)
u,inv,cnt=np.unique(np.round(P,6),axis=0,return_inverse=True,return_counts=True); inv=inv.ravel()
maxd=np.zeros(V.shape[1])
for g in np.where(cnt>1)[0]:
    s=inv==g; maxd=np.maximum(maxd,V[s].max(0)-V[s].min(0))
print('boundary nodes',len(P),'unique',len(u),'shared',(cnt>1).sum())
print('max difference between block copies per variable:',maxd)
