import numpy as np, os
def rec(f):
    n=np.frombuffer(f.read(4),'<i4')[0]
    d=f.read(n); n2=np.frombuffer(f.read(4),'<i4')[0]; assert n==n2,(n,n2)
    return d
def read_xyz_dbl(fn, nvar=3):
    f=open(fn,'rb'); nb=np.frombuffer(rec(f),'<i4')[0]
    dims=np.frombuffer(rec(f),'<i4').reshape(nb,-1)
    blocks=[]
    for b in range(nb):
        ni,nj,nk=dims[b,:3]; nv=dims[b,3] if dims.shape[1]>3 else nvar
        a=np.frombuffer(rec(f),'<f8').reshape(nv,nk,nj,ni)  # Fortran order: i fastest
        blocks.append(a.transpose(0,3,2,1))  # (nv,ni,nj,nk)
    return blocks
def read_xyzi_sgl(fn):
    f=open(fn,'rb'); nb=np.frombuffer(rec(f),'<i4')[0]
    dims=np.frombuffer(rec(f),'<i4').reshape(nb,3)
    blocks=[]
    for b in range(nb):
        ni,nj,nk=dims[b]; n=ni*nj*nk
        d=rec(f)
        xyz=np.frombuffer(d[:12*n],'<f4').reshape(3,nk,nj,ni).transpose(0,3,2,1)
        ib=np.frombuffer(d[12*n:],'<i4').reshape(nk,nj,ni).transpose(2,1,0)
        blocks.append((xyz,ib))
    return blocks
def read_puvwf_sgl(fn):
    f=open(fn,'rb'); nb=np.frombuffer(rec(f),'<i4')[0]
    dims=np.frombuffer(rec(f),'<i4').reshape(nb,3)
    blocks=[]
    for b in range(nb):
        ni,nj,nk=dims[b]; n=ni*nj*nk
        h=np.frombuffer(rec(f),'<f8'); d=rec(f)
        q=np.frombuffer(d,'<f4').reshape(5,nk,nj,ni).transpose(0,3,2,1)
        blocks.append(q)
    return blocks
