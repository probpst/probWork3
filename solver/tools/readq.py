import numpy as np
def rec(f):
    n=np.frombuffer(f.read(4),'<i4')[0]; d=f.read(n); f.read(4); return d
def read_x(fn):
    f=open(fn,'rb'); nb=np.frombuffer(rec(f),'<i4')[0]; dims=np.frombuffer(rec(f),'<i4').reshape(nb,3)
    out=[]
    for b in range(nb):
        ni,nj,nk=dims[b]; a=np.frombuffer(rec(f),'<f4').reshape(3,nk,nj,ni).transpose(0,3,2,1); out.append(a)
    return out
def read_q(fn):
    f=open(fn,'rb'); nb=np.frombuffer(rec(f),'<i4')[0]; dims=np.frombuffer(rec(f),'<i4').reshape(nb,4)
    out=[]
    for b in range(nb):
        ni,nj,nk,nv=dims[b]; a=np.frombuffer(rec(f),'<f4').reshape(nv,nk,nj,ni).transpose(0,3,2,1); out.append(a)
    return out
