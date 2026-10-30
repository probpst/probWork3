import numpy as np, sys
def rw(fn):
    with open(fn,'rb') as f: raw=f.read()
    p=0
    def rec():
        nonlocal p
        n=int(np.frombuffer(raw,'<i4',1,p)[0]); p+=4; d=raw[p:p+n]; p+=n+4; return d
    nb=int(np.frombuffer(rec(),'<i4')[0]); dims=np.frombuffer(rec(),'<i4').reshape(nb,4); out=[]
    for b in range(nb):
        ni,nj,nk,nv=dims[b]; a=np.frombuffer(rec(),'<f4').reshape(nv,nk,nj,ni).copy(); out.append(a)
    return nb,dims,out
def write(fn,nb,dims,blocks):
    with open(fn,'wb') as f:
        def rec(b): f.write(np.int32(len(b)).tobytes()); f.write(b); f.write(np.int32(len(b)).tobytes())
        rec(np.int32(nb).tobytes()); rec(dims.astype('<i4').tobytes())
        for a in blocks: rec(a.astype('<f4').tobytes())
def apply(a):
    ib=np.rint(a[0]).astype(int); bv=ib.astype('f4')
    act=(ib==1); fr=(ib==0); nb=np.zeros_like(act)
    nb[1:,:,:]|=act[:-1,:,:]; nb[:-1,:,:]|=act[1:,:,:]; nb[:,1:,:]|=act[:,:-1,:]; nb[:,:-1,:]|=act[:,1:,:]; nb[:,:,1:]|=act[:,:,:-1]; nb[:,:,:-1]|=act[:,:,1:]
    bv[fr&nb]=0.5
    h=(bv==0.5); z=(bv==0.0); nb2=np.zeros_like(h)
    nb2[1:,:,:]|=h[:-1,:,:]; nb2[:-1,:,:]|=h[1:,:,:]; nb2[:,1:,:]|=h[:,:-1,:]; nb2[:,:-1,:]|=h[:,1:,:]; nb2[:,:,1:]|=h[:,:,:-1]; nb2[:,:,:-1]|=h[:,:,1:]
    bv[z&nb2]=0.25
    h=(bv==0.25); z=(bv==0.0); nb3=np.zeros_like(h)
    nb3[1:,:,:]|=h[:-1,:,:]; nb3[:-1,:,:]|=h[1:,:,:]; nb3[:,1:,:]|=h[:,:-1,:]; nb3[:,:-1,:]|=h[:,1:,:]; nb3[:,:,1:]|=h[:,:,:-1]; nb3[:,:,:-1]|=h[:,:,1:]
    bv[z&nb3]=0.125; a[0]=bv; return a
if __name__=='__main__':
    fn=sys.argv[1]; nb,dims,bl=rw(fn); bl=[apply(a) for a in bl]; write(fn,nb,dims,bl); print(fn,'updated; 0.5 nodes:',sum(int((a[0]==0.5).sum()) for a in bl),' 0.25 nodes:',sum(int((a[0]==0.25).sum()) for a in bl),' 0.125 nodes:',sum(int((a[0]==0.125).sum()) for a in bl))
