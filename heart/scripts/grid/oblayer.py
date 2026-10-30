"""Global O-grid layer: every exterior patch of the core skeleton is extruded to the implicit wall.
Shared core nodes -> identical extrusion columns -> all BL/BL and BL/core interfaces point-matched."""
import numpy as np
from collections import defaultdict
import mb
from mb import Block, FACES, fax, rects, nkey

NU=3.30e-6
def h_first(U,D,yplus=25.0):
    Re=U*D/NU; Cf=0.079*Re**-0.25; ut=U*np.sqrt(Cf/2); return 2*yplus*NU/ut*1000.0   # mm, cell height (centre at y+)

def geo_s(L,h,m,gmin=0.80,gmax=1.25,iso=False):
    if iso:
        L=np.asarray(L); s=np.linspace(0,1,m+1)
        return np.broadcast_to(s,L.shape+(m+1,)).copy(),np.ones(L.shape)
    """node fractions s[0..m] (0=core,1=wall), first wall cell ~h, geometric; vectorised over nodes"""
    L=np.asarray(L); h=np.minimum(np.asarray(h),0.6*L)
    lo=np.full(L.shape,gmin); hi=np.full(L.shape,gmax)
    def tot(g): return np.where(abs(g-1)<1e-9,h*m,h*(1-g**m)/(1-g))
    for _ in range(60):
        g=0.5*(lo+hi); big=tot(g)>L; hi=np.where(big,g,hi); lo=np.where(big,lo,g)
    g=0.5*(lo+hi)
    c=g[...,None]**np.arange(m)      # wall->core relative sizes
    c=c/c.sum(-1,keepdims=True)
    fromwall=np.concatenate([np.zeros(L.shape+(1,)),np.cumsum(c,-1)],-1)   # 0..1 from wall
    return (1-fromwall)[...,::-1], g

def exterior_patches(blocks,part):
    P=[]
    for b,B in enumerate(blocks):
        for f in FACES:
            free=part[(b,f)]<0
            for t in B.tags:
                if t['face']==f: free[t['pr'][0]:t['pr'][1]+1,t['qr'][0]:t['qr'][1]+1]=False
            for pr,qr in rects(free): P.append((b,f,pr,qr))
    return P

import time as _time
_T0=[_time.time()]
def _tk(msg):
    print(f'      [{_time.time()-_T0[0]:6.1f}s] {msg}',flush=True)
def build(blocks,F,inout,m,h_of_comp,comp_names,clust=1.0,n_smooth=12,cosmin=0.5,verbose=True,cap_fac=2.5,n_tsmooth=20,n_untangle=150,jmin=0.02,dbg=None,n_wsmooth=0,walpha=0.05,rfloor=0.0,iso=True,curved=True,kk0=0.35):
    _tk('ifc,part=mb.connectivity(')
    ifc,part=mb.connectivity(blocks)
    patches=exterior_patches(blocks,part)
    _tk('# ---- surface node table')
    # ---- surface node table
    key2id={}; pos=[]; nrm=[]; quads=[]; pnodes=[]
    for (b,f,pr,qr) in patches:
        B=blocks[b]; a,hi=fax(f); idx=[slice(None)]*3
        idx[a]=-1 if hi else 0; Fc=B.X[tuple(idx)][pr[0]:pr[1]+2,qr[0]:qr[1]+2]
        idx[a]=-2 if hi else 1; G=B.X[tuple(idx)][pr[0]:pr[1]+2,qr[0]:qr[1]+2]
        dp=Fc[1:,:-1]-Fc[:-1,:-1]; dq=Fc[:-1,1:]-Fc[:-1,:-1]; cr=np.cross(dp,dq)
        if np.sum(cr*(G-Fc)[:-1,:-1])>0: Fc=Fc[:,::-1]           # make (p,q,out) right handed
        K=nkey(Fc); ids=np.zeros(Fc.shape[:2],int)
        for p in range(Fc.shape[0]):
            for q in range(Fc.shape[1]):
                k=tuple(K[p,q])
                if k not in key2id: key2id[k]=len(pos); pos.append(Fc[p,q]); nrm.append(np.zeros(3))
                ids[p,q]=key2id[k]
        pnodes.append((Fc,ids))
        # quad normals (outward): diag cross
        d1=Fc[1:,1:]-Fc[:-1,:-1]; d2=Fc[:-1,1:]-Fc[1:,:-1]; qn=np.cross(d1,d2)*0.5
        for p in range(ids.shape[0]-1):
            for q in range(ids.shape[1]-1):
                c=[ids[p,q],ids[p+1,q],ids[p+1,q+1],ids[p,q+1]]; quads.append(c)
                for n_ in c: nrm[n_]+=qn[p,q]
    pos=np.array(pos); nrm=np.array(nrm); N=len(pos)
    d0=nrm/np.linalg.norm(nrm,axis=1,keepdims=True)
    # ---- inlet/outlet rim constraints
    plane_n=np.zeros((N,3)); onplane=np.zeros(N,bool); rimR=np.zeros(N); rimC=np.zeros((N,3))
    for io in inout:
        dist=(pos-io['P'])@io['n']; rad=np.linalg.norm((pos-io['P'])-dist[:,None]*io['n'],axis=1)
        sel=(abs(dist)<1e-6)&(rad<io['r']*1.2)
        onplane|=sel; plane_n[sel]=io['n']; rimR[sel]=io['r']; rimC[sel]=io['P']
    def proj(d):
        d=d.copy(); d[onplane]-=np.einsum('ij,ij->i',d[onplane],plane_n[onplane])[:,None]*plane_n[onplane]
        return d/np.linalg.norm(d,axis=1,keepdims=True)
    # adjacency
    nb=defaultdict(set)
    for c in quads:
        for i in range(4): nb[c[i]].add(c[(i+1)%4]); nb[c[i]].add(c[(i-1)%4])
    nbl=[np.array(sorted(nb[i])) for i in range(N)]
    import scipy.sparse as sp
    rows=np.concatenate([[i]*len(l) for i,l in enumerate(nbl)]); cols=np.concatenate(nbl)
    Aop=sp.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(N,N)); Aop=sp.diags(1/np.asarray(Aop.sum(1)).ravel())@Aop
    _tk('d=proj(d0)')
    d=proj(d0)
    for it in range(n_smooth):
        avg=Aop@d; d=proj(0.5*d+0.5*avg)
        bad=np.einsum('ij,ij->i',d,d0)<cosmin
        if bad.any():
            dn=d.copy(); dn[bad]=0.5*d[bad]+0.5*d0[bad]; d=proj(dn)
    _tk('# ---- ray cast to wall')
    # ---- ray cast to wall
    t=np.zeros(N); x=pos.copy(); f=F(x)
    assert (f<0).all(), f'core node outside wall: {f.max()}'
    for it in range(400):
        act=f<-1e-4
        if not act.any(): break
        st=np.maximum(0.8*np.abs(f[act]),0.02); t[act]+=st; x[act]=pos[act]+t[act,None]*d[act]; f[act]=F(x[act])
    # bisection back to surface where overshot
    lo=np.maximum(t-np.maximum(0.8*np.abs(f),0.5),0); hi=t.copy()
    for it in range(40):
        mid=0.5*(lo+hi); fm=F(pos+mid[:,None]*d); out=fm>0; hi=np.where(out,mid,hi); lo=np.where(out,lo,mid)
    Wp=pos+hi[:,None]*d
    _tk('# ---- thickness field')
    # ---- thickness field: capped ray hit, min-smoothed, then untangled
    import scipy.sparse as sp
    rows=np.concatenate([[i]*len(l) for i,l in enumerate(nbl)]); cols=np.concatenate(nbl)
    A=sp.csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(N,N)); A=sp.diags(1/np.asarray(A.sum(1)).ravel())@A
    dist0=np.abs(F(pos)); T0=np.minimum(hi,cap_fac*dist0+1.0)
    # exact circles at inlets/outlets
    if onplane.any():
        v=pos[onplane]-rimC[onplane]; ax_=np.einsum('ij,ij->i',d[onplane],v)
        rr=np.linalg.norm(v,axis=1); T0[onplane]=np.sqrt(np.maximum(rimR[onplane]**2-rr**2+ax_**2,0))-ax_
    # ---- landing-point smoothing on the anatomical surface
    if n_wsmooth>0:
        def project(X,it=3,e=0.05):
            for _ in range(it):
                f=F(X); g=np.zeros_like(X)
                for c in range(3):
                    dx=np.zeros(3); dx[c]=e; g[:,c]=(F(X+dx)-f)/e
                X=X-(f/np.maximum((g*g).sum(1),1e-6))[:,None]*g
            return X
        def rims(X):
            if onplane.any():
                v=X[onplane]-rimC[onplane]; v-=np.einsum('ij,ij->i',v,plane_n[onplane])[:,None]*plane_n[onplane]
                X[onplane]=rimC[onplane]+rimR[onplane,None]*v/np.linalg.norm(v,axis=1,keepdims=True)
            return X
        W=pos+T0[:,None]*d; capped=hi>cap_fac*dist0+1.0
        W[capped]=project(W[capped],6); W=rims(W); anchor=W.copy()
        for it in range(n_wsmooth):
            Wn=(1-walpha)*(0.5*W+0.5*(A@W))+walpha*anchor
            Wn=rims(project(Wn,1))
            # visibility: segment p->W must stay inside
            ok=np.ones(N,bool)
            for t_ in (0.5,0.85):
                ok&=F(pos+t_*(Wn-pos))<0.0
            W=np.where(ok[:,None],Wn,W)
        W=rims(project(W,3)); vec=W-pos; T0=np.linalg.norm(vec,axis=1); d=vec/T0[:,None]
    _tk('rcap=np.ones(N); rho=np.o')
    nw=None
    if curved:
        W0=pos+T0[:,None]*d; e=0.05; f0=F(W0); gg=np.zeros_like(W0)
        for c in range(3):
            dx=np.zeros(3); dx[c]=e; gg[:,c]=(F(W0+dx)-f0)/e
        nw=gg/np.maximum(np.linalg.norm(gg,axis=1,keepdims=True),1e-9)
        nw=np.where((np.einsum('ij,ij->i',nw,d)>0.3)[:,None],nw,d)
        if onplane.any():   # keep inlet/outlet rim columns inside their plane
            v=nw[onplane]-np.einsum('ij,ij->i',nw[onplane],plane_n[onplane])[:,None]*plane_n[onplane]
            n2=np.linalg.norm(v,axis=1,keepdims=True)
            nw[onplane]=np.where(n2>1e-6,v/np.maximum(n2,1e-9),d[onplane])
    kkn=np.full(N,kk0); rcap=np.ones(N); rho=np.ones(N); T=T0.copy(); info_neg=-1
    _,which0=F(pos+T0[:,None]*d,which=True)
    # near-uniform layers: first cell = clust * (L/m)  (clust=1 -> uniform, <1 -> mild wall clustering)
    h0=clust*T0/m
    pid_of={id(ids):f'{blocks[b].name}.{f}{pr}{qr}' for ((b,f,pr,qr),(Fc,ids)) in zip(patches,pnodes)}
    def column(Fc,ids,Wp,S,nw,kk=None):
        """node positions along each column: straight, or arc-length-parameterised Hermite
        (tangent = core normal at the wall-layer start, wall normal at the wall)"""
        W=Wp[ids]; V=W-Fc
        if nw is None: return Fc[:,:,None,:]+S[...,None]*V[:,:,None,:]
        Tl=(np.linalg.norm(V,axis=-1)*(kkn[ids] if kk is None else kk))[...,None]
        P0=Fc[...,None,:]; P1=W[...,None,:]; M0=(Tl*d[ids])[...,None,:]; M1=(Tl*nw[ids])[...,None,:]
        t=np.linspace(0,1,41).reshape(1,1,41,1)
        h00=2*t**3-3*t**2+1; h10=t**3-2*t**2+t; h01=-2*t**3+3*t**2; h11=t**3-t**2
        dense=h00*P0+h10*M0+h01*P1+h11*M1                      # (p,q,41,3)
        seg=np.linalg.norm(np.diff(dense,axis=2),axis=-1)
        cum=np.concatenate([np.zeros(seg.shape[:2]+(1,)),np.cumsum(seg,axis=2)],axis=2)
        tot=cum[...,-1:]; target=S*np.maximum(tot,1e-12)        # (p,q,m+1)
        idx=np.clip(np.searchsorted(np.arange(1),0)+ (cum[...,None,:]<=target[...,:,None]).sum(-1)-1,0,39)
        c0=np.take_along_axis(cum,idx,axis=2); c1=np.take_along_axis(cum,idx+1,axis=2)
        w=np.where(c1>c0,(target-c0)/np.maximum(c1-c0,1e-12),0.0)[...,None]
        p0=np.take_along_axis(dense,idx[...,None].repeat(3,-1),axis=2)
        p1=np.take_along_axis(dense,(idx+1)[...,None].repeat(3,-1),axis=2)
        return p0+(p1-p0)*w
    patch_bad=[None]*len(pnodes); Tprev=np.full(N,np.nan); kprev=np.full(N,np.nan)
    for it in range(n_untangle):
        Wp=pos+T[:,None]*d; L=T
        changed=(~(np.abs(T-Tprev)<1e-12))|(~(np.abs(kkn-kprev)<1e-12)); Tprev=T.copy(); kprev=kkn.copy()
        bad=np.zeros(N,bool); nneg=0
        for pi,(Fc,ids) in enumerate(pnodes):
            if patch_bad[pi] is not None and not changed[ids].any():
                c,nn=patch_bad[pi]
            else:
                S,g=geo_s(L[ids],h0[ids],m,iso=iso)
                X=column(Fc,ids,Wp,S,nw)
                J=mb.hex_jac(X); c=np.argwhere(J<=jmin); nn=int((J<=0).sum()); patch_bad[pi]=(c,nn)
                if it==0 and dbg is not None and nn: dbg.append((pid_of[id(ids)],nn,J.shape))
            if len(c):
                nneg+=nn
                for dp_,dq_ in [(0,0),(1,0),(0,1),(1,1)]: bad[ids[c[:,0]+dp_,c[:,1]+dq_]]=True
        if verbose and it%10==0: print(f'    untangle it {it}: neg cells {nneg}, flagged nodes {bad.sum()}',flush=True)
        info_neg=nneg
        if not bad.any(): break
        kkn[bad]*=0.55
        rcap[bad]=np.maximum(rcap[bad]*0.94,rfloor)
        rho=np.minimum(rho,rcap)
        for k in range(4): rho=np.minimum(0.5*rho+0.5*(A@rho),rcap)
        T=rho*T0
    info_capped=int((hi>cap_fac*dist0+1.0).sum())
    dev=F(Wp)
    L=np.linalg.norm(Wp-pos,axis=1)
    _,which=F(Wp,which=True)
    h=clust*L/m
    S,g=geo_s(L,h,m,iso=iso)
    info=dict(neg=info_neg,rho=T/np.maximum(T0,1e-9),L=L,h=np.minimum(h,0.6*L),g=g,which=which,Wp=Wp,pos=pos,d=d,n_patch=len(patches),dev=dev,T0=T0)
    _tk('# ---- BL blocks')
    # ---- BL blocks
    out=[]
    for n,((b,f,pr,qr),(Fc,ids)) in enumerate(zip(patches,pnodes)):
        X=column(Fc,ids,Wp,S[ids],nw)
        B=Block(f'{blocks[b].name}.{f}.bl{n}',X,blocks[b].dom,'bl',blocks[b].comp)
        out.append(B)
    if verbose: print(f'  patches {len(patches)}, surface nodes {N}, capped rays {info_capped}, L {L.min():.2f}-{L.max():.2f} mm, wall dev from anatomy: mean {np.abs(dev).mean():.2f} p95 {np.percentile(np.abs(dev),95):.2f} max {np.abs(dev).max():.2f} mm')
    return out,info
