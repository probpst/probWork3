"""Core skeleton: chamber cores + swept tube cores, all point-matched by construction."""
import numpy as np, params as P
from paths import Path, unit
from mb import Block, fax, inplane

def smooth(t): t=np.clip(t,0,1); return t*t*(3-2*t)

def face_normal(B,f,pr=None,qr=None):
    """mean outward unit normal of (sub)face; pr,qr node ranges"""
    a,hi=fax(f); idx=[slice(None)]*3
    idx[a]=-1 if hi else 0; F=B.X[tuple(idx)]
    idx[a]=-2 if hi else 1; G=B.X[tuple(idx)]
    if pr: F=F[pr[0]:pr[1]+1,qr[0]:qr[1]+1]; G=G[pr[0]:pr[1]+1,qr[0]:qr[1]+1]
    dp=(F[1:,:-1]-F[:-1,:-1]); dq=(F[:-1,1:]-F[:-1,:-1]); n=np.cross(dp,dq)
    sgn=np.sign(np.einsum('pqi,pqi->pq',n,(G-F)[:-1,:-1]))
    n=(n*(-sgn)[...,None]).sum((0,1))
    return unit(n)

def sub(B,f,pr,qr):
    return B.face(f)[pr[0]:pr[1]+1,qr[0]:qr[1]+1].copy()

def orient_rh(S,T):
    """flip q so (dp x dq).T>0"""
    dp=S[-1,S.shape[1]//2]-S[0,S.shape[1]//2]; dq=S[S.shape[0]//2,-1]-S[S.shape[0]//2,0]
    return S if np.cross(dp,dq)@T>0 else S[:,::-1]

def rot(v,axis,th):
    return v*np.cos(th)+np.cross(axis,v)*np.sin(th)+axis*(axis@v)*(1-np.cos(th))

def stations(path,n1,n2,cf,ds_ref,lo=0.35,hi=2.0):
    """arc-length stations whose spacing follows the local cross-section spacing (near-cubic cells)"""
    s=path.s; tgt=np.clip(2*cf*path.r/max(n1-1,n2-1)*1.0,lo*ds_ref,hi*ds_ref)
    u=np.concatenate([[0],np.cumsum(np.diff(s)/(0.5*(tgt[1:]+tgt[:-1])))])
    n=max(int(round(u[-1])),4)
    return np.interp(np.linspace(0,u[-1],n+1),u,s)

SQ=0.30
def tube(name,path,S0,dom,comp,cf=0.5,ds=1.0,S1=None,twist=None,nst=None):
    """sweep a rectangular-section core along path. S0 (np,nq,3) start nodes; S1 optional end nodes.
    returns Block with i=path, j,k = section."""
    L=path.L
    if nst: s=np.linspace(0,L,nst+1)
    else: s=stations(path,S0.shape[0],S0.shape[1],cf,ds)
    n=len(s)-1
    C,T,r=path.at(s)
    S0=orient_rh(S0,T[0]); npp,nq=S0.shape[:2]
    mid0=nq//2; f1=S0[-1,mid0]-S0[0,mid0]; f1=unit(f1-(f1@T[0])*T[0])
    F1=[f1]
    for i in range(1,n+1):
        v=F1[-1]-(F1[-1]@T[i])*T[i]; F1.append(unit(v))
    F1=np.array(F1)
    th=np.zeros(n+1)
    if twist is not None: th+=twist(s)
    a0=0.5*np.linalg.norm(S0[-1,mid0]-S0[0,mid0]); b0=0.5*np.linalg.norm(S0[npp//2,-1]-S0[npp//2,0])
    lam=min(1.3*r[0],0.45*L); w0=1-smooth(s/lam); w1=np.zeros_like(s)
    if S1 is not None:
        assert S1.shape[:2]==(npp,nq), (name,S1.shape,S0.shape)
        cands=[S1,S1[::-1,::-1],S1[::-1,:],S1[:,::-1]]
        cands=[c for c in cands if (lambda d:np.cross(d[0],d[1])@T[-1]>0)((c[-1,nq//2]-c[0,nq//2],c[npp//2,-1]-c[npp//2,0]))]
        fe=F1[-1]; best=None
        for c in cands:
            g=c[-1,nq//2]-c[0,nq//2]; g=unit(g-(g@T[-1])*T[-1])
            ang=np.arctan2(np.cross(fe,g)@T[-1],fe@g)
            if best is None or abs(ang)<abs(best[1]): best=(c,ang)
        S1,dth=best; th+=dth*smooth(s/L)
        a1=0.5*np.linalg.norm(S1[-1,nq//2]-S1[0,nq//2]); b1=0.5*np.linalg.norm(S1[npp//2,-1]-S1[npp//2,0])
        lam1=min(1.6*r[-1],0.45*L); w1=smooth((s-(L-lam1))/lam1); w0=1-smooth(s/min(lam,0.45*L))
    F1=np.array([rot(F1[i],T[i],th[i]) for i in range(n+1)]); F2=np.cross(T,F1)
    am=cf*r; a=w0*a0+(1-w0)*am; b=w0*b0+(1-w0)*am
    if S1 is not None: a=a*(1-w1)+w1*a1; b=b*(1-w1)+w1*b1
    u=np.linspace(-1,1,npp); v=np.linspace(-1,1,nq)
    U,V=np.meshgrid(u,v,indexing='ij')
    # rounded-square (squircle) core section: lowers skew on the butterfly diagonals
    SU=U*np.sqrt(np.clip(1-V*V/2,0,None)); SV=V*np.sqrt(np.clip(1-U*U/2,0,None))
    U=(1-SQ)*U+SQ*SU; V=(1-SQ)*V+SQ*SV
    I=C[:,None,None,:]+(a[:,None,None,None]*U[None,...,None])*F1[:,None,None,:]+(b[:,None,None,None]*V[None,...,None])*F2[:,None,None,:]
    X=I+w0[:,None,None,None]*(S0-I[0])[None]
    if S1 is not None: X=X+w1[:,None,None,None]*(S1-I[-1])[None]
    X[0]=S0
    if S1 is not None: X[-1]=S1
    smooth2d(X)
    B=Block(name,X,dom,'core',comp); B.path=path; B.s=s; B.F1=F1; B.F2=F2; B.hs=(a,b)
    return B

def box(name,c,E,h,n,dom,comp):
    """E columns e1,e2,e3; h half sizes; n cells"""
    u=[np.linspace(-1,1,k+1) for k in n]
    X=c+u[0][:,None,None,None]*h[0]*E[:,0]+u[1][None,:,None,None]*h[1]*E[:,1]+u[2][None,None,:,None]*h[2]*E[:,2]
    B=Block(name,X,dom,'core',comp); B.c=c; B.E=E; B.h=h; return B

def ext(Rv,e): return np.sqrt(((Rv*e)**2).sum())   # ellipsoid support extent (axis aligned)

def pick_range(B,f,P0,halfw,other=None):
    """node range on face f of box-like block around point P0 with physical half width"""
    F=B.face(f); c=F.reshape(-1,3); d=np.linalg.norm(c-P0,axis=1); p,q=np.unravel_index(d.argmin(),F.shape[:2])
    sp=np.linalg.norm(np.diff(F[:,q],axis=0),axis=1).mean(); sq=np.linalg.norm(np.diff(F[p,:],axis=0),axis=1).mean()
    hp=max(2,int(round(halfw/sp))); hq=max(2,int(round(halfw/sq)))
    p=int(np.clip(p,hp+1,F.shape[0]-hp-2)); q=int(np.clip(q,hq+1,F.shape[1]-hq-2))
    return (p-hp,p+hp),(q-hq,q+hq)

def smooth2d(X,it=40,keep_ends=True):
    """Laplacian smoothing of each i-section's interior nodes; section boundary and the
    first/last sections (which are point-matched to neighbouring blocks) stay fixed"""
    n=X.shape[0]; w=np.ones(n)
    if keep_ends:
        r=min(4,max(1,n//6)); ramp=np.linspace(0,1,r+2)[1:-1]
        w[0]=0; w[-1]=0; w[1:1+len(ramp)]=ramp; w[-1-len(ramp):-1]=ramp[::-1]
    w=w[:,None,None,None]
    for _ in range(it):
        avg=0.25*(X[:,:-2,1:-1]+X[:,2:,1:-1]+X[:,1:-1,:-2]+X[:,1:-1,2:])
        X[:,1:-1,1:-1]=(1-w)*X[:,1:-1,1:-1]+w*avg
    return X

def smooth3d(X,it=60):
    for _ in range(it):
        X[1:-1,1:-1,1:-1]=(X[:-2,1:-1,1:-1]+X[2:,1:-1,1:-1]+X[1:-1,:-2,1:-1]+X[1:-1,2:,1:-1]+X[1:-1,1:-1,:-2]+X[1:-1,1:-1,2:])/6
    return X

def sbox(name,c,E,R,n,dom,comp,beta=0.62,gamma=0.7):
    """spherified box core inside an ellipsoid with frame-axis extents R (3,)"""
    u,v,w=np.meshgrid(*[np.linspace(-1,1,k+1) for k in n],indexing='ij')
    su=u*np.sqrt(np.clip(1-v*v/2-w*w/2+v*v*w*w/3,0,None))
    sv=v*np.sqrt(np.clip(1-u*u/2-w*w/2+u*u*w*w/3,0,None))
    sw=w*np.sqrt(np.clip(1-u*u/2-v*v/2+u*u*v*v/3,0,None))
    x=(1-gamma)*u+gamma*su; y=(1-gamma)*v+gamma*sv; z=(1-gamma)*w+gamma*sw
    X=c+(beta*R[0]*x)[...,None]*E[:,0]+(beta*R[1]*y)[...,None]*E[:,1]+(beta*R[2]*z)[...,None]*E[:,2]
    smooth3d(X)
    B=Block(name,X,dom,'core',comp); B.c=c; B.E=E; B.h=beta*np.asarray(R); return B
