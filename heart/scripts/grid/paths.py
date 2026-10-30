"""Centerline paths (dense polylines with radius) and implicit (signed-distance-like) wall functions."""
import numpy as np
from scipy.spatial import cKDTree

def unit(v): v=np.asarray(v,float); return v/np.linalg.norm(v)

def catmull(ctrl, n_per=40):
    """centripetal Catmull-Rom through ctrl points -> dense points and param of ctrl pts"""
    P=np.asarray(ctrl,float)
    Q=np.vstack([2*P[0]-P[1],P,2*P[-1]-P[-2]])
    out=[]
    for i in range(1,len(Q)-2):
        p0,p1,p2,p3=Q[i-1],Q[i],Q[i+1],Q[i+2]
        t0=0;t1=t0+np.linalg.norm(p1-p0)**0.5+1e-9;t2=t1+np.linalg.norm(p2-p1)**0.5+1e-9;t3=t2+np.linalg.norm(p3-p2)**0.5+1e-9
        t=np.linspace(t1,t2,n_per,endpoint=False)[:,None]
        A1=(t1-t)/(t1-t0)*p0+(t-t0)/(t1-t0)*p1;A2=(t2-t)/(t2-t1)*p1+(t-t1)/(t2-t1)*p2;A3=(t3-t)/(t3-t2)*p2+(t-t2)/(t3-t2)*p3
        B1=(t2-t)/(t2-t0)*A1+(t-t0)/(t2-t0)*A2;B2=(t3-t)/(t3-t1)*A2+(t-t1)/(t3-t1)*A3
        out.append((t2-t)/(t2-t1)*B1+(t-t1)/(t2-t1)*B2)
    out.append(P[-1][None]); return np.vstack(out)

class Path:
    """dense polyline, arclength param, radius, tangent"""
    def __init__(self, ctrl, radii, ds=0.25):
        D=catmull(ctrl, 60)
        seg=np.linalg.norm(np.diff(D,axis=0),axis=1); s=np.r_[0,np.cumsum(seg)]
        # radius at ctrl points -> interpolate by nearest arclength of ctrl pts
        cs=[s[np.argmin(np.linalg.norm(D-c,axis=1))] for c in np.asarray(ctrl,float)]
        L=s[-1]; n=max(int(L/ds),8); self.s=np.linspace(0,L,n+1)
        self.P=np.stack([np.interp(self.s,s,D[:,c]) for c in range(3)],1)
        self.r=np.interp(self.s,cs,radii); self.L=L
        T=np.gradient(self.P,axis=0); self.T=T/np.linalg.norm(T,axis=1)[:,None]
        self.tree=cKDTree(self.P)
    def at(self,s):
        s=np.clip(s,0,self.L)
        P=np.stack([np.interp(s,self.s,self.P[:,c]) for c in range(3)],-1)
        T=np.stack([np.interp(s,self.s,self.T[:,c]) for c in range(3)],-1)
        T=T/np.linalg.norm(T,axis=-1,keepdims=True)
        return P,T,np.interp(s,self.s,self.r)
    def proj(self,x):
        d,i=self.tree.query(x); return self.s[i]
    def _box(self):
        if not hasattr(self,'_lo') or self._n!=len(self.P):
            R=self.r.max(); self._lo=self.P.min(0)-R; self._hi=self.P.max(0)+R; self._n=len(self.P)
        return self._lo,self._hi
    def sdf(self,X,margin=20.0):
        X=np.asarray(X); lo,hi=self._box()
        db=np.linalg.norm(np.maximum(np.maximum(lo-X,X-hi),0),axis=-1)
        out=db.copy(); m=db<margin
        if m.any():
            d,i=self.tree.query(X[m]); out[m]=d-self.r[i]
        return out
    def dr(self,X):
        d,i=self.tree.query(X); return d,self.r[i]
    def nsdf(self,X):
        d,i=self.tree.query(X); return d/self.r[i]-1

# ---------------- implicit primitives (negative inside) -----------------
def smin(a,b,k):
    h=np.clip(0.5+0.5*(b-a)/k,0,1); return b*(1-h)+a*h-k*h*(1-h)
def smax(a,b,k): return -smin(-a,-b,k)

def ellipsoid(X,c,R,frame=np.eye(3)):
    q=((X-c)@frame)/R; return (np.linalg.norm(q,axis=-1)-1)*R.min()

def prolate(X,O,D,R,C):
    v=X-O; u=v@D; rad=np.linalg.norm(v-u[...,None]*D,axis=-1)
    return (np.sqrt((rad/R)**2+(u/C)**2)-1)*R

def plane(X,P0,n):  # negative on the side opposite n
    return (X-P0)@n

class Implicit:
    """F = smooth union of parts; each part callable X->f"""
    def __init__(self,k=3.0): self.parts=[]; self.k=k; self.names=[]; self.cuts=[]
    def add(self,name,f): self.parts.append(f); self.names.append(name)
    def __call__(self,X,which=False):
        X=np.asarray(X,float); F=None
        vals=[f(X) for f in self.parts]
        for v in vals: F=v if F is None else smin(F,v,self.k)
        for c in self.cuts: F=smax(F,c(X),0.5)
        if which: return F,np.argmin(np.stack(vals),axis=0)
        return F

def hermite_chain(pts,T0,T1=None,n_per=80):
    pts=[np.asarray(p,float) for p in pts]; n=len(pts)
    tans=[]
    for i in range(n):
        if i==0: t=unit(T0)
        elif i==n-1: t=unit(T1) if T1 is not None else unit(pts[-1]-pts[-2])
        else: t=unit(unit(pts[i+1]-pts[i])+unit(pts[i]-pts[i-1]))
        tans.append(t)
    out=[]
    for i in range(n-1):
        L=np.linalg.norm(pts[i+1]-pts[i]); u=np.linspace(0,1,n_per,endpoint=False)[:,None]
        h00=2*u**3-3*u**2+1;h10=u**3-2*u**2+u;h01=-2*u**3+3*u**2;h11=u**3-u**2
        out.append(h00*pts[i]+h10*L*tans[i]+h01*pts[i+1]+h11*L*tans[i+1])
    out.append(pts[-1][None]); return np.vstack(out)

def tilt(n,d,maxdeg=40):
    """direction between face normal n and desired d, limited to maxdeg from n"""
    n=unit(n); d=unit(d); ang=np.arccos(np.clip(n@d,-1,1)); m=np.radians(maxdeg)
    if ang<=m: return d
    ax=unit(np.cross(n,d)); return n*np.cos(m)+np.cross(ax,n)*np.sin(m)

class DPath(Path):
    """path from explicit dense polyline with radii at key points (by arclength fraction)"""
    def __init__(self,dense,key_pts,key_r,ds=0.25):
        seg=np.linalg.norm(np.diff(dense,axis=0),axis=1); s=np.r_[0,np.cumsum(seg)]
        cs=[s[np.argmin(np.linalg.norm(dense-c,axis=1))] for c in np.asarray(key_pts,float)]
        L=s[-1]; n=max(int(L/ds),8); self.s=np.linspace(0,L,n+1)
        self.P=np.stack([np.interp(self.s,s,dense[:,c]) for c in range(3)],1)
        self.r=np.interp(self.s,cs,key_r); self.L=L
        T=np.gradient(self.P,axis=0); self.T=T/np.linalg.norm(T,axis=1)[:,None]
        self.tree=cKDTree(self.P)

def mkpath(pts,radii,T0,T1=None):
    return DPath(hermite_chain(pts,T0,T1),pts,radii)


def frozen(F):
    """snapshot of an Implicit (parts and cuts copied) so later cuts don't recurse"""
    G=Implicit(k=F.k); G.parts=list(F.parts); G.names=list(F.names); G.cuts=list(F.cuts); return G
