"""Boundary motion v3 = surface-harmonic extension.
Driven vertices: ventricular bodies, upper atrial bodies (analytic chamber motion), fixed caps (0),
moving coronary caps (rigid). All other wall vertices (great vessels, valve funnels, atrial floor,
junctions) get the harmonic interpolant of the surface graph Laplacian -> stretch is spread over
the whole length of each structure. Linear in the driven data -> as smooth in time as the curves."""
import numpy as np,sys; sys.path.insert(0,'/home/claude/cyc'); sys.path.insert(0,'/home/claude/hg')
import scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.spatial import cKDTree
import params as P
from cycle import Motion
from deform_v1 import LV_V,RV_V,AO_D,PA_D,sstep
class BMotion2:
    def __init__(self,bm,inout,wdrive=0.97,u_atr=-58.0,wexp=2.0,far_mm=60.0,crotch_pts=None,crotch_r=14.0):
        self.crotch_pts=crotch_pts
        self.wexp=wexp
        X=bm['pos']; T=bm['tris']; self.X=X; self.T=T; self.dom=bm['dom']; N=len(X)
        fixed=[(x['P'],x['n'],x['r']) for x in inout if 'coronary' not in x['label']]
        moving=[(x['P'],x['n'],x['r']) for x in inout if 'coronary' in x['label']]
        self.mo=Motion(X,fixed,moving_caps=moving,dom=self.dom,moving_caps_dom=['L']*len(moving))
        mo=self.mo; u=(X-P.O_LV)@P.D_LV
        drive=np.zeros(N,bool)
        drive|=mo.w_lv>wdrive; drive|=mo.w_rv>wdrive
        drive|=(mo.w_la>wdrive)&(u<u_atr); drive|=(mo.w_ra>wdrive)&(u<u_atr)
        io=bm['io_vertex']>=0
        # junction crotch lines (LV base mitral|LVOT split, RV base split, PA end split): keep the
        # surrounding wall free so its motion is interpolated smoothly instead of prescribed
        crot=getattr(self,'crotch_pts',None)
        if crot is not None and len(crot):
            from scipy.spatial import cKDTree as _T0
            if isinstance(crot,tuple):                      # (points, domain labels): act on own domain only
                cp,cd=crot
                for dd in ('L','R'):
                    if (cd==dd).any():
                        m=self.dom==dd; near=np.zeros(N,bool); near[m]=_T0(cp[cd==dd]).query(X[m])[0]<crotch_r; drive&=~near
            else: drive&=~(_T0(crot).query(X)[0]<crotch_r)
        # walls far from every beating chamber follow the analytic field (~0): no RBF wiggle there
        from scipy.spatial import cKDTree as _T
        if drive.any():
            dfar=_T(X[drive]).query(X)[0]; drive|=dfar>far_mm
        capzone=mo.w_fix<1e-6
        mcz=np.zeros(N,bool)
        for c,w in mo.mcaps: mcz|=w>0.999
        # tiny / skinny triangles are never driven: their vertices follow the harmonic interpolant
        a=np.linalg.norm(np.cross(X[T[:,1]]-X[T[:,0]],X[T[:,2]]-X[T[:,0]]),axis=1)/2
        lmax=np.max(np.stack([np.linalg.norm(X[T[:,i]]-X[T[:,(i+1)%3]],axis=1) for i in range(3)]),0)
        skinny=(a<0.12)|(a/np.maximum(lmax**2,1e-12)<0.08)
        weak=np.zeros(N,bool); weak[T[skinny].ravel()]=True
        drive&=~weak
        self.drive=drive|io|capzone|mcz
        # graph Laplacian (uniform weights), per connected mesh; Dirichlet on driven set
        r=np.concatenate([T[:,0],T[:,1],T[:,2],T[:,1],T[:,2],T[:,0]]); c=np.concatenate([T[:,1],T[:,2],T[:,0],T[:,0],T[:,1],T[:,2]])
        ln=np.linalg.norm(X[r]-X[c],axis=1); w=1/np.maximum(ln,1e-3)**self.wexp
        W=sp.coo_matrix((w,(r,c)),shape=(N,N)).tocsr(); W.sum_duplicates()
        Lap=sp.diags(np.asarray(W.sum(1)).ravel())-W
        self.f=np.where(~self.drive)[0]; self.d=np.where(self.drive)[0]
        self.Lff=Lap[self.f][:,self.f].tocsc(); self.Lfd=Lap[self.f][:,self.d]
        self.lu=spl.splu(self.Lff)
        # implicit smoothing operator for the (small) pulse-distension field: (I + mu*Ln) x = b
        Dg=np.asarray(W.sum(1)).ravel(); Ln=sp.eye(N)-sp.diags(1/Dg)@W
        self.smooth_lu=spl.splu((sp.eye(N)+20.0*Ln).tocsc())
        self.io_or_cap=io|capzone
        self.Ltri=self.dom[T[:,0]]=='L'
        # contact coupling (pairs of close L/R walls move with the pair average)
        L=self.dom=='L'; dist=np.empty(N); part=np.empty(N,int)
        dd,ii=cKDTree(X[~L]).query(X[L]); dist[L]=dd; part[L]=np.where(~L)[0][ii]
        dd,ii=cKDTree(X[L]).query(X[~L]); dist[~L]=dd; part[~L]=np.where(L)[0][ii]
        self.partner=part; self.beta=1-sstep((dist-0.6)/5.0)
        mo.vfac=sstep((dist-0.5)/4.0)
        print(f'  driven vertices {len(self.d)} / free (harmonic) {len(self.f)}')
    def disp(self,s,parts=('lv','rv','atria','vessels')):
        Da=self.mo.disp(s,tuple(p for p in parts if p!='vessels'))
        D=Da.copy()
        rhs=-(self.Lfd@Da[self.d]); D[self.f]=np.stack([self.lu.solve(rhs[:,k]) for k in range(3)],1)
        if 'vessels' in parts:                       # pulse distension on top, smoothed on the surface
            Dv=self.mo.disp(s,('vessels',))
            Dv=np.stack([self.smooth_lu.solve(Dv[:,k]) for k in range(3)],1); Dv[self.io_or_cap]=0
            D+=Dv
        D=D+(0.5*self.beta)[:,None]*(D[self.partner]-D)
        return D
    def vol(self,V,sel):
        t=self.T[sel]; a,b,c=V[t[:,0]],V[t[:,1]],V[t[:,2]]; return abs(np.einsum('ij,ij->i',np.cross(a,b),c).sum()/6)/1000
    def calibrate(self,verbose=True):
        ss=np.linspace(0,1,500)
        for tag,curve,V0c,sel in (('lv',LV_V,self.mo.LV0,self.Ltri),('rv',RV_V,self.mo.RV0,~self.Ltri)):
            s_es=float(ss[np.argmin(curve(ss))]); tgt=float(curve(s_es))-V0c; Vr=self.vol(self.X,sel); lo,hi=0.5,2.5
            for _ in range(24):
                mid=0.5*(lo+hi); setattr(self.mo,'a_'+tag,mid); dv=self.vol(self.X+self.disp(s_es,(tag,)),sel)-Vr
                if dv>tgt: lo=mid
                else: hi=mid
            setattr(self.mo,'a_'+tag,0.5*(lo+hi)); dv=self.vol(self.X+self.disp(s_es,(tag,)),sel)-Vr
            if verbose: print(f'  calibrate {tag}: amplitude x{getattr(self.mo,"a_"+tag):.3f} target {tgt:+.1f} mL achieved {dv:+.1f} mL')
        return dict(a_lv=self.mo.a_lv,a_rv=self.mo.a_rv)
