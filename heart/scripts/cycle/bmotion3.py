"""Boundary motion v4: driven vertices (ventricular bodies, upper atria, inlets/outlets) follow the
analytic chamber motion; every other wall vertex follows a 3-D RBF (r^3 + linear) fitted to driven
control points, so walls that are close in space move consistently (volumetrically smooth field);
the small driven-set residual is blended in with a surface-harmonic extension."""
import numpy as np,sys; sys.path.insert(0,'/home/claude/cyc'); sys.path.insert(0,'/home/claude/hg')
import scipy.linalg as sla
from bmotion2 import BMotion2
def rbf_eval(C,coef,Y,chunk=6000):
    out=np.empty((len(Y),3)); w=coef[:-4]; lin=coef[-4:]
    for i in range(0,len(Y),chunk):
        y=Y[i:i+chunk]; d=np.sqrt(((y[:,None,:]-C[None,:,:])**2).sum(-1)); out[i:i+chunk]=(d*d*d)@w+y@lin[1:]+lin[0]
    return out

class BMotion3(BMotion2):
    def __init__(self,bm,inout,nctrl=(2600,1600),seed=0,rigid_groups=None,r_rig=10.0,r_blend=8.0,septum_mm=0.0,fade_far=(25.0,60.0),**kw):
        super().__init__(bm,inout,**kw)
        self.rigid=[]
        if fade_far is not None:
            from scipy.spatial import cKDTree as _Tf
            from deform_v1 import sstep as _sf
            mo=self.mo; core_drive=(mo.w_lv>0.97)|(mo.w_rv>0.97)|(mo.w_la>0.97)|(mo.w_ra>0.97)
            dist=_Tf(self.X[core_drive]).query(self.X)[0]; lo,hi=fade_far
            self.wfar=1-_sf((dist-lo)/(hi-lo))
        # thin septa (same domain): facing wall pairs closer than ~3 mm move together (off by default)
        if septum_mm>0:
          from scipy.spatial import cKDTree as _T
          from deform_v1 import sstep as _ss
          X=self.X; T=self.T; tn=np.cross(X[T[:,1]]-X[T[:,0]],X[T[:,2]]-X[T[:,0]]); vn=np.zeros_like(X)
          for k in range(3): np.add.at(vn,T[:,k],tn)
          vn/=np.linalg.norm(vn,axis=1,keepdims=True)
          tr=_T(X); p2=np.arange(len(X)); b2=np.zeros(len(X))
          for v,lst in enumerate(tr.query_ball_point(X,septum_mm)):
              best=None
              for w in lst:
                  if w==v or self.dom[w]!=self.dom[v]: continue
                  dv=X[w]-X[v]; dist=np.linalg.norm(dv)
                  if dist<1e-9 or vn[w]@vn[v]>-0.7: continue
                  if dv@vn[v]<0.7*dist: continue          # partner is across the solid (along the outward normal)
                  if best is None or dist<best[1]: best=(w,dist)
              if best is not None: p2[v]=best[0]; b2[v]=1-_ss((best[1]-0.5*septum_mm)/(0.5*septum_mm))
          self.p2=p2; self.b2=b2*(self.io_or_cap==False)
          print(f'  septum pairs: {(self.b2>0.01).sum()} vertices')
        if rigid_groups is not None:
            from scipy.spatial import cKDTree as _T
            from deform_v1 import sstep
            for pts in rigid_groups:
                d=_T(pts).query(self.X)[0]; w=1-sstep((d-r_rig)/r_blend)
                w[self.dom!=self.dom[_T(self.X).query(pts[len(pts)//2])[1]]]=0
                idx=np.where(w>1e-4)[0]; ww=np.zeros(len(self.X)); ww[idx]=w[idx]; self.rigid.append((idx,ww))
        X=self.X; rng=np.random.default_rng(seed); self.rbf=[]
        io=bm['io_vertex']>=0
        for dd,nc in zip(('L','R'),nctrl):
            dv=np.where(self.drive&(self.dom==dd))[0]
            capv=np.where(io&(self.dom==dd))[0]
            cand=dv[~np.isin(dv,capv)]
            # farthest-point sampling on the driven (non-cap) vertices, seeded with a cap subsample
            capsel=capv[::max(1,len(capv)//250)]
            sel=list(capsel); dmin=np.full(len(cand),1e9)
            for cc in capsel: dmin=np.minimum(dmin,np.linalg.norm(X[cand]-X[cc],axis=1))
            for _ in range(nc-len(sel)):
                k=int(np.argmax(dmin)); sel.append(cand[k]); dmin=np.minimum(dmin,np.linalg.norm(X[cand]-X[cand[k]],axis=1))
            sel=np.array(sel); C=X[sel]; n=len(C)
            M=np.block([[np.linalg.norm(C[:,None]-C[None],axis=2)**3,np.hstack([np.ones((n,1)),C])],[np.hstack([np.ones((n,1)),C]).T,np.zeros((4,4))]])
            vid=np.where(self.dom==dd)[0]; Y=X[vid]
            E=np.empty((len(vid),n+4),np.float32)
            for i in range(0,len(vid),4000):
                y=Y[i:i+4000]; d=np.sqrt(((y[:,None,:]-C[None,:,:])**2).sum(-1))
                E[i:i+4000,:n]=d*d*d; E[i:i+4000,n]=1; E[i:i+4000,n+1:]=y
            self.rbf.append(dict(sel=sel,C=C,lu=sla.lu_factor(M),vid=vid,E=E))
    def disp(self,s,parts=('lv','rv','atria','vessels')):
        Da=self.mo.disp(s,tuple(p for p in parts if p!='vessels'))
        D=np.zeros_like(Da)
        for r in self.rbf:
            coef=sla.lu_solve(r['lu'],np.vstack([Da[r['sel']],np.zeros((4,3))]))
            D[r['vid']]=r['E']@coef.astype(np.float32)
        res=np.zeros_like(Da); res[self.d]=Da[self.d]-D[self.d]
        rhs=-(self.Lfd@res[self.d]); res[self.f]=np.stack([self.lu.solve(rhs[:,k]) for k in range(3)],1)
        D=D+res                                      # exact on driven vertices, smooth elsewhere
        if hasattr(self,'b2'): D=D+(0.5*self.b2)[:,None]*(D[self.p2]-D)
        if getattr(self,'wfar',None) is not None: D=D*self.wfar[:,None]      # smooth fade far from the chambers
        # junction crotches move as rigid bodies (best-fit rotation+translation of the local wall motion)
        for idx,w in getattr(self,'rigid',[]):
            Xr=self.X[idx]; Yr=Xr+D[idx]; wr=w[idx]
            cx=(Xr*wr[:,None]).sum(0)/wr.sum(); cy=(Yr*wr[:,None]).sum(0)/wr.sum()
            H=((Xr-cx)*wr[:,None]).T@(Yr-cy); U,S_,Vt=np.linalg.svd(H); dd=np.sign(np.linalg.det(Vt.T@U.T))
            R=Vt.T@np.diag([1,1,dd])@U.T
            Drig=(Xr-cx)@R.T+cy-Xr
            D[idx]=(1-wr)[:,None]*D[idx]+wr[:,None]*Drig
        if 'vessels' in parts:
            Dv=self.mo.disp(s,('vessels',))
            Dv=np.stack([self.smooth_lu.solve(Dv[:,k]) for k in range(3)],1); Dv[self.io_or_cap]=0
            D+=Dv
        D=D+(0.5*self.beta)[:,None]*(D[self.partner]-D)
        return D
