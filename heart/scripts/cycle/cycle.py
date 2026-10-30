"""Cardiac-cycle wall motion on the v3 grid boundary (morphing guide).
Analytic, spatially smooth displacement field D(x,t) (usable for ANY point, incl. interior grid nodes)."""
import sys,numpy as np; sys.path.insert(0,'/home/claude/hg'); sys.path.insert(0,'/home/claude/cyc')
import params as P
import deform_v1 as DV
from deform_v1 import LV_V,RV_V,LA_V,RA_V,AO_D,PA_D,sstep,prolate_eta,ellipsoid_eta,polyline_radial,phase_name,PHASES

class Motion:
    def __init__(self,X,caps,fade=3.0,moving_caps=(),dom=None,moving_caps_dom=None):
        """X (n,3) points; caps: (centre, normal, radius) of inlets/outlets held fixed;
        moving_caps: caps that translate rigidly with the field at their centre (stay planar)"""
        X=np.asarray(X,float); self.X=X; d=P.D_LV
        eta,u,perp=prolate_eta(X,P.O_LV,d,P.LV_R,P.LV_C); self.w_lv=sstep((1.44-eta)/0.42); self.u_lv=u; self.perp_lv=perp
        eta_r,u_r,perp_r=prolate_eta(X,P.O_RV,d,P.RV_R,P.RV_C); self.w_rv=sstep((1.32-eta_r)/0.30)*(1-self.w_lv); self.u_rv=u_r; self.perp_rv=perp_r
        self.w_la=sstep((1.34-ellipsoid_eta(X,P.LA_C,P.LA_R))/0.34); self.w_ra=sstep((1.34-ellipsoid_eta(X,P.RA_C,P.RA_R))/0.34)
        dist,rad,vec=polyline_radial(X,P.AORTA_PTS,P.AORTA_R); self.w_ao=sstep((1.55-dist/np.maximum(rad,1e-6))/0.55); self.vec_ao=vec
        ws=[];vs=[]
        for pts,rr in ((P.PA_PTS,P.PA_R),(P.RPA_PTS,P.RPA_R),(P.LPA_PTS,P.LPA_R)):
            d1,r1,v1=polyline_radial(X,pts,rr); ws.append(sstep((1.55-d1/np.maximum(r1,1e-6))/0.55)); vs.append(v1)
        ws=np.stack(ws); sel=np.argmax(ws,0); self.w_pa=ws.max(0); self.vec_pa=np.stack(vs)[sel,np.arange(len(X))]
        # hold inlets/outlets: smooth fade to zero within `fade` radii of every cap centre
        self.w_fix=np.ones(len(X))
        for c,n,r in caps:
            dd=np.linalg.norm(X-c,axis=1); self.w_fix=np.minimum(self.w_fix,sstep((dd-1.0*r)/(fade*r)))
        if dom is not None:     # left-heart motions act only on the left domain, right-heart on the right
            L=(np.asarray(dom)=='L').astype(float); R=1-L
            self.w_lv=self.w_lv*L; self.w_la=self.w_la*L; self.w_ao=self.w_ao*L
            self.w_rv=sstep((1.32-eta_r)/0.30)*R; self.w_ra=self.w_ra*R; self.w_pa=self.w_pa*R
        self._dom=dom
        self.mcaps=[]
        for c,n,r in moving_caps:
            dd=np.linalg.norm(X-c,axis=1); w=1-sstep((dd-1.05*r)/(fade*r))
            self.mcaps.append((np.asarray(c,float),w))
        self._mcdom=list(moving_caps_dom) if moving_caps_dom is not None else [None]*len(self.mcaps)
        self.a_lv=self.a_rv=1.0; self.LV0=float(LV_V(0.0)); self.RV0=float(RV_V(0.0))
        ss=np.linspace(0,1,500); self.LAmax=float(LA_V(ss).mean()); self.RAmax=float(RA_V(ss).mean())   # reference = cycle-mean atrial volume
        self.vfac=np.ones(len(X))   # optional local limiter for vessel distension
    @staticmethod
    def _vent(ratio,k_long):
        ratio=max(ratio,0.12); fL=1-k_long*(1-ratio); return fL,float(np.sqrt(max(ratio/fL,1e-6)))
    def disp(self,s,parts=('lv','rv','atria','vessels')):
        X=self.X; D=np.zeros_like(X); d=P.D_LV
        if 'lv' in parts:
            ratio=1-self.a_lv*(1-float(LV_V(s))/self.LV0); fL,fR=self._vent(ratio,getattr(self,'k_lv',0.24))
            zeta=P.LV_C-self.u_lv; amp=np.clip((1-ratio)/0.624,0,1.6)
            phi=np.radians(getattr(self,'twist_deg',6.0))*amp*(1-2*zeta/(P.LV_C+P.LV_BASE_OFF))
            e1=np.cross(d,[0,0,1.]); e1/=np.linalg.norm(e1); e2=np.cross(d,e1)
            a=self.perp_lv@e1; b=self.perp_lv@e2; ca,sa=np.cos(phi),np.sin(phi)
            newp=(fR*(a*ca-b*sa))[:,None]*e1+(fR*(a*sa+b*ca))[:,None]*e2
            D+=self.w_lv[:,None]*((1-fL)*zeta[:,None]*d+(newp-self.perp_lv))
        if 'rv' in parts:
            ratio=1-self.a_rv*(1-float(RV_V(s))/self.RV0); fL,fR=self._vent(ratio,getattr(self,'k_rv',0.20))
            zeta=P.RV_C-self.u_rv; D+=self.w_rv[:,None]*((1-fL)*zeta[:,None]*d+(fR-1)*self.perp_rv)
        if 'atria' in parts:
            D+=self.w_la[:,None]*(((float(LA_V(s))/self.LAmax)**(1/3))-1)*(X-P.LA_C)
            D+=self.w_ra[:,None]*(((float(RA_V(s))/self.RAmax)**(1/3))-1)*(X-P.RA_C)
        if 'vessels' in parts:
            D+=(self.vfac*self.w_ao)[:,None]*(float(AO_D(s))-1)*self.vec_ao
            D+=(self.vfac*self.w_pa)[:,None]*(float(PA_D(s))-1)*self.vec_pa
        D=D*self.w_fix[:,None]
        if self.mcaps and not getattr(self,'_inner',False):
            for (c,w),cd in zip(self.mcaps,self._mcdom):
                if not np.any(w>0): continue
                sub=Motion.__new__(Motion); sub.__dict__.update(self._point_state(c,cd))
                Dc=sub.disp(s,parts)[0]
                D=D*(1-w)[:,None]+w[:,None]*Dc
        return D
    def _point_state(self,c,cd=None):
        if not hasattr(self,'_pcache'): self._pcache={}
        k=tuple(np.round(c,6))
        if k not in self._pcache:
            m=Motion(np.array([c]),[],dom=None if cd is None else np.array([cd])); m.a_lv=self.a_lv; m.a_rv=self.a_rv; m._inner=True
            self._pcache[k]=m.__dict__
        st=dict(self._pcache[k]); st['a_lv']=self.a_lv; st['a_rv']=self.a_rv; return st
