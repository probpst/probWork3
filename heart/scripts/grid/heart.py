"""Heart core skeleton + implicit walls (left & right circulations are separate domains)."""
import numpy as np, params as P
from paths import Path, unit, Implicit, ellipsoid, prolate, smax, mkpath, tilt, frozen
from core import tube, box, sbox, sub, face_normal, pick_range, ext, smooth
from mb import Block

D=P.D_LV; W=P.W_RV; V=np.cross(D,W)
def cyl_dir(phi): return np.cos(phi)*W+np.sin(phi)*V
RES=dict(ds=1.0,Nm1=20,Na1=16,N2=20,Nla3=20,Nlv3=46,Nt_arc=18,Nrvot=14,Nth=14,Nra1=18,Nrv3=48,beta_a=0.62,kblend=6.0,lean=0.7,lean_rv=0.4,Ngap=3,Ngap_rv=3,Ngap_pa=2,cedge=0.0,isthmus=11.0,vein_tilt=20.0,cclear=3.0,ctau=1.25,ck=1.5)

def wrap90(a): return (a+np.pi/4)%(np.pi/2)-np.pi/4
def mean90(angs):  # circular mean with 90deg period
    z=np.mean(np.exp(4j*np.asarray(angs))); return np.angle(z)/4

def curtain(pA,pB,O,Dax,u_edge,tau=0.75,reach=2.0,k=0.8,rs=14.0,gap=0.33,cores=None,clear=2.0):
    """thin septum on the radius-normalised midline of tubes A,B, only in their overlap zone, above u_edge (toward base)"""
    from paths import smin
    def C(X):
        X=np.asarray(X); sa=pA.sdf(X); sb=pB.sdf(X)
        c=reach-np.maximum(sa,sb)
        m=c>-3*k
        if m.any():
            Xm=X[m]; da,ra=pA.dr(Xm); db,rb=pB.dr(Xm); a=da/ra-1; b=db/rb-1; u=(Xm-O)@Dax
            cm=smin(tau-rs*np.abs(a-b), u_edge-u, k); cm=smin(cm,c[m],k)
            if cores is not None:
                mm=cm>-3*k
                if mm.any(): cm[mm]=smin(cm[mm], cores.query(Xm[mm])[0]-clear, k)
            c[m]=cm
        return c
    return C

def flare(path,Fparent,r_base_fn=None,fac=0.8):
    """ostial funnel: extra radius over chimney height H (distance from start to parent wall)"""
    P0,T0,_=path.at(0.0); H=0.0
    for t in np.arange(0.2,60,0.2):
        if Fparent(np.array([P0+t*T0]))[0]>0: H=t; break
    r0=path.r[0]; Lf=H+r0
    path.r=path.r+fac*H*np.clip(1-path.s/Lf,0,1)**2
    return H

def wall_dist(P0,n0,Fp,tmax=60):
    for t in np.arange(0.2,tmax,0.2):
        if Fp(np.array([P0+t*n0]))[0]>0: return t
    return tmax

def ostial_path(P0,n0,Fp,pts,rads,r0,run=0.8,Rfac=2.2,maxtilt=0.0):
    """leave the parent face perpendicular, straight through the parent wall, then a controlled arc
    (radius Rfac*r0) toward the anatomical centreline"""
    if maxtilt>0:
        tg=[p for p in pts if np.linalg.norm(p-P0)>3*r0]
        if tg: n0=tilt(n0,tg[0]-P0,maxtilt)
    H=wall_dist(P0,n0,Fp); Pm=P0+n0*H; P1=P0+n0*(H+run*r0)
    tgt=[(p,r) for p,r in zip(pts,rads) if np.linalg.norm(p-P1)>3*r0]
    if not tgt: tgt=[(pts[-1],rads[-1])]
    dt=unit(tgt[0][0]-P1); th=np.arccos(np.clip(n0@dt,-1,1)); ctrl=[P0,Pm,P1]; rr=[r0,r0,r0]
    if th>np.radians(10):
        b=unit(dt-(dt@n0)*n0); Rt=Rfac*r0; th=min(th,np.radians(80))
        for f in (0.5,1.0):
            a=th*f; ctrl.append(P1+Rt*(np.sin(a)*n0+(1-np.cos(a))*b)); rr.append(r0)
        P2=ctrl[-1]; d2=np.cos(th)*n0+np.sin(th)*b
        tgt=[(p,r) for p,r in tgt if (p-P2)@d2>2*r0] or [(P2+d2*4*r0,rads[-1])]
    return mkpath(ctrl+[p for p,_ in tgt],rr+[r for _,r in tgt],n0)

class Build:
    def __init__(s,**kw):
        s.R=dict(RES); sc=kw.pop('scale',1.0); s.R.update(kw)
        if sc!=1.0:
            for k in ('Nm1','Na1','N2','Nla3','Nlv3','Nt_arc','Nrvot','Nth','Nra1','Nrv3'):
                s.R[k]=max(4,int(round(s.R[k]*sc)))
        s.blocks=[]; s.inout=[]   # inlet/outlet specs
        s.FL=Implicit(k=s.R['kblend']); s.FR=Implicit(k=s.R['kblend']); s.comp_uD={}; s.paths={}
    def add(s,B): s.blocks.append(B); return B
    def exclude(s,g=0.8):
        """mutual exclusion: each circulation stays >= g away from the other's anatomy"""
        L0=frozen(s.FL); R0=frozen(s.FR)
        s.FL.cuts.append(lambda X:g-R0(X)); s.FR.cuts.append(lambda X:g-L0(X))
    def endcap(s,B,typ,label):
        B.tag('imax',typ,label)
        P0,T,r=B.path.at(B.path.L); s.inout.append(dict(dom=B.dom,P=P0,n=T,r=r,type=typ,label=label,block=B.name))
    def tube_part(s,F,name,path):
        F.add(name,path.sdf)

    # ------------------------------------------------------------ LEFT
    def left(s):
        R=s.R; ds=R['ds']
        e1=cyl_dir(np.radians(-13)); e2=np.cross(D,e1)
        G=R['Ngap']; N1=R['Nm1']+G+R['Na1']; N2=R['N2']   # LV base: [mitral | gap strip | LVOT]
        # LV frustum core
        ub,ua=-26.0,30.0; k=np.linspace(ub,ua,R['Nlv3']+1)
        prof=np.sqrt(np.clip(1-(k/P.LV_C)**2,0.05,1))*P.LV_R
        x=np.linspace(-1,1,N1+1); y=np.linspace(-1,1,N2+1)
        X=(P.O_LV+k[None,None,:,None]*D+(0.52*prof)[None,None,:,None]*x[:,None,None,None]*e1
           +(0.45*prof)[None,None,:,None]*y[None,:,None,None]*e2)
        LV=s.add(Block('LV',X,'L','core','LV'))
        # LA box
        MA=np.array(P.MITRAL_PTS[2]); e3a=unit(P.LA_C-MA); e1a=unit(e1-(e1@e3a)*e3a); e2a=np.cross(e3a,e1a)
        E=np.stack([e1a,e2a,e3a],1); Rx=np.array([ext(P.LA_R,e) for e in E.T])
        LA=s.add(sbox('LA',P.LA_C,E,Rx,(R['Nm1'],N2,R['Nla3']),'L','LA',beta=s.R['beta_a']))
        # mitral tube: LA kmin (full) -> LV kmin [0..Nm1]
        S0=LA.face('kmin'); P0=S0.mean((0,1))
        S1=sub(LV,'kmin',(0,R['Nm1']),(0,N2)); P1=S1.mean((0,1))
        T0=tilt(-e3a,MA-P0); Te=tilt(D+0.4*e1,P1-MA)
        T0=tilt(-e3a,P1-P0); Te=unit(D+s.R['lean']*e1)
        pm=mkpath([P0,P1],[18,17],T0,Te)
        MV=s.add(tube('MV',pm,S0,'L','mitral',cf=0.5,ds=ds,S1=S1))
        # aorta: LV kmin [Nm1..N1] -> ... outlet
        S0=sub(LV,'kmin',(R['Nm1']+G,N1),(0,N2)); P0=S0.mean((0,1))
        ctrl=[P0]+[np.array(p) for p in P.AORTA_PTS[3:]]
        T0=unit(-D+s.R['lean']*e1)
        pa=mkpath(ctrl,[12]+P.AORTA_R[3:],T0)
        # --- twist design (coronaries at root, arch branches) ---
        AO0=tube('tmp',pa,S0,'L','aorta',cf=0.5,ds=ds)
        def req(B,pt,d):
            sb=B.path.proj(pt); i=np.argmin(abs(B.s-sb)); T=B.path.at(B.s[i])[1]
            d=d-(d@T)*T; return B.s[i],np.arctan2(d@B.F2[i],d@B.F1[i])
        cor={n:req(AO0,np.array(p[0]),np.array(p[1])-np.array(p[0])) for n,(p,r) in P.CORONARIES.items()}
        br={n:req(AO0,np.array(p[0]),np.array(p[1])-np.array(p[0])) for n,(p,r) in P.BRANCHES.items()}
        th_r=mean90([wrap90(a) for sb,a in cor.values()]); th_a=mean90([wrap90(a) for sb,a in br.values()])
        th_a=th_r+wrap90(th_a-th_r)
        s_r=np.mean([v[0] for v in cor.values()]); s_a0=min(v[0] for v in br.values())
        def tw(sv):
            return np.where(sv<s_r,th_r*smooth(sv/(s_r-6)),th_r+(th_a-th_r)*smooth((sv-s_r-8)/(s_a0-s_r-24)))
        AO=s.add(tube('AO',pa,S0,'L','aorta',cf=0.5,ds=ds,twist=tw))
        s.endcap(AO,'outlet','outlet_descending_aorta')
        s.FL.add('LA',lambda X:ellipsoid(X,P.LA_C,P.LA_R))
        s.FL.add('LV',lambda X:smax(prolate(X,P.O_LV,D,P.LV_R,P.LV_C),-((X-P.O_LV)@D+32.0),2.0))
        s.FL.add('mitral',pm.sdf); s.FL.add('aorta',pa.sdf)
        from scipy.spatial import cKDTree
        tr=cKDTree(np.vstack([B.X.reshape(-1,3) for B in (LV,MV,AO)]))
        s.FL.cuts.append(curtain(pm,pa,P.O_LV,D,ub-s.R['cedge'],cores=tr,clear=s.R['cclear'],tau=s.R['ctau'],k=s.R['ck']))
        s.side_branches(AO,P.CORONARIES,'L',cf=0.55,dsb=0.6,label='outlet_',run=0.3,Rfac=1.5,search=(-3,3),maxtilt=15)
        s.side_branches(AO,P.BRANCHES,'L',cf=0.5,dsb=ds,label='outlet_',run=0.3,maxtilt=30)
        s.chamber_veins(LA,P.PULM_VEINS,['imin','imax'],'L',P.LA_C,P.LA_R,label='inlet_',keep={'imin':(1,1,4,1),'imax':(1,1,4,1)})
        s.comp_uD.update({'LA':(0.4,.045),'LV':(0.6,.05),'mitral':(0.9,.032),'aorta':(1.1,.027),
            'brachiocephalic_trunk':(0.9,.012),'left_common_carotid':(0.9,.012),'left_subclavian':(0.9,.012),
            'left_coronary_ostium':(0.25,.0045),'right_coronary_ostium':(0.25,.0045)})
        for n in P.PULM_VEINS: s.comp_uD[n]=(0.5,.013)

    def side_branches(s,AO,spec,dom,cf,dsb,label,run=0.8,Rfac=2.2,search=(-6,14),maxtilt=0.0):
        F=s.FL if dom=='L' else s.FR
        if not hasattr(s,'_used_ao'): s._used_ao={}
        for name,(pts,rad) in spec.items():
            pts=[np.array(p,float) for p in pts]; d=unit(pts[1]-pts[0])
            sb0=AO.path.proj(pts[0]); opts=[]
            for sh in np.arange(search[0],search[1]+0.1,1.0):
                i=int(np.argmin(abs(AO.s-(sb0+sh))))
                cand={'jmax':AO.F1[i],'jmin':-AO.F1[i],'kmax':AO.F2[i],'kmin':-AO.F2[i]}
                f=max(cand,key=lambda k:cand[k]@d); opts.append((cand[f]@d,-abs(sh),i,f))
            amax=max(o[0] for o in opts)
            opts=sorted([o for o in opts if o[0]>=amax-0.06],key=lambda o:o[1],reverse=True)+sorted(opts,reverse=True)
            dsa=AO.s[1]-AO.s[0]; hp=max(2,int(round(cf*rad[0]/dsa)))
            for _,_,i,f in opts:
                if all(abs(i-jj)>=2*hp+4 for jj in s._used_ao.get(f,[])): break
            s._used_ao.setdefault(f,[]).append(i)
            ni,nj,nk=AO.dims; across=(nk if f[0]=='j' else nj)
            a,b=AO.hs; hsz=(b[i] if f[0]=='j' else a[i]); sp=2*hsz/(across-1)
            hw=cf*rad[0]; hq=max(2,int(round(hw/sp))); hq=min(hq,across//2-2)
            dsa=AO.s[1]-AO.s[0]; hp=max(2,int(round(hw/dsa)))
            c=across//2; pr=(i-hp,i+hp); qr=(c-hq,c+hq)
            S0=sub(AO,f,pr,qr); n0=face_normal(AO,f,pr,qr); P0=S0.mean((0,1))
            r0=max(rad[0],0.9*max(hq*sp,hp*dsa)/cf)
            path=ostial_path(P0,n0,AO.path.sdf,pts[1:],rad[1:],r0,run=run,Rfac=Rfac,maxtilt=maxtilt)
            B=s.add(tube(name,path,S0,dom,name,cf=cf,ds=dsb))
            s.endcap(B,'outlet' if label=='outlet_' else 'inlet',label+name); F.add(name,path.sdf)

    def _separate(s,CH,plan,keep,gap):
        for it in range(20):
            moved=False
            for i in range(len(plan)):
                for j in range(i+1,len(plan)):
                    A,B=plan[i],plan[j]
                    if A[3]!=B[3]: continue
                    f=A[3]; npn,nqn=CH.face(f).shape[:2]; kp=keep.get(f,(1,1,1,1)) if keep else (1,1,1,1)
                    gp=max(B[4][0]-A[4][1],A[4][0]-B[4][1]); gq=max(B[5][0]-A[5][1],A[5][0]-B[5][1])
                    if gp>=gap or gq>=gap: continue
                    ax=4 if abs((A[4][0]+A[4][1])-(B[4][0]+B[4][1]))>=abs((A[5][0]+A[5][1])-(B[5][0]+B[5][1])) else 5
                    lo,hi=(kp[0],npn-1-kp[1]) if ax==4 else (kp[2],nqn-1-kp[3])
                    first,second=(A,B) if A[ax][0]<B[ax][0] else (B,A)
                    need=gap-(second[ax][0]-first[ax][1])
                    d1=min(need//2+need%2, first[ax][0]-lo); d2=min(need-d1, hi-second[ax][1]); d1=min(need-d2,first[ax][0]-lo)
                    first[ax]=[first[ax][0]-d1,first[ax][1]-d1]; second[ax]=[second[ax][0]+d2,second[ax][1]+d2]; moved=True
                    while second[ax][0]-first[ax][1]<gap:   # no room to shift: shrink instead
                        w1=first[ax][1]-first[ax][0]; w2=second[ax][1]-second[ax][0]
                        if max(w1,w2)<=2: raise RuntimeError(f'cannot separate {A[0]} {B[0]}')
                        if w1>=w2 and w1>2: first[ax][1]-=1
                        elif w2>2: second[ax][0]+=1
                        else: raise RuntimeError(f'cannot separate {A[0]} {B[0]}')
            if not moved: return

    def chamber_veins(s,CH,spec,faces,dom,cen,Rad,label,cf=0.5,dsv=None,keep=None):
        F=s.FL if dom=='L' else s.FR; plan=[]
        for name,(pts,rad) in spec.items():
            pts=[np.array(p,float) for p in pts]; dsv_=dsv or s.R['ds']
            best=None
            for f in faces:
                n=face_normal(CH,f); Fc=CH.face(f).reshape(-1,3).mean(0)
                dvec=unit(pts[-1]-pts[0]); sc=n@dvec
                if best is None or sc>best[0]: best=(sc,f,n,Fc)
            _,f,n,Fc=best
            # intersection of anatomical line with face plane
            dense=Path(pts,rad).P; t=(dense-Fc)@n; j=int(np.argmin(abs(t)))
            pr,qr=pick_range(CH,f,dense[j],cf*rad[0])
            if keep and f in keep:
                npn,nqn=CH.face(f).shape[:2]; (plo,phi,qlo,qhi)=keep[f]
                sh=max(0,plo-pr[0])-max(0,pr[1]-(npn-1-phi)); pr=(pr[0]+sh,pr[1]+sh)
                sh=max(0,qlo-qr[0])-max(0,qr[1]-(nqn-1-qhi)); qr=(qr[0]+sh,qr[1]+sh)
            plan.append([name,pts,rad,f,list(pr),list(qr)])
        # ---- separate ranges sharing a face (>= gap cells in p or q)
        s._separate(CH,plan,keep,gap=3)
        for name,pts,rad,f,pr,qr in plan:
            pr=tuple(pr); qr=tuple(qr); dsv_=dsv or s.R['ds']
            S0=sub(CH,f,pr,qr); P0=S0.mean((0,1)); n0=face_normal(CH,f,pr,qr)
            hd=0.5*max(np.linalg.norm(S0[-1,-1]-S0[0,0]),np.linalg.norm(S0[-1,0]-S0[0,-1]))
            path=ostial_path(P0,n0,lambda X:ellipsoid(X,cen,Rad),pts,rad,max(rad[0]*1.05,hd/0.8),maxtilt=s.R['vein_tilt'])
            s.paths[name]=path
            B=s.add(tube(name,path,S0,dom,name,cf=cf,ds=dsv_))
            s.endcap(B,'inlet',label+name); F.add(name,path.sdf)

    # ------------------------------------------------------------ RIGHT
    def right(s):
        R=s.R; ds=R['ds']; Gr=R['Ngap_rv']; Na=R['Nrvot']+Gr+R['Nt_arc']; Nth=R['Nth']   # RV base: [RVOT | gap | tricuspid]
        # RV crescent core
        def ro(phi,u):
            rr=P.RV_R*np.sqrt(np.clip(1-(u/P.RV_C)**2,0,1)); S=P.RV_SHIFT
            return S*np.cos(phi)+np.sqrt(np.clip(rr**2-(S*np.sin(phi))**2,0,None))
        def re(u): return (P.LV_R+P.LV_WALL)*np.sqrt(np.clip(1-(u/(P.LV_C+P.LV_WALL))**2,0,1))
        def phi_edge(u):
            ph=np.radians(np.linspace(0,120,1201)); ok=ph[ro(ph,u)>re(u)+0.5]; return ok.max()
        ub,ua=-24.0,38.0; uk=np.linspace(ub,ua,R['Nrv3']+1)
        t=np.linspace(0,1,Nth+1); X=np.zeros((Na+1,Nth+1,len(uk),3))
        for kk,u in enumerate(uk):
            pm=0.5*phi_edge(u); ph=np.linspace(pm,-pm,Na+1)
            r0=re(u); r1=ro(ph,u)
            rr=r0+(r1-r0)[:,None]*(0.30+0.42*t[None,:])
            X[:,:,kk]=P.O_LV+u*D+rr[...,None]*(np.cos(ph)[:,None,None]*W+np.sin(ph)[:,None,None]*V)
        from core import smooth3d; smooth3d(X)
        RV=s.add(Block('RV',X,'R','core','RV'))
        # RA box
        TA=np.array(P.TRICUSPID_PTS[2]); e3=np.array([0,0,1.]); hdir=TA-P.RA_C; hdir[2]=0; e1=unit(hdir); e2=np.cross(e3,e1)
        E=np.stack([e1,e2,e3],1); Rx=np.array([ext(P.RA_R,e) for e in E.T])
        RA=s.add(sbox('RA',P.RA_C,E,Rx,(R['Nra1'],R['Nt_arc'],Nth),'R','RA',beta=s.R['beta_a']))
        # tricuspid: RA imax full -> RV kmin [Nrvot..Na]
        S0=RA.face('imax'); P0=S0.mean((0,1))
        S1=sub(RV,'kmin',(R['Nrvot']+Gr,Na),(0,Nth)); P1=S1.mean((0,1))
        split=RV.face('kmin')[R['Nrvot']+Gr//2].mean(0); lean=unit(P1-split-((P1-split)@D)*D)
        Te=tilt(D-s.R['lean_rv']*lean,P1-TA); T0=tilt(e1,TA-P0)
        pt=mkpath([P0,TA,P1],[14,P.TV_ANNULUS_R,18],T0,Te)
        TV=s.add(tube('TV',pt,S0,'R','tricuspid',cf=0.5,ds=ds,S1=S1))
        # PA
        S0=sub(RV,'kmin',(0,R['Nrvot']),(0,Nth)); P0=S0.mean((0,1)); lean=unit(P0-split-((P0-split)@D)*D)
        ctrl=[P0]+[np.array(p) for p in P.PA_PTS[2:]]
        T0=tilt(-D+s.R['lean_rv']*lean,ctrl[1]-P0,35)
        pp=mkpath(ctrl,[19]+P.PA_R[2:],T0)
        PA0=tube('tmp',pp,S0,'R','pa',cf=0.5,ds=ds)
        dRL=np.array(P.RPA_PTS[1])-np.array(P.LPA_PTS[1]); T=PA0.path.at(PA0.path.L)[1]; dRL=unit(dRL-(dRL@T)*T)
        ang=np.arctan2(dRL@PA0.F2[-1],dRL@PA0.F1[-1]); th=wrap90(ang)
        L=pp.L; tw=lambda sv:th*smooth((sv-0.35*L)/(0.55*L))
        PA=s.add(tube('PA',pp,S0,'R','pa',cf=0.5,ds=ds,twist=tw))
        # which axis now aligns with dRL
        ax='j' if abs(dRL@PA.F1[-1])>abs(dRL@PA.F2[-1]) else 'k'
        vec=PA.F1[-1] if ax=='j' else PA.F2[-1]; sgn=np.sign(dRL@vec)   # +: RPA on high side
        ni,nj,nk=PA.dims; n_ax=nj if ax=='j' else nk; m=(n_ax-1)//2
        gp=R['Ngap_pa']; lo=(0,m-gp//2); hi=(m+gp-gp//2,n_ax-1); full=(0,(nk if ax=='j' else nj)-1)
        rng={'RPA':hi if sgn>0 else lo,'LPA':lo if sgn>0 else hi}
        Tpa=T; Pbif=PA.face('imax').mean((0,1))
        for nm,pts,rad in [('RPA',P.RPA_PTS,P.RPA_R),('LPA',P.LPA_PTS,P.LPA_R)]:
            pr,qr=(rng[nm],full) if ax=='j' else (full,rng[nm])
            S0=sub(PA,'imax',pr,qr); P0=S0.mean((0,1)); ln=unit(P0-Pbif-((P0-Pbif)@Tpa)*Tpa)
            T0=unit(Tpa+0.45*ln); far=[np.array(p) for p in pts[1:] if np.linalg.norm(np.array(p)-P0)>2.5*rad[0]]
            path=mkpath([P0]+far,[rad[0]*0.95]+rad[len(rad)-len(far):],tilt(T0,far[0]-P0,35))
            B=s.add(tube(nm,path,S0,'R',nm,cf=0.5,ds=ds)); s.endcap(B,'outlet','outlet_'+('right' if nm=='RPA' else 'left')+'_pulmonary_artery')
            s.FR.add(nm,path.sdf)
        s.FR.add('RA',lambda X:ellipsoid(X,P.RA_C,P.RA_R))
        def rvf(X):
            v=X-P.O_LV; u=v@D
            f=smax(prolate(X,P.O_LV+P.RV_SHIFT*W,D,P.RV_R,P.RV_C),-prolate(X,P.O_LV,D,P.LV_R+P.LV_WALL,P.LV_C+P.LV_WALL),2.5)
            return smax(f,-(u+P.RV_BASE_OFF),2.0)
        s.FR.add('RV',rvf); s.FR.add('tricuspid',pt.sdf); s.FR.add('pa',pp.sdf)
        from scipy.spatial import cKDTree
        tr=cKDTree(np.vstack([B.X.reshape(-1,3) for B in (RV,TV,PA)]))
        s.FR.cuts.append(curtain(pt,pp,P.O_LV,D,ub-s.R['cedge'],cores=tr,clear=s.R['cclear'],tau=s.R['ctau'],k=s.R['ck']))
        s.chamber_veins(RA,P.CAVAE,['kmin','kmax'],'R',P.RA_C,P.RA_R,label='inlet_',keep={'kmin':(1,4,1,1),'kmax':(1,4,1,1)})
        if s.R['isthmus']:
            ivc=[b for b in s.blocks if b.name=='ivc'][0]
            tr2=cKDTree(np.vstack([B.X.reshape(-1,3) for B in (RA,TV,ivc)]))
            s.FR.cuts.append(curtain(pt,s.paths['ivc'],P.O_LV,D,1e9,cores=tr2,clear=s.R['cclear'],reach=s.R['isthmus'],tau=s.R['ctau'],k=s.R['ck']))
        s.comp_uD.update({'RA':(0.3,.045),'RV':(0.5,.045),'tricuspid':(0.6,.035),'pa':(0.9,.027),
                          'RPA':(0.8,.021),'LPA':(0.8,.021),'svc':(0.25,.019),'ivc':(0.30,.023)})
