"""Local mesh optimisation: maximise the worst scaled Jacobian around each free node (boundary fixed,
shared interface nodes numbered once -> point matching preserved). Independent sets via greedy colouring."""
import pickle,numpy as np,sys,time
import scipy.sparse as sp
from mb import fax
A=pickle.load(open(sys.argv[1],'rb')); B=A['blocks']; bcs=A['bcs']
import os
jt=float(sys.argv[3]); sweeps=int(sys.argv[4]); tlimit=float(sys.argv[5]) if len(sys.argv)>5 else 240
STATES=[]
STEPS=tuple(float(x) for x in os.environ.get("STEPS","0.6,0.3,0.15,0.07,0.03,0.01").split(","))
for fn in sys.argv[6:]:
    if fn.endswith('.npy'): STATES.append(np.load(fn).astype(np.float64)); continue
    Z=np.load(fn); STATES+=[Z[k].astype(np.float64) for k in Z.files if k!='P']
XA=np.vstack([b.X.reshape(-1,3) for b in B]); off=np.cumsum([0]+[b.X.reshape(-1,3).shape[0] for b in B])
import os
from numbering import get as _get
gid,inv=_get()
P=XA[gid].copy(); N=len(P)
cells=[]
for n,b in enumerate(B):
    ni,nj,nk=b.dims; I=inv[off[n]:off[n+1]].reshape(ni,nj,nk)
    c=[I[:-1,:-1,:-1],I[1:,:-1,:-1],I[1:,1:,:-1],I[:-1,1:,:-1],I[:-1,:-1,1:],I[1:,:-1,1:],I[1:,1:,1:],I[:-1,1:,1:]]
    cells.append(np.stack([x.ravel() for x in c],1))
cells=np.vstack(cells); nc=len(cells)
fixed=np.zeros(N,bool)
for bi,S in bcs.items():
    b=B[bi]; ni,nj,nk=b.dims; I=inv[off[bi]:off[bi+1]].reshape(ni,nj,nk)
    for s in S:
        a,hi=fax(s['f']); idx=[slice(None)]*3; idx[a]=-1 if hi else 0
        fixed[I[tuple(idx)][s['pr'][0]:s['pr'][1]+2,s['qr'][0]:s['qr'][1]+2].ravel()]=True
# ---- wall nodes that may slide on the original wall surface (caps/rims stay fixed)
from scipy.spatial import cKDTree
BM=np.load('/home/claude/cyc/bmesh.npz'); Xs=BM['pos']; Ts=BM['tris']; iov=BM['io_vertex']
PR=XA[gid]
fidx=np.where(fixed)[0]; dd,vi=cKDTree(Xs).query(PR[fidx]); assert dd.max()<1e-5
wallv=(iov[vi]<0)
# vertices touching a cap triangle are rims -> fixed
capv=np.zeros(len(Xs),bool); capv[Ts[BM['tri_surf']>=0][np.isin(BM['tri_surf'][BM['tri_surf']>=0],np.where(~np.char.startswith(BM['surf_names'],'wall'))[0])].ravel()]=True
wallv&=~capv[vi]
slide=np.zeros(N,bool); slide[fidx[wallv]]=True
sv=np.full(N,-1); sv[fidx]=vi                       # global node -> STL vertex
# 1-ring triangle fans (original surface) for projection, plus vertex normals
vt=[[] for _ in range(len(Xs))]
for t,(a,b_,c_) in enumerate(Ts): vt[a].append(t); vt[b_].append(t); vt[c_].append(t)
KF=max(len(x) for x in vt)
fan=np.full((len(Xs),KF),-1)
for v,l in enumerate(vt): fan[v,:len(l)]=l
tn=np.cross(Xs[Ts[:,1]]-Xs[Ts[:,0]],Xs[Ts[:,2]]-Xs[Ts[:,0]]); vn=np.zeros_like(Xs)
for k in range(3): np.add.at(vn,Ts[:,k],tn)
vn/=np.linalg.norm(vn,axis=1,keepdims=True)
def closest_on_tri(p,a,b,c):
    ab=b-a; ac=c-a; ap=p-a; d1=(ab*ap).sum(-1); d2=(ac*ap).sum(-1)
    bp=p-b; d3=(ab*bp).sum(-1); d4=(ac*bp).sum(-1); cp=p-c; d5=(ab*cp).sum(-1); d6=(ac*cp).sum(-1)
    va=d3*d6-d5*d4; vb=d5*d2-d1*d6; vc=d1*d4-d3*d2
    den=np.where(np.abs(va+vb+vc)<1e-30,1e-30,va+vb+vc); v=vb/den; w=vc/den
    res=a+ab*v[...,None]+ac*w[...,None]
    # edge / vertex regions
    m=(d1<=0)&(d2<=0); res=np.where(m[...,None],a,res)
    m=(d3>=0)&(d4<=d3); res=np.where(m[...,None],b,res)
    m=(d6>=0)&(d5<=d6); res=np.where(m[...,None],c,res)
    m=(vc<=0)&(d1>=0)&(d3<=0); t=d1/np.where(np.abs(d1-d3)<1e-30,1e-30,d1-d3); res=np.where(m[...,None],a+ab*t[...,None],res)
    m=(vb<=0)&(d2>=0)&(d6<=0); t=d2/np.where(np.abs(d2-d6)<1e-30,1e-30,d2-d6); res=np.where(m[...,None],a+ac*t[...,None],res)
    m=(va<=0)&((d4-d3)>=0)&((d5-d6)>=0); t=(d4-d3)/np.where(np.abs((d4-d3)+(d5-d6))<1e-30,1e-30,(d4-d3)+(d5-d6)); res=np.where(m[...,None],b+(c-b)*t[...,None],res)
    return res
def project(nodes,pts):
    v=sv[nodes]; F=fan[v]; ok=F>=0; Fc=np.where(ok,F,0)
    a=Xs[Ts[Fc,0]]; b=Xs[Ts[Fc,1]]; c=Xs[Ts[Fc,2]]
    q=closest_on_tri(pts[:,None,:],a,b,c); d=np.linalg.norm(q-pts[:,None,:],axis=-1); d[~ok]=np.inf
    return q[np.arange(len(nodes)),d.argmin(1)]
print('sliding wall nodes:',slide.sum(),flush=True)
NBc=[(0,1,3,4),(1,2,0,5),(2,3,1,6),(3,0,2,7),(4,7,5,0),(5,4,6,1),(6,5,7,2),(7,6,4,3)]
def cjac(P,C):
    X=P[C]; J=np.full(len(C),np.inf)
    for o,a,b_,cc in NBc:
        u=X[:,a]-X[:,o]; v=X[:,b_]-X[:,o]; w=X[:,cc]-X[:,o]
        J=np.minimum(J,np.einsum('ij,ij->i',np.cross(u,v),w)/(np.linalg.norm(u,axis=1)*np.linalg.norm(v,axis=1)*np.linalg.norm(w,axis=1)+1e-30))
    return J
rowsC=np.repeat(np.arange(nc),8); NC=sp.csr_matrix((np.ones(nc*8),(cells.ravel(),rowsC)),shape=(N,nc))
Dz=[np.zeros(1)]+STATES
def cjacS(P,C):
    Jm=None
    for Dk in Dz:
        Jk=cjac(P if Dk.size==1 else P+Dk,C); Jm=Jk if Jm is None else np.minimum(Jm,Jk)
    return Jm
J=cjacS(P,cells); J0=J.copy()
print(f'before: minJ {J.min():.4f}  p1 {np.percentile(J,1):.3f}  <0.1 {(J<0.1).sum()}  <0.2 {(J<0.2).sum()}  <{jt} {(J<jt).sum()}',flush=True)
t=time.time()
for sw in range(sweeps):
    bad=np.where(J<jt)[0]
    act=np.zeros(N,bool); act[cells[bad].ravel()]=True; act&=(~fixed)|slide
    ids=np.where(act)[0]
    if not len(ids): break
    # incident cells (padded to 8)
    sub=NC[ids]; cnt=np.diff(sub.indptr); KM=int(cnt.max()); inc=np.full((len(ids),KM),-1)
    for k in range(KM):
        m=cnt>k; inc[m,k]=sub.indices[sub.indptr[:-1][m]+k]
    # greedy colouring: nodes sharing a cell get different colours
    colour=np.full(len(ids),-1); pos=np.full(N,-1); pos[ids]=np.arange(len(ids))
    CC=cells[np.unique(inc[inc>=0])]
    for i in np.argsort(J[np.where(inc>=0,inc,0)].min(1)*(inc>=0).any(1)):   # worst first
        nbc=inc[i][inc[i]>=0]; nn=pos[cells[nbc].ravel()]; nn=nn[nn>=0]
        used=set(colour[nn].tolist()); c=0
        while c in used: c+=1
        colour[i]=c
    improved=0
    for c in range(colour.max()+1):
        sel=np.where(colour==c)[0]; nid=ids[sel]; ic=inc[sel]; valid=ic>=0; icc=np.where(valid,ic,0)
        TAU=0.06
        def fobj(Pt):
            Cc=cells[icc.ravel()]
            Jr=cjac(Pt,Cc).reshape(icc.shape); Jr[~valid]=np.inf               # reference grid
            Jc=np.minimum(Jr,cjacS(Pt,Cc).reshape(icc.shape)); Jc[~valid]=np.inf
            sc=-np.where(valid,np.exp(-Jc/TAU)+2*np.exp(-Jr/TAU),0).sum(1)
            return sc, Jc.min(1), Jr.min(1)
        def fmin(Pt): return fobj(Pt)[0]
        f0,m0,r0=fobj(P)
        J0ref=None
        L=np.array([np.linalg.norm(np.ptp(P[cells[i[i>=0]]].reshape(-1,3),axis=0)) for i in ic])*0.25
        g=np.zeros((len(nid),3)); hh=(0.01*L)[:,None]
        for d in range(3):
            Pp=P.copy(); Pp[nid,d]+=hh[:,0]; fp=fmin(Pp)
            Pm=P.copy(); Pm[nid,d]-=hh[:,0]; fm=fmin(Pm); g[:,d]=(fp-fm)/(2*hh[:,0])
        isl=slide[nid]
        if isl.any():
            nrm=vn[sv[nid[isl]]]; g[isl]-=(g[isl]*nrm).sum(1)[:,None]*nrm
        gn=np.linalg.norm(g,axis=1); dirn=g/np.maximum(gn,1e-12)[:,None]
        best=f0.copy(); bestP=P[nid].copy()
        for step in STEPS:
            Pt=P.copy(); newp=P[nid]+dirn*(step*L)[:,None]
            if isl.any(): newp[isl]=project(nid[isl],newp[isl])
            Pt[nid]=newp; ft,mt,rt=fobj(Pt)
            upd=(ft>best+1e-12)&(mt>=np.minimum(m0,0.3)-1e-12)&(rt>=np.minimum(r0,float(os.environ.get('REFMIN','0.05')))-1e-12); best[upd]=ft[upd]; bestP[upd]=Pt[nid][upd]
        chg=best>f0+1e-9; P[nid[chg]]=bestP[chg]; improved+=int(chg.sum())
    aff=np.unique(NC[ids].indices); J[aff]=cjacS(P,cells[aff])
    print(f'  sweep {sw:2d}: active {len(ids):6d} colours {colour.max()+1:3d} improved {improved:6d}  minJ {J.min():.4f}  <0.1 {(J<0.1).sum()}  <0.2 {(J<0.2).sum()}  neg {(J<=0).sum()}  ({time.time()-t:.0f}s)',flush=True)
    if time.time()-t>tlimit: break
XN=P[inv]
for n,b in enumerate(B): b.X=np.ascontiguousarray(XN[off[n]:off[n+1]].reshape(b.X.shape))
Jref=cjac(P,cells); print(f'reference grid: minJ {Jref.min():.4f} neg {(Jref<=0).sum()} <0.05 {(Jref<0.05).sum()}')
mv=np.linalg.norm(P[fixed]-XA[gid][fixed],axis=1)
print(f'after: minJ {J.min():.4f}  p1 {np.percentile(J,1):.3f}  <0.1 {(J<0.1).sum()}  <0.2 {(J<0.2).sum()}  neg {(J<=0).sum()}  wall nodes slid: {(mv>1e-9).sum()} (max {mv.max():.3f} mm); caps moved {np.abs(P[fixed&~slide]-XA[gid][fixed&~slide]).max():.1e}')
pickle.dump(A,open(sys.argv[2],'wb'))
low=np.where(J<0.1)[0]; nf=fixed[cells[low]].sum(1)
print('cells J<0.1 by number of fixed (boundary) corners:',np.bincount(nf,minlength=9).tolist())
wf=low[nf>=4]
if len(wf):
    # for wall-attached bad cells, is the boundary quad itself distorted?  measure min corner angle quality of the fixed face
    print('   example worst cells (fixed corners, J):',[(int(fixed[cells[c]].sum()),round(float(J[c]),3)) for c in low[np.argsort(J[low])[:8]]])
