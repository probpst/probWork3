"""Local mesh optimisation: maximise the worst scaled Jacobian around each free node (boundary fixed,
shared interface nodes numbered once -> point matching preserved). Independent sets via greedy colouring."""
import pickle,numpy as np,sys,time
import scipy.sparse as sp
from mb import fax
A=pickle.load(open(sys.argv[1],'rb')); B=A['blocks']; bcs=A['bcs']
jt=float(sys.argv[3]); sweeps=int(sys.argv[4]); tlimit=float(sys.argv[5]) if len(sys.argv)>5 else 240
XA=np.vstack([b.X.reshape(-1,3) for b in B]); off=np.cumsum([0]+[b.X.reshape(-1,3).shape[0] for b in B])
key=np.round(XA/1e-6).astype(np.int64); _,gid,inv=np.unique(key,axis=0,return_index=True,return_inverse=True); inv=inv.ravel()
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
NBc=[(0,1,3,4),(1,2,0,5),(2,3,1,6),(3,0,2,7),(4,7,5,0),(5,4,6,1),(6,5,7,2),(7,6,4,3)]
def cjac(P,C):
    X=P[C]; J=np.full(len(C),np.inf)
    for o,a,b_,cc in NBc:
        u=X[:,a]-X[:,o]; v=X[:,b_]-X[:,o]; w=X[:,cc]-X[:,o]
        J=np.minimum(J,np.einsum('ij,ij->i',np.cross(u,v),w)/(np.linalg.norm(u,axis=1)*np.linalg.norm(v,axis=1)*np.linalg.norm(w,axis=1)+1e-30))
    return J
rowsC=np.repeat(np.arange(nc),8); NC=sp.csr_matrix((np.ones(nc*8),(cells.ravel(),rowsC)),shape=(N,nc))
J=cjac(P,cells); J0=J.copy()
print(f'before: minJ {J.min():.4f}  p1 {np.percentile(J,1):.3f}  <0.1 {(J<0.1).sum()}  <0.2 {(J<0.2).sum()}  <{jt} {(J<jt).sum()}',flush=True)
t=time.time()
for sw in range(sweeps):
    bad=np.where(J<jt)[0]
    act=np.zeros(N,bool); act[cells[bad].ravel()]=True; act&=~fixed
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
            Jc=cjac(Pt,cells[icc.ravel()]).reshape(icc.shape); Jc[~valid]=np.inf
            return -np.where(valid,np.exp(-Jc/TAU),0).sum(1), Jc.min(1)
        def fmin(Pt): return fobj(Pt)[0]
        f0,m0=fobj(P)
        L=np.array([np.linalg.norm(np.ptp(P[cells[i[i>=0]]].reshape(-1,3),axis=0)) for i in ic])*0.25
        g=np.zeros((len(nid),3)); hh=(0.01*L)[:,None]
        for d in range(3):
            Pp=P.copy(); Pp[nid,d]+=hh[:,0]; fp=fmin(Pp)
            Pm=P.copy(); Pm[nid,d]-=hh[:,0]; fm=fmin(Pm); g[:,d]=(fp-fm)/(2*hh[:,0])
        gn=np.linalg.norm(g,axis=1); dirn=g/np.maximum(gn,1e-12)[:,None]
        best=f0.copy(); bestP=P[nid].copy()
        for step in (0.6,0.3,0.15,0.07,0.03,0.01):
            Pt=P.copy(); Pt[nid]=P[nid]+dirn*(step*L)[:,None]; ft,mt=fobj(Pt)
            upd=(ft>best+1e-12)&(mt>0.002)&(mt>=np.minimum(0.9*m0,0.3)); best[upd]=ft[upd]; bestP[upd]=Pt[nid][upd]
        chg=best>f0+1e-9; P[nid[chg]]=bestP[chg]; improved+=int(chg.sum())
    aff=np.unique(NC[ids].indices); J[aff]=cjac(P,cells[aff])
    print(f'  sweep {sw:2d}: active {len(ids):6d} colours {colour.max()+1:3d} improved {improved:6d}  minJ {J.min():.4f}  <0.1 {(J<0.1).sum()}  <0.2 {(J<0.2).sum()}  neg {(J<=0).sum()}  ({time.time()-t:.0f}s)',flush=True)
    if time.time()-t>tlimit: break
XN=P[inv]
for n,b in enumerate(B): b.X=np.ascontiguousarray(XN[off[n]:off[n+1]].reshape(b.X.shape))
print(f'after: minJ {J.min():.4f}  p1 {np.percentile(J,1):.3f}  <0.1 {(J<0.1).sum()}  <0.2 {(J<0.2).sum()}  neg {(J<=0).sum()}  boundary moved {np.abs(P[fixed]-XA[gid][fixed]).max():.1e}')
pickle.dump(A,open(sys.argv[2],'wb'))
low=np.where(J<0.1)[0]; nf=fixed[cells[low]].sum(1)
print('cells J<0.1 by number of fixed (boundary) corners:',np.bincount(nf,minlength=9).tolist())
wf=low[nf>=4]
if len(wf):
    # for wall-attached bad cells, is the boundary quad itself distorted?  measure min corner angle quality of the fixed face
    print('   example worst cells (fixed corners, J):',[(int(fixed[cells[c]].sum()),round(float(J[c]),3)) for c in low[np.argsort(J[low])[:8]]])
