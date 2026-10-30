"""Diffusion (Laplacian) mesh morphing on the real multiblock grid - the standard CFD morpher.
Edge weight = 1/l^2 * (d_ref/d_wall)^alpha  (cells near walls are stiff -> move with the wall).
Global node numbering (shared interface nodes once) -> block point matching preserved."""
import numpy as np,sys,pickle,time
import scipy.sparse as sp, pyamg
from scipy.spatial import cKDTree
sys.path.insert(0,'/home/claude/hg'); from mb import fax
class DiffMorph:
    def __init__(self,A,Xb,alpha=1.5):
        B=A['blocks']; bcs=A['bcs']; self.B=B
        XA=np.vstack([b.X.reshape(-1,3) for b in B]); self.off=np.cumsum([0]+[b.X.reshape(-1,3).shape[0] for b in B])
        from numbering import get as _get
        gid,inv=_get(); self.inv=inv; P=XA[gid]; self.P=P; N=len(P)
        cells=[]
        for n,b in enumerate(B):
            ni,nj,nk=b.dims; I=self.inv[self.off[n]:self.off[n+1]].reshape(ni,nj,nk)
            c=[I[:-1,:-1,:-1],I[1:,:-1,:-1],I[1:,1:,:-1],I[:-1,1:,:-1],I[:-1,:-1,1:],I[1:,:-1,1:],I[1:,1:,1:],I[:-1,1:,1:]]
            cells.append(np.stack([x.ravel() for x in c],1))
        self.cells=np.vstack(cells)
        fixed=np.zeros(N,bool)
        for bi,S in bcs.items():
            b=B[bi]; ni,nj,nk=b.dims; I=self.inv[self.off[bi]:self.off[bi+1]].reshape(ni,nj,nk)
            for s in S:
                a,hi=fax(s['f']); idx=[slice(None)]*3; idx[a]=-1 if hi else 0
                fixed[I[tuple(idx)][s['pr'][0]:s['pr'][1]+2,s['qr'][0]:s['qr'][1]+2].ravel()]=True
        self.bnd=np.where(fixed)[0]; self.free=np.where(~fixed)[0]
        # map boundary grid nodes -> STL vertices
        d,i=cKDTree(Xb).query(P[self.bnd]); assert d.max()<1e-6, d.max(); self.b2v=i
        E=[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
        ea=np.concatenate([self.cells[:,a] for a,b_ in E]); eb=np.concatenate([self.cells[:,b_] for a,b_ in E])
        key=np.minimum(ea,eb).astype(np.int64)*N+np.maximum(ea,eb); _,u=np.unique(key,return_index=True); ea,eb=ea[u],eb[u]
        l=np.linalg.norm(P[ea]-P[eb],axis=1); mid=0.5*(P[ea]+P[eb])
        dw=cKDTree(P[self.bnd]).query(mid)[0]
        w=1/np.maximum(l,1e-4)**2*(1.0/np.maximum(dw,0.05))**alpha
        W=sp.coo_matrix((np.r_[w,w],(np.r_[ea,eb],np.r_[eb,ea])),shape=(N,N)).tocsr()
        Lap=(sp.diags(np.asarray(W.sum(1)).ravel())-W).tocsr()
        self.Lff=Lap[self.free][:,self.free].tocsr(); self.Lfb=Lap[self.free][:,self.bnd].tocsr()
        t=time.time(); self.ml=pyamg.smoothed_aggregation_solver(self.Lff,symmetry='symmetric'); print('  AMG setup %.0fs, free nodes %d'%(time.time()-t,len(self.free)),flush=True)
        self.x0=[None]*3
    def morph(self,Db,tol=1e-8):
        D=np.zeros((len(self.P),3)); D[self.bnd]=Db[self.b2v]
        for k in range(3):
            rhs=-(self.Lfb@D[self.bnd,k]); D[self.free,k]=self.ml.solve(rhs,x0=self.x0[k],tol=tol,accel='cg',maxiter=200); self.x0[k]=D[self.free,k]
        return D
    def jac(self,Pn):
        X=Pn[self.cells]; J=np.full(len(self.cells),np.inf)
        for o,a,b_,cc in [(0,1,3,4),(1,2,0,5),(2,3,1,6),(3,0,2,7),(4,7,5,0),(5,4,6,1),(6,5,7,2),(7,6,4,3)]:
            u=X[:,a]-X[:,o]; v=X[:,b_]-X[:,o]; w=X[:,cc]-X[:,o]
            J=np.minimum(J,np.einsum('ij,ij->i',np.cross(u,v),w)/(np.linalg.norm(u,axis=1)*np.linalg.norm(v,axis=1)*np.linalg.norm(w,axis=1)+1e-30))
        return J
