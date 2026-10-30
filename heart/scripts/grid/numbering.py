"""Stable global node numbering for every grid derived from assembled_o (same block structure):
ids are fixed once from the original coordinates, so all later grids / morph states share them."""
import numpy as np,pickle,os
F='/home/claude/hg/node_numbering.npz'
def get():
    if not os.path.exists(F):
        B=pickle.load(open('/home/claude/hg/assembled_base.pkl','rb'))['blocks']
        XA=np.vstack([b.X.reshape(-1,3) for b in B]); key=np.round(XA/1e-6).astype(np.int64)
        _,gid,inv=np.unique(key,axis=0,return_index=True,return_inverse=True)
        np.savez(F,gid=gid,inv=inv.ravel())
    Z=np.load(F); return Z['gid'],Z['inv']
