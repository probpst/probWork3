"""Boundary triangulation of the v3 PLOT3D grid with a vertex -> grid-node map.
Every STL vertex is a unique grid boundary node; triangle order is fixed."""
import pickle,numpy as np,sys; sys.path.insert(0,'/home/claude/hg')
from mb import fax, nkey
import os; A=pickle.load(open(os.environ.get('GRID','/home/claude/hg/assembled.pkl'),'rb')); blocks=A['blocks']; bcs=A['bcs']
key2id={}; pos=[]; dom=[]; tris=[]; tri_surf=[]; surf_names=[]; nodes_of={}
h_first={}   # vertex -> min wall-normal first-cell height
for bi,S in bcs.items():
    B=blocks[bi]
    for s in S:
        name='wall_'+('left' if B.dom=='L' else 'right') if s['type']=='wall' else s['label']
        if name not in surf_names: surf_names.append(name)
        sid=surf_names.index(name)
        a,hi=fax(s['f']); idx=[slice(None)]*3; idx[a]=-1 if hi else 0
        F=B.X[tuple(idx)][s['pr'][0]:s['pr'][1]+2,s['qr'][0]:s['qr'][1]+2]
        idx[a]=-2 if hi else 1
        G=B.X[tuple(idx)][s['pr'][0]:s['pr'][1]+2,s['qr'][0]:s['qr'][1]+2]
        # grid indices of these face nodes
        ip=[x for x in range(3) if x!=a]
        K=nkey(F); ids=np.zeros(F.shape[:2],int)
        for p in range(F.shape[0]):
            for q in range(F.shape[1]):
                k=tuple(K[p,q])
                if k not in key2id: key2id[k]=len(pos); pos.append(F[p,q]); dom.append(B.dom)
                v=key2id[k]; ids[p,q]=v
                ijk=[0,0,0]; ijk[a]=B.dims[a] if hi else 1; ijk[ip[0]]=s['pr'][0]+p+1; ijk[ip[1]]=s['qr'][0]+q+1
                nodes_of.setdefault(v,set()).add((bi+1,*ijk))
                hh=float(np.linalg.norm(G[p,q]-F[p,q])); h_first[v]=min(h_first.get(v,1e9),hh)
        flip=np.einsum('pqi,pqi->pq',np.cross(F[1:,:-1]-F[:-1,:-1],F[:-1,1:]-F[:-1,:-1]),(G-F)[:-1,:-1])>0
        use1=np.linalg.norm(F[1:,1:]-F[:-1,:-1],axis=-1)<=np.linalg.norm(F[1:,:-1]-F[:-1,1:],axis=-1)
        for p in range(F.shape[0]-1):
            for q in range(F.shape[1]-1):
                a0,a1,a2,a3=ids[p,q],ids[p+1,q],ids[p+1,q+1],ids[p,q+1]
                tt=[(a0,a1,a2),(a0,a2,a3)] if use1[p,q] else [(a0,a1,a3),(a1,a2,a3)]
                for t in tt: tris.append(t[::-1] if flip[p,q] else t); tri_surf.append(sid)
pos=np.array(pos); tris=np.array(tris,np.int32); dom=np.array(dom)
h=np.array([h_first[i] for i in range(len(pos))])
surf_of_vertex=np.full(len(pos),-1)
for t,sidx in zip(tris,tri_surf):
    for v in t:
        if not surf_names[sidx].startswith('wall'): surf_of_vertex[v]=sidx
np.savez(os.environ.get('OUT','/home/claude/cyc/bmesh.npz'),pos=pos,tris=tris,tri_surf=np.array(tri_surf,np.int16),dom=dom,h=h,
         surf_names=np.array(surf_names),io_vertex=surf_of_vertex)
with open(os.environ.get('MAP','/home/claude/cyc/vertex_to_grid_nodes.txt'),'w') as f:
    f.write('# STL vertex index (0-based)  ->  PLOT3D boundary node(s) as block i j k (1-based); several entries = shared node\n')
    for v in range(len(pos)): f.write(f'{v} '+' '.join(f'{b}:{i},{j},{k}' for b,i,j,k in sorted(nodes_of[v]))+'\n')
print('vertices',len(pos),'triangles',len(tris),'surfaces',len(surf_names),'L verts',(dom=="L").sum(),'R verts',(dom=="R").sum())
print('first-cell height mm: min %.3f p5 %.3f median %.3f'%(h.min(),np.percentile(h,5),np.median(h)))
