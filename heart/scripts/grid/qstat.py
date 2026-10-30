import pickle,numpy as np,sys; sys.path.insert(0,'/home/claude/hg'); import mb
d=pickle.load(open(sys.argv[1],'rb')); B=d['blocks']
J=np.concatenate([mb.hex_jac(b.X).ravel() for b in B])
NO=np.concatenate([mb.face_nonortho(b.X).ravel() for b in B])
def ar(X):
    e=[np.linalg.norm(np.diff(X,axis=a),axis=-1) for a in range(3)]
    e=[x[:,:-1,:-1] if a==0 else x[:-1,:,:-1] if a==1 else x[:-1,:-1,:] for a,x in enumerate(e)]
    m=np.stack(e); return m.max(0)/np.maximum(m.min(0),1e-12)
AR=np.concatenate([ar(b.X).ravel() for b in B])
def gr(X):
    g=[]
    for a in range(3):
        e=np.linalg.norm(np.diff(X,axis=a),axis=-1)
        if e.shape[a]<2: continue
        r=np.take(e,range(1,e.shape[a]),axis=a)/np.maximum(np.take(e,range(0,e.shape[a]-1),axis=a),1e-12)
        g.append(np.maximum(r,1/np.maximum(r,1e-12)).ravel())
    return np.concatenate(g)
G=np.concatenate([gr(b.X) for b in B])
print(f'cells {sum(int(np.prod(np.array(b.dims)-1)) for b in B):,} neg {int((J<=0).sum())}')
print(f'  scaled Jacobian: min {J.min():.3f}  p1 {np.percentile(J,1):.3f}  median {np.median(J):.3f}')
print(f'  non-orthogonality deg: median {np.median(NO):.1f}  p95 {np.percentile(NO,95):.1f}  p99.9 {np.percentile(NO,99.9):.1f}  max {NO.max():.1f}')
print(f'  aspect ratio: median {np.median(AR):.2f}  p95 {np.percentile(AR,95):.2f}  p99 {np.percentile(AR,99):.2f}  max {AR.max():.1f}')
print(f'  neighbour size ratio: p95 {np.percentile(G,95):.3f}  p99 {np.percentile(G,99):.3f}  max {G.max():.2f}')
