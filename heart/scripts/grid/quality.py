import pickle,numpy as np,sys,json; sys.path.insert(0,'/home/claude/hg')
import mb, heart, oblayer
from scipy.spatial import cKDTree
A=pickle.load(open('assembled.pkl','rb')); blocks=A['blocks']; bcs=A['bcs']
H=heart.Build(ds=1.0,scale=0.85); H.left(); H.right(); H.exclude(0.45)
# wall nodes per domain
W={'L':[],'R':[]}
for b,S in bcs.items():
    B=blocks[b]
    for s in S:
        if s['type']!='wall': continue
        F=B.face(s['f'])[s['pr'][0]:s['pr'][1]+2,s['qr'][0]:s['qr'][1]+2]; W[B.dom].append(F.reshape(-1,3))
WL=np.unique(np.vstack(W['L']).round(9),axis=0); WR=np.unique(np.vstack(W['R']).round(9),axis=0)
dLR=cKDTree(WR).query(WL)[0]
# penetration: any node of one domain inside the other's wall surface?  use grid nodes vs other implicit (conservative) + nearest distance
allL=np.vstack([B.X.reshape(-1,3) for B in blocks if B.dom=='L']); allR=np.vstack([B.X.reshape(-1,3) for B in blocks if B.dom=='R'])
inR=(H.FR(allL)<0).sum(); inL=(H.FL(allR)<0).sum()
print(f'L-R wall min distance {dLR.min():.2f} mm (1st pct {np.percentile(dLR,1):.2f}); L nodes inside R-anatomy {inR}, R nodes inside L-anatomy {inL}')
# quality
J=np.concatenate([mb.hex_jac(B.X).ravel() for B in blocks])
def ar(X):
    e=[np.linalg.norm(np.diff(X,axis=a),axis=-1) for a in range(3)]
    e=[x[:, :-1, :-1] if a==0 else x[:-1, :, :-1] if a==1 else x[:-1,:-1,:] for a,x in enumerate(e)]
    m=np.stack(e); return m.max(0)/np.maximum(m.min(0),1e-12)
AR=np.concatenate([ar(B.X).ravel() for B in blocks])
def growth(X):
    g=[]
    for a in range(3):
        e=np.linalg.norm(np.diff(X,axis=a),axis=-1)
        if e.shape[a]<2: continue
        r=np.take(e,range(1,e.shape[a]),axis=a)/np.maximum(np.take(e,range(0,e.shape[a]-1),axis=a),1e-12)
        g.append(np.maximum(r,1/np.maximum(r,1e-12)).ravel())
    return np.concatenate(g)
G=np.concatenate([growth(B.X) for B in blocks])
q=dict(min_scaled_jacobian=float(J.min()),jac_p1=float(np.percentile(J,1)),jac_median=float(np.median(J)),
       frac_jac_below_0p2=float((J<0.2).mean()),aspect_max=float(AR.max()),aspect_p99=float(np.percentile(AR,99)),
       growth_max=float(G.max()),growth_p99=float(np.percentile(G,99)),LR_min_dist_mm=float(dLR.min()))
print(json.dumps(q,indent=1))
json.dump(q,open('quality.json','w'),indent=1)
