"""Multiblock infrastructure: blocks, face access, automatic point-matched connectivity,
BC-surface decomposition, reference-format .bcs / PLOT3D writers, validation."""
import numpy as np
from collections import defaultdict

FACES=['imin','imax','jmin','jmax','kmin','kmax']
def fax(f): return 'ijk'.index(f[0]), f[1:]=='max'
def inplane(a): return [x for x in range(3) if x!=a]

class Block:
    def __init__(self,name,X,dom,kind='core',comp=None):
        self.name=name; self.X=np.ascontiguousarray(X,float); self.dom=dom; self.kind=kind; self.comp=comp
        self.tags=[]   # dicts: face, (p0,p1),(q0,q1) cell ranges (0-based, inclusive), type, label
    @property
    def dims(self): return self.X.shape[:3]
    def face(self,f):
        a,hi=fax(f); idx=[slice(None)]*3; idx[a]=-1 if hi else 0
        return self.X[tuple(idx)]           # (n_p, n_q, 3), p,q = in-plane axes in order
    def tag(self,f,typ,label,pr=None,qr=None):
        F=self.face(f); pr=pr or (0,F.shape[0]-2); qr=qr or (0,F.shape[1]-2)
        self.tags.append(dict(face=f,pr=pr,qr=qr,type=typ,label=label))

def nkey(P,tol=1e-6):
    return np.round(np.asarray(P)/tol).astype(np.int64)

def cell_keys(F):
    K=nkey(F); n0,n1=K.shape[0]-1,K.shape[1]-1
    c=np.stack([K[:-1,:-1],K[1:,:-1],K[1:,1:],K[:-1,1:]],2)       # n0,n1,4,3
    # order-independent key: sort the 4 node keys lexicographically
    c=c.reshape(-1,4,3); v=c[:,:,0]*1_000_003**2+c[:,:,1]*1_000_003+c[:,:,2]
    v=np.sort(v,axis=1)
    return [tuple(r) for r in v], (n0,n1)

def rects(mask):
    """greedy decomposition of boolean 2-D mask into rectangles (inclusive cell ranges)"""
    m=mask.copy(); out=[]
    while m.any():
        p,q=np.argwhere(m)[0]
        q1=q
        while q1+1<m.shape[1] and m[p,q1+1]: q1+=1
        p1=p
        while p1+1<m.shape[0] and m[p1+1,q:q1+1].all(): p1+=1
        m[p:p1+1,q:q1+1]=False; out.append(((p,p1),(q,q1)))
    return out

def connectivity(blocks):
    """find every point-coincident block-to-block interface. returns list of interfaces
    and per (block,face) cell ownership arrays."""
    table=defaultdict(list); shapes={}
    for b,B in enumerate(blocks):
        for f in FACES:
            ks,sh=cell_keys(B.face(f)); shapes[(b,f)]=sh
            for n,k in enumerate(ks): table[k].append((b,f,n))
    part={}   # (b,f) -> array of partner id (-1 none)
    pairs={}
    for (b,f),sh in shapes.items(): part[(b,f)]=np.full(sh,-1,int)
    plist=[]
    bad=[]
    for k,L in table.items():
        if len(L)==1: continue
        if len(L)>2: bad.append(L); continue
        (b1,f1,n1),(b2,f2,n2)=L
        for (ba,fa,na),(bb,fb,nb) in [((b1,f1,n1),(b2,f2,n2)),((b2,f2,n2),(b1,f1,n1))]:
            key=(ba,fa,bb,fb)
            if key not in pairs: pairs[key]=len(plist); plist.append(key)
            sh=shapes[(ba,fa)]; part[(ba,fa)][np.unravel_index(na,sh)]=pairs[key]
    if bad: raise RuntimeError(f'{len(bad)} cell faces shared by >2 blocks')
    ifaces=[]
    for pid,(ba,fa,bb,fb) in enumerate(plist):
        M=part[(ba,fa)]==pid
        for (p0,p1),(q0,q1) in rects(M):
            ifaces.append(dict(b=ba,f=fa,pr=(p0,p1),qr=(q0,q1),db=bb,df=fb))
    for I in ifaces: I.update(donor_range(blocks,I))
    return ifaces,part

def donor_range(blocks,I):
    import itertools
    A=blocks[I['b']].face(I['f']); D=blocks[I['db']].face(I['df'])
    (p0,p1),(q0,q1)=I['pr'],I['qr']
    corners=[(p0,q0),(p1+1,q0),(p0,q1+1)]
    Dk=nkey(D).reshape(-1,3); lut=defaultdict(list)
    for n,k in enumerate(map(tuple,Dk)): lut[k].append(n)
    cand=[[np.unravel_index(n,D.shape[:2]) for n in lut[tuple(nkey(A[c]))]] for c in corners]
    Asub=A[p0:p1+2,q0:q1+2]; best=None
    for (a0,b0),(a1,b1),(a2,b2) in itertools.product(*cand):
        if a1!=a0 and b1==b0 and a2==a0 and b2!=b0: order='same'; dp=(a0,a1); dq=(b0,b2)
        elif b1!=b0 and a1==a0 and b2==b0 and a2!=a0: order='transposed'; dp=(a0,a2); dq=(b0,b1)
        else: continue
        sp=1 if dp[1]>=dp[0] else -1; sq=1 if dq[1]>=dq[0] else -1
        Dsub=D[dp[0]:(dp[1]+sp if dp[1]+sp>=0 else None):sp, dq[0]:(dq[1]+sq if dq[1]+sq>=0 else None):sq]
        if order=='transposed': Dsub=Dsub.transpose(1,0,2)
        if Dsub.shape!=Asub.shape: continue
        err=np.abs(Asub-Dsub).max()
        if best is None or err<best['err']: best=dict(dpr=dp,dqr=dq,order=order,err=err)
    if best is None: raise RuntimeError('no consistent donor mapping')
    return best

def bc_surfaces(blocks,part):
    """all non-interface cell faces -> rectangles with physical type"""
    out=defaultdict(list)
    for b,B in enumerate(blocks):
        for f in FACES:
            free=part[(b,f)]<0
            if not free.any(): continue
            typ=np.full(free.shape,'wall',object); lab=np.full(free.shape,'wall',object)
            for t in B.tags:
                if t['face']!=f: continue
                sl=(slice(t['pr'][0],t['pr'][1]+1),slice(t['qr'][0],t['qr'][1]+1))
                typ[sl]=t['type']; lab[sl]=t['label']
            code=np.full(free.shape,'',object); code[free]=(typ[free]+'|'+lab[free])
            for c in set(code[free]):
                for (p0,p1),(q0,q1) in rects(code==c):
                    t,l=c.split('|'); out[b].append(dict(f=f,pr=(p0,p1),qr=(q0,q1),type=t,label=l))
    return out

# ----------------------------- quality ---------------------------------
def hex_jac(X):
    """min corner Jacobian (scaled by edge lengths) and volume per cell"""
    c=[X[:-1,:-1,:-1],X[1:,:-1,:-1],X[1:,1:,:-1],X[:-1,1:,:-1],X[:-1,:-1,1:],X[1:,:-1,1:],X[1:,1:,1:],X[:-1,1:,1:]]
    # corner: (origin, e_i, e_j, e_k) neighbours
    nb=[(0,1,3,4),(1,2,0,5),(2,3,1,6),(3,0,2,7),(4,7,5,0),(5,4,6,1),(6,5,7,2),(7,6,4,3)]
    J=[]
    for o,a,b,cc in nb:
        u=c[a]-c[o];v=c[b]-c[o];w=c[cc]-c[o]
        det=np.einsum('...i,...i',np.cross(u,v),w)
        J.append(det/(np.linalg.norm(u,axis=-1)*np.linalg.norm(v,axis=-1)*np.linalg.norm(w,axis=-1)+1e-30))
    J=np.stack(J); return J.min(0)

def hex_vol(X):
    # 5-tet decomposition-free: use divergence-free approximation via 2 prisms (sufficient for totals)
    a=X[:-1,:-1,:-1];b=X[1:,:-1,:-1];c=X[1:,1:,:-1];d=X[:-1,1:,:-1]
    e=X[:-1,:-1,1:];f=X[1:,:-1,1:];g=X[1:,1:,1:];h=X[:-1,1:,1:]
    def tet(p,q,r,s): return np.einsum('...i,...i',np.cross(q-p,r-p),s-p)/6
    return tet(a,b,d,e)+tet(b,c,d,g)+tet(b,d,e,g)+tet(b,f,g,e)+tet(d,g,h,e)

# ----------------------------- writers ---------------------------------
def _fortran_E(v):
    s='%.15E'%v
    if v==0: return '0.0000000000000000E+00'
    neg=s[0]=='-'; s=s.lstrip('-'); m,e=s.split('E'); d=m.replace('.','')
    return ('-' if neg else '')+'0.'+d+'E%+03d'%(int(e)+1)

def write_grd_ascii(path,blocks):
    with open(path,'w') as fh:
        fh.write('%12d\n'%len(blocks))
        for B in blocks: fh.write('%12d%12d%12d\n'%B.dims)
        for B in blocks:
            for c in range(3):
                a=B.X[...,c].transpose(2,1,0).ravel()
                strs=['%28s'%_fortran_E(float(v)) for v in a]
                for i in range(0,len(strs),4): fh.write(''.join(strs[i:i+4])+'\n')

def write_grd_binary(path,blocks):
    import struct
    def rec(fh,b): fh.write(struct.pack('<i',len(b))); fh.write(b); fh.write(struct.pack('<i',len(b)))
    with open(path,'wb') as fh:
        rec(fh,struct.pack('<i',len(blocks)))
        rec(fh,np.array([B.dims for B in blocks],'<i4').tobytes())
        for B in blocks:
            rec(fh,np.concatenate([B.X[...,c].transpose(2,1,0).ravel() for c in range(3)]).astype('<f8').tobytes())

BCCODE={'wall':(20,'no-slip wall'),'inlet':(9,'Inlet'),'outlet':(12,'Outlet')}

def face_ranges(B,f,pr,qr):
    """cell ranges on face -> node index ranges (1-based) for i,j,k"""
    a,hi=fax(f); ip=inplane(a); dims=B.dims; R=[None]*3
    R[a]=(dims[a],dims[a]) if hi else (1,1)
    R[ip[0]]=(pr[0]+1,pr[1]+2); R[ip[1]]=(qr[0]+1,qr[1]+2)
    return R

def write_bcs(path,blocks,ifaces,bcs):
    L=[]
    def ln(vals,cmt):
        s=''.join(['%8d'%vals[0]]+['%5d'%v for v in vals[1:]])
        L.append(s.ljust(24)+(cmt if cmt.startswith('!') else '! '+cmt))
    byb=defaultdict(list)
    for I in ifaces: byb[I['b']].append(('if',I))
    for b,S in bcs.items():
        for s in S: byb[b].append(('bc',s))
    for b,B in enumerate(blocks):
        items=sorted(byb[b],key=lambda t:(FACES.index(t[1]['f']),t[1]['pr'],t[1]['qr']))
        ln([len(items)],'number of boundary surfaces')
        for kind,s in items:
            a,hi=fax(s['f']); d=(a+1)*(-1 if hi else 1)
            R=face_ranges(B,s['f'],s['pr'],s['qr'])
            if kind=='bc':
                code,nm=BCCODE[s['type']]
                ln([code],f"!*** {s['label']} = type #{code} {nm}")
                ln([d],'coordinate direction normal to surface')
                for c,ax in zip(R,'ijk'): ln(list(c),f'start, end in {ax}-direction')
                ln([0],'flag for order of first derivative bc')
            else:
                D=blocks[s['db']]; da,dhi=fax(s['df']); dd=(da+1)*(-1 if dhi else 1)
                ln([92],f"!*** patched MB #{s['db']+1}= type #92")
                ln([d],'coordinate direction normal to surface')
                for c,ax in zip(R,'ijk'): ln(list(c),f'start, end in {ax}-direction')
                ln([0],'flag for order of first derivative bc')
                ln([s['db']+1],'donor mesh')
                ln([dd],'coordinate direction normal to surface')
                DR=[None]*3; DR[da]=(D.dims[da],)*2 if dhi else (1,1)
                ip=inplane(da)
                DR[ip[0]]=(s['dpr'][0]+1,s['dpr'][1]+1); DR[ip[1]]=(s['dqr'][0]+1,s['dqr'][1]+1)
                for c,ax in zip(DR,'ijk'): ln(list(c),f'donor mesh start, end in {ax}-direction')
    open(path,'w',newline='').write('\r\n'.join(L)+'\r\n')

def face_nonortho(X):
    """max face non-orthogonality angle (deg) per cell: angle between face normal and
    the line joining the centroids of the two cells sharing that face (internal faces of a block)"""
    c=0.125*sum(X[i:i+X.shape[0]-1 or None,j:j+X.shape[1]-1 or None,k:k+X.shape[2]-1 or None]
                for i in (0,1) for j in (0,1) for k in (0,1))
    out=np.zeros(c.shape[:3])
    for a in range(3):
        if c.shape[a]<2: continue
        sl=[slice(None)]*3
        sl[a]=slice(0,c.shape[a]-1); c0=c[tuple(sl)]
        sl[a]=slice(1,c.shape[a]); c1=c[tuple(sl)]
        v=c1-c0; v/=np.maximum(np.linalg.norm(v,axis=-1,keepdims=True),1e-30)
        # face between the cells: quad of nodes at index a+1 in node space
        idx=[slice(None)]*3
        idx[a]=slice(1,X.shape[a]-1)
        for b in range(3):
            if b!=a: idx[b]=slice(None)
        Fn=X[tuple(idx)]
        p=[slice(None)]*3; q=[slice(None)]*3
        o1=[x for x in range(3) if x!=a][0]; o2=[x for x in range(3) if x!=a][1]
        sA=[slice(None)]*3; sB=[slice(None)]*3
        sA[o1]=slice(1,None); sB[o1]=slice(0,-1)
        d1=Fn[tuple(sA)][...,:]-Fn[tuple(sB)][...,:]
        sA=[slice(None)]*3; sB=[slice(None)]*3
        sA[o2]=slice(1,None); sB[o2]=slice(0,-1)
        d2=Fn[tuple(sA)]-Fn[tuple(sB)]
        s1=[slice(None)]*3; s1[o2]=slice(0,-1); s2=[slice(None)]*3; s2[o1]=slice(0,-1)
        n=np.cross(d1[tuple(s1)],d2[tuple(s2)])
        n/=np.maximum(np.linalg.norm(n,axis=-1,keepdims=True),1e-30)
        ang=np.degrees(np.arccos(np.clip(np.abs(np.einsum('...i,...i',n,v)),0,1)))
        sl=[slice(None)]*3; sl[a]=slice(0,c.shape[a]-1); o=out[tuple(sl)]
        np.maximum(o,ang,out=o)
        sl[a]=slice(1,c.shape[a]); o=out[tuple(sl)]
        np.maximum(o,ang,out=o)
    return out
