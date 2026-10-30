"""assembled.pkl -> heart_grid.grd (ASCII, reference format) / .g (binary) / .bcs / block_list.txt / STL
usage: python3 export.py OUTDIR"""
import pickle,numpy as np,sys,os,json
import mb
from mb import fax
import trimesh
out=sys.argv[1] if len(sys.argv)>1 else 'out'; os.makedirs(out+'/stl/boundaries',exist_ok=True)
import os; A=pickle.load(open(os.environ.get('GRID','assembled.pkl'),'rb')); B=A['blocks']; bcs=A['bcs']
mb.write_bcs(f'{out}/heart_grid.bcs',B,A['ifc'],bcs)
mb.write_grd_binary(f'{out}/heart_grid.g',B)
with open(f'{out}/block_list.txt','w') as f:
    f.write('# block  domain  kind  ni nj nk  name\n')
    for n,b in enumerate(B): f.write(f'{n+1:4d}  {b.dom}  {b.kind:4s}  {b.dims[0]:4d} {b.dims[1]:4d} {b.dims[2]:4d}  {b.name}\n')
def fE(a):
    a=np.asarray(a,float); m=np.char.replace(np.char.lstrip(np.char.mod('%.15E',a),'-'),'.','')
    mant=np.array([x.split('E')[0] for x in m.tolist()]); ex=np.array([int(x.split('E')[1])+1 for x in m.tolist()])
    o=np.char.add(np.char.add(np.where(a<0,'-0.','0.'),mant),np.char.mod('E%+03d',ex)); o[a==0]='0.0000000000000000E+00'
    return np.char.rjust(o,28)
with open(f'{out}/heart_grid.grd','w') as fh:
    fh.write('%12d\n'%len(B))
    for b in B: fh.write('%12d%12d%12d\n'%b.dims)
    for b in B:
        for c in range(3):
            v=fE(b.X[...,c].transpose(2,1,0).ravel()); n=len(v); full=n//4*4
            if full: fh.write('\n'.join(np.char.add(np.char.add(v[0:full:4],v[1:full:4]),np.char.add(v[2:full:4],v[3:full:4])))+'\n')
            if n>full: fh.write(''.join(v[full:])+'\n')
# ---- STL straight from the grid boundary faces (quads -> 2 triangles, outward normals)
groups={}
for bi,S in bcs.items():
    b=B[bi]
    for s in S:
        F=b.face(s['f'])[s['pr'][0]:s['pr'][1]+2,s['qr'][0]:s['qr'][1]+2]
        a,hi=fax(s['f']); idx=[slice(None)]*3; idx[a]=-2 if hi else 1
        G=b.X[tuple(idx)][s['pr'][0]:s['pr'][1]+2,s['qr'][0]:s['qr'][1]+2]
        p00=F[:-1,:-1];p10=F[1:,:-1];p11=F[1:,1:];p01=F[:-1,1:]
        flip=np.einsum('pqi,pqi->pq',np.cross(p10-p00,p01-p00),(G-F)[:-1,:-1])>0
        use1=np.linalg.norm(p11-p00,axis=-1)<=np.linalg.norm(p10-p01,axis=-1)
        T=[]
        for x0,x1,x2,x3,fl,u in zip(p00.reshape(-1,3),p10.reshape(-1,3),p11.reshape(-1,3),p01.reshape(-1,3),flip.ravel(),use1.ravel()):
            for t in ([(x0,x1,x2),(x0,x2,x3)] if u else [(x0,x1,x3),(x1,x2,x3)]): T.append(t[::-1] if fl else t)
        groups.setdefault((b.dom,'wall' if s['type']=='wall' else s['label']),[]).append(np.array(T))
def mesh(tris):
    V=np.vstack(tris).reshape(-1,3); m=trimesh.Trimesh(vertices=V,faces=np.arange(len(V)).reshape(-1,3),process=True); m.merge_vertices(digits_vertex=7); return m
rep={}
for dom,nm in (('L','left_heart'),('R','right_heart')):
    mesh(groups[(dom,'wall')]).export(f'{out}/stl/{nm}_wall.stl')
    for k in [k for k in groups if k[0]==dom and k[1]!='wall']: mesh(groups[k]).export(f'{out}/stl/boundaries/{k[1]}.stl')
    full=mesh([t for k in groups if k[0]==dom for t in groups[k]]); full.export(f'{out}/stl/{nm}_fluid.stl')
    rep[nm]=dict(triangles=int(len(full.faces)),watertight=bool(full.is_watertight),volume_mL=round(float(full.volume)/1000,1))
trimesh.util.concatenate([trimesh.load(f'{out}/stl/left_heart_fluid.stl'),trimesh.load(f'{out}/stl/right_heart_fluid.stl')]).export(f'{out}/stl/heart_full_merged.stl')
json.dump(rep,open(f'{out}/stl/stl_report.json','w'),indent=1); print(rep)
