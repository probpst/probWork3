import numpy as np,json,sys,csv; sys.path.insert(0,'/home/claude/cyc'); sys.path.insert(0,'/home/claude/hg')
from scipy.spatial import cKDTree
from deform_v1 import LV_V,RV_V,LA_V,RA_V,AO_D,PA_D,phase_name
bm=np.load('bmesh.npz'); X=bm['pos']; T=bm['tris']; dom=bm['dom']; h=bm['h']
D=np.load('Dsamples.npy',mmap_mode='r'); NS=D.shape[0]; fr=json.load(open('frame_sets.json'))['0.2_0.15_10']
Tc=60/70.0; out='/home/claude/out4/cycle_stl/'
dt=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')])
n0=np.cross(X[T[:,1]]-X[T[:,0]],X[T[:,2]]-X[T[:,0]]); a0=np.linalg.norm(n0,axis=1)
Lt=dom[T[:,0]]=='L'; iL=np.where(dom=='L')[0]; iR=np.where(dom=='R')[0]
def vol(Y,sel):
    t=T[sel]; return abs(np.einsum('ij,ij->i',np.cross(Y[t[:,0]],Y[t[:,1]]),Y[t[:,2]]).sum()/6)/1000
rows=[]; prev=None; worst=dict(flips=0,area=9,lr=9)
V0L,V0R=vol(X,Lt),vol(X,~Lt)
for n,i in enumerate(fr):
    Y=(X+D[i]).astype(np.float64); s=i/NS
    n1=np.cross(Y[T[:,1]]-Y[T[:,0]],Y[T[:,2]]-Y[T[:,0]]); nn=np.linalg.norm(n1,axis=1)
    fl=int(((n0*n1).sum(1)<=0).sum()); ar=float((nn/a0).min()); lr=float(cKDTree(Y[iR]).query(Y[iL])[0].min())
    worst=dict(flips=max(worst['flips'],fl),area=min(worst['area'],ar),lr=min(worst['lr'],lr))
    rec=np.zeros(len(T),dt); rec['n']=n1/np.maximum(nn,1e-30)[:,None]; rec['v']=Y[T].astype(np.float32)
    fn=f'heart_cycle_{n:03d}.stl'
    with open(out+fn,'wb') as f:
        f.write(f'heart cycle frame {n} t={s*Tc:.5f}s'.ljust(80).encode()[:80]); f.write(np.uint32(len(T)).tobytes()); f.write(rec.tobytes())
    step=0.0 if prev is None else float(np.linalg.norm(D[i]-D[prev],axis=1).max())
    rows.append(dict(frame=n,file=fn,time_s=round(s*Tc,6),phase_fraction=round(s,6),phase=phase_name(s),
        LV_target_mL=round(float(LV_V(s)),2),RV_target_mL=round(float(RV_V(s)),2),
        left_domain_volume_mL=round(vol(Y,Lt),2),right_domain_volume_mL=round(vol(Y,~Lt),2),
        max_wall_displacement_from_reference_mm=round(float(np.linalg.norm(D[i],axis=1).max()),3),
        max_wall_step_from_previous_frame_mm=round(step,3)))
    prev=i
with open('/home/claude/out4/cycle_frames.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print('frames',len(rows),'| surface checks over all frames: flips',worst['flips'],'min area ratio %.3f'%worst['area'],'min L-R gap %.2f mm'%worst['lr'])
print('left-domain volume range %.1f-%.1f mL, right %.1f-%.1f mL'%(min(r['left_domain_volume_mL'] for r in rows),max(r['left_domain_volume_mL'] for r in rows),min(r['right_domain_volume_mL'] for r in rows),max(r['right_domain_volume_mL'] for r in rows)))
