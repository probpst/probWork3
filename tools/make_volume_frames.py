"""Volume morph of the heart grid at the 97 cycle frames (diffusion morph, validated tool).
Writes /home/claude/frames/disp_NNN.npy  (global-node displacement, float32, N x 3, mm)
and   /home/claude/frames/frames_log.csv (frame, s, inverted cells, min scaled Jacobian, cells<0.05).
Resumable: frames already written are skipped.  Also checks the linear midpoint between
consecutive frames (the solver interpolates in time between frames)."""
import sys,os,csv,time,pickle,numpy as np
sys.path.insert(0,'/home/claude/hg'); sys.path.insert(0,'/home/claude/cyc')
os.chdir('/home/claude/cyc')
from diffmorph import DiffMorph
A=pickle.load(open('/home/claude/hg/assembled_s5.pkl','rb')); bm=dict(np.load('bmesh.npz'))
D=np.load('Dsamples.npy',mmap_mode='r'); NS=D.shape[0]
fr=[int(round(float(r['phase_fraction'])*NS)) for r in csv.DictReader(open('/home/claude/egg_out/cycle_frames.csv'))]
M=DiffMorph(A,bm['pos'],alpha=1.0)
np.save('/home/claude/frames/P0.npy',M.P)
if not os.path.exists('/home/claude/frames/cells.npy'): np.save('/home/claude/frames/cells.npy',M.cells)
out='/home/claude/frames/'; log=out+'frames_log.csv'
done=set()
if os.path.exists(log):
    done={int(r['frame']) for r in csv.DictReader(open(log))}
else:
    open(log,'w').write('frame,s,i_sample,neg,minJ,below005,mid_neg,mid_minJ,sec\n')
prev=None
for n,i in enumerate(fr):
    f=out+f'disp_{n:03d}.npy'
    if n in done and os.path.exists(f):
        prev=np.load(f).astype(np.float64); continue
    t=time.time(); Dn=M.morph(np.asarray(D[i],dtype=np.float64)); J=M.jac(M.P+Dn)
    mneg=-1; mmin=np.nan
    if prev is not None:
        Jm=M.jac(M.P+0.5*(prev+Dn)); mneg=int((Jm<=0).sum()); mmin=float(Jm.min())
    np.save(f,Dn.astype(np.float32))
    with open(log,'a') as fo: fo.write(f'{n},{i/NS:.6f},{i},{int((J<=0).sum())},{J.min():.5f},{int((J<0.05).sum())},{mneg},{mmin:.5f},{time.time()-t:.0f}\n')
    print(f'frame {n:03d} s={i/NS:.4f} neg {(J<=0).sum()} minJ {J.min():.4f} mid_neg {mneg} ({time.time()-t:.0f}s)',flush=True)
    prev=Dn
print('ALL FRAMES DONE')
