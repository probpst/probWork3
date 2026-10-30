import numpy as np,sys,pickle,time,os; sys.path.insert(0,'/home/claude/cyc')
from diffmorph import DiffMorph
from bmotion3 import BMotion3
bm=dict(np.load('bmesh.npz')); A=pickle.load(open(sys.argv[1],'rb'))
from crotch import crotch_points; B=BMotion3(bm,A['inout'],crotch_pts=crotch_points(A)); cal=pickle.load(open('calib_r.pkl','rb')); B.mo.a_lv=cal['a_lv']; B.mo.a_rv=cal['a_rv']
M=DiffMorph(A,B.X,alpha=1.0)
out={}
for s in [float(x) for x in sys.argv[3:]]:
    t=time.time(); D=M.morph(B.disp(s)); J=M.jac(M.P+D); out[s]=D.astype(np.float32)
    np.save(sys.argv[2].replace('.npz','')+f'_s{int(round(s*1000)):04d}.npy',D.astype(np.float32))
    print(f'  s={s:.2f} neg {(J<=0).sum()} minJ {J.min():.3f} ({time.time()-t:.0f}s)',flush=True)
np.savez(sys.argv[2],P=M.P,**{f's{int(round(s*1000)):04d}':v for s,v in out.items()})
