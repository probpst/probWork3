import numpy as np,sys,pickle,os,json; sys.path.insert(0,'/home/claude/cyc')
from diffmorph import DiffMorph
from bmotion3 import BMotion3
from crotch import crotch_points
G=os.environ.get('GRID','/home/claude/hg/assembled_s4.pkl'); A=pickle.load(open(G,'rb')); bm=dict(np.load('bmesh.npz'))
B=BMotion3(bm,A['inout'],crotch_pts=crotch_points(A)); cal=pickle.load(open('calib_final.pkl','rb')) if os.path.exists('calib_final.pkl') else B.calibrate(False)
B.mo.a_lv=cal['a_lv']; B.mo.a_rv=cal['a_rv']; pickle.dump(cal,open('calib_final.pkl','wb'))
M=DiffMorph(A,B.X,alpha=1.0)
res=json.load(open('fullcycle.json')) if os.path.exists('fullcycle.json') else {}
for s in [float(x) for x in sys.argv[1:]]:
    D=M.morph(B.disp(s)); J=M.jac(M.P+D); res[f'{s:.3f}']=dict(neg=int((J<=0).sum()),minJ=float(J.min()),below005=int((J<0.05).sum()))
    print(f'  s={s:.3f} neg {(J<=0).sum()} minJ {J.min():.4f} <0.05 {(J<0.05).sum()}',flush=True)
    json.dump(res,open('fullcycle.json','w'),indent=1)
