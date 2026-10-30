import numpy as np,sys,pickle,time; sys.path.insert(0,'/home/claude/cyc')
from bmotion3 import BMotion3
from crotch import crotch_points
A=pickle.load(open('/home/claude/hg/assembled_s5.pkl','rb')); bm=dict(np.load('bmesh.npz'))
B=BMotion3(bm,A['inout'],crotch_pts=crotch_points(A)); cal=pickle.load(open('calib_final.pkl','rb')); B.mo.a_lv=cal['a_lv']; B.mo.a_rv=cal['a_rv']
NS=int(sys.argv[1]); i0,i1=int(sys.argv[2]),int(sys.argv[3])
Dm=np.lib.format.open_memmap('Dsamples.npy',mode='r+' if i0>0 else 'w+',dtype=np.float32,shape=(NS,len(B.X),3))
t=time.time()
for i in range(i0,min(i1,NS)):
    Dm[i]=B.disp(i/NS)
Dm.flush(); print('samples',i0,'-',min(i1,NS),'%.0fs'%(time.time()-t))
