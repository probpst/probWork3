set -e
G=$1; T=$2
GRID=/home/claude/hg/$G.pkl OUT=/home/claude/cyc/bmesh_$T.npz MAP=/home/claude/cyc/map_$T.txt python3 bmesh.py | head -1
cp bmesh_$T.npz bmesh.npz; cp map_$T.txt vertex_to_grid_nodes.txt
python3 - <<PY
import numpy as np,sys,pickle; sys.path.insert(0,'/home/claude/cyc')
from bmotion3 import BMotion3
bm=dict(np.load('bmesh.npz')); A=pickle.load(open('/home/claude/hg/$G.pkl','rb'))
from crotch import crotch_points; B=BMotion3(bm,A['inout'],crotch_pts=crotch_points(A)); cal=B.calibrate(False); pickle.dump(cal,open('calib_r.pkl','wb')); print('calib',cal)
PY
