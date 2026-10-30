"""Diffusion-morph check of the grid at given cycle phases, with the blocks that hold inverted cells.
usage: GRID=/home/claude/hg/assembled_s5.pkl python3 check_morph.py 0.33 0.45 ...
(calibration read from calib_final.pkl; rebuild bmesh.npz with outer1.sh first if GRID changed)"""
import os, sys, pickle, collections, time
import numpy as np
sys.path.insert(0, '/home/claude/cyc'); sys.path.insert(0, '/home/claude/hg')
from diffmorph import DiffMorph
from bmotion3 import BMotion3
from crotch import crotch_points
A = pickle.load(open(os.environ.get('GRID', '/home/claude/hg/assembled_s5.pkl'), 'rb')); bm = dict(np.load('/home/claude/cyc/bmesh.npz'))
B = BMotion3(bm, A['inout'], crotch_pts=crotch_points(A)); c = pickle.load(open('/home/claude/cyc/calib_final.pkl', 'rb'))
B.mo.a_lv = c['a_lv']; B.mo.a_rv = c['a_rv']
M = DiffMorph(A, B.X, alpha=1.0); names = np.concatenate([[b.name] * int(np.prod(np.array(b.dims) - 1)) for b in M.B])
for s in [float(x) for x in sys.argv[1:]]:
    t = time.time(); D = M.morph(B.disp(s)); J = M.jac(M.P + D); neg = J <= 0
    print(f's={s:.3f} inverted {neg.sum()}  min scaled Jacobian {J.min():.4f}  cells<0.05 {(J < 0.05).sum()}  '
          f'{collections.Counter(names[neg]).most_common(5)}  ({time.time() - t:.0f}s)', flush=True)
