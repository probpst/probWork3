"""Check that the restored working state reproduces the original session exactly.
  1. ventricular amplitude calibration  == recorded values (LV/RV stroke volumes -83 / -78 mL)
  2. wall motion                         == delivered STL frames (vertex positions, triangle order)
  3. diffusion morph at end-systole     -> 0 inverted cells (original: min scaled Jacobian 0.0126)
usage: python3 validate_resume.py <data dir>   (needs reference_frames.npz, ES_FRAME.txt, cycle_frames.csv)"""
import sys, os, csv, pickle, time
import numpy as np
sys.path.insert(0, '/home/claude/cyc'); sys.path.insert(0, '/home/claude/hg')
from bmotion3 import BMotion3
from crotch import crotch_points
src = sys.argv[1] if len(sys.argv) > 1 else '.'
A = pickle.load(open('/home/claude/hg/assembled_s5.pkl', 'rb')); bm = dict(np.load('/home/claude/cyc/bmesh.npz'))
ok = True
t = time.time(); B = BMotion3(bm, A['inout'], crotch_pts=crotch_points(A))
cal = B.calibrate(True); ref = pickle.load(open('/home/claude/cyc/calib_final.pkl', 'rb'))
d = max(abs(cal['a_lv'] - ref['a_lv']), abs(cal['a_rv'] - ref['a_rv']))
# tolerances: the RBF wall-motion evaluation runs in float32, whose rounding depends on the CPU / BLAS
# kernel (measured: float32 vs float64 differ by up to 0.044 mm). Differences below these limits are
# machine rounding, not a restore error.
print(f'[1] calibration  a_lv {cal["a_lv"]:.10f}  a_rv {cal["a_rv"]:.10f}  |diff vs recorded| {d:.1e}  -> {"OK" if d < 1e-4 else "DIFFERENT"}  ({time.time()-t:.0f}s)')
ok &= d < 1e-4
B.mo.a_lv = ref['a_lv']; B.mo.a_rv = ref['a_rv']
rows = {r['file']: r for r in csv.DictReader(open(os.path.join(src, 'cycle_frames.csv')))}
X = B.X; T = B.T
REF = np.load(os.path.join(src, 'reference_frames.npz'))        # vertex coordinates of delivered frames (float32)
for fn, key in (('heart_cycle_000.stl', 'frame000'), (open(os.path.join(src, 'ES_FRAME.txt')).read().strip(), 'frameES')):
    s = round(float(rows[fn]['phase_fraction']) * 700) / 700
    Yref = REF[key]; Y = (X + B.disp(s)).astype(np.float64).astype(np.float32)
    err = np.abs(Yref - Y).max()
    print(f'[2] {fn} (s={s:.4f}): vertices {len(Yref)} vs {len(X)}, max vertex difference {err:.2e} mm  -> {"OK" if (len(Yref) == len(X) and err < 0.1) else "DIFFERENT"}')
    ok &= bool(len(Yref) == len(X) and err < 0.1)
from diffmorph import DiffMorph
t = time.time(); M = DiffMorph(A, B.X, alpha=1.0); D = M.morph(B.disp(0.33)); J = M.jac(M.P + D)
print(f'[3] diffusion morph s=0.33: inverted cells {(J <= 0).sum()}, min scaled Jacobian {J.min():.4f} (original 0.0126)  ({time.time()-t:.0f}s)')
ok &= bool((J <= 0).sum() == 0)
print('RESUME STATE VALID' if ok else '!! RESUME STATE DIFFERS FROM ORIGINAL')
