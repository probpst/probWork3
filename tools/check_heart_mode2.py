#!/usr/bin/env python3
"""check_heart_mode2.py CASEDIR [--csv cycle_frames.csv --block-list block_list.txt --dt 5e-4]
For every grid written by the Mode-2 (or Mode-0) run of the heart case: min / p1 scaled Jacobian, inverted
cells, and the fluid volume of the left and right domains against the frame table's
left_domain_volume_mL / right_domain_volume_mL interpolated at the grid's time (periodic)."""
import sys, os, glob, csv, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('case'); ap.add_argument('--csv', required=True); ap.add_argument('--block-list', required=True)
ap.add_argument('--dt', type=float, default=5e-4); ap.add_argument('--period', type=float, default=60/70.0)
a = ap.parse_args()
def rec(f):
    n = np.frombuffer(f.read(4), '<i4')[0]; b = f.read(n); f.read(4); return b
def read_grd(fn):
    f = open(fn, 'rb'); nb = np.frombuffer(rec(f), '<i4')[0]; dims = np.frombuffer(rec(f), '<i4').reshape(nb, 3)
    return [np.frombuffer(rec(f), '<f4').reshape(3, nk, nj, ni).transpose(0, 3, 2, 1).astype(np.float64) for ni, nj, nk in dims]
def cell_stats(g):
    X = g.transpose(1, 2, 3, 0)
    c = [X[:-1, :-1, :-1], X[1:, :-1, :-1], X[1:, 1:, :-1], X[:-1, 1:, :-1], X[:-1, :-1, 1:], X[1:, :-1, 1:], X[1:, 1:, 1:], X[:-1, 1:, 1:]]
    J = np.full(c[0].shape[:3], np.inf)
    for o, p, q, r in [(0,1,3,4),(1,2,0,5),(2,3,1,6),(3,0,2,7),(4,7,5,0),(5,4,6,1),(6,5,7,2),(7,6,4,3)]:
        u = c[p] - c[o]; v = c[q] - c[o]; w = c[r] - c[o]
        J = np.minimum(J, np.einsum('...i,...i->...', np.cross(u, v), w) / (np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1) * np.linalg.norm(w, axis=-1) + 1e-30))
    tot = 0.0
    for (p, q, r, s_) in [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]:
        n = np.cross(c[r] - c[p], c[s_] - c[q]) / 2.0; xc = (c[p] + c[q] + c[r] + c[s_]) / 4
        tot = tot + np.einsum('...i,...i->...', xc, n)
    return J.ravel(), float(np.abs(tot.sum()) / 3.0)
dom = [l.split()[1] for l in open(a.block_list) if l.strip() and not l.startswith('#')]
rows = list(csv.DictReader(open(a.csv))); tf = np.array([float(r['time_s']) for r in rows])
VL = np.array([float(r['left_domain_volume_mL']) for r in rows]); VR = np.array([float(r['right_domain_volume_mL']) for r in rows])
def interp(t, y):
    s = t % a.period; return np.interp(s, np.r_[tf, a.period], np.r_[y, y[0]])
files = sorted(glob.glob(os.path.join(a.case, 'Check_grid_movement', '*.grd'))) or sorted(glob.glob(os.path.join(a.case, 'Volume_solutions_grids', '*_????????.grd')))
print(f'{"step":>6} {"t[s]":>8} {"minJ":>8} {"p1 J":>7} {"inverted":>8} {"V_left[mL]":>11} {"table":>8} {"err%":>6} {"V_right":>9} {"table":>8} {"err%":>6}')
for fn in files:
    step = int(os.path.basename(fn).split('_')[-1][:8]); t = step * a.dt
    G = read_grd(fn); Js = []; vl = vr = 0.0
    for g, d in zip(G, dom):
        J, vol = cell_stats(g); Js.append(J)
        if d == 'L': vl += vol
        else: vr += vol
    J = np.concatenate(Js); vl *= 1e6; vr *= 1e6; el = interp(t, VL); er = interp(t, VR)
    print(f'{step:6d} {t:8.4f} {J.min():8.4f} {np.percentile(J, 1):7.3f} {(J <= 0).sum():8d} {vl:11.2f} {el:8.2f} {100*(vl-el)/el:6.2f} {vr:9.2f} {er:8.2f} {100*(vr-er)/er:6.2f}')
