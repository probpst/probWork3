#!/usr/bin/env python3
"""make_valve_zones.py - valve zone files (valves.F90) for the heart grid: the cells of the valve
tube blocks whose centre lies within +-half_thickness of the annulus plane.

Valve planes from heart/scripts/grid/params.py (mm): mitral MITRAL_PTS[MV_ANNULUS_IDX], tricuspid
TRICUSPID_PTS[TV_ANNULUS_IDX], aortic AORTA_PTS[AO_ANNULUS_IDX], pulmonary PA_PTS[PV_ANNULUS_IDX];
the plane normal = the centre-line tangent there.  Zone blocks = the tube core block and its
boundary-layer blocks (block names 'MV', 'MV.*', 'TV', 'TV.*', 'AO', 'AO.*', 'PA', 'PA.*').
Weights taper linearly over the outer 'taper' of the half thickness.
usage: make_valve_zones.py --pkl assembled_s5.pkl --params-dir heart/scripts/grid OUTDIR [--half 1.5 --taper 0.5]
       (or --grid heart_grid.g --block-list block_list.txt instead of --pkl)
"""
import sys, os, argparse, struct, numpy as np

def read_g(fn):
    with open(fn, 'rb') as f:
        def rec():
            n = struct.unpack('<i', f.read(4))[0]; d = f.read(n); f.read(4); return d
        nb = struct.unpack('<i', rec())[0]; dims = np.frombuffer(rec(), '<i4').reshape(nb, 3)
        return [np.ascontiguousarray(np.frombuffer(rec(), '<f8').reshape(3, nk, nj, ni).transpose(0, 3, 2, 1)) for ni, nj, nk in dims]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir'); ap.add_argument('--pkl'); ap.add_argument('--grid'); ap.add_argument('--block-list')
    ap.add_argument('--params-dir', required=True); ap.add_argument('--half', type=float, default=1.5, help='half thickness [mm]')
    ap.add_argument('--taper', type=float, default=0.5, help='linear taper width [mm] inside the half thickness')
    ap.add_argument('--rfac', type=float, default=1.8, help='in-plane radius limit = rfac x annulus radius (the long curved '
                    'aorta / pulmonary tubes re-cross the annulus plane far away)')
    a = ap.parse_args()
    sys.path.insert(0, a.params_dir); import params as P
    if a.pkl:
        import pickle; sys.path.insert(0, os.path.dirname(os.path.abspath(a.pkl)))
        A = pickle.load(open(a.pkl, 'rb')); X = [b.X for b in A['blocks']]; names = [b.name for b in A['blocks']]
    else:
        G = read_g(a.grid); X = [np.moveaxis(g, 0, -1) for g in G]
        names = [l.split()[6] for l in open(a.block_list) if l.strip() and not l.startswith('#')]
    def tangent(pts, i):
        p = np.array(pts, float); t = p[min(i + 1, len(p) - 1)] - p[max(i - 1, 0)]; return t / np.linalg.norm(t)
    valves = {
        'mitral':    (np.array(P.MITRAL_PTS[P.MV_ANNULUS_IDX], float), tangent(P.MITRAL_PTS, P.MV_ANNULUS_IDX), 'MV', P.MITRAL_R[P.MV_ANNULUS_IDX]),
        'tricuspid': (np.array(P.TRICUSPID_PTS[P.TV_ANNULUS_IDX], float), tangent(P.TRICUSPID_PTS, P.TV_ANNULUS_IDX), 'TV', P.TRICUSPID_R[P.TV_ANNULUS_IDX]),
        'aortic':    (np.array(P.AORTA_PTS[P.AO_ANNULUS_IDX], float), tangent(P.AORTA_PTS, P.AO_ANNULUS_IDX), 'AO', P.AORTA_R[P.AO_ANNULUS_IDX]),
        'pulmonary': (np.array(P.PA_PTS[P.PV_ANNULUS_IDX], float), tangent(P.PA_PTS, P.PV_ANNULUS_IDX), 'PA', P.PA_R[P.PV_ANNULUS_IDX]),
    }
    os.makedirs(a.outdir, exist_ok=True)
    for vn, (Pc, n, pre, R) in valves.items():
        rows = []; area_cells = 0
        for b, (x, nm) in enumerate(zip(X, names)):
            if not (nm == pre or nm.startswith(pre + '.')): continue
            c = 0.125 * (x[:-1, :-1, :-1] + x[1:, :-1, :-1] + x[:-1, 1:, :-1] + x[:-1, :-1, 1:] +
                         x[1:, 1:, :-1] + x[1:, :-1, 1:] + x[:-1, 1:, 1:] + x[1:, 1:, 1:])
            dn = (c - Pc) @ n
            rad = np.linalg.norm((c - Pc) - dn[..., None] * n, axis=-1)
            w = np.clip((a.half - np.abs(dn)) / max(a.taper, 1e-9), 0, 1) * (rad <= a.rfac * R)
            idx = np.argwhere(w > 0)
            for i, j, k in idx:
                rows.append(f'{b + 1} {i + 1} {j + 1} {k + 1} {w[i, j, k]:.4f}')
            area_cells += len(idx)
        with open(os.path.join(a.outdir, f'valve_{vn}.dat'), 'w') as f:
            f.write(f'# {vn} valve zone: plane point {Pc.tolist()} mm, normal {np.round(n, 4).tolist()}, half thickness {a.half} mm, '
                    f'blocks {pre}*; columns block i j k weight (1-based cells)\n')
            f.write('\n'.join(rows) + '\n')
        print(f'{vn:10s}: {len(rows)} cells in blocks {pre}*  (plane {Pc}, n {np.round(n,3)}, R {R} mm, r <= {a.rfac*R:.1f} mm)')

if __name__ == '__main__':
    main()
