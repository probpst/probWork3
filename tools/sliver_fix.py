#!/usr/bin/env python3
"""sliver_fix.py - rigidify the sliver cells of the heart grid in the morph frames.

The reference grid has a few dozen sliver cells at the branch inlets (edges of 0.02 - 0.1 mm next to
1.5 mm edges).  The wall motion there has a small kink (about 0.03 mm across the sliver), which is a
100 % strain for such a cell and flips it in some frames (frame 49: LPA inlet).  Fix: the nodes of every
connected group of sliver cells get the group's MEAN displacement (a rigid translation of the group; the
wall changes by at most the kink size, ~0.03 - 0.09 mm, below the frame tolerance of 0.15 mm).
usage: sliver_fix.py FRAMEDIR [--edge 0.1] [--out FRAMEDIR_fixed]    (FRAMEDIR holds P0.npy, cells.npy, disp_NNN.npy)
As a module: groups = sliver_groups(P, cells, edge); D = apply_fix(D, groups)
"""
import sys, os, glob, argparse, numpy as np

EDGES = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]

def jac(Pn, cells):
    X = Pn[cells]; J = np.full(len(cells), np.inf)
    for o, a, b, c in [(0,1,3,4),(1,2,0,5),(2,3,1,6),(3,0,2,7),(4,7,5,0),(5,4,6,1),(6,5,7,2),(7,6,4,3)]:
        u = X[:, a] - X[:, o]; v = X[:, b] - X[:, o]; w = X[:, c] - X[:, o]
        J = np.minimum(J, np.einsum('ij,ij->i', np.cross(u, v), w) / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1) * np.linalg.norm(w, axis=1) + 1e-30))
    return J

def sliver_groups(P, cells, edge=0.05):
    """node groups connected by SHORT edges (< edge mm): the displacement is averaged within a group, which
    removes the relative motion across the thin direction of a sliver cell and leaves the long directions free"""
    ea = np.concatenate([cells[:, a] for a, b in EDGES]); eb = np.concatenate([cells[:, b] for a, b in EDGES])
    l = np.linalg.norm(P[ea] - P[eb], axis=1); sel = l < edge
    ea, eb = ea[sel], eb[sel]
    if len(ea) == 0: return []
    nodes = np.unique(np.concatenate([ea, eb])); idx = {int(n): i for i, n in enumerate(nodes)}
    parent = np.arange(len(nodes))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b in zip(ea, eb):
        ra, rb = find(idx[int(a)]), find(idx[int(b)])
        if ra != rb: parent[ra] = rb
    roots = np.array([find(i) for i in range(len(nodes))])
    groups = {}
    for n, r in zip(nodes, roots): groups.setdefault(int(r), []).append(int(n))
    return [np.array(g) for g in groups.values() if len(g) > 1]

def apply_fix(D, groups):
    D = D.copy()
    for g in groups:
        D[g] = D[g].mean(0)
    return D

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('framedir'); ap.add_argument('--edge', type=float, default=0.05); ap.add_argument('--out')
    a = ap.parse_args()
    out = a.out or a.framedir.rstrip('/') + '_fixed'; os.makedirs(out, exist_ok=True)
    P = np.load(os.path.join(a.framedir, 'P0.npy')); cells = np.load(os.path.join(a.framedir, 'cells.npy'))
    groups = sliver_groups(P, cells, a.edge)
    print(f'short edges (< {a.edge} mm): {len(groups)} node groups, sizes max {max(len(g) for g in groups)} mean {np.mean([len(g) for g in groups]):.1f}, nodes {sum(len(g) for g in groups)}')
    np.save(os.path.join(out, 'P0.npy'), P); np.save(os.path.join(out, 'cells.npy'), cells)
    files = sorted(glob.glob(os.path.join(a.framedir, 'disp_*.npy'))); prev = None; worst = (1, 0); wmid = 1; maxchg = 0
    for f in files:
        D = np.load(f).astype(np.float64); Df = apply_fix(D, groups); maxchg = max(maxchg, np.abs(Df - D).max())
        J = jac(P + Df, cells); neg = int((J <= 0).sum())
        if prev is not None:
            Jm = jac(P + 0.5 * (prev + Df), cells); wmid = min(wmid, Jm.min())
        np.save(os.path.join(out, os.path.basename(f)), Df.astype(np.float32))
        worst = min(worst, (J.min(), -neg)); prev = Df
        print(f'{os.path.basename(f)}: neg {neg} minJ {J.min():.4f}', flush=True)
    print(f'done: min scaled Jacobian over all frames {worst[0]:.4f} (inverted max {-worst[1]}), midpoints min {wmid:.4f}, '
          f'max node change by the fix {maxchg:.4f} mm')

if __name__ == '__main__':
    main()
