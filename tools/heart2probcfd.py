#!/usr/bin/env python3
"""heart2probcfd.py - convert the heart grid (PLOT3D .g binary or .grd ASCII, mm, CFDShip-style
.bcs with types 20/9/12/92) into a probCFD case:
    <out>/<case>.grd                  ASCII multi-block PLOT3D in METRES
    <out>/<case>.bcs                  probCFD boundary conditions (type 92 omitted: auto-detected)
    <out>/Input_polyMesh_probCFD.csv  per-block MPI splitting for N ranks (+ Object_ID row)
    <out>/<case>_faces.txt            summary of the physical faces (label, block, face, nodes)
BC mapping (all patches of the heart are either no-slip walls or pressure boundaries):
    20 wall          -> u, v, w Dirichlet 0, p zero gradient
    9 inlet / 12 outlet -> u, v, w zero gradient, p Dirichlet [Pa]  (pressure boundary, flow in or out)
Pressure values: --p default=<Pa> and --p <label-substring>=<Pa> (later rules override earlier ones).
usage: heart2probcfd.py grid.g|grid.grd heart.bcs outdir case [--nranks N] [--p rule ...] [--object-id k]
"""
import sys, os, struct, argparse, numpy as np

def read_g(fn):
    with open(fn, 'rb') as f:
        def rec():
            n = struct.unpack('<i', f.read(4))[0]; d = f.read(n); f.read(4); return d
        nb = struct.unpack('<i', rec())[0]
        dims = np.frombuffer(rec(), '<i4').reshape(nb, 3)
        blocks = []
        for ni, nj, nk in dims:
            a = np.frombuffer(rec(), '<f8').reshape(3, nk, nj, ni).transpose(0, 3, 2, 1)
            blocks.append(np.ascontiguousarray(a))
    return blocks

def read_grd(fn):
    toks = iter(open(fn).read().split())
    nb = int(next(toks)); dims = [(int(next(toks)), int(next(toks)), int(next(toks))) for _ in range(nb)]
    blocks = []
    for ni, nj, nk in dims:
        n = ni * nj * nk
        a = np.array([float(next(toks)) for _ in range(3 * n)]).reshape(3, nk, nj, ni).transpose(0, 3, 2, 1)
        blocks.append(np.ascontiguousarray(a))
    return blocks

def read_cfdship_bcs(fn):
    """returns per block a list of (type, dir, ir, jr, kr, label); type-92 records skipped"""
    raw = open(fn, encoding='utf-8', errors='replace').read().replace('\r', '').split('\n')
    lines = [(l.split('!')[0].split(), l.split('!', 1)[1].strip() if '!' in l else '') for l in raw]
    lines = [(t, c) for t, c in lines if t]
    p = 0; blocks = []
    while p < len(lines):
        nsurf = int(lines[p][0][0]); p += 1; surfs = []
        for _ in range(nsurf):
            typ = int(lines[p][0][0]); cmt = lines[p][1]
            d = int(lines[p + 1][0][0])
            ir = tuple(map(int, lines[p + 2][0][:2])); jr = tuple(map(int, lines[p + 3][0][:2])); kr = tuple(map(int, lines[p + 4][0][:2]))
            p += 6
            if typ == 92:
                p += 5; continue
            label = cmt.split('=')[0].replace('***', '').strip() if cmt else f'type{typ}'
            surfs.append((typ, d, ir, jr, kr, label))
        blocks.append(surfs)
    return blocks

def split_plan(dims, target):
    """split counts (si, sj, sk) so that every part has about `target` cells (>= 2 cells per direction)"""
    nc = [d - 1 for d in dims]; sp = [1, 1, 1]
    while np.prod(nc) / np.prod(sp) > 1.5 * target:
        # split the direction with the most cells per part
        a = int(np.argmax([nc[q] / sp[q] for q in range(3)]))
        if nc[a] // (sp[a] + 1) < 2: break
        sp[a] += 1
    return sp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('grid'); ap.add_argument('bcs'); ap.add_argument('outdir'); ap.add_argument('case')
    ap.add_argument('--nranks', type=int, default=20)
    ap.add_argument('--p', action='append', default=[], help='pressure rule label-substring=Pa (or default=Pa)')
    ap.add_argument('--object-id', type=int, default=0, help='Object_ID row value for every block (0 = static background)')
    ap.add_argument('--scale', type=float, default=1e-3, help='grid units -> metres (mm: 1e-3)')
    ap.add_argument('--domains-per-rank', type=float, default=3.0, help='target sub-domains per rank (round-robin balance)')
    a = ap.parse_args()
    G = read_g(a.grid) if a.grid.endswith('.g') else read_grd(a.grid)
    B = read_cfdship_bcs(a.bcs)
    assert len(G) == len(B), (len(G), len(B))
    os.makedirs(a.outdir, exist_ok=True)
    nb = len(G)
    # ---- grid (metres)
    with open(os.path.join(a.outdir, a.case + '.grd'), 'w') as f:
        f.write(f'{nb}\n')
        for g in G: f.write(f'{g.shape[1]} {g.shape[2]} {g.shape[3]}\n')
        for g in G:
            for c in range(3):
                v = (g[c] * a.scale).transpose(2, 1, 0).ravel()
                n4 = (len(v) // 4) * 4
                if n4: np.savetxt(f, v[:n4].reshape(-1, 4), fmt='%.16g')
                if n4 < len(v): np.savetxt(f, v[n4:][None, :], fmt='%.16g')
    # ---- pressure rules
    rules = []
    for r in a.p:
        k, v = r.split('='); rules.append((k.strip(), float(v)))
    def pval(label):
        p = 0.0
        for k, v in rules:
            if k == 'default' or k in label: p = v
        return p
    # ---- bcs
    L = [f'! probCFD boundary conditions for case {a.case}: converted from {os.path.basename(a.bcs)}',
         '! walls (CFDShip type 20): u, v, w Dirichlet 0, p zero gradient;',
         '! inlets (9) / outlets (12): pressure boundaries, u, v, w zero gradient, p Dirichlet [Pa].',
         '! type-92 block-to-block faces omitted (probCFD detects them itself).',
         f'{nb:8d}                ! number of grid blocks']
    S = ['# label  block  face  i-range  j-range  k-range  nodes  p[Pa]']
    nwall = nio = 0
    for b, surfs in enumerate(B):
        L.append(f'{len(surfs):8d}                ! block {b + 1}: {len(surfs)} boundary surfaces')
        for typ, d, ir, jr, kr, label in surfs:
            L.append(f'{d:8d}                ! {label} (type {typ})')
            for nm, r in zip('ijk', (ir, jr, kr)):
                L.append(f'{min(r):8d} {max(r):4d}           ! {nm}-range (nodes)')
            if typ == 20:
                L += [f'{1:8d}   0.0          ! u: Dirichlet 0 (no-slip wall)', f'{1:8d}   0.0          ! v', f'{1:8d}   0.0          ! w',
                      f'{0:8d}   0.0          ! p: zero gradient']
                nwall += 1; pp = ''
            else:
                pp = pval(label)
                L += [f'{0:8d}   0.0          ! u: zero gradient (pressure boundary)', f'{0:8d}   0.0          ! v', f'{0:8d}   0.0          ! w',
                      f'{1:8d}   {pp:.6g}          ! p: Dirichlet [Pa]']
                nio += 1
            nn = (abs(ir[1] - ir[0]) + 1) * (abs(jr[1] - jr[0]) + 1) * (abs(kr[1] - kr[0]) + 1)
            fname = {1: 'i-min', -1: 'i-max', 2: 'j-min', -2: 'j-max', 3: 'k-min', -3: 'k-max'}[d]
            S.append(f'{label:32s} {b + 1:4d}  {fname:5s}  {min(ir):3d}-{max(ir):<3d}  {min(jr):3d}-{max(jr):<3d}  {min(kr):3d}-{max(kr):<3d}  {nn:6d}  {pp}')
    open(os.path.join(a.outdir, a.case + '.bcs'), 'w').write('\n'.join(L) + '\n')
    open(os.path.join(a.outdir, a.case + '_faces.txt'), 'w').write('\n'.join(S) + '\n')
    # ---- csv (MPI splitting)
    cells = np.array([np.prod(np.array(g.shape[1:]) - 1) for g in G]); tot = int(cells.sum())
    target = tot / (a.nranks * a.domains_per_rank)
    sp = [split_plan(g.shape[1:], target) for g in G]
    ndom = sum(int(np.prod(s)) for s in sp)
    load = np.zeros(a.nranks); d = 0
    for b, s in enumerate(sp):
        for _ in range(int(np.prod(s))):
            load[d % a.nranks] += cells[b] / np.prod(s); d += 1
    rows = [['probCFD polyMesh Input'] + [f'block-{b + 1}' for b in range(nb)],
            ['i_split'] + [str(s[0]) for s in sp], ['j_split'] + [str(s[1]) for s in sp], ['k_split'] + [str(s[2]) for s in sp],
            ['n_of_BLK'] + [str(int(np.prod(s))) for s in sp], ['Object_ID'] + [str(a.object_id)] * nb]
    open(os.path.join(a.outdir, 'Input_polyMesh_probCFD.csv'), 'w').write('\n'.join(','.join(r) for r in rows) + '\n')
    print(f'blocks {nb}, cells {tot:,}, walls {nwall}, pressure faces {nio}')
    print(f'sub-domains {ndom} for {a.nranks} ranks (round-robin): cells per rank min {load.min():,.0f} max {load.max():,.0f} '
          f'mean {load.mean():,.0f} (imbalance {load.max() / load.mean():.2f})')

if __name__ == '__main__':
    main()
