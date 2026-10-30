"""Rebuild the working state of the v4 heart grid from the delivered files.

inputs : heart_grid.g (binary PLOT3D), heart_grid.bcs, block_list.txt
output : /home/claude/hg/assembled_s5.pkl  (+ assembled_base.pkl, assembled.pkl)
         dict(blocks=[mb.Block...], ifc=[...], bcs={b:[...]}, inout=[...])

The pickle is exactly what the grid / cycle scripts expect (same structure the original
session produced with assemble.py + optimize_motion.py). Self-checks:
  * every type-92 entry of the .bcs file is re-found by the point-matching search
  * the physical-BC rectangles re-derived by mb.bc_surfaces equal the .bcs entries
  * every inlet/outlet rim is an exact circle (fit residual reported)
usage: python3 restore_state.py <dir with heart_grid.g, heart_grid.bcs, block_list.txt>
"""
import sys, os, struct, pickle, re
import numpy as np
sys.path.insert(0, '/home/claude/hg')
import mb, params as P
from mb import fax, FACES

src = sys.argv[1] if len(sys.argv) > 1 else '.'
HG = '/home/claude/hg'

# ---------------------------------------------------------------- grid (binary PLOT3D)
raw = open(os.path.join(src, 'heart_grid.g'), 'rb').read(); p = 0
def rec():
    global p
    n = struct.unpack('<i', raw[p:p + 4])[0]; d = raw[p + 4:p + 4 + n]
    assert struct.unpack('<i', raw[p + 4 + n:p + 8 + n])[0] == n; p += 8 + n; return d
nb = struct.unpack('<i', rec())[0]
dims = np.frombuffer(rec(), '<i4').reshape(nb, 3)
meta = [l.split() for l in open(os.path.join(src, 'block_list.txt')) if not l.startswith('#') and l.strip()]
assert len(meta) == nb, (len(meta), nb)
blocks = []
for b, (ni, nj, nk) in enumerate(dims):
    a = np.frombuffer(rec(), '<f8'); n = ni * nj * nk
    X = np.stack([a[c * n:(c + 1) * n].reshape(nk, nj, ni).transpose(2, 1, 0) for c in range(3)], -1)
    num, dom, kind, mi, mj, mk, name = meta[b]
    assert int(num) == b + 1 and (int(mi), int(mj), int(mk)) == (ni, nj, nk), meta[b]
    blocks.append(mb.Block(name, X, dom, kind, name.split('.')[0]))
assert p == len(raw)
print(f'grid: {nb} blocks, {sum(int(np.prod(d - 1)) for d in dims):,} cells')

# ---------------------------------------------------------------- .bcs (with labels from comments)
L = []
for line in open(os.path.join(src, 'heart_grid.bcs')).read().replace('\r', '').split('\n'):
    if not line.strip(): continue
    num, _, cmt = line.partition('!')
    L.append(([int(x) for x in num.split()], cmt))
k = 0; file_bc = []; file_92 = []
for b in range(nb):
    ns = L[k][0][0]; k += 1
    for _ in range(ns):
        typ = L[k][0][0]; cmt = L[k][1]; d = L[k + 1][0][0]
        rng = [L[k + 2][0], L[k + 3][0], L[k + 4][0]]; k += 6
        a = abs(d) - 1; f = ('ijk'[a]) + ('min' if d > 0 else 'max'); ip = [x for x in range(3) if x != a]
        pr = (rng[ip[0]][0] - 1, rng[ip[0]][1] - 2); qr = (rng[ip[1]][0] - 1, rng[ip[1]][1] - 2)
        if typ == 92:
            donor = L[k][0][0]; k += 5; file_92.append((b, f, pr, qr, donor - 1))
        else:
            m = re.search(r'\*\*\*\s*(\S+)\s*=\s*type', cmt); label = m.group(1)
            t = {20: 'wall', 9: 'inlet', 12: 'outlet'}[typ]
            file_bc.append((b, f, pr, qr, t, label))
            if t != 'wall': blocks[b].tag(f, t, label, pr, qr)
assert k == len(L)
print(f'.bcs: {len(file_92)} block-to-block entries, {len(file_bc)} physical-BC entries')

# ---------------------------------------------------------------- connectivity + BC surfaces
ifc, part = mb.connectivity(blocks)
bcs = mb.bc_surfaces(blocks, part)
got92 = sorted((I['b'], I['f'], tuple(I['pr']), tuple(I['qr']), I['db']) for I in ifc)
exp92 = sorted((b, f, tuple(pr), tuple(qr), d) for b, f, pr, qr, d in file_92)
assert got92 == exp92, 'block-to-block entries differ from the .bcs file'
assert max(I['err'] for I in ifc) < 1e-9
gotbc = sorted((b, s['f'], tuple(s['pr']), tuple(s['qr']), s['type'], s['label']) for b, S in bcs.items() for s in S)
assert gotbc == sorted(file_bc), 'physical-BC rectangles differ from the .bcs file'
print('check: interfaces and BC rectangles identical to the .bcs file; max point mismatch %.1e mm' % max(I['err'] for I in ifc))

# ---------------------------------------------------------------- inlet / outlet discs
def face_nodes(b, s):
    B = blocks[b]; a, hi = fax(s['f']); idx = [slice(None)] * 3; idx[a] = -1 if hi else 0
    F = B.X[tuple(idx)][s['pr'][0]:s['pr'][1] + 2, s['qr'][0]:s['qr'][1] + 2]
    idx[a] = -2 if hi else 1
    G = B.X[tuple(idx)][s['pr'][0]:s['pr'][1] + 2, s['qr'][0]:s['qr'][1] + 2]
    return F, G
wallkeys = set()
for b, S in bcs.items():
    for s in S:
        if s['type'] == 'wall': wallkeys |= set(map(tuple, mb.nkey(face_nodes(b, s)[0]).reshape(-1, 3)))
caps = {}
for b, S in bcs.items():
    for s in S:
        if s['type'] == 'wall': continue
        F, G = face_nodes(b, s); c = caps.setdefault(s['label'], dict(F=[], G=[], type=s['type'], dom=blocks[b].dom, core=None))
        c['F'].append(F.reshape(-1, 3)); c['G'].append((G - F).reshape(-1, 3))
        if blocks[b].kind == 'core': c['core'] = blocks[b].name
inout = {}
for label, c in caps.items():
    F = np.vstack(c['F']); Pm = F.mean(0); U, S_, Vt = np.linalg.svd(F - Pm); n = Vt[2]
    if (np.vstack(c['G']) @ n).mean() > 0: n = -n                      # outward (away from the fluid)
    rim = np.array([x for x in F if tuple(mb.nkey(x)) in wallkeys]); rim = np.unique(rim.round(9), axis=0)
    e1 = Vt[0]; e2 = Vt[1]; u = (rim - Pm) @ e1; v = (rim - Pm) @ e2
    Am = np.c_[2 * u, 2 * v, np.ones(len(u))]; sol = np.linalg.lstsq(Am, u * u + v * v, rcond=None)[0]
    C = Pm + sol[0] * e1 + sol[1] * e2; d = np.linalg.norm(rim - C, axis=1)
    inout[label] = dict(dom=c['dom'], P=C, n=n, r=float(d.mean()), type=c['type'], label=label, block=c['core'],
                        rim_err=float(np.abs(d - d.mean()).max()), planarity=float(np.abs((F - C) @ n).max()))
order = (['outlet_descending_aorta'] + ['outlet_' + k for k in P.CORONARIES] + ['outlet_' + k for k in P.BRANCHES]
         + ['inlet_' + k for k in P.PULM_VEINS] + ['outlet_right_pulmonary_artery', 'outlet_left_pulmonary_artery']
         + ['inlet_' + k for k in P.CAVAE])
assert sorted(order) == sorted(inout), (sorted(order), sorted(inout))
inout = [inout[k] for k in order]
print('inlet/outlet discs: %d, worst rim-circle error %.1e mm, worst planarity %.1e mm'
      % (len(inout), max(x['rim_err'] for x in inout), max(x['planarity'] for x in inout)))

A = dict(blocks=blocks, ifc=ifc, bcs=bcs, inout=inout)
for fn in ('assembled_s5.pkl', 'assembled_base.pkl', 'assembled.pkl'):
    pickle.dump(A, open(os.path.join(HG, fn), 'wb'))
nn = os.path.join(HG, 'node_numbering.npz')
if os.path.exists(nn): os.remove(nn)
print('written:', ', '.join(os.path.join(HG, f) for f in ('assembled_s5.pkl', 'assembled_base.pkl', 'assembled.pkl')))
