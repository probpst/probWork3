#!/usr/bin/env python3
"""Sliced-plane plot of the overset blank status ([Blank_var]) of a case.

Reads Volume_solutions_grids/<case>_00000000.grd and <case>_00000000_iblank.p3d
written by a grid check (Ovst_only = 1) - or the last solution .p3d, whose
last variable is the node blank - plus the exact cell status
<case>_NNNNNNNN_cblank.p3d, cuts the plane axis = value through every
block and draws two figures:

  <out>/<case>_<axis><value>_status.png   blank status of every block
                                           (1 active, 0.5/0.25/0.125 fringe
                                           layers, 0 fringe, -1 hole); with a
                                           <case>_cblank.q present the exact
                                           CELL status decides active / hole
  <out>/<case>_<axis><value>_active.png   the ACTIVE cells only (everything
                                           else blanked): every fluid point
                                           of the plane must be covered by
                                           the active cells of some block,
                                           the body interiors must be empty

Three panels: background block only (the last block), body / refinement
blocks (all others), all blocks.

usage: python3 tools/blank_slice.py <case_dir> <x|y|z> <value>
                 [--zoom x0,x1,y0,y1] [--out dir] [--panels "label:1,2;label:3"]
needs numpy, pyvista, matplotlib (pip install pyvista matplotlib).
(The pre-update-13 layout <case>.x / <case>_iblank.q / <case>_cblank.q in the case folder is read too.)
"""
import sys, os, glob, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from readq import read_q, read_x
import numpy as np, pyvista as pv, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

ap = argparse.ArgumentParser()
ap.add_argument('case_dir'); ap.add_argument('axis', choices='xyz'); ap.add_argument('value', type=float)
ap.add_argument('--zoom', default=None, help='x0,x1,y0,y1 of the plot window (plane coordinates)')
ap.add_argument('--out', default=None, help='output folder (default: <case_dir>/blank_check)')
ap.add_argument('--panels', default=None, help='"label:blocks;label:blocks" (1-based block lists)')
ap.add_argument('--title', default=None)
a = ap.parse_args()

case_dir = a.case_dir.rstrip('/'); name = os.path.basename(os.path.abspath(case_dir))
# update 13 layout: Volume_solutions_grids/<case>_00000000.grd with
# <case>_00000000_iblank.p3d (grid check) or <case>_NNNNNNNN.p3d (solution,
# blank = last variable; since update 16b the exact cell status) and, when present, the
# older <case>_NNNNNNNN_cblank.p3d; the old layout
# (<case>.x, <case>_iblank.q, <case>_cblank.q in the case folder) still works
vd = os.path.join(case_dir, 'Volume_solutions_grids')
xf = sorted(glob.glob(os.path.join(vd, '*.grd'))) if os.path.isdir(vd) else []
if xf:
    xf = xf[0]; base = xf[:-4]
    qf = base + '_iblank.p3d'
    if not os.path.exists(qf):
        sols = sorted(f for f in glob.glob(os.path.join(vd, '*.p3d')) if 'blank' not in os.path.basename(f))
        if not sols: sys.exit(f'{qf} not found: run the case with Ovst_only = 1 first')
        qf = sols[-1]; base = qf[:-4]
    cf = base + '_cblank.p3d'
else:
    xf = sorted(glob.glob(os.path.join(case_dir, '*.x')))
    xf = [f for f in xf if not os.path.basename(f)[:-2].split('_')[-1].isdigit()] or xf   # <case>.x before <case>_NNNNNNNN.x
    if not xf: sys.exit('no Volume_solutions_grids/*.grd and no *.x found: run the case with Ovst_only = 1 first')
    xf = xf[0]; qf = xf[:-2] + '_iblank.q'
    if not os.path.exists(qf):
        sys.exit(f'{qf} not found: run the case with Ovst_only = 1 first')
    cf = xf[:-2] + '_cblank.q'
X = read_x(xf); Q = read_q(qf); nb = len(X)
# exact cell status (cell-centred, written by the solver beside every
# iblank / solution file): used when present - the node blank cannot show
# a single hole cell (a background cell crossed by a thin blade tip) nor
# keep the solved neighbours of a hole whole
C = read_q(cf) if os.path.exists(cf) else None
ax_n = 'xyz'.index(a.axis); normal = [0, 0, 0]; normal[ax_n] = 1; origin = [0, 0, 0]; origin[ax_n] = a.value
a1, a2 = [i for i in range(3) if i != ax_n]

sl = []
for ib, (xx, qq) in enumerate(zip(X, Q)):
    g = pv.StructuredGrid(np.ascontiguousarray(xx[0]), np.ascontiguousarray(xx[1]), np.ascontiguousarray(xx[2]))
    g.point_data['b'] = qq[-1].ravel(order='F')
    if C is not None:
        g.cell_data['c'] = C[ib][-1].ravel(order='F')
    s = g.slice(normal=normal, origin=origin)
    sl.append(s.triangulate() if s.n_points > 0 else None)

levels = [-1, 0, 0.125, 0.25, 0.5, 1]
cols = ['#7f7f7f', '#e6550d', '#fd8d3c', '#fdae6b', '#fdd0a2', '#3182bd']
cmap = ListedColormap(cols); bnd = [-1.5, -0.5, 0.06, 0.19, 0.37, 0.75, 1.5]; norm = BoundaryNorm(bnd, cmap.N)
def snap(v): return np.array(levels)[np.abs(np.array(levels)[None, :] - v[:, None]).argmin(1)]

def draw(ax, blocks, active_only):
    for b in blocks:
        s = sl[b]
        if s is None: continue
        pts = s.points; tri = s.faces.reshape(-1, 4)[:, 1:]; v = s.point_data['b']
        if 'c' in s.cell_data:                       # exact: colour / keep by the parent cell's status
            cst = np.asarray(s.cell_data['c'])
            fc = np.where(cst > 0.5, 1.0, np.where(cst < -0.5, -1.0, snap(v[tri].mean(1))))
            fc = np.where((cst > -0.5) & (cst < 0.5) & (fc == 1.0), 0.0, fc)
            keep = (cst > 0.5) if active_only else np.ones(len(tri), bool)
        else:
            fc = snap(v[tri].mean(1))
            keep = (v[tri].min(1) > 0.999) if active_only else np.ones(len(tri), bool)
        if keep.any():
            ax.tripcolor(pts[:, a1], pts[:, a2], tri[keep], facecolors=fc[keep], cmap=cmap, norm=norm,
                         edgecolors='k', linewidth=0.05)
    ax.set_aspect('equal')

if a.panels:
    pan = [(p.split(':')[0], [int(b) - 1 for b in p.split(':')[1].split(',')]) for p in a.panels.split(';')]
else:
    pan = [('Background only', [nb - 1]), ('Body / refinement grids', list(range(nb - 1))), ('All blocks', list(range(nb)))]
    if nb == 1: pan = [('All blocks', [0])]
if a.zoom:
    z = list(map(float, a.zoom.split(','))); lim = [(z[0], z[1]), (z[2], z[3])]
else:
    lo = np.min([[x[c].min() for c in range(3)] for x in X[:-1]] if nb > 1 else [[x[c].min() for c in range(3)] for x in X], 0)
    hi = np.max([[x[c].max() for c in range(3)] for x in X[:-1]] if nb > 1 else [[x[c].max() for c in range(3)] for x in X], 0)
    m = 0.05 * (hi - lo); lim = [(lo[a1] - m[a1], hi[a1] + m[a1]), (lo[a2] - m[a2], hi[a2] + m[a2])]
out = a.out or os.path.join(case_dir, 'blank_check'); os.makedirs(out, exist_ok=True)
title = a.title or name
pref = os.path.join(out, f'{name}_{a.axis}{a.value:g}' + ('_zoom' if a.zoom else ''))
for tag, act in [('status', False), ('active', True)]:
    n = len(pan); fig, axs = plt.subplots(1, n, figsize=(5.2 * n, 5.4)); axs = np.atleast_1d(axs)
    for ax, (lab, bl) in zip(axs, pan):
        draw(ax, bl, act); ax.set_title(lab, fontsize=10)
        ax.set_xlabel('xyz'[a1] + ' [m]'); ax.set_ylabel('xyz'[a2] + ' [m]')
        ax.set_xlim(*lim[0]); ax.set_ylim(*lim[1])
    if act: h = [Patch(color=cols[-1], label='active cells only (everything else blanked)')]
    else: h = [Patch(color=cols[5], label='active 1'), Patch(color=cols[4], label='0.5'), Patch(color=cols[3], label='0.25'),
               Patch(color=cols[2], label='0.125'), Patch(color=cols[1], label='fringe 0'), Patch(color=cols[0], label='hole -1')]
    fig.legend(handles=h, loc='lower center', ncol=len(h), frameon=False)
    fig.suptitle(f'{title} - {a.axis} = {a.value:g} plane' + (' (active cells only)' if act else ' (blank status)'))
    fig.tight_layout(rect=[0, 0.06, 1, 0.95]); fig.savefig(f'{pref}_{tag}.png', dpi=150); plt.close(fig)
    print('wrote', f'{pref}_{tag}.png')
