#!/usr/bin/env python3
"""Verification cases of the fluid-structure interaction (update 24).

usage: python3 tools/fsi_cases.py <case> <outdir> [options]
  pfs2a  : FSI-PfS-2a (Kalmbach & Breuer 2013, J. Fluids Struct. 42; Fabbri et al. 2023,
           Comput. Fluids 265): fixed cylinder D = 22 mm with a flexible para-rubber plate
           50 x 2 mm (E 4.1 MPa, nu 0.48, rho 1090) and a steel rear mass 10 x 2 mm bonded to
           its end, in a water channel 0.338 x 0.240 m (slip walls), u = 1.385 m/s, Re_D = 30470.
           Quasi-2-D (span 0.4 mm, 2 cells, slip walls; option span=<m>, nz=<cells>): 3-block fully
           staggered OVERSET grid - ONE elliptic O-grid around cylinder + plate (24c; its outer
           boundary the convex hull of circles 12 mm off the body, the seam at the cylinder front),
           a refinement block and the channel background.  Body 1 = cylinder + plate, SW_FSI = 1,
           the cylinder surface clamped, the steel Submat region at the rear.  one=0: the two-body
           grid of update 24 (cylinder O-grid + plate O-grid as a child body).  Options: dt=<s>
           end=<s> start=<s> coarse=0/1 (default 1) span=<m> nz=<n> one=0/1.  The span must be
           thin: with 10 mm (2 cells of 5 mm) a spanwise sawtooth mode of the wake grew from the
           cylinder's separation line.
  kp505  : cases/KP505_FSI - the KP505 overset propeller case with the hub as body 1 (imposed
           roll) and every blade a flexible child body (SW_FSI = 1, clamped inside the hub, the
           tip-cap block a passive patch).  Options: E=<Pa> nu=<-> rho=<kg/m3> end=<s>
"""
import sys, os, re, shutil
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bgr_cases import write_multiblock_grid
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')

def write_grid_xyz(path, blocks):
    """blocks: list of (X, Y, Z) arrays (ni, nj, nk) -> formatted PLOT3D multi-block"""
    with open(path, 'w') as f:
        f.write("%8d\n" % len(blocks))
        for (X, Y, Z) in blocks:
            f.write("%8d%8d%8d\n" % X.shape)
        for (X, Y, Z) in blocks:
            for a in (X, Y, Z):
                v = a.transpose(2, 1, 0).ravel(); n4 = (len(v)//4)*4
                if n4: np.savetxt(f, v[:n4].reshape(-1, 4), fmt='%22.15e')
                if n4 < len(v): np.savetxt(f, v[n4:][None, :], fmt='%22.15e')

def geom_spacing(d1, n, L):
    """n cells from a wall, first cell d1, total length L: geometric ratio by bisection -> distances of the n+1 nodes"""
    lo, hi = 1.0, 3.0
    for _ in range(100):
        r = 0.5*(lo + hi)
        tot = d1*(r**n - 1)/(r - 1)
        if tot > L: hi = r
        else: lo = r
    r = 0.5*(lo + hi)
    d = d1*r**np.arange(n)
    return np.concatenate([[0.0], np.cumsum(d)])*L/np.sum(d)

def rounded_rect(x0, x1, h, ds_flat, ds_arc):
    """closed curve of a plate from x0 to x1 (end-circle centres at x0 + r, x1 - r), thickness h, counter-clockwise
    starting at the trailing end mid-point; returns points (n, 2), the last = the first, and the outward normals"""
    r = h/2
    cx0, cx1 = x0 + r, x1 - r
    pts = []
    # trailing arc (right) from angle -90 to 90 (ccw), x = cx1
    na = max(4, int(round(np.pi*r/ds_arc)))
    nf = max(4, int(round((cx1 - cx0)/ds_flat)))
    for a in np.linspace(-0.5*np.pi, 0.5*np.pi, na + 1)[:-1]:
        pts.append((cx1 + r*np.cos(a), r*np.sin(a)))
    for x in np.linspace(cx1, cx0, nf + 1)[:-1]:          # top flat, towards the leading end
        pts.append((x, r))
    for a in np.linspace(0.5*np.pi, 1.5*np.pi, na + 1)[:-1]:   # leading arc
        pts.append((cx0 + r*np.cos(a), r*np.sin(a)))
    for x in np.linspace(cx0, cx1, nf + 1)[:-1]:          # bottom flat
        pts.append((x, -r))
    pts.append(pts[0])
    P = np.array(pts)[::-1]                             # clockwise: right-handed (i, j outward, k = +z) grid
    # outward normals from the tangent (clockwise curve: outward = (-ty, tx))
    T = np.gradient(P, axis=0)
    T[0] = P[1] - P[-2]; T[-1] = T[0]
    N = np.column_stack([-T[:,1], T[:,0]])
    N /= np.linalg.norm(N, axis=1)[:, None]
    return P, N

def ogrid_from_curve(P, N, dist, z):
    """O-grid: layer j at distance dist[j] along the outward normal; returns X, Y, Z (ni, nj, nk)"""
    ni, nj, nk = len(P), len(dist), len(z)
    X = np.zeros((ni, nj, nk)); Y = np.zeros((ni, nj, nk)); Z = np.zeros((ni, nj, nk))
    for j, d in enumerate(dist):
        Q = P + d*N                 # the offset of a convex curve never folds
        for k, zz in enumerate(z):
            X[:, j, k] = Q[:, 0]; Y[:, j, k] = Q[:, 1]; Z[:, j, k] = zz
    return X, Y, Z

def cyl_ogrid(R, nth, dist, z, xc=0.0, yc=0.0):
    th = -np.linspace(0.0, 2*np.pi, nth + 1)       # clockwise (right-handed with j outward, k = +z); node 0 = node nth
    ni, nj, nk = nth + 1, len(dist), len(z)
    X = np.zeros((ni, nj, nk)); Y = np.zeros((ni, nj, nk)); Z = np.zeros((ni, nj, nk))
    for j, d in enumerate(dist):
        for k, zz in enumerate(z):
            X[:, j, k] = xc + (R + d)*np.cos(th); Y[:, j, k] = yc + (R + d)*np.sin(th); Z[:, j, k] = zz
    return X, Y, Z

def cart(x, y, z):
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    return X, Y, Z

def rec(d, ir, jr, kr, lines, title):
    n = {1:'i-min',-1:'i-max',2:'j-min',-2:'j-max',3:'k-min',-3:'k-max'}[d]
    L = [f"  ! --- {n}: {title}", f"      {d:2d}                ! face direction ({n})",
         f"    {ir[0]:4d} {ir[1]:4d}           ! i-range (nodes)",
         f"    {jr[0]:4d} {jr[1]:4d}           ! j-range (nodes)",
         f"    {kr[0]:4d} {kr[1]:4d}           ! k-range (nodes)"]
    for fl, val, c in lines: L.append(f"       {fl}   {val:<12}! {c}")
    return L
WALL = [('1','0.000','u: Dirichlet 0 (no-slip wall)'),('1','0.000','v: Dirichlet 0'),('1','0.000','w: Dirichlet 0'),('0','0.0','p: zero gradient')]
def slip(ax, c='free slip'):
    L = [('0','0.0',f'u: zero gradient ({c})'),('0','0.0',f'v: zero gradient ({c})'),('0','0.0',f'w: zero gradient ({c})'),('0','0.0','p: zero gradient')]
    L[ax] = ('1','0.000',['u','v','w'][ax]+': Dirichlet 0 (impermeable)')
    return L
def inlet(u): return [('1','%.3f' % u, 'u: Dirichlet %.3f m/s (inflow)' % u),('1','0.000','v: Dirichlet 0'),('1','0.000','w: Dirichlet 0'),('0','0.0','p: zero gradient')]
def outlet(): return [('0','0.0','u: zero gradient'),('0','0.0','v: zero gradient'),('0','0.0','w: zero gradient'),('1','0.000','p: Dirichlet 0 Pa')]

def body_outline(R, h, x_tip, nint_):
    """closed outline of the cylinder (R at the origin) + the plate (thickness h, rounded end at x_tip) as ONE
    curve, CLOCKWISE (right-handed grid: i around, j outward, k = +z), starting at the cylinder front (x = -R),
    nint_ intervals of nearly equal length; returns pts (nint_+1, 2) (last = first) and the segment label of
    every node: 0 cylinder, 1 upper flat, 2 rounded end, 3 lower flat"""
    r = h/2; xj = np.sqrt(R*R - r*r); thj = np.arcsin(r/R); cx = x_tip - r
    La = R*(np.pi - thj); Lf = cx - xj; Lt = np.pi*r
    Ltot = 2*La + 2*Lf + Lt
    n_t = max(4, int(round(nint_*Lt/Ltot)))                     # at least 4 intervals on the rounded end
    n_a = int(round((nint_ - n_t)*La/(2*La + 2*Lf))); n_f = (nint_ - n_t - 2*n_a)//2; n_a = (nint_ - n_t - 2*n_f)//2
    n_t = nint_ - 2*n_a - 2*n_f
    pts = []; seg = []
    for a in np.linspace(np.pi, thj, n_a + 1)[:-1]: pts.append((R*np.cos(a), R*np.sin(a))); seg.append(0)
    for x in np.linspace(xj, cx, n_f + 1)[:-1]: pts.append((x, r)); seg.append(1)
    for a in np.linspace(0.5*np.pi, -0.5*np.pi, n_t + 1)[:-1]: pts.append((cx + r*np.cos(a), r*np.sin(a))); seg.append(2)
    # (24c: the lower flat starts at the end of the arc, already written - no duplicate node; a duplicate
    # made a zero-width sliver column, the root of first-row and fringe velocity spikes)
    for x in np.linspace(cx, xj, n_f + 1)[:-1]: pts.append((x, -r)); seg.append(3)
    for a in np.linspace(-thj, -np.pi, n_a + 1)[:-1]: pts.append((R*np.cos(a), R*np.sin(a))); seg.append(0)
    pts.append(pts[0])
    P = np.array(pts)
    dsw = np.linalg.norm(np.diff(P, axis=0), axis=1)
    assert dsw.min() > 1e-3*dsw.max(), 'outline: a (near-)duplicate node'
    return P, np.array(seg)

def hull_curve(R1, R2, x2, m=4000):
    """convex hull of the circles (R1 at the origin) and (R2 at (x2, 0)), clockwise from the front (x = -R1),
    as a dense polyline (m points, closed)"""
    phi = np.arcsin((R1 - R2)/x2)
    a_u = 0.5*np.pi - phi            # tangent points: angle of the upper normal
    pts = []
    # big arc clockwise from 180 deg to a_u
    for a in np.linspace(np.pi, a_u, 400)[:-1]: pts.append((R1*np.cos(a), R1*np.sin(a)))
    p1 = np.array([R1*np.cos(a_u), R1*np.sin(a_u)]); p2 = np.array([x2 + R2*np.cos(a_u), R2*np.sin(a_u)])
    for f in np.linspace(0, 1, 400)[:-1]: pts.append(tuple(p1 + f*(p2 - p1)))
    for a in np.linspace(a_u, -a_u, 400)[:-1]: pts.append((x2 + R2*np.cos(a), R2*np.sin(a)))
    p3 = np.array([x2 + R2*np.cos(-a_u), R2*np.sin(-a_u)]); p4 = np.array([R1*np.cos(-a_u), R1*np.sin(-a_u)])
    for f in np.linspace(0, 1, 400)[:-1]: pts.append(tuple(p3 + f*(p4 - p3)))
    for a in np.linspace(-a_u, -np.pi, 400): pts.append((R1*np.cos(a), R1*np.sin(a)))
    return np.array(pts)

def resample(C, f):
    """points at the arc-length fractions f along the polyline C"""
    L = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))]); L /= L[-1]
    return np.column_stack([np.interp(f, L, C[:, 0]), np.interp(f, L, C[:, 1])])

def elliptic_ogrid(P, Q, nj, d1, nfine=81, iters=6000):
    """O-grid between the outline P and the outer curve Q (same i-count, both closed): the Winslow (inverse
    Laplace) equations on a uniform nfine-level computational grid give smooth, nearly orthogonal i-lines;
    every i-line is then re-sampled with a geometric distribution of first cell d1 (the wall spacing is the
    same on every line) and nj nodes"""
    ni = len(P) - 1                                   # periodic
    s = np.linspace(0, 1, nfine)
    X = P[:-1, None, :]*(1 - s[None, :, None]) + Q[:-1, None, :]*s[None, :, None]
    for it in range(iters):
        xp = np.roll(X, -1, axis=0); xm = np.roll(X, 1, axis=0)
        x_xi = 0.5*(xp - xm)
        x_et = np.zeros_like(X); x_et[:, 1:-1] = 0.5*(X[:, 2:] - X[:, :-2])
        a = (x_et**2).sum(-1); g = (x_xi**2).sum(-1); b = (x_xi*x_et).sum(-1)
        cross = np.zeros_like(X)
        cross[:, 1:-1] = 0.25*(np.roll(X, -1, axis=0)[:, 2:] - np.roll(X, 1, axis=0)[:, 2:]
                               - np.roll(X, -1, axis=0)[:, :-2] + np.roll(X, 1, axis=0)[:, :-2])
        Xn = X.copy()
        den = 2*(a + g)
        Xn[:, 1:-1] = (a[:, 1:-1, None]*(xp[:, 1:-1] + xm[:, 1:-1]) + g[:, 1:-1, None]*(X[:, 2:] + X[:, :-2])
                       - 2*b[:, 1:-1, None]*cross[:, 1:-1])/den[:, 1:-1, None]
        dX = np.abs(Xn - X).max()
        X = X + 1.0*(Xn - X)
        if dX < 1e-10: break
    # re-sample every i-line: geometric distribution, first cell d1
    G = np.zeros((ni, nj, 2))
    for i in range(ni):
        C = X[i]
        L = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))])
        lo, hi = 1.0, 3.0
        for _ in range(100):
            r = 0.5*(lo + hi)
            if d1*(r**(nj-1) - 1)/(r - 1) > L[-1]: hi = r
            else: lo = r
        r = 0.5*(lo + hi)
        sd = np.concatenate([[0], np.cumsum(d1*r**np.arange(nj-1))]); sd *= L[-1]/sd[-1]
        G[i, :, 0] = np.interp(sd, L, C[:, 0]); G[i, :, 1] = np.interp(sd, L, C[:, 1])
    G = np.concatenate([G, G[:1]], axis=0)
    return G, it, dX

def pfs2a(out, opts):
    name = os.path.basename(out.rstrip('/'))
    os.makedirs(out, exist_ok=True)
    dt = float(opts.get('dt', 5.0e-4)); end_t = float(opts.get('end', 1.0)); start_t = float(opts.get('start', 0.15))
    coarse = int(opts.get('coarse', 1))
    one = int(opts.get('one', 1))
    D = 0.022; R = D/2; U = 1.385
    h = 0.002; x_root = 0.004; x_rubber_end = 0.061; x_tip = 0.071      # plate from the cylinder to the rear-mass end
    W = float(opts.get('span', 0.0004)); nzc = int(opts.get('nz', 2))  # quasi-2-D: nz cells (default 2) across the thin span
    assert nzc >= 2, 'nz >= 2: the solver needs at least 3 nodes in every direction'
    z = np.linspace(-W/2, W/2, nzc + 1)
    if one:
        # ---- ONE BODY (24c): an elliptic O-grid around cylinder + plate as a single wall patch - no
        #      overset junction between a cylinder grid and a plate grid (with two body grids the
        #      staircase cut of the plate root carried spurious velocities).  The outer boundary is the
        #      CONVEX hull of the circles R + 12 mm and h/2 + 12 mm (a concave outer boundary left a
        #      fringe cell without a valid donor - 190 m/s); the seam (i = 1 = ni) at the cylinder
        #      front (next to the holes of the plate the seam's fringe spiked); the first cell 0.08 mm
        #      on every line.
        nint_ = 200 if coarse else 300
        P, seg = body_outline(R, h, x_tip, nint_)
        H = hull_curve(R + 0.012, h/2 + 0.012, x_tip - h/2)
        Lp = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
        Q = resample(H, Lp/Lp[-1])
        G, _, _ = elliptic_ogrid(P, Q, 21 if coarse else 31, 0.00008 if coarse else 0.00005)
        ni_, nj_ = G.shape[0], G.shape[1]
        b1 = (np.repeat(G[:, :, None, 0], len(z), axis=2), np.repeat(G[:, :, None, 1], len(z), axis=2),
              np.broadcast_to(z[None, None, :], (ni_, nj_, len(z))).copy())
        hr = 0.0015 if coarse else 0.0010
        b2 = cart(np.arange(-0.030, 0.130 + 1e-9, hr), np.arange(-0.045, 0.045 + 1e-9, hr), z)
        hb = 0.0035 if coarse else 0.0025
        nxb = int(round(0.338/hb)); nyb = int(round(0.240/hb))
        b3 = cart(np.linspace(-0.077, 0.261, nxb + 1), np.linspace(-0.120, 0.120, nyb + 1), z)
        blocks = [b1, b2, b3]; pb = 1                                   # the flexible body's block
        xw = P[:, 0]; yw = P[:, 1]
        i_clamp = np.where(np.hypot(xw, yw) <= R + 1e-9)[0]            # the cylinder surface (incl. the junction nodes)
    else:
        # ---- the two-body grid of update 24 (cylinder O-grid + plate O-grid, the plate root inside the cylinder)
        nth = 72 if coarse else 96
        d1 = 0.00012 if coarse else 0.0001
        dist = geom_spacing(d1, 16 if coarse else 22, 0.013)
        b1 = cyl_ogrid(R, nth, dist, z)
        P, N = rounded_rect(x_root, x_tip, h, 0.0010 if coarse else 0.0007, 0.00045 if coarse else 0.0003)
        distp = geom_spacing(0.00008 if coarse else 0.00005, 14 if coarse else 18, 0.0075)
        b2 = ogrid_from_curve(P, N, distp, z)
        hr = 0.0015 if coarse else 0.0010
        b3 = cart(np.arange(-0.030, 0.130 + 1e-9, hr), np.arange(-0.045, 0.045 + 1e-9, hr), z)
        hb = 0.0035 if coarse else 0.0025
        nxb = int(round(0.338/hb)); nyb = int(round(0.240/hb))
        b4 = cart(np.linspace(-0.077, 0.261, nxb + 1), np.linspace(-0.120, 0.120, nyb + 1), z)
        blocks = [b1, b2, b3, b4]; pb = 2
        xw = P[:, 0]; yw = P[:, 1]
        i_clamp = np.where(np.hypot(xw, yw) <= R + 0.0005)[0]
    write_grid_xyz(os.path.join(out, name + '.grd'), blocks)
    n = [b[0].shape for b in blocks]
    i_steel = np.where(xw >= x_rubber_end - 1e-9)[0]
    def ranges(idx, ni):
        idx = sorted(set(int(i) + 1 for i in idx))          # 1-based; the seam node ni = node 1
        if ni in idx and 1 not in idx: idx.append(1); idx = sorted(idx)
        out_ = []; s_ = idx[0]; p_ = idx[0]
        for i in idx[1:]:
            if i == p_ + 1: p_ = i
            else: out_.append((s_, p_)); s_ = i; p_ = i
        out_.append((s_, p_))
        return out_
    rc = ranges(i_clamp, n[pb-1][0]); rs = ranges(i_steel, n[pb-1][0])
    # ---- bcs ----
    nb = len(blocks)
    L = [f"! probCFD boundary conditions for case '{name}' (grid {name}.grd)",
         "! FSI-PfS-2a (Kalmbach & Breuer 2013): cylinder D = 22 mm at the origin, flexible plate behind it,",
         f"! water channel x -0.077 .. 0.261, y -0.120 .. 0.120 m (slip walls), quasi-2-D span z {-W/2:.4f} .. {W/2:.4f} ({nzc} cells).",
         (f"! block 1 = ONE elliptic O-grid around cylinder + plate: j = 1 = no-slip wall, j-max = overset boundary (not listed),"
          if one else "! blocks 1 (cylinder O-grid) and 2 (plate O-grid): j = 1 = no-slip wall, j-max = overset boundary (not listed),"),
         f"! i = 1 / i = ni coincide (periodic block-to-block face, detected); block {pb+1} = refinement box (overset),",
         f"! block {pb+2} = channel background.  k-min / k-max of every block = the span walls (slip).",
         "! Record layout per face: direction (|d| = 1 i, 2 j, 3 k; + min, - max), i-, j-, k-ranges (NODE indices),",
         "!   then one 'flag value' line per variable u, v, w, p: 1 Dirichlet, 0 zero gradient.",
         f"       {nb}                ! number of grid blocks"]
    def blockrecs(bi, title, recs):
        L.append(f"! ==================== block {bi}: {title}")
        L.append(f"       {len(recs)}                ! number of boundary surfaces of block {bi}")
        for r in recs: L.extend(r)
    if not one:
        ni, nj, nk = n[0]
        blockrecs(1, f"cylinder O-grid, {ni} x {nj} x {nk} nodes", [rec(2,(1,ni),(1,1),(1,nk),WALL,'cylinder surface, no-slip wall'),
                  rec(3,(1,ni),(1,nj),(1,1),slip(2),'span wall z = -W/2'), rec(-3,(1,ni),(1,nj),(nk,nk),slip(2),'span wall z = +W/2')])
    ni, nj, nk = n[pb-1]
    blockrecs(pb, f"{'cylinder + plate' if one else 'plate'} O-grid, {ni} x {nj} x {nk} nodes", [rec(2,(1,ni),(1,1),(1,nk),WALL,'body surface, no-slip wall (flexible plate, SW_FSI; the cylinder part clamped)'),
              rec(3,(1,ni),(1,nj),(1,1),slip(2),'span wall'), rec(-3,(1,ni),(1,nj),(nk,nk),slip(2),'span wall')])
    ni, nj, nk = n[pb]
    blockrecs(pb+1, f"refinement box {hr*1000:.1f} mm, {ni} x {nj} x {nk} nodes (i/j faces overset)", [rec(3,(1,ni),(1,nj),(1,1),slip(2),'span wall'),
              rec(-3,(1,ni),(1,nj),(nk,nk),slip(2),'span wall')])
    ni, nj, nk = n[pb+1]
    blockrecs(pb+2, f"channel background {hb*1000:.1f} mm, {ni} x {nj} x {nk} nodes", [rec(1,(1,1),(1,nj),(1,nk),inlet(U),'inlet x = -0.077 m'),
              rec(-1,(ni,ni),(1,nj),(1,nk),outlet(),'outlet x = 0.261 m'), rec(2,(1,ni),(1,1),(1,nk),slip(1),'channel wall y = -0.12 m (slip)'),
              rec(-2,(1,ni),(nj,nj),(1,nk),slip(1),'channel wall y = +0.12 m (slip)'), rec(3,(1,ni),(1,nj),(1,1),slip(2),'span wall'),
              rec(-3,(1,ni),(1,nj),(nk,nk),slip(2),'span wall')])
    open(os.path.join(out, name + '.bcs'), 'w').write("\n".join(L) + "\n")
    if one:
        csv = ("probCFD polyMesh Input,block-1,block-2,block-3\ni_split,1,2,2\nj_split,1,1,1\nk_split,1,1,1\nn_of_BLK,1,2,2\n"
               "Object_ID,1,0,0\n")
    else:
        csv = ("probCFD polyMesh Input,block-1,block-2,block-3,block-4\ni_split,1,1,2,2\nj_split,1,1,1,1\nk_split,1,1,1,1\nn_of_BLK,1,1,2,2\n"
               "Object_ID,1,1,0,0\nChild_Object_ID,0,1,0,0\n")
    open(os.path.join(out, 'Input_polyMesh_probCFD.csv'), 'w').write(csv)
    # ---- namelist ----
    nofsi = "\n".join(f" NoFSI_block_num{q+1}      = {pb}\n NoFSI_i_range{q+1}        = {a}, {b}\n NoFSI_j_range{q+1}        = 1\n NoFSI_k_range{q+1}        = 1, {n[pb-1][2]}"
                      for q, (a, b) in enumerate(rc))
    submat = "\n".join(f" Submat_block_num{q+1}     = {pb}\n Submat_i_range{q+1}       = {a}, {b}\n Submat_j_range{q+1}       = 1\n Submat_k_range{q+1}       = 1, {n[pb-1][2]}\n"
                       f" Submat_Density{q+1}       = 7850.0\n Submat_Youngs_modulus{q+1} = 2.1e11\n Submat_Poisson_ratio{q+1} = 0.3"
                       for q, (a, b) in enumerate(rs))
    # monitor node: the upper wall node 2 mm before the end of the rubber part (x = 0.059) at mid-span
    iu = [i for i in range(len(P) - 1) if P[i,1] > 0 and abs(P[i,0] - 0.059) < 0.0006]
    imon = (iu[0] + 1) if iu else 1
    itip = int(np.argmax(P[:-1, 0])) + 1
    kmid = nzc//2 + 1
    if one:
        six = f"""$SixDOF
 Body                  = 1              !cylinder + plate: ONE flexible body (24c), fixed in space - the cylinder part clamped
 Mass                  = 0.0            !computed from Density and the enclosed volume
 CG                    = 0.0375 0.0 0.0
 Rot_Center            = 0.0  0.0  0.0
 Motion_DOF            = 0 0 0 0 0 0
 Inertia               = 1e-6 1e-6 1e-6
"""
        end_body = " Body_End              = -1\n/end\n"
    else:
        six = """$SixDOF
 Body                  = 1              !the cylinder: fixed
 Mass                  = 1.0
 CG                    = 0.0  0.0  0.0
 Rot_Center            = 0.0  0.0  0.0
 Motion_DOF            = 0 0 0 0 0 0
 Inertia               = 1 1 1
 Child_body            = 1              !the plate with its rear mass: rigidly attached to the cylinder, FLEXIBLE (update 24)
 Mother_body           = 1
 Mass                  = 0.0            !computed from Density and the enclosed volume of the plate surface
 CG                    = 0.0375 0.0 0.0
 Rot_Center            = 0.011  0.0 0.0
 Motion_DOF            = 0 0 0 0 0 0
 Inertia               = 1e-6 1e-6 1e-6
"""
        end_body = " Child_body_End        = -1\n Body_End              = -1\n/end\n"
    six += f"""!---- fluid-structure interaction (update 24) ----
 SW_FSI                = 1              !0 rigid (default) / 1 flexible
 Youngs_modulus        = 4.1e6          !para-rubber [Pa]
 Poisson_ratio         = 0.48
 Density               = 1090.0         !kg/m3 (Mass is then computed = Density x enclosed volume)
 Num_NoFSI_Group       = {len(rc)}              !clamped node groups: {'the cylinder surface' if one else 'the plate nodes inside the cylinder'} (block {pb}, j = 1 wall)
{nofsi}
 Num_FSI_Submat_regions = {len(rs)}             !the steel rear mass (x >= 0.061 m)
{submat}
 FSI_Start_T           = {start_t}            !the structure released after the flow has developed
 Damping_ratio         = 0.0
 FSI_Min_Thickness     = 0.01
 FSI_Max_Iter          = 6              !sub-iterations per step: 6 - 9 reach FSI_Tol after the release of the pluck (strong added mass: the plate moves 8 x its own mass of water), 6 bounds the cost
 FSI_Tol               = 5.0e-4         !as Fabbri et al. (2023), Eq. (23)
 FSI_Relax             = 0.4
 FSI_Mode              = 1
 FSI_Plane_Strain      = 3              !quasi-2-D: cylindrical bending only (no spanwise displacement, rotation about z only)
 FSI_Test_Motion       = 0.011, 11.25, 2, 0   !PLUCK (start-up aid, 24c): the static tip-load shape of the structure (inextensible), 11 mm at the tip, imposed at 11.25 Hz (ramped over one period) up to its maximum at 2.25 periods (t = start + 0.2 s, zero velocity), then released and solved; from the round-off asymmetry of the wake the oscillation does not build up within the affordable time
 FSI_Monitor_node1     = {pb}, {imon}, 1, {kmid}    !upper wall node 2 mm before the rubber end (x = 0.059), mid-span: the paper's probe
 FSI_Monitor_node2     = {pb}, {itip}, 1, {kmid}    !the trailing end of the rear mass, mid-span
""" + end_body
    nml = f"""!=====================================================================
! probCFD_부산 name list - FSI-PfS-2a: fixed cylinder with a flexible plate and rear mass (update 24)
!   Kalmbach & Breuer (2013) J. Fluids Struct. 42, 369-387; Fabbri et al. (2023) Comput. Fluids 265, 105963
!   D = 22 mm, u = 1.385 m/s (Re_D = 30470), para-rubber plate 50 x 2 mm (E 4.1 MPa, nu 0.48, rho 1090 kg/m3),
!   steel rear mass 10 x 2 mm; measured: second swiveling mode, f = 11.25 Hz, tip amplitude dy/D = 0.65.
!   Quasi-2-D sandbox grid (span {W*1000:.1f} mm, {nzc} cells: thin, so that no spanwise mode can grow); {len(blocks)} overset blocks.
!=====================================================================

$Case_title
  Mode                 = 0              !0 standard run, 1 restart from the latest <case>_NNNNNNNN.rsto
  Case_name            = "{name}"
/end

$Control
 Char_L                = 0.022          !cylinder diameter [m]
 Char_U                = 1.385          !inflow [m/s]
 Del_T                 = {dt:.1e}        !time step [s]
 End_T                 = {end_t}
 FreeSurf_Z            = 0.0
 Grvty                 = 0, 0, 0        !single phase, horizontal plane: no gravity
/end

$Fluid_Properties
 Num_of_Phase          = 1
 Lou_w                 = 1000.0         !water
 mu_w                  = 0.001
 Lou_a                 = 1.2
 mu_a                  = 0.000018
/end

$Reference_Pressure
 SW_Ref_Prs            = 0              !the outlet has p = 0
 Loc_X                 = 0.0
 Loc_Y                 = 0.0
 Loc_Z                 = 0.0
 RefPrs                = 0.0
/end

$Cell_Initial_Condition
 u_cel_initial         = 1.385          !channel flow from the start (the inflow ramps over t_ramp from this value: a running start)
 v_cel_initial         = 0.0
 w_cel_initial         = 0.0
 a_cel_initial         = 0.0
 p_cel_initial         = 0.0
/end

$Solver_Iteration
 Iner_Itn              = 1              !the FSI sub-iterations repeat the fractional step (FSI_Max_Iter)
 Time_scheme           = 1
 Max_Itr_Vel           = 200
 Min_Itr_Vel           = 1
 Rsd_Lmt_Vel           = 1.0e-9
 SOR_fct_Vel           = 0.8
 Max_Itr_Prs           = 20000
 Min_Itr_Prs           = 1
 Trl_Lmt_Prs           = 1.0e-7
 Rsd_Lmt_Prs           = 1.0e-4        !the plate is very compliant (0.3 mm per Pa statically): a tighter tolerance only adds sub-iterations; 1e-5 costs 2 x the CG iterations on the thin-span grid for the same coupling residual
 SOR_fct_Prs           = 1.0
 Max_Itr_VoF           = 200
 Min_Itr_VoF           = 1
 Rsd_Lmt_VoF           = 1.0e-9
 SOR_fct_VoF           = 1.0
/end

$Export_Interval
 Output_Intvl          = {int(round(0.02/dt))}            !volume + surface files every 20 ms
 SW_Prs_Frc_surf       = 1              !wall surface files (pressure / friction + Disp_x/y/z, vonMises of the flexible plate)
 Frame_Prs_Frc_surf    = 0
 Restart_Intvl         = 20             !restart file every 20 steps (the sandbox runs the case in chunks)
 Data_Intvl            = 1
 Log_Intvl             = 20
/end

$Numerical_Scheme
 Cnvt_scheme           = 2              !2nd-order upwind
 curvi_form            = 1
 Cnvt_mode             = 3
 Dffs_mode             = 3
 SW_Quasi_2D           = 3              !strictly two-dimensional (update 24c): spanwise velocity 0, fields averaged across the span
                                        !  (a laminar thin-span run still grew a spanwise mode of w from the cylinder's separation line)
 cfl_switch            = 0.3
 vn_switch             = 0.5
 poi_solver            = 2
 poi_precond           = 2
 vof_scheme            = 2
 vof_cfl_max           = 0.5
 t_ramp                = 0.0
 b2b_tol               = 1.0e-8
/end
 TURB_model            = 0              !laminar NS on the quasi-2-D grid (the paper: LES on 40 M cells; SST on the thin
                                        !  span grew a spanwise mode at the cylinder's separation and doubled the cost)
 IDDES_model           = 0
 Max_Itr_Turb          = 20
 Min_Itr_Turb          = 1
 Rsd_Lmt_Turb          = 1.0e-4
 SOR_fct_Turb          = 1.0
 Cnvt_scheme_Turb      = 1
 Cnvt_mode_Turb        = 2
 Dffs_mode_Turb        = 2
 cfl_switch_Turb       = 0.5
 vn_switch_Turb        = 0.5

$Overset
 SW_Ovst               = -1
 Surfovl_tol           = 1.0e-5
 Ovst_fringe           = 2
 Ovst_vol_tol          = 1.0e-6
 Ovst_Reassemble_Tol   = 0.0002         !re-assembly when the plate has moved 0.2 mm since the last assembly (update 20 / 24)
 Ovst_only             = 0
 Ovst_inner_max        = 150
 Ovst_inner_tol        = 1.0e-1         !thin-span (quasi-2-D) grid: the inner coupled solve does not reach 1e-2 within Ovst_inner_max,
                                        !  and 1e-1 changes the forces by < 0.1 %
 poi_verbose           = 0
/end

""" + six
    open(os.path.join(out, 'probCFD.nml'), 'w').write(nml)
    print(name, 'written: blocks', n, ' clamped i-ranges', rc, ' steel i-ranges', rs, ' monitor i', imon, itip)
    tot = sum((a-1)*(b-1)*(c-1) for a, b, c in n)
    print('  cells', tot)

def kp505(out, opts):
    """cases/KP505_FSI: the KP505 overset case with the hub as body 1 (imposed roll) and every blade a
    flexible child body (SW_FSI = 1), clamped where its root block lies inside the hub"""
    name = os.path.basename(out.rstrip('/'))
    os.makedirs(out, exist_ok=True)
    src = os.path.join(root, 'cases', 'KP505'); base = 'KP505'
    for f in (base + '.grd', base + '.bcs', base + '.bcs_b2b', base + '_body1_imposed.dat'):
        shutil.copy(os.path.join(src, f), os.path.join(out, f.replace(base, name)))
    b = open(os.path.join(out, name + '.bcs')).read().replace("case '%s' (grid %s.grd)" % (base, base), "case '%s' (grid %s.grd)" % (name, name))
    open(os.path.join(out, name + '.bcs'), 'w').write(b)
    rows = open(os.path.join(src, 'Input_polyMesh_probCFD.csv')).read().strip().split('\n')
    rows.append(','.join(['Child_Object_ID', '0'] + [str((bl - 2)//2 + 1) for bl in range(2, 12)] + ['0', '0']))
    open(os.path.join(out, 'Input_polyMesh_probCFD.csv'), 'w').write('\n'.join(rows) + '\n')
    s = open(os.path.join(src, 'probCFD.nml')).read().replace('"%s"' % base, '"%s"' % name).replace(base + '_body1_imposed.dat', name + '_body1_imposed.dat')
    s = re.sub(r'^ Restart_Intvl .*$', ' Restart_Intvl         = 100          !restart file <case>_NNNNNNNN.rsto every 100 steps (the structural states are included)', s, flags=re.M)
    if 'end' in opts:
        s = re.sub(r'^ End_T .*$', ' End_T                 = %s' % opts['end'], s, flags=re.M)
    # impulsive start: the towing-tank ramp of the wall velocity (t_ramp = 0.02 of cases/KP505) lets fluid
    # pass through the walls while it lasts; with the 0.1-mm cells at the blade roots the implicit upwind
    # of that inflow (mirrored ghost value) makes the momentum sweeps diverge at the first step
    s = re.sub(r'^ t_ramp .*$', ' t_ramp                = 0.0          !impulsive start (the towing-tank ramp of cases/KP505 diverges with the tiny root cells, see README)', s, flags=re.M)
    E = float(opts.get('E', 3.0e9)); nu = float(opts.get('nu', 0.35)); rho = float(opts.get('rho', 1300.0))
    head = s[:s.index('$SixDOF')]
    six = """$SixDOF
!  SNAME / Fossen convention: body axes x forward (= -x grid), y starboard, z down.
!  Body 1 = HUB (block 1, csv row 6) with the IMPOSED ROLL about the propeller axis (x)
!  from %s_body1_imposed.dat (n = 41.7 rps, J = 0.7 at V = 1 m/s, ramped over
!  t_ramp = 0.02 s; see cases/KP505).  Child bodies 1-5 = the five BLADES (blocks 2-3,
!  4-5, 6-7, 8-9, 10-11: csv row 7), rigidly attached to the hub (Motion_DOF = 0: they
!  turn with it) and FLEXIBLE (SW_FSI = 1, update 24): the root block's j = 1 wall patch
!  (57 x 42 nodes, O-grid around the sections, k = span) is the sheet of the structural
!  mid-surface; the tip-cap block (35 x 46) overlaps it and is attached to that mid-surface
!  (the code reports it as a PASSIVE patch).  Clamped: the root stations inside the hub,
!  k = 1 .. 13 of the root block (r < 3.2 mm = hub radius).  Material: a stiff engineering
!  plastic (3-D printed model blade, E = %.1e Pa) so that the deformation is visible on the
!  scale of the model (aluminium: E = 70e9, nu = 0.33, rho = 2700).  Blocks 12 (refinement)
!  and 13 (background) are static.
 Body                  = 1
 Mass                  = 1.0            !not used (imposed / fixed DOFs)
 CG                    = 0.0  0.0  -1.0 !on the propeller axis
 Rot_Center            = 0.0  0.0  -1.0 !propeller axis, mid hub
 Motion_DOF            = 0 0 0 1 0 0    !surge sway heave roll pitch yaw: roll imposed
 Inertia               = 1.0 1.0 1.0
""" % (name, E)
    for c in range(1, 6):
        rb = 2*c
        six += """
!---- child body %d: blade %d (blocks %d-%d), flexible --------------------------------
 Child_body            = %d
 Mother_body           = 1
 Mass                  = 0.0            !computed = Density x enclosed volume of the blade surface
 CG                    = 0.0  0.0  -1.0 !overwritten by the shell's mass data (printed in the log)
 Rot_Center            = 0.0  0.0  -1.0 !the blade turns about the hub's axis with the hub
 Motion_DOF            = 0 0 0 0 0 0    !rigidly attached to the hub (no relative DOF)
 Inertia               = 1e-6 1e-6 1e-6
 SW_FSI                = 1
 Youngs_modulus        = %.1e          ![Pa]
 Poisson_ratio         = %.2f
 Density               = %.1f         ![kg/m3]
 Num_NoFSI_Group       = 1              !clamped: the root stations inside the hub
 NoFSI_block_num1      = %d
 NoFSI_i_range1        = 1, 57
 NoFSI_j_range1        = 1
 NoFSI_k_range1        = 1, 13
 Num_FSI_Submat_regions = 0
 FSI_Start_T           = 0.02           !released after the rotation ramp
 Damping_ratio         = 0.02
 FSI_Min_Thickness     = 0.01
 FSI_Max_Iter          = 8
 FSI_Tol               = 1.0e-4
 FSI_Relax             = 0.5
 FSI_Mode              = 1
 FSI_Monitor_node1     = %d, 15, 1, 42    !mid-chord upper node at the tip station of the root block
 Child_body_End        = -1
""" % (c, c, rb, rb + 1, c, E, nu, rho, rb, rb)
    six += """ Body_End              = -1
/end
"""
    open(os.path.join(out, 'probCFD.nml'), 'w').write(head + six)
    print(name, 'written')

if __name__ == '__main__':
    case, out = sys.argv[1], sys.argv[2]
    opts = dict(a.split('=') for a in sys.argv[3:])
    if case == 'pfs2a': pfs2a(out, opts)
    elif case == 'kp505': kp505(out, opts)
    else: raise SystemExit('unknown case')
