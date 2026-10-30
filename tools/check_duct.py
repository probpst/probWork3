#!/usr/bin/env python3
"""check_duct.py CASEDIR - verify a Duct squeeze run: for every written volume solution
(Volume_solutions_grids/Duct_NNNNNNNN.grd + .p3d) compare
  * the grid volume with the analytic V(t)                       (grid motion / frame interpolation)
  * the outflow through the open end with -dV/dt                 (ALE continuity: GCL + flux re-projection)
  * the maximum |div| from the monitor file                      (discrete continuity)
Also reads Check_grid_movement/ when a Mode = 2 run exists."""
import sys, os, json, glob, numpy as np
d = sys.argv[1]; P = json.load(open(os.path.join(d, 'duct_params.json')))
L, h0, A, T, dt = P['L'], P['h0'], P['A'], P['T'], P['dt']
def rec(f):
    n = np.frombuffer(f.read(4), '<i4')[0]; b = f.read(n); f.read(4); return b
def read_grd(fn):
    f = open(fn, 'rb'); nb = np.frombuffer(rec(f), '<i4')[0]; dims = np.frombuffer(rec(f), '<i4').reshape(nb, 3)
    out = []
    for ni, nj, nk in dims:
        a = np.frombuffer(rec(f), '<f4').reshape(3, nk, nj, ni).transpose(0, 3, 2, 1).astype(np.float64); out.append(a)
    return out
def read_p3d(fn):
    f = open(fn, 'rb'); nb = np.frombuffer(rec(f), '<i4')[0]; dims = np.frombuffer(rec(f), '<i4').reshape(nb, 4)
    out = []
    for ni, nj, nk, nv in dims:
        a = np.frombuffer(rec(f), '<f4').reshape(nv, nk, nj, ni).transpose(0, 3, 2, 1).astype(np.float64); out.append(a)
    return out
def hex_volume(g):
    X = g.transpose(1, 2, 3, 0)
    c = [X[:-1, :-1, :-1], X[1:, :-1, :-1], X[1:, 1:, :-1], X[:-1, 1:, :-1], X[:-1, :-1, 1:], X[1:, :-1, 1:], X[1:, 1:, 1:], X[:-1, 1:, 1:]]
    # volume by the divergence theorem over the 6 faces (bilinear faces -> exact for the two-triangle split average)
    def face(a, b, c_, d_):
        return (np.cross(c_ - a, d_ - b)) / 2.0
    tot = 0.0
    for (a, b, c_, d_) in [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]:
        n = face(c[a], c[b], c[c_], c[d_]); xc = (c[a] + c[b] + c[c_] + c[d_]) / 4
        tot = tot + np.einsum('...i,...i->...', xc, n)
    return float(tot.sum() / 3.0)
def V(t):
    st = np.sin(np.pi * t / T) ** 2
    return h0 ** 2 * (L - A * st * L + A ** 2 * st ** 2 * 3 * L / 8)
def dVdt(t):
    st = np.sin(np.pi * t / T) ** 2; dst = 2 * np.sin(np.pi * t / T) * np.cos(np.pi * t / T) * np.pi / T
    return h0 ** 2 * (-A * L * dst + A ** 2 * 2 * st * dst * 3 * L / 8)
vd = os.path.join(d, 'Volume_solutions_grids')
rows = []
for gf in sorted(glob.glob(os.path.join(vd, 'Duct_????????.grd'))):
    step = int(os.path.basename(gf)[5:13]); t = step * dt
    G = read_grd(gf); vol = sum(hex_volume(g) for g in G)
    pf = gf.replace('.grd', '.p3d'); q = None
    if os.path.exists(pf):
        Q = read_p3d(pf); g = G[-1]; u = Q[-1][0]            # block 2, u at the nodes; outlet = i-max face
        # outflow through the i-max face: trapezoidal (node-average) integration of u over the face quads
        y = g[1, -1]; z = g[2, -1]; uu = u[-1]
        dy = y[1:, :-1] - y[:-1, :-1]; dz = z[:-1, 1:] - z[:-1, :-1]
        uf = 0.25 * (uu[:-1, :-1] + uu[1:, :-1] + uu[:-1, 1:] + uu[1:, 1:])
        q = float((uf * dy * dz).sum())
    rows.append((step, t, vol, V(t), q, -dVdt(t)))
print(f'{"step":>6} {"t[s]":>8} {"V_grid[mL]":>11} {"V_exact":>9} {"rel.err":>9} {"Q_out[mL/s]":>12} {"-dV/dt":>9} {"rel.err":>9}')
worst_v = 0; worst_q = 0
for step, t, vol, ve, q, qe in rows:
    ev = (vol - ve) / ve; worst_v = max(worst_v, abs(ev))
    if q is not None and abs(qe) > 0.05 * abs(dVdt(T / 4)):
        eq = (q - qe) / qe; worst_q = max(worst_q, abs(eq)); sq = f'{q*1e6:12.4f} {qe*1e6:9.4f} {eq:9.2e}'
    elif q is not None:
        sq = f'{q*1e6:12.4f} {qe*1e6:9.4f}      (small)'
    else:
        sq = ''
    print(f'{step:6d} {t:8.4f} {vol*1e6:11.5f} {ve*1e6:9.5f} {ev:9.2e} {sq}')
print(f'worst grid-volume error {worst_v:.2e}, worst outflow error (where |dV/dt| > 5% of peak) {worst_q:.2e}')
mon = os.path.join(d, 'Duct_monitor.dat')
if os.path.exists(mon):
    M = np.loadtxt(mon, comments='#')
    if M.ndim == 1: M = M[None, :]
    print(f'monitor: {len(M)} steps, max CFL {M[:,2].max():.3f}, max |div|max {M[:,7].max():.3e}, max itm {int(M[:,9].max())}, '
          f'max res(u*) {M[:,10].max():.2e}, max itp {int(M[:,11].max())}, max res(p) {M[:,12].max():.2e}')
pf = os.path.join(d, 'Output_data', 'Duct_patch_flux.dat')
if os.path.exists(pf):
    # solver-integrated outlet flux (patch_flux.F90): exact face fluxes, every Data_Intvl steps
    F = np.loadtxt(pf, comments='#')
    if F.ndim == 1: F = F[None, :]
    tt = F[:, 1]; Q = F[:, -1]                                   # last column = sum over the groups = outflow
    qe = -dVdt(tt); peak = abs(qe).max()
    sel = abs(qe) > 0.05 * peak
    rel = np.abs(Q[sel] - qe[sel]) / np.abs(qe[sel])
    print(f'patch flux (solver): {len(F)} records, outflow vs -dV/dt where |dV/dt| > 5% of peak: '
          f'max rel.err {rel.max():.3e}, mean {rel.mean():.3e}; cumulative volume out {np.trapezoid(Q, tt)*1e6:.4f} mL '
          f'vs exact {(V(tt[0]) - V(tt[-1]))*1e6:.4f} mL')
    for i in range(0, len(F), max(1, len(F) // 10)):
        print(f'  t {tt[i]:8.4f}  Q_out {Q[i]*1e6:10.4f} mL/s  -dV/dt {qe[i]*1e6:10.4f}  rel.err {((Q[i]-qe[i])/qe[i] if abs(qe[i])>0.05*peak else float("nan")):9.2e}')
    if P.get('valve') and F.shape[1] >= 5:
        # columns: step, time, group 1 (upstream end), group 2 (outlet), sum.  Valve closed for phase 0..0.5:
        hdr = [l for l in open(pf) if l.startswith('#   step')][0].split()
        Qa = F[:, 2]; Qb = F[:, 3]; ph = (tt / T) % 1.0
        closed = (ph > 0.03) & (ph < 0.47); opn = (ph > 0.53) & (ph < 0.97)
        # the segment upstream of the valve (x < 0.25 L) also squeezes: its own -dV/dt leaves through the upstream end
        xq = np.linspace(0, 0.25 * L, 2001)
        def Vup(t):
            st = np.sin(np.pi * t / T) ** 2
            return h0 ** 2 * np.trapezoid((1 - A * np.sin(np.pi * xq / L) ** 2 * st) ** 2, xq)
        dVup = np.gradient(np.array([Vup(x) for x in tt]), tt)
        leak = Qa - (-dVup)
        print(f'valve test: groups {hdr[3]} (upstream end, behind the valve) and {hdr[4]} (outlet)')
        print(f'  closed phase (0.03-0.47): max |Q_upstream - (-dV_upstream/dt)| {np.abs(leak[closed]).max()*1e6:.4f} mL/s '
              f'(leak through the valve) vs max |Q_outlet| {np.abs(Qb[closed]).max()*1e6:.4f} mL/s '
              f'-> leak ratio {np.abs(leak[closed]).max()/max(np.abs(Qb[closed]).max(),1e-30):.3e}')
        print(f'  open phase   (0.53-0.97): max |Q_upstream| {np.abs(Qa[opn]).max()*1e6:.4f} mL/s, max |Q_outlet| {np.abs(Qb[opn]).max()*1e6:.4f} mL/s')
        vf = os.path.join(d, 'Output_data', 'Duct_valves.dat')
        if os.path.exists(vf):
            Vv = np.loadtxt(vf, comments='#'); print(f'  closure c: min {Vv[:,3].min():.3f} max {Vv[:,3].max():.3f}')
cd = os.path.join(d, 'Check_grid_movement')
if os.path.isdir(cd):
    print('Mode 2 grids:')
    for gf in sorted(glob.glob(os.path.join(cd, 'Duct_????????.grd')))[:200]:
        step = int(os.path.basename(gf)[5:13]); t = step * dt
        vol = sum(hex_volume(g) for g in read_grd(gf))
        print(f'  step {step:6d} t {t:8.4f} V_grid {vol*1e6:10.5f} mL  V_exact {V(t)*1e6:10.5f}  rel.err {(vol-V(t))/V(t):9.2e}')
