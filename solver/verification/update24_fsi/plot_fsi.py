#!/usr/bin/env python3
"""Plots of the FSI verification (update 24).

usage: python3 plot_fsi.py grid  <case dir> <step> [out.png]   deformed grid of the plate / cylinder blocks (z = 0 plane)
       python3 plot_fsi.py tip   <case dir> [out.png]           monitor-node histories dy(t), the spectrum, mean / amplitude / frequency
       python3 plot_fsi.py shell <case dir> [out.png]           the mid-surface mesh with the thickness and the clamped / steel nodes
       python3 plot_fsi.py mode  <case dir> <step1> [step2 ...] [out=file.png]   deformed centre line of the plate (Disp of the
                                                                 wall surface file, upper side, mid-span) at the given steps
"""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def rec(f):
    n = np.frombuffer(f.read(4), '<i4')[0]; d = f.read(n); f.read(4); return d

def read_grd(fn):
    with open(fn, 'rb') as f:
        nb = np.frombuffer(rec(f), '<i4')[0]; dims = np.frombuffer(rec(f), '<i4').reshape(nb, 3); out = []
        for b in range(nb):
            ni, nj, nk = dims[b]
            a = np.frombuffer(rec(f), '<f4').reshape(3, nk, nj, ni)
            out.append(a)
    return out

def read_p3d(fn):
    with open(fn, 'rb') as f:
        nb = np.frombuffer(rec(f), '<i4')[0]; dims = np.frombuffer(rec(f), '<i4').reshape(nb, 4); out = []
        for b in range(nb):
            ni, nj, nk, nv = dims[b]
            a = np.frombuffer(rec(f), '<f4').reshape(nv, nk, nj, ni)
            out.append(a)
    return out

def case_name(d):
    for l in open(os.path.join(d, 'probCFD.nml')):
        if 'Case_name' in l: return l.split('=')[1].split('!')[0].strip().strip('"')

def plot_grid(d, step, out):
    case = case_name(d)
    g = read_grd(os.path.join(d, 'Volume_solutions_grids', '%s_%08d.grd' % (case, step)))
    fig, ax = plt.subplots(figsize=(11, 5.5))
    cols = ['tab:blue', 'tab:red', 'tab:gray', 'tab:green']
    for b, a in enumerate(g):
        x = a[0, 1]; y = a[1, 1]                      # k = 1 (mid-span) plane
        ni = x.shape[1]; nj = x.shape[0]
        if b >= 2:
            step_i = 1 if b == 2 else 1
            for j in range(0, nj, 2): ax.plot(x[j, :], y[j, :], color=cols[b], lw=0.25, alpha=0.6)
            for i in range(0, ni, 2): ax.plot(x[:, i], y[:, i], color=cols[b], lw=0.25, alpha=0.6)
        else:
            for j in range(nj): ax.plot(x[j, :], y[j, :], color=cols[b], lw=0.4)
            for i in range(ni): ax.plot(x[:, i], y[:, i], color=cols[b], lw=0.4)
    ax.set_aspect('equal'); ax.set_xlim(-0.03, 0.11); ax.set_ylim(-0.04, 0.04)
    ax.set_title('%s step %d: cylinder (blue) and plate (red) O-grids, refinement (gray), background (green), z = 0' % (case, step))
    ax.set_xlabel('x [m]'); ax.set_ylabel('y [m]')
    fig.tight_layout(); fig.savefig(out, dpi=160); print('wrote', out)

def plot_tip(d, out, D=0.022, tstart=None):
    case = case_name(d)
    files = sorted(glob.glob(os.path.join(d, 'Output_data', '%s_fsi_*_node*.dat' % case)))
    fig, axs = plt.subplots(2, 1, figsize=(10, 7))
    res = []
    for fn in files:
        a = np.loadtxt(fn)
        t = a[:, 1]; dy = a[:, 3]
        lab = os.path.basename(fn).replace(case + '_fsi_', '').replace('.dat', '')
        axs[0].plot(t, dy/D, lw=0.9, label=lab)
        # statistics over the last 3 periods (from the spectrum) of the active part
        m = t > (t[-1] - 0.5) if tstart is None else t > tstart
        if m.sum() > 50:
            tt = t[m]; yy = dy[m]
            dt = np.median(np.diff(tt))
            yy0 = yy - yy.mean()
            sp = np.abs(np.fft.rfft(yy0*np.hanning(len(yy0))))
            fr = np.fft.rfftfreq(len(yy0), dt)
            f0 = fr[np.argmax(sp[1:]) + 1]
            axs[1].plot(fr, sp/sp.max(), lw=0.9, label='%s: f = %.2f Hz' % (lab, f0))
            # peaks of the last periods
            n = int(round(3.0/f0/dt)); yl = yy[-n:]
            res.append((lab, f0, yl.max()/D, yl.min()/D, (yl.max() - yl.min())/2/D, yl.mean()/D))
    axs[0].set_xlabel('t [s]'); axs[0].set_ylabel('d_y / D'); axs[0].grid(alpha=0.3); axs[0].legend(fontsize=8)
    axs[0].axhline(0.648, color='k', ls='--', lw=0.6); axs[0].axhline(-0.648, color='k', ls='--', lw=0.6)
    axs[1].set_xlim(0, 40); axs[1].set_xlabel('f [Hz]'); axs[1].set_ylabel('spectrum (normalised)'); axs[1].grid(alpha=0.3)
    axs[1].axvline(11.25, color='k', ls='--', lw=0.6, label='experiment 11.25 Hz'); axs[1].legend(fontsize=8)
    fig.suptitle('%s: monitor nodes (experiment: f = 11.25 Hz, d_y/D amplitude 0.648 = dashed)' % case)
    fig.tight_layout(); fig.savefig(out, dpi=150); print('wrote', out)
    for r in res:
        print('%-22s f = %6.2f Hz   dy/D max %6.3f  min %6.3f  amplitude %6.3f  mean %6.3f (last 3 periods)' % r)

def plot_shell(d, out):
    case = case_name(d)
    fn = glob.glob(os.path.join(d, 'Output_data', '%s_fsi_*_shell.grd' % case))[0]
    g = read_grd(fn); q = read_p3d(fn.replace('.grd', '.p3d'))
    # the two coordinate axes of the largest extent of the mid-surface
    allc = np.concatenate([a[:, 0].reshape(3, -1) for a in g], axis=1)
    ext = allc.max(1) - allc.min(1); ax1, ax2 = sorted(np.argsort(ext)[-2:])
    fig, ax = plt.subplots(figsize=(11, 3.5))
    for b, a in enumerate(g):
        x = a[ax1, 0]; y = a[ax2, 0]; t = q[b][0, 0]; cl = q[b][2, 0]; mt = q[b][1, 0]
        sc = ax.scatter(x, y, c=t*1000, s=12, cmap='viridis')
        ax.scatter(x[cl > 0.5], y[cl > 0.5], marker='x', c='r', s=18, label='clamped' if b == 0 else None)
        ax.scatter(x[mt > 0.5], y[mt > 0.5], marker='s', facecolors='none', edgecolors='k', s=30, label='Submat (steel)' if b == 0 else None)
    plt.colorbar(sc, ax=ax, label='thickness [mm]')
    ax.set_aspect('equal'); ax.legend(fontsize=8); ax.set_title('%s: mid-surface nodes from the ray pairing' % case)
    ax.set_xlabel('xyz'[ax1] + ' [m]'); ax.set_ylabel('xyz'[ax2] + ' [m]')
    fig.tight_layout(); fig.savefig(out, dpi=150); print('wrote', out)

def plot_mode(d, steps, out, D=0.022):
    """the plate's deformed centre line: wall nodes of the plate patch (zone 2 of the surface file) at mid-span,
    x = undeformed x of the upper side, y = Disp_y; the second swiveling mode shows a node along the plate"""
    case = case_name(d)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for st in steps:
        g = read_grd(os.path.join(d, 'Volume_solutions_grids', '%s_surf_%08d.grd' % (case, st)))
        q = read_p3d(os.path.join(d, 'Volume_solutions_grids', '%s_surf_%08d.p3d' % (case, st)))
        gb = g[1]; qb = q[1]                                   # zone 2 = the plate patch (ni x 1 x nk)
        nk = gb.shape[1]; km = nk//2
        x = gb[0, km, 0, :]; y = gb[1, km, 0, :]; dy = qb[7, km, 0, :]; dx = qb[6, km, 0, :]
        up = y - dy > 0                                       # upper side: undeformed y > 0
        xs = (x - dx)[up]; o = np.argsort(xs)
        t = None
        try:
            mon = np.loadtxt(glob.glob(os.path.join(d, 'Output_data', '%s_fsi_*_node1.dat' % case))[0])
            t = mon[mon[:, 0] == st, 1][0]
        except Exception:
            pass
        ax.plot(xs[o]*1e3, dy[up][o]/D, lw=1.2, label='step %d' % st + ('' if t is None else ' (t = %.3f s)' % t))
    ax.axhline(0, color='k', lw=0.5); ax.set_xlabel('x [mm] (undeformed, from the cylinder centre)'); ax.set_ylabel('d_y / D')
    ax.set_title('%s: deformed centre line of the plate (upper side, mid-span)' % case); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=150); print('wrote', out)

if __name__ == '__main__':
    mode = sys.argv[1]; d = sys.argv[2]
    if mode == 'grid': plot_grid(d, int(sys.argv[3]), sys.argv[4] if len(sys.argv) > 4 else 'fsi_grid.png')
    elif mode == 'tip': plot_tip(d, sys.argv[3] if len(sys.argv) > 3 else 'fsi_tip.png')
    elif mode == 'shell': plot_shell(d, sys.argv[3] if len(sys.argv) > 3 else 'fsi_shell.png')
    elif mode == 'mode':
        args = [a for a in sys.argv[3:] if not a.startswith('out=')]; outs = [a[4:] for a in sys.argv[3:] if a.startswith('out=')]
        plot_mode(d, [int(a) for a in args], outs[0] if outs else 'fsi_mode.png')
