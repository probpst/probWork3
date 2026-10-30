#!/usr/bin/env python3
"""make_duct_case.py - verification case of the prescribed grid deformation (morph_frames.F90):
a closed-end square duct whose side walls squeeze in and out periodically.  The fluid is pushed
out of the open end (pressure outlet): the outflow must equal -dV/dt of the analytic volume.

  domain  x in [0, L], y, z in [-h0/2, h0/2];  2 blocks split at x = L/2 (block-to-block interface)
  walls   x = 0 (closed end), the four side walls (no-slip, moving);  x = L pressure outlet p = 0
  motion  y' = y f(x,t), z' = z f(x,t),  f = 1 - A sin^2(pi x / L) sin^2(pi t / T)   (outlet face rigid)
  volume  V(t) = h0^2 int_0^L f^2 dx  (closed form below);  Q_out(t) = -dV/dt

usage: make_duct_case.py OUTDIR [--nx 40 --ny 20 --L 0.1 --h0 0.02 --A 0.3 --T 0.5 --nframes 50 --dt 1e-3 --nranks 2]
Writes Duct.grd, Duct.bcs, Input_polyMesh_probCFD.csv, Duct_frames.bin, probCFD.nml, duct_params.json.
"""
import sys, os, json, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from write_frames_bin import write_header, write_block

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir'); ap.add_argument('--nx', type=int, default=40); ap.add_argument('--ny', type=int, default=20)
    ap.add_argument('--L', type=float, default=0.1); ap.add_argument('--h0', type=float, default=0.02)
    ap.add_argument('--A', type=float, default=0.3); ap.add_argument('--T', type=float, default=0.5)
    ap.add_argument('--nframes', type=int, default=50); ap.add_argument('--dt', type=float, default=1e-3)
    ap.add_argument('--nranks', type=int, default=2); ap.add_argument('--end', type=float, default=0.0)
    ap.add_argument('--mode', type=int, default=0)
    ap.add_argument('--valve', action='store_true', help='valve test: BOTH ends pressure outlets, a porous valve zone at '
                    'x = 0.2-0.3 L closed during the first half period (squeeze) and open during the second (expansion)')
    ap.add_argument('--valve-K', type=float, default=1e5)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    L, h0, A, T = a.L, a.h0, a.A, a.T
    nxb = a.nx // 2                                   # cells per block in x
    x1 = np.linspace(0, L / 2, nxb + 1); x2 = np.linspace(L / 2, L, nxb + 1)
    y = np.linspace(-h0 / 2, h0 / 2, a.ny + 1)
    blocks = []
    for xs in (x1, x2):
        X, Y, Z = np.meshgrid(xs, y, y, indexing='ij')
        blocks.append(np.stack([X, Y, Z]))            # (3, ni, nj, nk)
    nb = len(blocks)
    # ---- grid
    with open(os.path.join(a.outdir, 'Duct.grd'), 'w') as f:
        f.write(f'{nb}\n')
        for g in blocks: f.write(f'{g.shape[1]} {g.shape[2]} {g.shape[3]}\n')
        for g in blocks:
            for c in range(3):
                v = g[c].transpose(2, 1, 0).ravel()
                np.savetxt(f, v.reshape(-1, 4) if len(v) % 4 == 0 else v[None, :], fmt='%.16g')
    # ---- bcs: block 1: i-min wall, j/k walls; block 2: j/k walls, i-max outlet (i-faces between blocks: b2b)
    ni, nj, nk = blocks[0].shape[1:]
    def wall(d, ir, jr, kr, name):
        return [f'{d}   ! {name}', f'{ir[0]} {ir[1]}', f'{jr[0]} {jr[1]}', f'{kr[0]} {kr[1]}',
                '1 0.0   ! u no-slip', '1 0.0   ! v', '1 0.0   ! w', '0 0.0   ! p zero gradient']
    def outlet(d, ir, jr, kr, name):
        return [f'{d}   ! {name}', f'{ir[0]} {ir[1]}', f'{jr[0]} {jr[1]}', f'{kr[0]} {kr[1]}',
                '0 0.0   ! u zero gradient', '0 0.0   ! v', '0 0.0   ! w', '1 0.0   ! p Dirichlet 0 Pa']
    L1 = ['! Duct squeeze test: closed end x = 0, moving side walls, pressure outlet x = L', f'{nb}']
    sides = lambda: (wall(2, (1, ni), (1, 1), (1, nk), 'j-min wall') + wall(-2, (1, ni), (nj, nj), (1, nk), 'j-max wall') +
                     wall(3, (1, ni), (1, nj), (1, 1), 'k-min wall') + wall(-3, (1, ni), (1, nj), (nk, nk), 'k-max wall'))
    if a.valve:
        L1 += ['5   ! block 1'] + outlet(1, (1, 1), (1, nj), (1, nk), 'upstream end') + sides()
    else:
        L1 += ['5   ! block 1'] + wall(1, (1, 1), (1, nj), (1, nk), 'closed end') + sides()
    L1 += ['5   ! block 2'] + sides() + outlet(-1, (ni, ni), (1, nj), (1, nk), 'pressure outlet')
    open(os.path.join(a.outdir, 'Duct.bcs'), 'w').write('\n'.join(L1) + '\n')
    # ---- csv: one rank per block (or split block i-direction)
    per = max(1, a.nranks // nb)
    rows = [['probCFD polyMesh Input'] + [f'block-{b + 1}' for b in range(nb)], ['i_split'] + [str(per)] * nb,
            ['j_split'] + ['1'] * nb, ['k_split'] + ['1'] * nb, ['n_of_BLK'] + [str(per)] * nb, ['Object_ID'] + ['1'] * nb]
    open(os.path.join(a.outdir, 'Input_polyMesh_probCFD.csv'), 'w').write('\n'.join(','.join(r) for r in rows) + '\n')
    # ---- frames: displacement of every node, nframes over one period (periodic file; the last frame < T)
    times = np.arange(a.nframes) * T / a.nframes
    dims = np.array([g.shape[1:] for g in blocks])
    with open(os.path.join(a.outdir, 'Duct_frames.bin'), 'wb') as f:
        write_header(f, dims, times, True, T, 1.0)
        for t in times:
            st = np.sin(np.pi * t / T) ** 2
            for g in blocks:
                fx = 1.0 - A * np.sin(np.pi * g[0] / L) ** 2 * st
                d = np.zeros_like(g); d[1] = g[1] * (fx - 1); d[2] = g[2] * (fx - 1)
                write_block(f, d)
    # ---- analytic volume: V(t) = h0^2 int_0^L (1 - A s(x) st)^2 dx, s = sin^2(pi x/L); int s = L/2, int s^2 = 3L/8
    def V(t):
        st = np.sin(np.pi * t / T) ** 2
        return h0 ** 2 * (L - 2 * A * st * L / 2 + A ** 2 * st ** 2 * 3 * L / 8)
    def dVdt(t):
        st = np.sin(np.pi * t / T) ** 2; dst = 2 * np.sin(np.pi * t / T) * np.cos(np.pi * t / T) * np.pi / T
        return h0 ** 2 * (-A * L * dst + A ** 2 * 2 * st * dst * 3 * L / 8)
    Umax = float(np.max([-dVdt(t) for t in np.linspace(0, T, 400)])) / h0 ** 2
    json.dump(dict(L=L, h0=h0, A=A, T=T, nx=a.nx, ny=a.ny, nframes=a.nframes, dt=a.dt, Umax_outlet=Umax,
                   V0=V(0.0), Vmin=float(min(V(t) for t in np.linspace(0, T, 400)))), open(os.path.join(a.outdir, 'duct_params.json'), 'w'), indent=1)
    end_t = a.end if a.end > 0 else T
    nml = f"""!=====================================================================
! probCFD name list - Duct squeeze: verification of the prescribed grid deformation (H1)
!=====================================================================
$Case_title
  Mode                 = {a.mode}            !0 run, 1 restart, 2 grid-movement check (frames) -> Check_grid_movement/
  Case_name            = "Duct"
/end
$Control
 Char_L                = {h0}
 Char_U                = {max(Umax, 0.05):.4g}
 Del_T                 = {a.dt}
 End_T                 = {end_t}
 FreeSurf_Z            = 0.0
 Grvty                 = 0, 0, 0
/end
$Fluid_Properties
 Num_of_Phase          = 1
 Lou_w                 = 1000.0
 mu_w                  = 0.001
 Lou_a                 = 1.2
 mu_a                  = 0.000018
/end
$Reference_Pressure
 SW_Ref_Prs            = 0
 Loc_X                 = 0.0
 Loc_Y                 = 0.0
 Loc_Z                 = 0.0
 RefPrs                = 0.0
/end
$Cell_Initial_Condition
 u_cel_initial         = 0.0
 v_cel_initial         = 0.0
 w_cel_initial         = 0.0
 a_cel_initial         = 0.0
 p_cel_initial         = 0.0
/end
$Solver_Iteration
 Iner_Itn              = 1
 Time_scheme           = 1
 Max_Itr_Vel           = 200
 Min_Itr_Vel           = 1
 Rsd_Lmt_Vel           = 1.0e-9
 SOR_fct_Vel           = 0.8
 Max_Itr_Prs           = 20000
 Min_Itr_Prs           = 1
 Trl_Lmt_Prs           = 1.0e-7
 Rsd_Lmt_Prs           = 1.0e-6
 SOR_fct_Prs           = 1.0
 Max_Itr_VoF           = 200
 Min_Itr_VoF           = 1
 Rsd_Lmt_VoF           = 1.0e-9
 SOR_fct_VoF           = 1.0
/end
$Export_Interval
 Output_Intvl          = 25
 SW_Prs_Frc_surf       = 0
 Frame_Prs_Frc_surf    = 0
 Restart_Intvl         = 0
 Data_Intvl            = 1
 Log_Intvl             = 25
/end
$Numerical_Scheme
 Cnvt_scheme           = 2
 curvi_form            = 1
 Cnvt_mode             = 3
 Dffs_mode             = 3
 cfl_switch            = 0.3
 vn_switch             = 0.5
 poi_solver            = 2
 poi_precond           = 2
 vof_scheme            = 2
 vof_cfl_max           = 0.5
 t_ramp                = 0.0
 b2b_tol               = 1.0e-8
 TURB_model            = 0
 IDDES_model           = 0
/end
$Overset
 SW_Ovst               = 0
/end
$SixDOF
 Body                  = 1              !the whole duct: one part with the prescribed deformation
 Mass                  = 1.0
 CG                    = {L/2} 0.0 0.0
 Rot_Center            = {L/2} 0.0 0.0
 Motion_DOF            = 0 0 0 0 0 0
 Inertia               = 1e-6 1e-6 1e-6
 SW_FSI                = 1
 FSI_Frames            = "Duct_frames.bin"
 FSI_Frames_Interp     = 3
 Body_End              = -1
/end
"""
    if a.valve:
        # zone: cells of block 1 with centre x in [0.2 L, 0.3 L] (weight 1)
        xc = 0.5 * (x1[1:] + x1[:-1]); sel = np.where((xc >= 0.2 * L) & (xc <= 0.3 * L))[0]
        with open(os.path.join(a.outdir, 'valve_zone.dat'), 'w') as f:
            f.write('# duct valve zone: block i j k weight (1-based cell indices)\n')
            for i in sel:
                for j in range(a.ny):
                    for k in range(a.ny):
                        f.write(f'1 {i + 1} {j + 1} {k + 1} 1.0\n')
        nml += f"""$Valves
 Num_Valves            = 1
 Valve_Period          = {T}
 Valve_name1           = "duct_valve"
 Valve_file1           = "valve_zone.dat"
 Valve_K1              = {a.valve_K:g}            !1/s: closed-valve resistance
 Valve_open1           = 0.5, 1.0        !open during the expansion half (inflow from both ends), closed during the squeeze
 Valve_ramp1           = 0.01
 Valve_mode1           = 1
/end
"""
        print(f'valve zone: {len(sel)} i-planes x {a.ny*a.ny} cells, closed for phase 0 .. 0.5')
    open(os.path.join(a.outdir, 'probCFD.nml'), 'w').write(nml)
    json.dump(dict(json.load(open(os.path.join(a.outdir, 'duct_params.json'))), valve=bool(a.valve)),
              open(os.path.join(a.outdir, 'duct_params.json'), 'w'), indent=1)
    print(f'wrote {a.outdir}: {nb} blocks of {ni}x{nj}x{nk} nodes, {2*nxb*a.ny*a.ny} cells, {a.nframes} frames, '
          f'U_outlet max {Umax:.3f} m/s, V0 {V(0)*1e6:.2f} mL, Vmin {min(V(t) for t in np.linspace(0,T,400))*1e6:.2f} mL')

if __name__ == '__main__':
    main()
