#!/usr/bin/env python3
"""Single-block wave-tank cases for the incident-wave model (update 18).

usage: python3 tools/wave_cases.py <case> <outdir> [options]
  regular  : regular wave (depth 1 m, kd = 6.3), lambda = 1 m, a = 0.01 m (ka = 0.063), flat start with a
             two-period ramp, relaxation zones at both ends, three wave probes (validation case)
  current  : the same wave on a mean current U = 0.5 m/s (encounter frequency check)
  init     : the same wave started from the exact field (Wave_Init = 1), no ramp
  jonswap  : irregular long-crested JONSWAP sea, Hs = 0.02 m, Tp = 0.8 s (stability / spectrum check)
options: --nx N (cells per wavelength, default 32)  --tend T (end time [s])  --split N (i-splits, default 2)
Grid: x = 0 .. 5 m (5 wavelengths), y two cells between symmetry planes, z = -1 .. 0.25 m with a
uniform 2 mm band around the still-water level (|z| <= 0.03 m) and tanh stretching to the bottom / top.
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from energy_cases import write_grid, face

g = 9.81

def zone_nodes(z0, z1, n, beta, cluster_at_end):
    """n cells from z0 to z1 clustered towards z1 (cluster_at_end) or z0 (one-sided tanh)."""
    s = np.linspace(0, 1, n + 1)
    if cluster_at_end:
        t = 1 - np.tanh(beta*(1 - s))/np.tanh(beta)
    else:
        t = np.tanh(beta*s)/np.tanh(beta)
    return z0 + (z1 - z0)*t

def nml(name, dt, tend, outi, umean, wave_block, split_note=""):
    return f"""!=====================================================================
! probCFD_부산 name list - {name} (incident waves, update 18)
!   wave tank 5 x 1.25 m (2-D: two cells in y), still water level z = 0, depth 1 m
{split_note}!=====================================================================
$Case_title
  Mode                 = 0
  Case_name            = "{name}"
/end

$Control
 Char_L                = 1.0
 Char_U                = 1.0
 Del_T                 = {dt}
 End_T                 = {tend}
 FreeSurf_Z            = 0.0
 Grvty                 = 0, 0, -9.81
/end

$Fluid_Properties
 Num_of_Phase          = 2
 Lou_w                 = 1000.0
 mu_w                  = 0.001
 Lou_a                 = 1.2
 mu_a                  = 1.8e-5
/end

$Reference_Pressure
 SW_Ref_Prs            = 0              !the top boundary has p = 0 (Dirichlet)
 Loc_X                 = 0.0
 Loc_Y                 = 0.0
 Loc_Z                 = 0.0
 RefPrs                = 0.0
/end

$Cell_Initial_Condition
 u_cel_initial         = {umean}
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
 Rsd_Lmt_Prs           = 1.0e-7
 SOR_fct_Prs           = 1.0
 Max_Itr_VoF           = 200
 Min_Itr_VoF           = 1
 Rsd_Lmt_VoF           = 1.0e-9
 SOR_fct_VoF           = 1.0
/end

$Export_Interval
 Output_Intvl          = {outi}
 Restart_Intvl         = 0
 Data_Intvl            = 1
 Log_Intvl             = 50
/end

$Numerical_Scheme
 Cnvt_scheme           = 2
 curvi_form            = 1
 Prs_form              = 2              !incremental pressure form (gravity; update 16)
 Cnvt_mode             = 3
 Dffs_mode             = 3
 cfl_switch            = 0.3
 vn_switch             = 0.5
 poi_solver            = 2
 poi_precond           = 2
 vof_scheme            = 3              !geometric PLIC-VoF
 vof_cfl_max           = 0.5
 t_ramp                = 0.0
 b2b_tol               = 1.0e-8
/end
 TURB_model            = 0
 IDDES_model           = 0

$Overset
 SW_Ovst               = 0
/end

$Cavitation
 SW_Cav                = 0
/end

$Energy
 SW_Energy             = 0
/end

{wave_block}"""

def wave_regular(a, T, ramp, zin, zout, probes, init=0, extra=""):
    px = ", ".join("%.3f" % p for p in probes)
    return f"""$Wave
 SW_Wave               = 1              !0 off, 1 regular wave, 2 irregular (spectrum / Wave_File)
 Wave_Amp              = {a}            !amplitude [m]
 Wave_Period           = {T:.6f}        !period [s] (wave number from the dispersion relation)
 Wave_Phase            = 0.0            !phase at X = 0, t = 0 [deg]
 Wave_Heading          = 0.0            !propagation direction from +x [deg]
 Water_Depth           = 1.0            !depth below FreeSurf_Z [m] (0 = deep water)
 Wave_Order            = 1              !1 Airy, 2 Stokes 2nd order
 Wave_Ramp_T           = {ramp}         !half-cosine amplitude ramp [s]
 Wave_Init             = {init}         !1: start from the exact wave field (no ramp)
 Relax_In_L            = {zin}          !relaxation zone from the x-min face [m] (0 = none)
 Relax_Out_L           = {zout}         !relaxation zone from the x-max face [m]
 Relax_Side_L          = 0.0            !relaxation zones from the y faces [m]
 Wave_Probe_N          = {len(probes)}
 Wave_Probe_X          = {px}
 Wave_Probe_Y          = {", ".join("0.025" for p in probes)}
{extra}/end
"""

def wave_jonswap(Hs, Tp, ramp, zin, zout, probes, ncomp=60, seed=7):
    px = ", ".join("%.3f" % p for p in probes)
    return f"""$Wave
 SW_Wave               = 2              !irregular sea from a spectrum
 Spectrum              = 2              !1 Bretschneider, 2 JONSWAP, 3 Pierson-Moskowitz
 Hs                    = {Hs}
 Tp                    = {Tp}
 Gamma_JONSWAP         = 3.3
 N_Wave_Comp           = {ncomp}
 N_Wave_Dir            = 1              !1 long-crested, > 1 cos^2 spreading
 Wave_Seed             = {seed}
 Wave_Heading          = 0.0
 Water_Depth           = 1.0
 Wave_Ramp_T           = {ramp}
 Relax_In_L            = {zin}
 Relax_Out_L           = {zout}
 Relax_Side_L          = 0.0
 Wave_Probe_N          = {len(probes)}
 Wave_Probe_X          = {px}
 Wave_Probe_Y          = {", ".join("0.025" for p in probes)}
/end
"""

def main():
    args = sys.argv[1:]
    case, out = args[0], args[1]
    nxl = 32; tend = None; split = 2
    i = 2
    while i < len(args):
        if args[i] == '--nx': nxl = int(args[i+1]); i += 2
        elif args[i] == '--tend': tend = float(args[i+1]); i += 2
        elif args[i] == '--split': split = int(args[i+1]); i += 2
        else: raise SystemExit('unknown option ' + args[i])
    os.makedirs(out, exist_ok=True)
    name = os.path.basename(out.rstrip('/'))
    lam = 1.0; d = 1.0; Lx = 5.0
    k = 2*np.pi/lam; om = np.sqrt(g*k*np.tanh(k*d)); T = 2*np.pi/om; c = om/k
    # grid
    nx = int(round(Lx/lam*nxl))
    x = np.linspace(0, Lx, nx + 1)
    y = np.linspace(0, 0.05, 3)
    zb = zone_nodes(-d, -0.03, 24, 2.2, True)
    zm = np.linspace(-0.03, 0.03, 31)
    zt = zone_nodes(0.03, 0.25, 10, 2.0, False)
    z = np.concatenate([zb[:-1], zm, zt[1:]])
    write_grid(os.path.join(out, name + '.grd'), x, y, z)
    nz = len(z) - 1
    umean = 0.5 if case == 'current' else 0.0
    wave = ["12  %.3f        ! u: incident wave (mean %.2f m/s + orbital velocity)" % (umean, umean),
            "12  0.0         ! v: incident wave", "12  0.0         ! w: incident wave",
            "12  0.0         ! p: incident wave (zero-gradient row, wave-pressure ghost)",
            "12  0.0         ! alpha: incident wave (fill below the exact elevation)"]
    side = ["0   0.0         ! u: symmetry", "1   0.0         ! v: symmetry (impermeable)", "0   0.0         ! w",
            "0   0.0         ! p: zero gradient", "0   0.0         ! alpha: zero gradient"]
    bottom = ["0   0.0         ! u: free slip", "0   0.0         ! v", "1   0.0         ! w: impermeable bottom",
              "0   0.0         ! p: zero gradient", "0   0.0         ! alpha: zero gradient"]
    top = ["0   0.0         ! u: open top (air)", "0   0.0         ! v", "0   0.0         ! w",
           "1   0.0         ! p: Dirichlet 0 Pa", "0   0.0         ! alpha: zero gradient"]
    bcs = ("! %s: wave tank (incident waves, update 18): x = 0 .. %.1f m, depth %.1f m, still water at z = 0\n"
           "!   flag 12 = incident wave on u, v, w, p and alpha together (the u/v/w values are the mean velocity)\n"
           % (name, Lx, d))
    bcs += "       1                ! number of grid blocks\n       6                ! number of boundary surfaces of block 1\n"
    bcs += face(1, (1, 1), (1, 3), (1, nz+1), wave, "x-min: incident wave (inlet)")
    bcs += face(-1, (nx+1, nx+1), (1, 3), (1, nz+1), wave, "x-max: incident wave (outlet, behind the relaxation zone)")
    bcs += face(2, (1, nx+1), (1, 1), (1, nz+1), side, "y-min (symmetry)")
    bcs += face(-2, (1, nx+1), (3, 3), (1, nz+1), side, "y-max (symmetry)")
    bcs += face(3, (1, nx+1), (1, 3), (1, 1), bottom, "z-min: bottom")
    bcs += face(-3, (1, nx+1), (1, 3), (nz+1, nz+1), top, "z-max: open top")
    open(os.path.join(out, name + '.bcs'), 'w').write(bcs)
    probes = [1.5, 2.0, 3.0]
    dt = 0.004
    if case in ('regular', 'current', 'init'):
        a = 0.01
        if case == 'init':
            wb = wave_regular(a, T, 0.0, 1.0, 1.5, probes, init=1)
            te = tend or 6*T
        else:
            wb = wave_regular(a, T, 2*T, 1.0, 1.5, probes)
            te = tend or 12*T
        s = nml(name, dt, "%.4f" % te, 200, umean, wb)
    elif case == 'jonswap':
        wb = wave_jonswap(0.02, 0.8, 1.6, 1.0, 1.5, probes)
        te = tend or 8.0
        s = nml(name, dt, "%.4f" % te, 250, 0.0, wb)
    else:
        raise SystemExit('unknown case')
    open(os.path.join(out, 'probCFD.nml'), 'w').write(s)
    open(os.path.join(out, 'Input_polyMesh_probCFD.csv'), 'w').write(
        "probCFD polyMesh Input,block-1\ni_split,%d\nj_split,1\nk_split,1\nn_of_BLK,%d\nObject_ID,0\n" % (split, split))
    print(f"{name}: lambda {lam} m, T {T:.4f} s, c {c:.4f} m/s, kd {k*d:.2f}; grid {nx} x 2 x {nz} cells "
          f"({nx*2*nz} cells, dx {x[1]-x[0]:.4f} m, dz_min {np.diff(z).min():.4f} m), dt {dt} s "
          f"({T/dt:.0f} steps/period), End_T {te:.3f} s ({te/dt:.0f} steps)")

if __name__ == '__main__':
    main()
