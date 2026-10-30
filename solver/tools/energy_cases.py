#!/usr/bin/env python3
"""Single-block Cartesian box cases for the energy equation (update 16):
verification tests and the natural-convection cavity benchmark.

usage: python3 tools/energy_cases.py <case> <outdir>
  slab_transient : 1-D transient conduction, T(0) = 310, T(1) = 290 K, T_init 300 (series solution)
  slab_steady    : 1-D steady conduction, heat flux 100 W/m2 at x = 0, convective h = 10, T_inf = 300 at x = 1
  step_advect    : uniform flow u = 1 m/s, inlet T = 310 K into T_init = 300 K (scheme boundedness)
  cavity         : natural convection in a square cavity, Ra = 1e5, Pr = 0.71 (de Vahl Davis 1983)
"""
import numpy as np, sys, os

def write_grid(path, x, y, z):
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    with open(path, 'w') as f:
        f.write("%8d\n" % 1); f.write("%8d%8d%8d\n" % (len(x), len(y), len(z)))
        for a in (X, Y, Z):
            v = a.transpose(2, 1, 0).ravel(); n4 = (len(v)//4)*4
            if n4: np.savetxt(f, v[:n4].reshape(-1, 4), fmt='%22.15e')
            if n4 < len(v): np.savetxt(f, v[n4:][None, :], fmt='%22.15e')

def face(d, ir, jr, kr, lines, comment):
    s = "  ! --- %s\n" % comment
    s += "      %2d                ! face direction\n" % d
    s += "    %4d %4d           ! i-range (nodes)\n" % ir
    s += "    %4d %4d           ! j-range (nodes)\n" % jr
    s += "    %4d %4d           ! k-range (nodes)\n" % kr
    for l in lines: s += "   " + l + "\n"
    return s

WALL = ["1   0.0         ! u: no-slip wall", "1   0.0         ! v", "1   0.0         ! w", "0   0.0         ! p: zero gradient"]
def SLIP(normal):
    return ["%d   0.0         ! u" % (1 if normal == 'x' else 0), "%d   0.0         ! v" % (1 if normal == 'y' else 0),
            "%d   0.0         ! w" % (1 if normal == 'z' else 0), "0   0.0         ! p: zero gradient"]

CSV = "probCFD polyMesh Input,block-1\ni_split,1\nj_split,1\nk_split,1\nn_of_BLK,1\nObject_ID,0\n"

def tanh_nodes(L, n, beta):
    s = np.linspace(-1, 1, n + 1)
    return L*0.5*(1 + np.tanh(beta*s)/np.tanh(beta))

def nml(name, L, U, dt, tend, rho, mu, cp, k, beta, grav, tinit, out, refprs, energy_block, extra="", prs_form=1):
    return f"""!=====================================================================
! probCFD_부산 name list - {name} (energy equation, update 16)
!=====================================================================
$Case_title
  Mode                 = 0
  Case_name            = "{name}"
/end

$Control
 Char_L                = {L}
 Char_U                = {U}
 Del_T                 = {dt}
 End_T                 = {tend}
 FreeSurf_Z            = 100.0
 Grvty                 = {grav}
/end

$Fluid_Properties
 Num_of_Phase          = 1
 Lou_w                 = {rho}
 mu_w                  = {mu}
 Lou_a                 = 1.2
 mu_a                  = 1.8e-5
/end

$Reference_Pressure
 SW_Ref_Prs            = {refprs}
 Loc_X                 = {L/2}
 Loc_Y                 = 0.0
 Loc_Z                 = {L/2}
 RefPrs                = 0.0
/end

$Cell_Initial_Condition
 u_cel_initial         = 0.0
 v_cel_initial         = 0.0
 w_cel_initial         = 0.0
 a_cel_initial         = 1.0
 p_cel_initial         = 0.0
/end

$Solver_Iteration
 Iner_Itn              = 1
 Time_scheme           = 1
 Max_Itr_Vel           = 200
 Min_Itr_Vel           = 1
 Rsd_Lmt_Vel           = 1.0e-9
 SOR_fct_Vel           = 0.8
 Max_Itr_Prs           = 5000
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
 Output_Intvl          = {out}
 Restart_Intvl         = 0
 Data_Intvl            = 1
 Log_Intvl             = 50
/end

$Numerical_Scheme
 Cnvt_scheme           = 2
 curvi_form            = 1
 Prs_form              = {prs_form}   !1 non-incremental (full pressure), 2 incremental (grad p^n in the predictor; use 2 for
                                      !  gravity / buoyancy-driven flows on wall-clustered grids, update 16)
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
/end
 TURB_model            = 0
 IDDES_model           = 0

$Overset
 SW_Ovst               = 0
/end

$Cavitation
 SW_Cav                = 0
/end

{energy_block}
{extra}"""

def energy_block(tinit, cp, k, beta, scheme=2, cmode=3, dmode=3, maxit=50, rsd=1e-6, extra=""):
    return f"""$Energy
!---------------------------------------------------------------------------------------------
! Heat / energy equation (update 16, energy.F90).  T [K] at the cells, 'Temperature' in the
! .p3d files, Output_data/<case>_heat.dat with the heat flow of every heat surface.  The .bcs
! file needs a T line on every surface: 0 zero gradient (adiabatic), 1 Dirichlet T [K],
! 3 heat flux q_w [W/m2] into the fluid, 4 convective "4  h [W/(m2 K)]  T_inf [K]".
!---------------------------------------------------------------------------------------------
 SW_Energy             = 1            !0 off, 1 on
 T_init                = {tinit}      !initial temperature [K] = Boussinesq reference temperature
 Cp_l                  = {cp}         !liquid (fluid "w") specific heat [J/(kg K)]
 k_l                   = {k}          !liquid thermal conductivity [W/(m K)]
 Cp_a                  = 1005.0       !air (two-phase only)
 k_a                   = 0.0257
 Cp_v                  = 1996.0       !vapour (cavitation only)
 k_v                   = 0.0184
 Beta_l                = {beta}       !thermal expansion [1/K]: 0 = no buoyancy (T passive scalar)
 Beta_a                = 0.0
 Pr_t                  = 0.85         !turbulent Prandtl number (optional)

 Max_Itr_Engy          = {maxit}      !SOR sweeps max / min, residual limit (relative to the first sweep), relaxation
 Min_Itr_Engy          = 1
 Rsd_Lmt_Engy          = {rsd}
 SOR_fct_Engy          = 1.0
 Cnvt_scheme_Engy      = {scheme}     !1 1st-order upwind, 2 2nd-order upwind, 3 QUICK, 4 2nd-order ENO (minmod)
 Cnvt_mode_Engy        = {cmode}      !convection: 1 explicit (sub-cycled to cfl_switch_Engy), 2 implicit, 3 hybrid per direction
 Dffs_mode_Engy        = {dmode}      !diffusion : 1 explicit (sub-cycled to vn_switch_Engy),  2 implicit, 3 hybrid per direction
 cfl_switch_Engy       = 0.3
 vn_switch_Engy        = 0.5
{extra}/end
"""

def main():
    case, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    name = os.path.basename(out.rstrip('/'))
    if case in ('slab_transient', 'slab_steady'):
        L = 1.0; nx = 40
        x = np.linspace(0, L, nx + 1); y = np.linspace(0, 0.1, 3); z = np.linspace(0, 0.1, 3)
        write_grid(os.path.join(out, name + '.grd'), x, y, z)
        if case == 'slab_transient':
            tl = ["1   310.0       ! T: Dirichlet 310 K"]; tr = ["1   290.0       ! T: Dirichlet 290 K"]
            dt, tend, outi = 0.001, 0.1, 50
        else:
            tl = ["3   100.0       ! T: heat flux 100 W/m2 into the fluid"]; tr = ["4   10.0  300.0 ! T: convective h = 10 W/(m2 K), T_inf = 300 K"]
            dt, tend, outi = 0.01, 4.0, 400
        bcs = "! %s: 1-D conduction test (energy equation, update 16); fluid rho = cp = k = 1 -> alpha = 1 m2/s\n" % name
        bcs += "       1                ! number of grid blocks\n       6                ! number of boundary surfaces of block 1\n"
        bcs += face(1, (1, 1), (1, 3), (1, 3), WALL + tl, "x-min")
        bcs += face(-1, (nx+1, nx+1), (1, 3), (1, 3), WALL + tr, "x-max")
        bcs += face(2, (1, nx+1), (1, 1), (1, 3), SLIP('y') + ["0   0.0         ! T: adiabatic"], "y-min (symmetry)")
        bcs += face(-2, (1, nx+1), (3, 3), (1, 3), SLIP('y') + ["0   0.0         ! T: adiabatic"], "y-max (symmetry)")
        bcs += face(3, (1, nx+1), (1, 3), (1, 1), SLIP('z') + ["0   0.0         ! T: adiabatic"], "z-min (symmetry)")
        bcs += face(-3, (1, nx+1), (1, 3), (3, 3), SLIP('z') + ["0   0.0         ! T: adiabatic"], "z-max (symmetry)")
        open(os.path.join(out, name + '.bcs'), 'w').write(bcs)
        eb = energy_block(300.0, 1.0, 1.0, 0.0, scheme=2, cmode=2, dmode=2, maxit=200, rsd=1e-8)
        open(os.path.join(out, 'probCFD.nml'), 'w').write(nml(name, 1.0, 1.0, dt, tend, 1.0, 1.0, 1.0, 1.0, 0.0, "0, 0, 0", 300.0, outi, 1, eb))
    elif case == 'step_advect':
        L = 1.0; nx = 100
        x = np.linspace(0, L, nx + 1); y = np.linspace(0, 0.1, 3); z = np.linspace(0, 0.1, 3)
        write_grid(os.path.join(out, name + '.grd'), x, y, z)
        inlet = ["1   1.0         ! u: inlet 1 m/s", "1   0.0         ! v", "1   0.0         ! w", "0   0.0         ! p: zero gradient",
                 "1   310.0       ! T: inlet 310 K"]
        outlet = ["0   0.0         ! u: zero gradient", "0   0.0         ! v", "0   0.0         ! w", "1   0.0         ! p: outlet p = 0",
                  "0   0.0         ! T: zero gradient"]
        bcs = "! %s: step advection test (energy equation, update 16)\n" % name
        bcs += "       1                ! number of grid blocks\n       6                ! number of boundary surfaces of block 1\n"
        bcs += face(1, (1, 1), (1, 3), (1, 3), inlet, "x-min inlet")
        bcs += face(-1, (nx+1, nx+1), (1, 3), (1, 3), outlet, "x-max outlet")
        bcs += face(2, (1, nx+1), (1, 1), (1, 3), SLIP('y') + ["0   0.0         ! T: adiabatic"], "y-min (symmetry)")
        bcs += face(-2, (1, nx+1), (3, 3), (1, 3), SLIP('y') + ["0   0.0         ! T: adiabatic"], "y-max (symmetry)")
        bcs += face(3, (1, nx+1), (1, 3), (1, 1), SLIP('z') + ["0   0.0         ! T: adiabatic"], "z-min (symmetry)")
        bcs += face(-3, (1, nx+1), (1, 3), (3, 3), SLIP('z') + ["0   0.0         ! T: adiabatic"], "z-max (symmetry)")
        open(os.path.join(out, name + '.bcs'), 'w').write(bcs)
        eb = energy_block(300.0, 1.0, 1.0e-6, 0.0, scheme=4, cmode=3, dmode=3, maxit=50, rsd=1e-8)
        s = nml(name, 1.0, 1.0, 0.002, 0.5, 1.0, 1.0e-6, 1.0, 1.0e-6, 0.0, "0, 0, 0", 300.0, 50, 0, eb)
        s = s.replace(" u_cel_initial         = 0.0", " u_cel_initial         = 1.0")
        open(os.path.join(out, 'probCFD.nml'), 'w').write(s)
    elif case == 'cavity':
        # air at 300 K; Ra = g beta dT L^3 / (nu alpha) = 1e5 exactly by the choice of beta
        L = 0.1; rho = 1.1614; mu = 1.846e-5; cp = 1007.0; k = 0.0263; dT = 1.0; g = 9.81
        nu = mu/rho; alpha = k/(rho*cp); Pr = nu/alpha
        beta = 1.0e5*nu*alpha/(g*dT*L**3)
        Ub = np.sqrt(g*beta*dT*L)                        # buoyancy velocity scale
        n = 64; beta_s = 1.6
        x = tanh_nodes(L, n, beta_s); z = tanh_nodes(L, n, beta_s); y = np.linspace(0, L/n, 3)
        write_grid(os.path.join(out, name + '.grd'), x, y, z)
        Th, Tc = 300.0 + 0.5*dT, 300.0 - 0.5*dT
        bcs = ("! %s: natural convection in a square cavity (de Vahl Davis 1983), Ra = 1e5, Pr = %.3f\n"
               "!   hot wall x = 0 (T = %.2f K), cold wall x = L (T = %.2f K), adiabatic top / bottom (z), gravity -z,\n"
               "!   y = spanwise symmetry planes (2 cells).  All pressure BCs zero gradient: SW_Ref_Prs = 1.\n" % (name, Pr, Th, Tc))
        bcs += "       1                ! number of grid blocks\n       6                ! number of boundary surfaces of block 1\n"
        bcs += face(1, (1, 1), (1, 3), (1, n+1), WALL + ["1   %.4f    ! T: hot wall" % Th], "x-min: hot wall")
        bcs += face(-1, (n+1, n+1), (1, 3), (1, n+1), WALL + ["1   %.4f    ! T: cold wall" % Tc], "x-max: cold wall")
        bcs += face(2, (1, n+1), (1, 1), (1, n+1), SLIP('y') + ["0   0.0         ! T: adiabatic"], "y-min (symmetry)")
        bcs += face(-2, (1, n+1), (3, 3), (1, n+1), SLIP('y') + ["0   0.0         ! T: adiabatic"], "y-max (symmetry)")
        bcs += face(3, (1, n+1), (1, 3), (1, 1), WALL + ["0   0.0         ! T: adiabatic wall"], "z-min: bottom wall (adiabatic)")
        bcs += face(-3, (1, n+1), (1, 3), (n+1, n+1), WALL + ["0   0.0         ! T: adiabatic wall"], "z-max: top wall (adiabatic)")
        open(os.path.join(out, name + '.bcs'), 'w').write(bcs)
        eb = energy_block(300.0, cp, k, "%.6e" % beta, scheme=2, cmode=3, dmode=3, maxit=50, rsd=1e-6)
        dt = 0.05; tend = 200.0
        s = nml(name, L, "%.5e" % Ub, dt, tend, rho, mu, cp, k, beta, "0, 0, -9.81", 300.0, 400, 1, eb, prs_form=2)
        # diffusion fully implicit: with the hybrid switch two explicit directions of VN < 0.5 each add up to
        # VN_x + VN_z > 0.5 in the wall-clustered air cells (explicit 2-D limit 0.5) -> odd-even instability
        s = s.replace(" Dffs_mode             = 3", " Dffs_mode             = 2            !implicit everywhere: VN_x + VN_z of the small air cells exceeds the explicit 2-D limit 0.5")
        s = s.replace(" Dffs_mode_Engy        = 3      !diffusion : 1 explicit (sub-cycled to vn_switch_Engy),  2 implicit, 3 hybrid per direction",
                      " Dffs_mode_Engy        = 2      !diffusion : 1 explicit (sub-cycled to vn_switch_Engy),  2 implicit, 3 hybrid per direction (2 here, see Dffs_mode)")
        open(os.path.join(out, 'probCFD.nml'), 'w').write(s)
        print(f"cavity: nu {nu:.4e} alpha {alpha:.4e} Pr {Pr:.4f} beta {beta:.6e} 1/K  Ra = {g*beta*dT*L**3/(nu*alpha):.1f}  "
              f"U_b {Ub:.4e} m/s  alpha/L {alpha/L:.4e} m/s  first cell {x[1]-x[0]:.3e} m  dt {dt} s ({tend/dt:.0f} steps)")
    else:
        raise SystemExit("unknown case")
    open(os.path.join(out, 'Input_polyMesh_probCFD.csv'), 'w').write(CSV)
    print("wrote", out)

if __name__ == '__main__' and sys.argv[1] not in ('comp_box', 'comp_wave'):
    main()

# ---------------------------------------------------------------------------
# update 17: compressible-mixture tests (python3 tools/energy_cases.py comp_box|comp_wave <outdir>)
def cav_block(extra=""):
    return f"""$Cavitation
 SW_Cav                = 1
 Cav_model             = 1
 Vapor_scheme          = 2
 Cav_Start_T           = 1.0e9        !no mass transfer: the mixture is a fixed vapour-liquid mixture (vap_cel_initial)
 Prs_Vapor             = 2340.0
 Lou_v                 = 0.0173
 mu_v                  = 9.8e-6
 Cav_compressible      = 1            !Wood / Wallis mixture compressibility in the pressure equation
 Snd_liquid            = 1480.0
 vap_cel_initial       = 0.3          !uniform vapour fraction 0.3 of the liquid cells (an isothermal gas, no mass transfer)
{extra}/end
"""

def comp_case(case, out):
    import numpy as np, os
    name = os.path.basename(out.rstrip('/'))
    L = 1.0 if case == 'comp_box' else 2.0
    nx = 50 if case == 'comp_box' else 200
    x = np.linspace(0, L, nx + 1); y = np.linspace(0, 0.1, 3); z = np.linspace(0, 0.1, 3)
    os.makedirs(out, exist_ok=True)
    write_grid(os.path.join(out, name + '.grd'), x, y, z)
    piston = ["1   0.01        ! u: piston / inlet 0.01 m/s", "1   0.0         ! v", "1   0.0         ! w", "0   0.0         ! p: zero gradient"]
    bcs = "! %s: compressible-mixture test (update 17): closed tube, piston u = 0.01 m/s at x = 0, vapour fraction 0.3\n" % name
    bcs += "       1                ! number of grid blocks\n       6                ! number of boundary surfaces of block 1\n"
    bcs += face(1, (1, 1), (1, 3), (1, 3), piston, "x-min: piston")
    bcs += face(-1, (nx+1, nx+1), (1, 3), (1, 3), WALL, "x-max: closed end")
    bcs += face(2, (1, nx+1), (1, 1), (1, 3), SLIP('y'), "y-min (symmetry)")
    bcs += face(-2, (1, nx+1), (3, 3), (1, 3), SLIP('y'), "y-max (symmetry)")
    bcs += face(3, (1, nx+1), (1, 3), (1, 1), SLIP('z'), "z-min (symmetry)")
    bcs += face(-3, (1, nx+1), (1, 3), (3, 3), SLIP('z'), "z-max (symmetry)")
    open(os.path.join(out, name + '.bcs'), 'w').write(bcs)
    dt = 0.005 if case == 'comp_box' else 0.0005
    tend = 1.0 if case == 'comp_box' else 0.08
    outi = 40 if case == 'comp_box' else 40
    s = nml(name, 1.0, 1.0, dt, tend, 1000.0, 1.0e-3, 4182.0, 0.6, 0.0, "0, 0, 0", 293.15, outi, 1,
            "$Energy\n SW_Energy = 0\n/end\n")
    s = s.replace("$Cavitation\n SW_Cav                = 0\n/end\n", cav_block())
    s = s.replace(" Iner_Itn              = 1", " Iner_Itn              = 2")
    s = s.replace(" Log_Intvl             = 50", " Log_Intvl             = 20")
    s = s.replace(" p_cel_initial         = 0.0", " p_cel_initial         = 1.0e5      !the mixture at 1 bar: vapour compressibility 1/p")
    open(os.path.join(out, 'probCFD.nml'), 'w').write(s)
    open(os.path.join(out, 'Input_polyMesh_probCFD.csv'), 'w').write(CSV)
    av = 0.3; inv = av/1.0e5 + (1 - av)/(1000.0*1480.0**2); rhom = (1 - av)*1000.0 + av*0.0173
    c = 1/np.sqrt(rhom*inv)
    print(f"{name}: Wood law av = {av}: rho_m c^2 = {1/inv:.2f} Pa, c_mix = {c:.4f} m/s; piston 0.01 m/s -> "
          f"dp/dt = {0.01/(L*inv):.3f} Pa/s (box), plateau p = rho c u = {rhom*c*0.01:.3f} Pa (wave)")

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] in ('comp_box', 'comp_wave'):
    comp_case(sys.argv[1], sys.argv[2])
