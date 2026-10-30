#!/usr/bin/env python3
"""Test cases for the force models (update 19).

usage: python3 tools/force_cases.py <case> <outdir>
  disk   : open-water actuator disk fixed in a uniform stream (single phase, box
           2 x 1 x 1 m, U = 1 m/s, D = 0.25 m, n = 8 1/s): writes probCFD.nml with the
           Hough-Ordway disk (Force_type 1, B4-70-like K_T / K_Q curves) and
           probCFD_BEM.nml with the BEM disk (Force_type 2, 4 blades, P/D = 1)
  spring : the NACA66mod_Overset foil (single phase, overset) with a free heave DOF
           and a point-force spring (Force_type 3): the foil oscillates about the
           spring equilibrium with the period 2 pi sqrt(m/k) = 0.1 s (a stiff spring:
           the overset foil costs 30-45 s per step, 100 steps = two periods)
"""
import numpy as np, sys, os, shutil, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from energy_cases import write_grid, face

def disk_nml(name, force_block):
    return f"""!=====================================================================
! probCFD_부산 name list - {name}: actuator disk in a uniform stream (force models, update 19)
!=====================================================================
$Case_title
  Mode                 = 0
  Case_name            = "{name}"
/end

$Control
 Char_L                = 1.0
 Char_U                = 1.0
 Del_T                 = 0.005
 End_T                 = 3.0
 FreeSurf_Z            = 100.0
 Grvty                 = 0, 0, 0
/end

$Fluid_Properties
 Num_of_Phase          = 1
 Lou_w                 = 1000.0
 mu_w                  = 0.001
/end

$Reference_Pressure
 SW_Ref_Prs            = 0
/end

$Cell_Initial_Condition
 u_cel_initial         = 1.0
 v_cel_initial         = 0.0
 w_cel_initial         = 0.0
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
/end

$Export_Interval
 Output_Intvl          = 200
 Restart_Intvl         = 0
 Data_Intvl            = 1
 Log_Intvl             = 50
/end

$Numerical_Scheme
 Cnvt_scheme           = 2
 curvi_form            = 1
 Prs_form              = 1
 Cnvt_mode             = 3
 Dffs_mode             = 3
 cfl_switch            = 0.3
 vn_switch             = 0.5
 poi_solver            = 2
 poi_precond           = 2
 t_ramp                = 0.0
 b2b_tol               = 1.0e-8
/end
 TURB_model            = 0
 IDDES_model           = 0

$Overset
 SW_Ovst               = 0
/end

$SixDOF
!---- the disk is fixed in space (Mother_body = 0): no 6-DoF body in this case ----------------
{force_block}/end
"""

HO = """ Force_model           = 1
 Force_type            = 1              !1 Hough-Ordway actuator disk
 Mother_body           = 0              !0 = no mother: fixed in space
 Prop_P1               = 0.75  0.0  0.0 !upstream disk face centre [m]
 Prop_P2               = 0.80  0.0  0.0 !downstream face centre [m] (P1->P2 = direction the water is pushed)
 Prop_D                = 0.25           !diameter [m]
 Prop_Hub              = 0.20           !hub ratio r_h/R
 Prop_RPS              = 8.0            !rev/s (+ right-handed about P1->P2)
 Prop_KT               = 0.45  -0.35  -0.10    !K_T = c0 + c1 J + c2 J^2  (B4-70-like, P/D = 1)
 Prop_KQ               = 0.062 -0.045  -0.012  !K_Q = c0 + c1 J + c2 J^2
 Force_model_End       = -1
"""
BEM = """ Force_model           = 1
 Force_type            = 2              !2 blade-element momentum disk (pure BEM)
 Mother_body           = 0
 Prop_P1               = 0.75  0.0  0.0
 Prop_P2               = 0.80  0.0  0.0
 Prop_D                = 0.25
 Prop_Hub              = 0.20
 Prop_RPS              = 8.0
 Prop_Blades           = 4
 Prop_Sections         = 9
 Prop_rR               = 0.20  0.30  0.40  0.50  0.60  0.70  0.80  0.90  1.00   !r/R
 Prop_cD               = 0.21  0.24  0.27  0.29  0.30  0.29  0.26  0.19  0.00   !chord / D (B4-70-like)
 Prop_PD               = 1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00  1.00   !pitch / D
 Prop_fc               = 0.020 0.022 0.024 0.024 0.022 0.020 0.016 0.010 0.0    !camber / chord
 Prop_CD0              = 0.012 0.010 0.009 0.008 0.008 0.008 0.009 0.010 0.012  !section minimum drag
 Force_model_End       = -1
"""

def main():
    case, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    name = os.path.basename(out.rstrip('/'))
    if case == 'disk':
        nx, ny, nz = 64, 32, 32
        x = np.linspace(0.0, 2.0, nx + 1); y = np.linspace(-0.5, 0.5, ny + 1); z = np.linspace(-0.5, 0.5, nz + 1)
        write_grid(os.path.join(out, name + '.grd'), x, y, z)
        inlet = ["1   1.0         ! u: inlet 1 m/s", "1   0.0         ! v", "1   0.0         ! w", "0   0.0         ! p: zero gradient"]
        outlet = ["0   0.0         ! u: zero gradient", "0   0.0         ! v", "0   0.0         ! w", "1   0.0         ! p: outlet p = 0"]
        def slip(normal):
            return ["%d   0.0         ! u" % (1 if normal == 'x' else 0), "%d   0.0         ! v" % (1 if normal == 'y' else 0),
                    "%d   0.0         ! w" % (1 if normal == 'z' else 0), "0   0.0         ! p: zero gradient"]
        bcs = "! %s: actuator disk in a uniform stream (force models, update 19)\n" % name
        bcs += "       1                ! number of grid blocks\n       6                ! number of boundary surfaces of block 1\n"
        bcs += face(1, (1, 1), (1, ny+1), (1, nz+1), inlet, "x-min inlet")
        bcs += face(-1, (nx+1, nx+1), (1, ny+1), (1, nz+1), outlet, "x-max outlet")
        bcs += face(2, (1, nx+1), (1, 1), (1, nz+1), slip('y'), "y-min (slip)")
        bcs += face(-2, (1, nx+1), (ny+1, ny+1), (1, nz+1), slip('y'), "y-max (slip)")
        bcs += face(3, (1, nx+1), (1, ny+1), (1, 1), slip('z'), "z-min (slip)")
        bcs += face(-3, (1, nx+1), (1, ny+1), (nz+1, nz+1), slip('z'), "z-max (slip)")
        open(os.path.join(out, name + '.bcs'), 'w').write(bcs)
        open(os.path.join(out, 'probCFD.nml'), 'w').write(disk_nml(name, HO))
        open(os.path.join(out, 'probCFD_BEM.nml'), 'w').write(disk_nml(name, BEM))
        open(os.path.join(out, 'Input_polyMesh_probCFD.csv'), 'w').write(
            "probCFD polyMesh Input,block-1\ni_split,2\nj_split,1\nk_split,1\nn_of_BLK,2\nObject_ID,0\n")
        print(f"{name}: {nx} x {ny} x {nz} cells, dx {x[1]-x[0]:.4f} m, disk D 0.25 m at x = 0.75..0.80 m")
    elif case == 'spring':
        src = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'NACA66mod_Overset')
        for f in ('NACA66mod_Overset.grd', 'NACA66mod_Overset.bcs', 'Input_polyMesh_probCFD.csv'):
            shutil.copy(os.path.join(src, f), os.path.join(out, f.replace('NACA66mod_Overset', name)))
        s = open(os.path.join(src, 'probCFD.nml')).read()
        s = s.replace('"NACA66mod_Overset"', '"%s"' % name)
        s = re.sub(r'^ End_T .*$', ' End_T                 = 0.2          !two spring periods (30-45 s per step: the overset foil)', s, flags=re.M)
        s = re.sub(r'^ Log_Intvl .*$', ' Log_Intvl             = 10', s, flags=re.M)
        s = re.sub(r'^ Output_Intvl .*$', ' Output_Intvl          = 100', s, flags=re.M)
        s = re.sub(r'^ Rsd_Lmt_Prs .*$', ' Rsd_Lmt_Prs           = 1.0e-5      !loose: the test is about the body dynamics', s, flags=re.M)
        i = s.index('$SixDOF'); j = s.index('/end', i)
        m = 1.0e4; T = 0.1; k = 4*np.pi**2*m/T**2; zeq = 0.005
        six = f"""$SixDOF
!---- point-force test (update 19): the foil (body 1) is free in heave with a spring force
!     F_z = -k (z - z_eq) - m g applied at the quarter chord by force model 1; m = {m:.0f} kg, k = {k:.4e} N/m:
!     the heave oscillates about z_eq = {zeq} m (SNAME z down) with the period 2 pi sqrt(m/k) = {T} s and
!     the amplitude {zeq} m (the fluid forces of the foil are < 0.1 % of the spring force)
 Motion_Start_T        = 0.0
 Body                  = 1
 Mass                  = {m:.1f}
 CG                    = 0.25  0.1  -0.97
 Rot_Center            = 0.25  0.1  -0.97
 Motion_DOF            = 0 0 2 0 0 0    !free heave only
 Inertia               = 1.0 1.0 1.0
 Prod_Inertia          = 0 0 0
 Added_Mass            = 0 0 0 0 0 0
 Motion_Damp           = 0 0 0 0 0 0
 Body_End              = -1

 Force_model           = 1
 Force_type            = 3              !point force with user formulas (mother's body frame, SNAME)
 Mother_body           = 1
 Force_Point           = 0.25  0.1  -0.97   !quarter chord = Rot_Center (no moment of the force)
 Force_Fz              = "-{k:.6e}*(z - {zeq}) - {m:.1f}*g"
 Force_model_End       = -1
"""
        s = s[:i] + six + s[j:]
        open(os.path.join(out, 'probCFD.nml'), 'w').write(s)
        print(f"{name}: m {m} kg, k {k:.4e} N/m, period {T} s, equilibrium z = {zeq} m (SNAME, down)")
    else:
        raise SystemExit('unknown case')

if __name__ == '__main__':
    main()
