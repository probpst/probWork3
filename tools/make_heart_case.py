#!/usr/bin/env python3
"""make_heart_case.py - build the probCFD heart case from the v4 grid and the 97 volume-morph frames.

  OUT/heart.grd        initial grid = reference grid + displacement of frame 0 (the cycle's t = 0 state), metres
  OUT/heart.bcs        walls no-slip, 14 inlet/outlet discs as pressure boundaries
  OUT/Input_polyMesh_probCFD.csv   MPI splitting for --nranks
  OUT/heart_frames.bin 97 frames relative to that initial grid (periodic, period 60/70 s), float32 mm, scale 1e-3
  OUT/probCFD.nml      name list (blood: rho 1060, mu 3.5e-3; Body 1 = whole heart with FSI_Frames)
usage: make_heart_case.py OUT --grid heart_grid.g --bcs heart_grid.bcs --frames FRAMEDIR --numbering node_numbering.npz
                          --csv cycle_frames.csv [--nranks 20] [--dt 5e-4] [--cycles 1] [--mode 0] [--p rule ...]
"""
import sys, os, argparse, subprocess, csv, numpy as np
here = os.path.dirname(os.path.abspath(__file__))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('--grid', required=True); ap.add_argument('--bcs', required=True)
    ap.add_argument('--frames', required=True); ap.add_argument('--numbering', required=True); ap.add_argument('--csv', required=True)
    ap.add_argument('--nranks', type=int, default=20); ap.add_argument('--dt', type=float, default=5e-4)
    ap.add_argument('--cycles', type=float, default=1.0); ap.add_argument('--mode', type=int, default=0)
    ap.add_argument('--bpm', type=float, default=70.0); ap.add_argument('--p', action='append', default=['default=0'])
    ap.add_argument('--char-u', type=float, default=1.0); ap.add_argument('--char-l', type=float, default=0.03)
    ap.add_argument('--output-intvl', type=int, default=50); ap.add_argument('--log-intvl', type=int, default=10)
    ap.add_argument('--turb', type=int, default=0)
    ap.add_argument('--valves', help='directory with valve_{mitral,tricuspid,aortic,pulmonary}.dat zone files: adds the $Valves '
                    'group with open windows derived from the LV / RV volume curves of the frame table')
    ap.add_argument('--valve-K', type=float, default=1.0e6, help='closed-valve resistance [1/s]')
    ap.add_argument('--valve-ramp', type=float, default=0.01, help='opening / closing ramp [s]')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    T = 60.0 / a.bpm
    d0 = os.path.join(a.frames, 'disp_000.npy')
    cmd = [sys.executable, os.path.join(here, 'heart2probcfd.py'), a.grid, a.bcs, a.out, 'heart', '--nranks', str(a.nranks),
           '--add-disp', d0, '--numbering', a.numbering] + sum([['--p', r] for r in a.p], [])
    print(' '.join(cmd)); subprocess.check_call(cmd)
    cmd = [sys.executable, os.path.join(here, 'write_frames_bin.py'), os.path.join(a.out, 'heart_frames.bin'), '--global', a.frames,
           '--grid', a.grid, '--numbering', a.numbering, '--csv', a.csv, '--period', repr(T), '--scale', '1e-3', '--zero-first']
    print(' '.join(cmd)); subprocess.check_call(cmd)
    rows = list(csv.DictReader(open(a.csv)))
    end_t = a.cycles * T
    nml = f"""!=====================================================================
! probCFD name list - HEART: blood flow in the beating heart, prescribed wall motion (frame file, H1)
!   grid heart.grd = v4 heart grid at the cycle's t = 0 state (143 blocks, 949,295 cells), metres
!   frames heart_frames.bin: {len(rows)} frames of one cardiac cycle ({a.bpm:.0f} bpm, T = {T:.6f} s), periodic
!=====================================================================
$Case_title
  Mode                 = {a.mode}            !0 run, 1 restart, 2 grid-movement check of the frame motion -> Check_grid_movement/
  Case_name            = "heart"
/end
$Control
 Char_L                = {a.char_l}           !reference length [m] (about the aortic root diameter)
 Char_U                = {a.char_u}           !reference velocity [m/s] (peak aortic velocity ~ 1 m/s)
 Del_T                 = {a.dt}
 End_T                 = {end_t:.6f}         !{a.cycles:g} cycle(s)
 FreeSurf_Z            = 0.0
 Grvty                 = 0, 0, 0
/end
$Fluid_Properties
 Num_of_Phase          = 1
 Lou_w                 = 1060.0         !blood density [kg/m3]
 mu_w                  = 0.0035         !blood viscosity [Pa s] (Newtonian)
 Lou_a                 = 1.2
 mu_a                  = 0.000018
/end
$Reference_Pressure
 SW_Ref_Prs            = 0              !the pressure level is set by the Dirichlet inlet / outlet discs
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
 Rsd_Lmt_Vel           = 1.0e-8
 SOR_fct_Vel           = 0.8
 Max_Itr_Prs           = 20000
 Min_Itr_Prs           = 1
 Trl_Lmt_Prs           = 1.0e-7
 Rsd_Lmt_Prs           = 1.0e-5
 SOR_fct_Prs           = 1.0
 Max_Itr_VoF           = 200
 Min_Itr_VoF           = 1
 Rsd_Lmt_VoF           = 1.0e-9
 SOR_fct_VoF           = 1.0
/end
$Export_Interval
 Output_Intvl          = {a.output_intvl}
 SW_Prs_Frc_surf       = 1
 Frame_Prs_Frc_surf    = 0
 Restart_Intvl         = 200
 Data_Intvl            = 1
 Log_Intvl             = {a.log_intvl}
/end
$Numerical_Scheme
 Cnvt_scheme           = 2              !2nd-order upwind
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
 TURB_model            = {a.turb}              !0 laminar; 4 = k-omega SST
 IDDES_model           = 0
/end
$Overset
 SW_Ovst               = 0
/end
$SixDOF
 Body                  = 1              !the whole heart: one part, fixed in space, PRESCRIBED DEFORMATION from the frame file
 Mass                  = 1.0
 CG                    = 0.0 0.0 1.25
 Rot_Center            = 0.0 0.0 1.25
 Motion_DOF            = 0 0 0 0 0 0
 Inertia               = 1e-6 1e-6 1e-6
 SW_FSI                = 1
 FSI_Frames            = "heart_frames.bin"
 FSI_Frames_Interp     = 3
 Body_End              = -1
/end
"""
    if a.valves:
        import shutil
        t = np.array([float(r['time_s']) for r in rows]); s = np.array([float(r['phase_fraction']) for r in rows])
        win = {}
        for ch, (av, mv) in (('LV_target_mL', ('aortic', 'mitral')), ('RV_target_mL', ('pulmonary', 'tricuspid'))):
            Vc = np.array([float(r[ch]) for r in rows])
            tt = np.r_[t, T]; VV = np.r_[Vc, Vc[0]]; dV = np.gradient(VV, tt)[:-1]; thr = 0.02 * np.abs(dV).max()
            ej = np.where(dV < -thr)[0]; fi = np.where(dV > thr)[0]
            win[av] = (max(s[ej[0]] - 0.005, 0.0), s[ej[-1]] + 0.005)          # ejection: the outflow valve is open
            win[mv] = (max(s[fi[0]] - 0.005, 0.0), 1.0)                         # filling to end-diastole: the inflow valve is open
        vn = ['mitral', 'tricuspid', 'aortic', 'pulmonary']
        for v in vn: shutil.copy(os.path.join(a.valves, f'valve_{v}.dat'), a.out)
        g = f'''$Valves
 Num_Valves            = 4
 Valve_Period          = {T:.9f}
'''
        for i, v in enumerate(vn, 1):
            g += f''' Valve_name{i}           = "{v}"
 Valve_file{i}           = "valve_{v}.dat"
 Valve_K{i}              = {a.valve_K:g}
 Valve_open{i}           = {win[v][0]:.3f}, {win[v][1]:.3f}      !phase fractions (from the {'ejection' if v in ('aortic','pulmonary') else 'filling'} phase of the volume curve)
 Valve_ramp{i}           = {a.valve_ramp}
 Valve_mode{i}           = 1
'''
        nml += g + '/end\n'
        print('valves:', {k: (round(v[0], 3), round(v[1], 3)) for k, v in win.items()})
    open(os.path.join(a.out, 'probCFD.nml'), 'w').write(nml)
    print('wrote', os.path.join(a.out, 'probCFD.nml'))

if __name__ == '__main__':
    main()
