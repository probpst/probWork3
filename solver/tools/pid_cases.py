#!/usr/bin/env python3
"""Test cases of the PID controllers (update 21).

usage: python3 tools/pid_cases.py <case> <outdir>
  steer    : NACA66mod_Overset foil (body 1, overset) steered by a PID controller about the span
             axis (y) through the quarter chord: the controlled variable is the foil's own pitch,
             targets +4 deg then -4 deg (tolerance 0.1 deg, hold 0.01 s), rate limit 300 deg/s,
             P = 0.2, I = 300 (a steering controller tracking an angle: the mechanics of the
             rotation, the target sequence and the rate limit; the lift of the force file changes
             sign with the pitch); Rsd_Lmt_Prs 1e-4 to keep the overset foil at ~30 s per step
  throttle : the Cylinder_MovingBGR cylinder (overset, 1.0 m/s stream, background fixed) spun by a
             throttle controller about the vertical axis (pointing down): target own yaw_rate =
             60 deg/s (= 1/6 rev/s), P = 0.001, I = 0.3 rev/s per deg/s - a Magnus lift appears
  direct   : the ActuatorDisk box (single block, no wall) declared body 1 (m = 100 kg, free in heave,
             damping 2 /s) with a point force (Force_type 3) whose Force_Fz is a direct-force PID
             controller: heave targets +5 mm then -5 mm (SNAME z down), P = 4e4 N/m, D = 2500 N s/m
             (omega = 20 rad/s, damping ratio 0.6), tolerance 0.3 mm, hold 0.05 s
  disk     : the same box (damping 20 /s) pushed down by a Hough-Ordway disk (D = 0.25 m, water
             pushed up: T = 1.76 n^2 N) whose rps is a throttle controller: heave targets +5 mm
             then +10 mm, rps in [0, 12] (one-sided actuator), P = 1000 rev/s per m, I = 1000, D = 200
  disksteer: the same box with the disk pushing water downstream (n = 10 rev/s fixed, T = 115 N)
             steered about -y by a steering controller: the vertical component T sin(delta) pushes
             the box to +5 mm then +10 mm, delta in [-40, 40] deg, rate <= 200 deg/s, P = 3000 deg/m,
             I = 3000, D = 300
"""
import numpy as np, sys, os, shutil, re
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')

def foil_base(out, name, end_t, rsd='1.0e-5'):
    src = os.path.join(root, 'NACA66mod_Overset')
    for f in ('NACA66mod_Overset.grd', 'NACA66mod_Overset.bcs', 'Input_polyMesh_probCFD.csv'):
        shutil.copy(os.path.join(src, f), os.path.join(out, f.replace('NACA66mod_Overset', name)))
    s = open(os.path.join(src, 'probCFD.nml')).read()
    s = s.replace('"NACA66mod_Overset"', '"%s"' % name)
    s = re.sub(r'^ End_T .*$', ' End_T                 = %s' % end_t, s, flags=re.M)
    s = re.sub(r'^ Log_Intvl .*$', ' Log_Intvl             = 10', s, flags=re.M)
    s = re.sub(r'^ Output_Intvl .*$', ' Output_Intvl          = 1000', s, flags=re.M)
    s = re.sub(r'^ Rsd_Lmt_Prs .*$', ' Rsd_Lmt_Prs           = %s      !loose: the tests are about the controller' % rsd, s, flags=re.M)
    i = s.index('$SixDOF'); j = s.index('/end', i) + 4
    return s[:i], s[j:]

def box_base(out, name, end_t):
    """the ActuatorDisk box (2 x 1 x 1 m, 64 x 32 x 32 cells, inflow 1 m/s, slip sides) declared body 1
    (Object_ID 1): a single-block moving grid without any wall, so the body feels no fluid force -
    the controller tests run on a clean 1-DOF mass with artificial damping (about 3 s per step)."""
    src = os.path.join(root, 'cases', 'ActuatorDisk')
    shutil.copy(os.path.join(src, 'ActuatorDisk.grd'), os.path.join(out, name + '.grd'))
    open(os.path.join(out, name + '.bcs'), 'w').write(open(os.path.join(src, 'ActuatorDisk.bcs')).read().replace('ActuatorDisk', name))
    open(os.path.join(out, 'Input_polyMesh_probCFD.csv'), 'w').write(
        open(os.path.join(src, 'Input_polyMesh_probCFD.csv')).read().replace('Object_ID,0', 'Object_ID,1'))
    s = open(os.path.join(src, 'probCFD.nml')).read().replace('"ActuatorDisk"', '"%s"' % name)
    s = s.replace('ActuatorDisk: actuator disk in a uniform stream (force models, update 19)',
                  name + ': PID controller test on the ActuatorDisk box declared body 1 (update 21)')
    s = re.sub(r'^ End_T .*$', ' End_T                 = %s' % end_t, s, flags=re.M)
    s = re.sub(r'^ Log_Intvl .*$', ' Log_Intvl             = 20', s, flags=re.M)
    s = re.sub(r'^ Output_Intvl .*$', ' Output_Intvl          = 1000', s, flags=re.M)
    s = re.sub(r'^ Rsd_Lmt_Prs .*$', ' Rsd_Lmt_Prs           = 1.0e-6', s, flags=re.M)
    i = s.index('$SixDOF'); j = s.index('/end', i) + 4
    return s[:i], s[j:]

def main():
    case, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    name = os.path.basename(out.rstrip('/'))
    if case == 'steer':
        head, tail = foil_base(out, name, 0.1, rsd='1.0e-4')
        six = f"""$SixDOF
!---- PID steering test (update 21): the foil (body 1) is rotated about its span axis (y, through the
!     quarter chord) by a steering controller whose controlled variable is the foil's own pitch:
!     targets +4 deg, then -4 deg; rate limit 300 deg/s (the first steps run at the rate limit)
 Body                  = 1
 Mass                  = 10.0
 CG                    = 0.5   0.1  -0.97
 Rot_Center            = 0.25  0.1  -0.97
 Motion_DOF            = 0 0 0 0 0 0    !ignored: the part is rotated by its PID controller
 Inertia               = 1.0 1.0 1.0
 PID_use               = 1
 PID_type              = 2              !steering: the output is the angle [deg] about PID_Axis
 PID_Axis              = 0 1 0          !span axis (grid frame): + = bow up (SNAME pitch)
 PID_P                 = 0.2            !the rotation of a grid part is kinematic (the angle IS the output):
 PID_I                 = 300.0          !  small P (P < 0.7 for the discrete loop), the integral does the work
 PID_D                 = 0.0
 PID_Target_body       = 1              !the foil itself
 PID_Target_type       = 1              !location (eta)
 PID_Target_var        = pitch
 PID_Target_value      = 4.0  -4.0      ![deg]
 PID_Target_tol        = 0.1
 PID_Target_hold_T     = 0.01
 PID_Out_init          = 0.0
 PID_Out_limit         = -10.0  10.0
 PID_Rate_limit        = 300.0          ![deg/s]: 0.6 deg per step
 Body_End              = -1
/end
"""
        open(os.path.join(out, 'probCFD.nml'), 'w').write(head + six + tail)
    elif case == 'throttle':
        src = os.path.join(root, 'cases', 'Cylinder_MovingBGR')
        shutil.copy(os.path.join(src, 'Cylinder_MovingBGR.grd'), os.path.join(out, name + '.grd'))
        shutil.copy(os.path.join(src, 'Input_polyMesh_probCFD.csv'), os.path.join(out, 'Input_polyMesh_probCFD.csv'))
        open(os.path.join(out, name + '.bcs'), 'w').write(open(os.path.join(src, 'Cylinder_MovingBGR.bcs')).read().replace('Cylinder_MovingBGR', name))
        s = open(os.path.join(src, 'probCFD.nml')).read().replace('"Cylinder_MovingBGR"', '"%s"' % name)
        s = re.sub(r'^  Moving_BGR .*$', '  Moving_BGR           = 0', s, flags=re.M)
        s = re.sub(r'^ End_T .*$', ' End_T                 = 0.2            !40 steps', s, flags=re.M)
        s = re.sub(r'^ u_cel_initial .*$', ' u_cel_initial         = 1.0', s, flags=re.M)
        s = re.sub(r'^ t_ramp .*$', ' t_ramp                = 0.0', s, flags=re.M)
        i = s.index('$SixDOF'); j = s.index('/end', i) + 4
        six = """$SixDOF
!---- PID throttle test (update 21): the cylinder (body 1) is spun about the vertical axis by a
!     throttle controller; target: its own yaw rate 60 deg/s (= 1/6 rev/s about -z grid)
 Body                  = 1
 Mass                  = 1000.0
 CG                    = 0.0  0.0  -0.25
 Rot_Center            = 0.0  0.0  -0.25
 Motion_DOF            = 0 0 0 0 0 0    !ignored: the part is rotated by its PID controller
 Inertia               = 100 100 100
 PID_use               = 1
 PID_type              = 1              !throttle: the output is the rotation rate [rev/s] about PID_Axis
 PID_Axis              = 0 0 -1         !vertical, pointing DOWN (grid -z = SNAME +z): a positive rev/s is a
                                        !  positive yaw rate, so the gains are positive
 PID_P                 = 0.001          ![rev/s per deg/s]: the rotation of a grid part is kinematic (the
 PID_I                 = 0.3            !  controller sets the rate), so the loop is algebraic - small P, large I
 PID_D                 = 0.0
 PID_Target_body       = 1
 PID_Target_type       = 2              !velocity (nu)
 PID_Target_var        = yaw_rate
 PID_Target_value      = 60.0           ![deg/s]
 PID_Target_tol        = 1.0
 PID_Target_hold_T     = 0.0
 PID_Out_init          = 0.0
 PID_Out_limit         = -2.0  2.0
 PID_Rate_limit        = 0.0
 Body_End              = -1
/end
"""
        open(os.path.join(out, 'probCFD.nml'), 'w').write(s[:i] + six + s[j:])
    elif case in ('direct', 'disk', 'disksteer'):
        head, tail = box_base(out, name, {'direct': 1.0, 'disk': 1.5, 'disksteer': 1.5}[case])
        body = """ Motion_Start_T        = 0.0
 Body                  = 1
 Mass                  = 100.0
 CG                    = 1.0  0.0  0.0
 Rot_Center            = 1.0  0.0  0.0
 Motion_DOF            = 0 0 2 0 0 0    !free heave only
 Motion_Damp           = 0 0 %s 0 0 0   !artificial heave damping [1/s] (the box has no wall: no fluid force)
 Inertia               = 10 10 10
 Body_End              = -1
"""
        if case == 'direct':
            note = ("direct-force PID test (update 21): the box (body 1, m = 100 kg, no wall) is free in heave\n"
                    "!     with a damping of 2 /s; a point force at its centre carries a direct-force controller on\n"
                    "!     Force_Fz that drives the heave to +5 mm and then -5 mm (SNAME z down):\n"
                    "!     P = 4e4 N/m (omega = 20 rad/s), D = 2500 N s/m (damping ratio 0.6); tolerance 0.3 mm, hold 0.05 s")
            fm = """ Force_model           = 1
 Force_type            = 3              !point force: the formula gives 0, the controller the rest
 Mother_body           = 1
 Force_Point           = 1.0  0.0  0.0
 PID_use               = 1
 PID_type              = 3              !direct force: the output is added to PID_Component
 PID_Component         = Force_Fz       !SNAME z down
 PID_P                 = 4.0e4          ![N/m]
 PID_I                 = 0.0
 PID_D                 = 2500.0         ![N s/m]
 PID_Target_body       = 1
 PID_Target_type       = 1              !location (eta)
 PID_Target_var        = z              !heave [m]
 PID_Target_value      = 0.005  -0.005
 PID_Target_tol        = 0.0003
 PID_Target_hold_T     = 0.05
 PID_Out_limit         = -1.0e4  1.0e4
 Force_model_End       = -1
"""
            damp = '2'
        else:
            damp = '20'
            if case == 'disk':
                note = ("force-model PID test (update 21): the box (body 1, m = 100 kg, no wall) is free in heave with\n"
                        "!     a damping of 20 /s and is pushed DOWN (SNAME +z) by a Hough-Ordway disk (D = 0.25 m, water\n"
                        "!     pushed up, J = 0 -> T = rho n^2 D^4 K_T(0) = 1.76 n^2 N) whose rps is a throttle controller:\n"
                        "!     heave targets +5 mm then +10 mm, rps in [0, 12] (one-sided actuator), P = 1000 rev/s per m,\n"
                        "!     I = 1000, D = 200")
                ctrl = """ Prop_P1               = 1.0  0.0  -0.125  !water pushed UP (P1 -> P2 = +z grid): the reaction pushes the box down
 Prop_P2               = 1.0  0.0  -0.075
 Prop_RPS              = 0.0            !initial rps (= the controller's initial output; PID_Out_init is ignored)
 PID_use               = 1
 PID_type              = 1              !throttle: the output replaces Prop_RPS
 PID_P                 = 1000.0         ![rev/s per m]
 PID_I                 = 1000.0
 PID_D                 = 200.0          ![rev/s per m/s]
 PID_Target_body       = 1
 PID_Target_type       = 1
 PID_Target_var        = z              !heave, SNAME z down
 PID_Target_value      = 0.005  0.010
 PID_Target_tol        = 0.0003
 PID_Target_hold_T     = 0.05
 PID_Out_limit         = 0.0  12.0      !rps >= 0
 PID_Rate_limit        = 0.0
"""
            else:
                note = ("force-model PID test (update 21): the box (body 1, m = 100 kg, no wall) is free in heave with\n"
                        "!     a damping of 20 /s; a Hough-Ordway disk (D = 0.25 m, water pushed downstream, n = 10 rev/s\n"
                        "!     fixed, J = 0.4 -> T = 115 N) is steered about -y (grid) through its centre by a steering\n"
                        "!     controller: the vertical component T sin(delta) pushes the box down to +5 mm then +10 mm,\n"
                        "!     delta in [-40, 40] deg, rate <= 200 deg/s, P = 3000 deg/m, I = 3000, D = 300")
                ctrl = """ Prop_P1               = 0.875 0.0  -0.1  !water pushed downstream (+x grid)
 Prop_P2               = 1.125 0.0  -0.1
 Prop_RPS              = 10.0           !fixed (a steering controller does not touch the rps)
 PID_use               = 1
 PID_type              = 2              !steering: the disk axis turns about PID_Axis through PID_Axis_Center
 PID_Axis              = 0 -1 0         !-y (grid): a positive angle tilts the push direction (+x) towards +z (up),
                                        !so the reaction on the box points down (SNAME +z)
 PID_Axis_Center       = 1.0  0.0  -0.1 !(default: the disk centre)
 PID_P                 = 3000.0         ![deg per m]
 PID_I                 = 3000.0
 PID_D                 = 300.0          ![deg per m/s]
 PID_Target_body       = 1
 PID_Target_type       = 1
 PID_Target_var        = z
 PID_Target_value      = 0.005  0.010
 PID_Target_tol        = 0.0003
 PID_Target_hold_T     = 0.05
 PID_Out_init          = 0.0
 PID_Out_limit         = -40.0  40.0    ![deg]
 PID_Rate_limit        = 200.0          ![deg/s]
"""
            fm = f""" Force_model           = 1
 Force_type            = 1              !Hough-Ordway disk
 Mother_body           = 1
{ctrl} Prop_D                = 0.25
 Prop_Hub              = 0.2
 Prop_KT               = 0.45  -0.35  -0.10
 Prop_KQ               = 0.062 -0.045  -0.012
 Force_model_End       = -1
"""
        six = "$SixDOF\n!---- " + note + "\n" + body % damp + "\n" + fm + "/end\n"
        open(os.path.join(out, 'probCFD.nml'), 'w').write(head + six + tail)
    else:
        raise SystemExit('unknown case')
    print(name, 'written')

if __name__ == '__main__':
    main()
