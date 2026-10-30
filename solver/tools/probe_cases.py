#!/usr/bin/env python3
"""Test cases of the probe points and the wall surface files (update 22).

usage: python3 tools/probe_cases.py <case> <outdir>
  cylinder : cases/Cylinder_MovingBGR (towed cylinder, the background following it, 20 steps) with
             two probe points following body 1 (0.1 m behind the cylinder and 0.05 m beside it),
             two background points (one inside the cylinder -> NoData, one 1.5 m behind) and the
             surface files (SW_Prs_Frc_surf = 1, earth-fixed frame)
  foil     : cases/PID_Steer (the foil steered to +4 deg, 5 steps) with two probe points following
             the foil (near the trailing edge and above the quarter chord), one background point in
             the wake and the surface files
"""
import sys, os, shutil, re
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')

def main():
    case, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    name = os.path.basename(out.rstrip('/'))
    if case == 'cylinder':
        src = os.path.join(root, 'cases', 'Cylinder_MovingBGR'); base = 'Cylinder_MovingBGR'
        for f in (base + '.grd', base + '.bcs', 'Input_polyMesh_probCFD.csv', base + '_body1_imposed.dat'):
            shutil.copy(os.path.join(src, f), os.path.join(out, f.replace(base, name)))
        s = open(os.path.join(src, 'probCFD.nml')).read().replace('"%s"' % base, '"%s"' % name)
        s = re.sub(r'^ End_T .*$', ' End_T                 = 0.1            !20 steps (sandbox test; 0.3 = the update-20 run)', s, flags=re.M)
        s = re.sub(r'^ Output_Intvl .*$', """ Output_Intvl          = 10
 SW_Prs_Frc_surf       = 1              !update 22: wall surface files <case>_surf_NNNNNNNN.grd / .p3d with every volume solution
 Frame_Prs_Frc_surf    = 0              !  0 earth-fixed (grid axes), 1 background frame""", s, flags=re.M)
        s = re.sub(r'^( Body\s*=\s*1[^\n]*\n)', r"""\1 Probe_point1          = 0.6   0.0  -0.25 !update 22: probe points following the cylinder (0.1 m behind it,
 Probe_point2          = 0.0   0.55 -0.25 !  0.05 m beside it) -> Output_data/<case>_body1_probe<n>.dat
""", s, count=1, flags=re.M)
        s += """
$Probe_background
!---- update 22: probe points following the background (which follows the cylinder here, Moving_BGR = 1)
 Probe_point1          = 2.0   0.0  -0.25 !1.5 m behind the cylinder -> Output_data/<case>_probe1.dat
 Probe_point2          = 0.0   0.0  -0.25 !the cylinder's centre: inside the solid -> "NoData" every step
/end
"""
        open(os.path.join(out, 'probCFD.nml'), 'w').write(s)
    elif case == 'foil':
        src = os.path.join(root, 'cases', 'PID_Steer'); base = 'PID_Steer'
        for f in (base + '.grd', base + '.bcs', 'Input_polyMesh_probCFD.csv'):
            shutil.copy(os.path.join(src, f), os.path.join(out, f.replace(base, name)))
        s = open(os.path.join(src, 'probCFD.nml')).read().replace('"%s"' % base, '"%s"' % name)
        s = re.sub(r'^ End_T .*$', ' End_T                 = 0.01           !5 steps (sandbox test; 0.1 = the update-21 run)', s, flags=re.M)
        s = re.sub(r'^ Output_Intvl .*$', """ Output_Intvl          = 5
 SW_Prs_Frc_surf       = 1              !update 22: wall surface files with every volume solution
 Frame_Prs_Frc_surf    = 0""", s, flags=re.M)
        s = re.sub(r'^( Body\s*=\s*1[^\n]*\n)', r"""\1 Probe_point1          = 1.05  0.1  -0.97 !update 22: probe points turning with the foil (5 cm behind the trailing
 Probe_point3          = 0.25  0.1  -0.90 !  edge; 7 cm above the quarter chord) -> Output_data/<case>_body1_probe<n>.dat
""", s, count=1, flags=re.M)
        s += """
$Probe_background
 Probe_point1          = 1.5   0.1  -0.97 !update 22: fixed point in the wake -> Output_data/<case>_probe1.dat
/end
"""
        open(os.path.join(out, 'probCFD.nml'), 'w').write(s)
    else:
        raise SystemExit('unknown case')
    print(name, 'written')

if __name__ == '__main__':
    main()
