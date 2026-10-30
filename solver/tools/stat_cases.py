#!/usr/bin/env python3
"""Test cases of the on-the-fly statistics (update 23).

usage: python3 tools/stat_cases.py <case> <outdir>
  disk     : cases/ActuatorDisk (static Cartesian box, a Hough-Ordway disk starting in a uniform stream)
             run for 10 steps with the solution written EVERY step (Output_Intvl = 1) and the statistics
             of steps 4-10 (Statistics_Start_T = 0.02) on the grid (type 1), a box, a cylinder (polar
             velocities) and a plane: check_stats.py compares the type-1 mean / rms / Reynolds stresses
             with the per-step volume files (exact to the file precision) and recomputes eps_res
  cylinder : cases/Probes_TowedCylinder (towed cylinder, moving background, body- and background-
             following probe points, 20 steps) with the statistics of steps 4-20: the grid blocks in
             the part frames, a box following the cylinder that holds the probe point behind it, a
             background box holding the background probe point, a background box across the cylinder
             (coverage 0 inside the solid) and a ring around the cylinder with polar velocities:
             check_stats.py compares the lattice points with the probe time histories (exact)
"""
import sys, os, shutil, re
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')

STAT_DISK = """
$Statistics
!---- update 23: statistics of steps 4-10 (t = 0.02 .. 0.05 s) on the grid, a box, a ring and a plane
 SW_statistics           = 1
 Statistics_type1        = 1
 Statistics_type2        = 1
 Statistics_type3        = 1
 Statistics_type4        = 1
 Statistics_first_group  = 1
 Statistics_second_group = 1
 Statistics_Start_T      = 0.02
 Statistics_End_T        = 0
 Statistics_Sample_Intvl = 1
 Statistics_Output_Intvl = 5
 Statistics_Type1_Frame  = 1
 Statistics_box          = 1              !box around the disk (x = 0.9 .. 1.5), earth-fixed
 Mother_body             = 0
 Statistics_box_Origin   = 0.9 -0.2 -0.2
 Statistics_box_Axis1    = 1 0 0
 Statistics_box_Axis2    = 0 1 0
 Statistics_box_Length   = 0.6 0.4 0.4
 Statistics_box_nPoints  = 13 9 9
 Statistics_box_End      = -1
 Statistics_cyl          = 1              !the slipstream: ring/disc lattice on the disk axis, polar velocities
 Statistics_cyl_Facecenter = 0.85 0.0 0.0
 Statistics_cyl_Axis     = 1 0 0
 Statistics_cyl_Ref      = 0 0 1
 Statistics_cyl_Radius   = 0.0 0.15
 Statistics_cyl_Theta    = 0 360
 Statistics_cyl_Length   = 0.5
 Statistics_cyl_nPoints  = 7 25 11
 Statistics_cyl_Polar_Vel = 1
 Statistics_cyl_End      = -1
 Statistics_plane        = 2              !a cross plane 0.4 m behind the disk
 Statistics_plane_Origin = 1.2 -0.3 -0.3
 Statistics_plane_Axis1  = 0 1 0
 Statistics_plane_Axis2  = 0 0 1
 Statistics_plane_Length = 0.6 0.6
 Statistics_plane_nPoints = 25 25
 Statistics_plane_End    = -1
/end
"""
STAT_CYL = """
$Statistics
!---- update 23: statistics of steps 4-20 (t = 0.02 .. 0.1 s)
 SW_statistics           = 1
 Statistics_type1        = 1              !the grid blocks in the part frames (the O-grid follows the cylinder)
 Statistics_type2        = 1
 Statistics_type3        = 1
 Statistics_type4        = 0
 Statistics_first_group  = 1
 Statistics_second_group = 1
 Statistics_Start_T      = 0.02
 Statistics_Output_Intvl = 10
 Statistics_Type1_Frame  = 1
 Statistics_box          = 1              !follows the cylinder: its point (2,2,2) = the probe point (0.6 0 -0.25)
 Mother_body             = 1
 Statistics_box_Origin   = 0.5 -0.1 -0.35
 Statistics_box_Axis1    = 1 0 0
 Statistics_box_Axis2    = 0 1 0
 Statistics_box_Length   = 0.4 0.2 0.2
 Statistics_box_nPoints  = 5 3 3
 Statistics_box_End      = -1
 Statistics_box          = 2              !follows the background: its point (2,2,2) = the background probe point (2 0 -0.25)
 Mother_body             = 0
 Statistics_box_Origin   = 1.8 -0.2 -0.45
 Statistics_box_Axis1    = 1 0 0
 Statistics_box_Axis2    = 0 1 0
 Statistics_box_Length   = 0.4 0.4 0.4
 Statistics_box_nPoints  = 3 3 3
 Statistics_box_End      = -1
 Statistics_box          = 3              !background plane across the cylinder: the points inside the solid get Coverage 0
 Mother_body             = 0
 Statistics_box_Origin   = -0.8 -0.8 -0.25
 Statistics_box_Axis1    = 1 0 0
 Statistics_box_Axis2    = 0 1 0
 Statistics_box_Length   = 1.6 1.6 0.0
 Statistics_box_nPoints  = 17 17 1
 Statistics_box_End      = -1
 Statistics_cyl          = 1              !ring around the cylinder, following it, polar velocities (u_tan = the flow around it)
 Mother_body             = 1
 Statistics_cyl_Facecenter = 0.0 0.0 -0.4
 Statistics_cyl_Axis     = 0 0 1
 Statistics_cyl_Ref      = 1 0 0
 Statistics_cyl_Radius   = 0.55 0.8
 Statistics_cyl_Theta    = 0 360
 Statistics_cyl_Length   = 0.3
 Statistics_cyl_nPoints  = 3 37 4
 Statistics_cyl_Polar_Vel = 1
 Statistics_cyl_End      = -1
/end
"""

def main():
    case, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    name = os.path.basename(out.rstrip('/'))
    if case == 'disk':
        src = os.path.join(root, 'cases', 'ActuatorDisk'); base = 'ActuatorDisk'
        for f in (base + '.grd', base + '.bcs', 'Input_polyMesh_probCFD.csv'):
            shutil.copy(os.path.join(src, f), os.path.join(out, f.replace(base, name)))
        s = open(os.path.join(src, 'probCFD.nml')).read().replace('"%s"' % base, '"%s"' % name)
        s = re.sub(r'^ End_T .*$', ' End_T                 = 0.05           !10 steps (the statistics check)', s, flags=re.M)
        s = re.sub(r'^ Output_Intvl .*$', ' Output_Intvl          = 1              !every step: the reference of check_stats.py', s, flags=re.M)
        s = re.sub(r'^ Log_Intvl .*$', ' Log_Intvl             = 5', s, flags=re.M)
        s += STAT_DISK
    elif case == 'cylinder':
        src = os.path.join(root, 'cases', 'Probes_TowedCylinder'); base = 'Probes_TowedCylinder'
        for f in (base + '.grd', base + '.bcs', 'Input_polyMesh_probCFD.csv', base + '_body1_imposed.dat'):
            shutil.copy(os.path.join(src, f), os.path.join(out, f.replace(base, name)))
        s = open(os.path.join(src, 'probCFD.nml')).read().replace('"%s"' % base, '"%s"' % name)
        s += STAT_CYL
    else:
        raise SystemExit('unknown case')
    open(os.path.join(out, 'probCFD.nml'), 'w').write(s)
    print(name, 'written')

if __name__ == '__main__':
    main()
