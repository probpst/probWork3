#!/usr/bin/env python3
"""Test cases for the moving background (update 20).

usage: python3 tools/bgr_cases.py <case> <outdir>
  cyl_static : Cylinder_Overset_1phase, fluid at rest, inlet ramped 0 -> 1.5 m/s over t_ramp =
               0.1 s (u_cel_initial = 0: the whole channel accelerates, no impulsive start, no
               flow through the wall) - the reference
  cyl_surge  : the same cylinder, fluid at rest, inlet ramped 0 -> 1.0 m/s, body 1 towed
               upstream with the imposed surge velocity 0.5 min(t/0.1, 1) m/s and the background
               following the surge (Moving_BGR = 1, DoF 1 0 0 0 0 1): the velocity of the body
               relative to the stream is the same -1.5 min(t/0.1, 1) m/s as in cyl_static.
               The two frames differ by a uniform acceleration during the ramp, so the flow
               fields are the same and the forces differ only by the Froude-Krylov term
               rho V (a_static - a_towed) = 1000 * 2.749 * (15 - 10) = 13.7 kN while t < 0.1 s;
               after the ramp the forces must match (Galilean invariance of the moving
               background).  (A towing-tank start, u_cel_initial = stream speed, is NOT
               frame-equivalent: its wall is permeable during the ramp with the flux
               u_init (1 - ramp), 1.5 against 1.0 in the two frames.)
  cyl_surge_tol : cyl_surge with Ovst_Reassemble_Tol = 0.02 m (the body does not move relative
               to the background: the assembly of step 1 is re-used in every later step)
  cyl_surge_nobg : cyl_surge with Moving_BGR = 0 (the body moves through the static background:
               the update-19 moving-body path, re-assembled every step) - must match cyl_surge
  wave_surge : the RegularWave_Tank with a dummy overset block (body 1, no walls) imposed to
               surge at 0.15 m/s (1.2 m in 10 periods) and the tank (background) following it:
               the earth-fixed probes (x = 2, 2.5, 3 m) must still see the Airy wave
  wave_yaw   : the same, tank 0.1 m wide, with an imposed yaw of the dummy body (2 deg over 2 s,
               then held): the tank yaws about the body; the probes must still see the Airy wave
"""
import numpy as np, sys, os, shutil, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from energy_cases import face
import wave_cases

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')

def write_multiblock_grid(path, blocks):
    with open(path, 'w') as f:
        f.write("%8d\n" % len(blocks))
        for (x, y, z) in blocks:
            f.write("%8d%8d%8d\n" % (len(x), len(y), len(z)))
        for (x, y, z) in blocks:
            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
            for a in (X, Y, Z):
                v = a.transpose(2, 1, 0).ravel(); n4 = (len(v)//4)*4
                if n4: np.savetxt(f, v[:n4].reshape(-1, 4), fmt='%22.15e')
                if n4 < len(v): np.savetxt(f, v[n4:][None, :], fmt='%22.15e')

def cyl(case, out):
    name = os.path.basename(out.rstrip('/'))
    src = os.path.join(root, 'Cylinder_Overset_1phase')
    shutil.copy(os.path.join(src, 'Cylinder_Overset_1phase.grd'), os.path.join(out, name + '.grd'))
    bcs = open(os.path.join(src, 'Cylinder_Overset_1phase.bcs')).read()
    s = open(os.path.join(src, 'probCFD.nml')).read()
    s = s.replace('"Cylinder_Overset_1phase"', '"%s"' % name)
    s = re.sub(r'^ End_T .*$', ' End_T                 = 0.3            !60 steps: 20 of ramp + 40', s, flags=re.M)
    s = re.sub(r'^ Log_Intvl .*$', ' Log_Intvl             = 10', s, flags=re.M)
    s = re.sub(r'^ Output_Intvl .*$', ' Output_Intvl          = 60', s, flags=re.M)
    s = re.sub(r'^ t_ramp .*$', ' t_ramp                = 0.1            !inlet ramped from rest over 0.1 s (same ramp as the imposed surge)', s, flags=re.M)
    s = re.sub(r'^ Data_Intvl .*$', ' Data_Intvl            = 1', s, flags=re.M)
    csv = "probCFD polyMesh Input,block-1,block-2,block-3\ni_split,1,1,2\nj_split,1,1,1\nk_split,2,2,2\nn_of_BLK,2,2,4\nObject_ID,%s\n"
    if case == 'cyl_static':
        bcs = bcs.replace('       1   1.000       ! u: Dirichlet 1 m/s (free stream along +x)', '       1   1.500       ! u: Dirichlet 1.5 m/s (free stream along +x)')
        s = re.sub(r'^ u_cel_initial .*$', ' u_cel_initial         = 0.0            !fluid at rest (ramped inflow)', s, flags=re.M)
        csv = csv % '1,1,0'
        six = """$SixDOF
 Body                  = 1              !the cylinder (static reference: all DOFs fixed)
 Mass                  = 1000.0
 CG                    = 0.0  0.0  -0.25
 Rot_Center            = 0.0  0.0  -0.25
 Motion_DOF            = 0 0 0 0 0 0
 Inertia               = 100 100 100
 Body_End              = -1
/end
"""
    else:
        s = re.sub(r'^ u_cel_initial .*$', ' u_cel_initial         = 0.0            !fluid at rest (ramped inflow)', s, flags=re.M)
        csv = csv % '1,1,0'
        tol = 0.02 if case == 'cyl_surge_tol' else 0.0
        mbg = 0 if case == 'cyl_surge_nobg' else 1
        s = s.replace('  Mode                 = 0', '  Mode                 = 0\n  Moving_BGR           = %d              !the background follows body 1 (update 20)\n'
                      '  Moving_BGR_DoF       = 1 0 0 0 0 1    !surge and yaw permitted\n  BGR_Following_Body   = 1\n  Moving_BGR_Filter_T  = 0.0' % mbg)
        s = s.replace(' Ovst_only             = 0', ' Ovst_Reassemble_Tol   = %.3f          !re-assembly tolerance [m] (update 20)\n Ovst_only             = 0' % tol)
        six = """$SixDOF
 Body                  = 1              !the cylinder, towed upstream: surge velocity 0.5 min(t/0.1, 1) m/s (imposed)
 Mass                  = 1000.0
 CG                    = 0.0  0.0  -0.25
 Rot_Center            = 0.0  0.0  -0.25
 Motion_DOF            = 1 0 0 0 0 0    !surge imposed from <case>_body1_imposed.dat (forward = -x grid)
 Inertia               = 100 100 100
 Body_End              = -1
/end
"""
        t = np.linspace(0, 1.0, 1001); tr = 0.1
        xs = np.where(t < tr, 0.5*t**2/(2*tr), 0.5*(t - tr/2))      # velocity 0.5 min(t/tr, 1): the same ramp as t_ramp
        np.savetxt(os.path.join(out, name + '_body1_imposed.dat'), np.column_stack([t, xs, 0*t, 0*t, 0*t, 0*t, 0*t]),
                   fmt='%12.8f', header='t x y z [m] phi theta psi [deg] (SNAME, relative to the start): surge velocity 0.5 min(t/0.1, 1)')
    i = s.index('$SixDOF') if '$SixDOF' in s else len(s)
    s = s[:i] + six if '$SixDOF' in s else s + six
    open(os.path.join(out, 'probCFD.nml'), 'w').write(s)
    open(os.path.join(out, name + '.bcs'), 'w').write(bcs.replace("Cylinder_Overset_1phase", name))
    open(os.path.join(out, 'Input_polyMesh_probCFD.csv'), 'w').write(csv)
    print(name, 'written')

def wave(case, out):
    name = os.path.basename(out.rstrip('/'))
    # the tank of wave_cases.regular (block 1) + a dummy overset block (block 2) inside the water;
    # wave_yaw: the tank is 0.1 m wide (4 cells) so that the earth-fixed probes stay inside the
    # yawed tank (2 deg about x = 2.5 m moves a probe 1 m away by 35 mm)
    lam = 1.0; d = 1.0; Lx = 5.0; nxl = 32
    W, ny = (0.05, 2) if case == 'wave_surge' else (0.10, 4)
    nx = int(round(Lx/lam*nxl))
    x = np.linspace(0, Lx, nx + 1); y = np.linspace(0, W, ny + 1)
    zb = wave_cases.zone_nodes(-d, -0.03, 24, 2.2, True); zm = np.linspace(-0.03, 0.03, 31); zt = wave_cases.zone_nodes(0.03, 0.25, 10, 2.0, False)
    z = np.concatenate([zb[:-1], zm, zt[1:]])
    x2 = np.linspace(2.2, 2.8, 13); y2 = np.linspace(0, W, ny + 1); z2 = np.linspace(-0.45, -0.15, 13)
    write_multiblock_grid(os.path.join(out, name + '.grd'), [(x, y, z), (x2, y2, z2)])
    nz = len(z) - 1
    # bcs of the tank (block 1) as wave_cases.regular; block 2 has no listed surface (all overset boundaries)
    wave_bc = ["12  0.0         ! u: incident wave (earth-frame value)", "12  0.0         ! v", "12  0.0         ! w",
               "12  0.0         ! p", "12  0.0         ! alpha"]
    side = ["0   0.0         ! u: symmetry", "1   0.0         ! v: symmetry (impermeable)", "0   0.0         ! w",
            "0   0.0         ! p: zero gradient", "0   0.0         ! alpha: zero gradient"]
    bottom = ["0   0.0         ! u: free slip", "0   0.0         ! v", "1   0.0         ! w: impermeable bottom",
              "0   0.0         ! p: zero gradient", "0   0.0         ! alpha: zero gradient"]
    top = ["0   0.0         ! u: open top (air)", "0   0.0         ! v", "0   0.0         ! w",
           "1   0.0         ! p: Dirichlet 0 Pa", "0   0.0         ! alpha: zero gradient"]
    bcs = "! %s: wave tank (block 1, background) + dummy overset block 2 (body 1, no walls) - moving background test (update 20)\n" % name
    bcs += "       2                ! number of grid blocks\n       6                ! number of boundary surfaces of block 1\n"
    bcs += face(1, (1, 1), (1, ny+1), (1, nz+1), wave_bc, "x-min: incident wave")
    bcs += face(-1, (nx+1, nx+1), (1, ny+1), (1, nz+1), wave_bc, "x-max: incident wave")
    bcs += face(2, (1, nx+1), (1, 1), (1, nz+1), side, "y-min (symmetry)")
    bcs += face(-2, (1, nx+1), (ny+1, ny+1), (1, nz+1), side, "y-max (symmetry)")
    bcs += face(3, (1, nx+1), (1, ny+1), (1, 1), bottom, "z-min: bottom")
    bcs += face(-3, (1, nx+1), (1, ny+1), (nz+1, nz+1), top, "z-max: open top")
    bcs += "       0                ! number of boundary surfaces of block 2 (all faces overset boundaries)\n"
    open(os.path.join(out, name + '.bcs'), 'w').write(bcs)
    g = 9.81; k = 2*np.pi/lam; om = np.sqrt(g*k*np.tanh(k*d)); T = 2*np.pi/om
    probes = [2.0, 2.5, 3.0]          # earth-fixed; inside the free region of the tank over the whole run
    wb = wave_cases.wave_regular(0.01, T, 2*T, 1.0, 1.5, probes)
    # the six-face zone lengths of update 20 (x-min 1.0 m, x-max 1.5 m) instead of the update-18 aliases
    wb = re.sub(r' Relax_In_L .*\n Relax_Out_L .*\n Relax_Side_L .*\n',
                ' Relax_L               = 1.0 1.5 0 0 0 0   !relaxation zones from the faces x-min x-max y-min y-max z-min z-max [m]\n'
                '                                        !  (attached to the grid faces: they move with the tank; the probes are earth-fixed)\n', wb)
    wb = wb.replace(' Wave_Probe_Y          = ' + ', '.join('0.025' for p in probes), ' Wave_Probe_Y          = ' + ', '.join('%.3f' % (W/2) for p in probes))
    s = wave_cases.nml(name, 0.004, "%.4f" % (10*T), 200, 0.0, wb)
    s = s.replace("  Mode                 = 0", "  Mode                 = 0\n  Moving_BGR           = 1              !the tank follows body 1 (update 20)\n"
                  "  Moving_BGR_DoF       = %s    !permitted DOFs\n  BGR_Following_Body   = 1\n  Moving_BGR_Filter_T  = 0.0"
                  % ("1 0 0 0 0 0" if case == 'wave_surge' else "0 0 0 0 0 1"))
    s = s.replace(" SW_Ovst               = 0", " SW_Ovst               = -1           !block 2 is an overset block")
    t = np.linspace(0, 10, 1001)
    if case == 'wave_surge':
        motion = "surge: x_SNAME = -0.15 t (the tank translates in +x grid at 0.15 m/s through the earth-fixed wave: 1.2 m in 8 s)"
        tab = np.column_stack([t, -0.15*t, 0*t, 0*t, 0*t, 0*t, 0*t])
    else:
        motion = "yaw: psi = 2 deg x min(t/2, 1) (the tank yaws about the body centre)"
        tab = np.column_stack([t, 0*t, 0*t, 0*t, 0*t, 0*t, 2.0*np.minimum(t/2.0, 1.0)])
    np.savetxt(os.path.join(out, name + '_body1_imposed.dat'), tab, fmt='%12.6f', header='t x y z [m] phi theta psi [deg]')
    s += """
$SixDOF
!---- body 1: dummy overset block 2 (no wall: it carries no fluid force), imposed %s
 Body                  = 1
 Mass                  = 1.0
 CG                    = 2.5  %.3f  -0.30
 Rot_Center            = 2.5  %.3f  -0.30
 Motion_DOF            = %s
 Inertia               = 1 1 1
 Body_End              = -1
/end
""" % (motion, W/2, W/2, "1 0 0 0 0 0" if case == 'wave_surge' else "0 0 0 0 0 1")
    open(os.path.join(out, 'probCFD.nml'), 'w').write(s)
    open(os.path.join(out, 'Input_polyMesh_probCFD.csv'), 'w').write(
        "probCFD polyMesh Input,block-1,block-2\ni_split,2,1\nj_split,1,1\nk_split,1,1\nn_of_BLK,2,1\nObject_ID,0,1\n")
    print(name, 'written: tank', nx, 'x', ny, 'x', nz, 'cells + dummy block 12 x', ny, 'x 12')

def main():
    case, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    if case.startswith('cyl'):
        cyl(case, out)
    elif case.startswith('wave'):
        wave(case, out)
    else:
        raise SystemExit('unknown case')

if __name__ == '__main__':
    main()
