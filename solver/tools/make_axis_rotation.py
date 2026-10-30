#!/usr/bin/env python3
"""Imposed-motion table of a constant-speed ROTATION ABOUT AN ARBITRARY AXIS
(a propeller on an inclined shaft, a rudder on a raked stock, ...) for the
probCFD 6-DoF module (update 14).

The imposed-motion file of a body / child / sub-child holds Euler angles
(Z-Y-X: R = Rz(psi) Ry(theta) Rx(phi)) of the attitude relative to the
mother frame (= the SNAME earth frame E for a body, or the mother's frame
for a child / sub-child; all frames are aligned with E at t = 0: x forward
= -x grid, y starboard = y grid, z down = -z grid).  A rotation by the
angle a(t) = 2 pi n t about a unit axis d through the centre of rotation is
  R(t) = I + sin(a) [d]x + (1 - cos(a)) [d]x^2        (Rodrigues)
and this script writes the Euler angles of R(t) (roll unwrapped, so that
the cubic Hermite interpolation of the solver never sees a jump) with a
fine time sampling.  The translation columns are zero.

usage: make_axis_rotation.py <out.dat> <n [rev/s]> <dx dy dz (GRID frame)>
                            <t_end [s]> [--dt 2.5e-4] [--t0 0.0] [--ramp 0.0]
                            [--stop_deg A]
  n > 0 : positive rotation about the axis by the right-hand rule, the axis
          given in GRID coordinates (it is converted to E internally); for a
          propeller shaft pointing forward (-x grid) a positive n is a
          right-handed propeller going ahead (clockwise seen from astern)
  --t0  : the rotation starts at t0 (the angle is 0 before)
  --ramp: the rate ramps linearly from 0 to n over this time after t0
  --stop_deg A : the rotation stops at the angle A [deg] and the part stays
          there (a rudder / stern drive turned to a fixed angle)
examples (GPPH): propeller 10 rev/s about the shaft axis tilted 2.27 deg,
  python3 tools/make_axis_rotation.py GPPH_body1_child1_sub1_imposed.dat 10 \
          -0.99921 0 0.03964 2.1
stern drive turning at 45 deg/s (0.125 rev/s) about its raked steering axis
up to +45 deg, then held:
  python3 tools/make_axis_rotation.py GPPH_body1_child1_imposed.dat 0.125 \
          -0.0998 0 0.9950 2.1 --stop_deg 45
"""
import sys, argparse, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('out'); ap.add_argument('n', type=float)
ap.add_argument('dx', type=float); ap.add_argument('dy', type=float); ap.add_argument('dz', type=float)
ap.add_argument('t_end', type=float)
ap.add_argument('--dt', type=float, default=2.5e-4)
ap.add_argument('--t0', type=float, default=0.0)
ap.add_argument('--ramp', type=float, default=0.0)
ap.add_argument('--stop_deg', type=float, default=None)
a = ap.parse_args()

Fm = np.diag([-1.0, 1.0, -1.0])                     # grid -> SNAME earth frame
d = Fm @ np.array([a.dx, a.dy, a.dz]); d /= np.linalg.norm(d)
K = np.array([[0, -d[2], d[1]], [d[2], 0, -d[0]], [-d[1], d[0], 0]])

def angle(t):                                        # rotation angle [rad] at time t
    tt = max(t - a.t0, 0.0)
    if a.ramp > 0 and tt < a.ramp:
        al = 2 * np.pi * a.n * 0.5 * tt * tt / a.ramp
    else:
        al = 2 * np.pi * a.n * (tt - 0.5 * a.ramp)
    if a.stop_deg is not None:                       # stop and hold at the given angle
        al = min(al, np.radians(abs(a.stop_deg))) if a.n >= 0 else max(al, -np.radians(abs(a.stop_deg)))
    return al

def euler(R):                                        # Z-Y-X Euler angles of R
    th = -np.arcsin(max(-1.0, min(1.0, R[2, 0])))
    ph = np.arctan2(R[2, 1], R[2, 2]); ps = np.arctan2(R[1, 0], R[0, 0])
    return ph, th, ps

t = np.arange(0.0, a.t_end + 0.5 * a.dt, a.dt)
rows = []; prev = None
for ti in t:
    al = angle(ti)
    R = np.eye(3) + np.sin(al) * K + (1 - np.cos(al)) * (K @ K)
    ph, th, ps = euler(R)
    e = np.array([ph, th, ps])
    if prev is not None:                             # unwrap every angle (continuity)
        e = prev + (e - prev + np.pi) % (2 * np.pi) - np.pi
    prev = e; rows.append((ti, *np.degrees(e)))
with open(a.out, 'w') as f:
    f.write(f'# probCFD imposed motion: rotation about the axis d = ({a.dx} {a.dy} {a.dz}) (grid frame) through the\n')
    f.write(f'#   centre of rotation at n = {a.n} rev/s (start t0 = {a.t0} s, ramp {a.ramp} s'
            + (f', stops at {a.stop_deg} deg' if a.stop_deg is not None else '') + '); written by tools/make_axis_rotation.py\n')
    f.write('# columns: t [s]  x y z [m] (translation, zero)  phi theta psi [deg] (Z-Y-X Euler angles of R(t), roll unwrapped)\n')
    for ti, ph, th, ps in rows:
        f.write(f'{ti:.6f}  0 0 0  {ph:.8f} {th:.8f} {ps:.8f}\n')
print(f'wrote {a.out}: {len(rows)} rows, axis in E = ({d[0]:.5f} {d[1]:.5f} {d[2]:.5f}), {a.n*360*a.dt:.3f} deg per row')
