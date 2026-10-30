#!/usr/bin/env python3
"""Natural-convection cavity (cases/NaturalConvection_Cavity): compare the
probCFD result with the de Vahl Davis (1983) benchmark at Ra = 1e5.

usage: python3 tools/cavity_benchmark.py <case folder> [step]
  reads Output_data/<case>_heat.dat (Nu of the hot / cold wall history) and the
  last (or the given) .p3d solution: u_max on the vertical mid-line x = L/2,
  w_max on the horizontal mid-line z = L/2, both in alpha/L units, with their
  positions, and the local Nusselt number along the hot wall.
Benchmark (Ra = 1e5, Pr = 0.71): Nu_avg 4.519, Nu_max 7.717 at z/L = 0.081,
  Nu_min 0.729 at z/L = 1, u_max 34.73 at z/L = 0.855, w_max 68.59 at x/L = 0.066.
"""
import numpy as np, sys, os, glob, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from readq import read_x, read_q

case = sys.argv[1].rstrip('/'); name = os.path.basename(case)
nml = open(os.path.join(case, 'probCFD.nml')).read()
def val(key):
    return float(re.search(r'^\s*%s\s*=\s*([^\s!]+)' % key, nml, re.M).group(1))
L = val('Char_L'); rho = val('Lou_w'); mu = val('mu_w'); cp = val('Cp_l'); k = val('k_l'); beta = val('Beta_l')
nu = mu/rho; alpha = k/(rho*cp); Pr = nu/alpha
bcs = open(os.path.join(case, name + '.bcs')).read()
Tw = [float(v) for v in re.findall(r'^\s*1\s+([0-9.]+)\s*!\s*T:', bcs, re.M)]
Th, Tc = max(Tw), min(Tw); dT = Th - Tc
Ra = 9.81*beta*dT*L**3/(nu*alpha)
print(f'{name}: L {L} m, Pr {Pr:.3f}, Ra {Ra:.4g}, T_hot {Th} K, T_cold {Tc} K')

# ---- Nu history from the heat file -------------------------------------
h = np.loadtxt(os.path.join(case, 'Output_data', name + '_heat.dat'))
grd = read_x(sorted(glob.glob(os.path.join(case, 'Volume_solutions_grids', name + '_*.grd')))[0])[0]
Ly = grd[1].max() - grd[1].min()                      # spanwise width
A = L*Ly
Nu_h = h[:, 5]*L/(k*dT*A); Nu_c = -h[:, 7]*L/(k*dT*A)
print(f'Nu (hot wall) history: t = {h[-1,1]:.1f} s  Nu_hot {Nu_h[-1]:.4f}  Nu_cold {Nu_c[-1]:.4f}   (benchmark 4.519)')
for frac in (0.25, 0.5, 0.75):
    i = int(frac*(len(h) - 1)); print(f'   t = {h[i,1]:7.1f} s  Nu_hot {Nu_h[i]:.4f}')

# ---- profiles from the solution file ------------------------------------
files = sorted(glob.glob(os.path.join(case, 'Volume_solutions_grids', name + '_*.p3d')))
f = files[-1] if len(sys.argv) < 3 else os.path.join(case, 'Volume_solutions_grids', '%s_%08d.p3d' % (name, int(sys.argv[2])))
q = read_q(f)[0]; u, w, T = q[0], q[2], q[4]
x = grd[0][:, 0, 0]; z = grd[2][0, 0, :]; jm = grd.shape[2]//2
im = np.argmin(abs(x - L/2)); km = np.argmin(abs(z - L/2))
uprof = u[im, jm, :]; wprof = w[:, jm, km]
iu = abs(uprof).argmax(); iw = abs(wprof).argmax()
print(f'solution {os.path.basename(f)}:')
print(f'   u_max on x = L/2 : {abs(uprof[iu])*L/alpha:.2f} (alpha/L) at z/L = {z[iu]/L:.3f}   (benchmark 34.73 at 0.855)')
print(f'   w_max on z = L/2 : {abs(wprof[iw])*L/alpha:.2f} (alpha/L) at x/L = {x[iw]/L:.3f}   (benchmark 68.59 at 0.066)')
# local Nu along the hot wall from the node temperatures: -dT/dx at x = 0 by a one-sided 2nd-order difference
Tw0 = T[0, jm, :]; T1 = T[1, jm, :]; T2 = T[2, jm, :]
d1 = x[1] - x[0]; d2 = x[2] - x[0]
dTdx = -(d2**2*(T1 - Tw0) - d1**2*(T2 - Tw0))/(d1*d2*(d2 - d1))   # = -(dT/dx)_wall
Nu_loc = dTdx*L/dT                                   # heat flux into the fluid: -k dT/dx > 0 at the hot wall
inner = slice(1, len(z) - 1)
print(f'   local Nu hot wall: max {Nu_loc[inner].max():.3f} at z/L = {z[inner][Nu_loc[inner].argmax()]/L:.3f} (benchmark 7.717 at 0.081), '
      f'min {Nu_loc[inner].min():.3f} at z/L = {z[inner][Nu_loc[inner].argmin()]/L:.3f} (benchmark 0.729 at 1.0)')
print(f'   mean of the local Nu (trapezoid over z): {np.trapezoid(Nu_loc, z)/L:.4f}')
