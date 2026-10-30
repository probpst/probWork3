# FSI_PfS2a — cylinder with a flexible rubber plate and a steel rear mass (update 24 / 24c)

The FSI-PfS-2a benchmark of Kalmbach & Breuer (J. Fluids Struct. 42, 2013) as
computed by Fabbri et al. (Comput. Fluids 265, 2023, 105963): a fixed
cylinder D = 22 mm with a para-rubber plate 50 x 2 mm (E = 4.1 MPa,
nu = 0.48, rho = 1090 kg/m3) and a steel rear mass 10 x 2 mm
(E = 2.1e11 Pa, nu = 0.3, rho = 7850) bonded to its end, in a water channel
0.338 x 0.240 m (slip walls) at u = 1.385 m/s (Re_D = 30 470).  The
experiment finds the second swiveling mode at f = 11.25 Hz with an amplitude
d_y/D = 0.648 (Fabbri et al.: 11.73 Hz, 0.629 on 40 M cells, LES).

Written by `python3 tools/fsi_cases.py pfs2a cases/FSI_PfS2a coarse=1 dt=5e-4
end=1.0 start=0.15` (24c generator; `one=0` gives the two-body grid of update 24,
`coarse=0` a finer grid).

* Grid (quasi-2-D, span 0.4 mm, 2 cells, slip span walls): ONE elliptic O-grid
  around cylinder + plate (201 x 21 x 3 nodes, first cell 0.08 mm, the outer
  boundary the convex hull of circles 12 mm off the body, the seam at the
  cylinder front), a 1.5-mm refinement box (107 x 61 x 3) and the channel
  background (98 x 70 x 3): 34 106 cells, overset.
* `$SixDOF`: Body 1 = cylinder + plate, fixed in space, `SW_FSI = 1`
  (`Youngs_modulus 4.1e6`, `Poisson_ratio 0.48`, `Density 1090`), the cylinder
  surface clamped (two `NoFSI` groups), the steel end a `Submat` region,
  `FSI_Start_T = 0.15`, `FSI_Max_Iter 6`, `FSI_Tol 5e-4`, `FSI_Relax 0.4`,
  `FSI_Plane_Strain = 3`, `FSI_Test_Motion = 0.011, 11.25, 2, 0` (the PLUCK:
  the static tip-load shape, 11 mm at the tip, imposed at 11.25 Hz for 2.25
  periods, then released and solved - without it the laminar quasi-2-D wake
  does not start the oscillation within an affordable time); monitor node 1 =
  the upper wall node 2 mm before the end of the rubber (x = 0.059 m, the
  paper's probe), monitor node 2 = the trailing end of the rear mass.
* Laminar Navier-Stokes with `SW_Quasi_2D = 3`; restart file every 20 steps.

## How the sandbox run was made (`run_sandbox.sh`)

Phase 1: `Del_T = 5e-4`, `End_T = 0.15` (300 steps, rigid, the flow develops).
Phase 2: `Mode = 1`, `Del_T = 1e-3`, `End_T = 1.0` - the solver counts the steps
as End_T / Del_T, so this continues to step 1000 = t 0.85 s (pluck 0.15 - 0.35
s, free from 0.35 s).  `run_sandbox.sh <dir> 2` makes both phases (and resumes
after an interruption).  About 15 s per rigid step and 80 - 120 s per free
step on 2 cores.

## Result (`sandbox_results/`, `../../verification/update24_fsi/fsi_pfs2a_result.png`)

Over t = 0.60 - 0.85 s: f = 11.7 Hz (experiment 11.25, LES 11.73); d_y/D
amplitude 0.54 (the six extremes 0.46 - 0.65) at node 1 (experiment 0.648, LES 0.629), 0.87 (0.78 - 1.01) at
the trailing end, which swivels 5 ms behind the rubber end.  Details:
`verification/update24_fsi/README.md`, section 2.

Evaluate with `python3 ../../verification/update24_fsi/plot_fsi.py tip .`
(histories, spectrum, mean / amplitude / frequency of the last periods),
`... grid . <step>` (the morphed overset grid) and `... shell .` (the mid-surface).
