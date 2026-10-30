# probCFD_부산 — HANDOVER (read this first in the new room)

Date: 2026-09-06 (update 3: OVERSET grids; the previous handover of
2026-09-04 is kept below unchanged where still valid).  This package is
the complete, validated state of the project.  Everything needed to continue is
inside: solver sources, four runnable cases, the reference files the
user supplied, small analysis tools and the full history.

## 0. Latest state (update 24, 2026-10-19)
Update 24 added the FLUID-STRUCTURE INTERACTION of flexible parts
(`fsi_shell.F90`: co-rotational MITC4 shell, generalized-alpha, banded
Cholesky, stand-alone tests `fsi_test.F90`; `fsi.F90`: the shell
mid-surface built from the part's wall surface grid by ray pairing -
thickness measured from the grid, no structural mesh - load / motion
transfer, mesh morphing of the part's blocks, ALE grid velocity and flux
re-projection, strongly coupled sub-iterations with interface
quasi-Newton relaxation inside the time step, monitor-node files, the
surface files' Disp / vonMises, restart).  `SW_FSI = 1` with
`Youngs_modulus`, `Poisson_ratio`, `Density` / `Mass`, `Num_NoFSI_Group`
+ `NoFSI_*` clamped ranges, `Num_FSI_Submat_regions` + `Submat_*`, and the
optional `FSI_*` keys inside a Body / Child_body / Sub_child_body group
(template).  Section 36; verification `verification/update24_fsi` (CSM3
structural test; the FSI-PfS-2a benchmark in `cases/FSI_PfS2a`, written
by `tools/fsi_cases.py`); run-through case `cases/KP505_FSI` (flexible
blades as child bodies of the hub); 32 source files.  Also fixed: an
index slip in the overset coincident-wall trimming, the position units of
the stop-criterion reports.
Update 24b/24c (2026-10-28, the FSI-PfS-2a verification): the coupling
restarts every sub-iteration from the level-n state (fluxes, Cartesian
field; the old metrics kept through the step; the full pressure kept as the
Poisson initial guess), the IQN works in a scaled space, morphing columns
move rigidly with their wall node, the structural step falls back to 2/4/8
sub-steps; new keys `SW_Quasi_2D` (span averaging + w = 0 on extruded
quasi-2-D grids), `FSI_Plane_Strain` (span-axis displacement / the two other
rotations of the structure held) and `FSI_Test_Motion` (imposed bending
mode, static deflection, or the PLUCK that starts the benchmark: a vibration
mode of the structure imposed for 2.25 periods, then released); slab end
caps in the hole cutting (bodies spanning two parallel planar patches - the
2-D extrusion - are closed by ear-clipped cap triangles; before, thin spans
left active cells inside the body); the single-body elliptic O-grid of the
benchmark with a convex outer boundary (`tools/fsi_cases.py pfs2a`); and the
fix of the old metrics of the flux re-projection (move_grids overwrote the
level-n metrics of the morphing blocks with the rigid rotation of the
undeformed ones, which injected fvw (m_def - m_undef) . u into the fluxes
every step and blew up the benchmark at large deflection).
`PROBCFD_DUMP_GRID=<step>` (environment) writes the grid of a step in double
precision (a debugging aid).  The morph diffuses the column rotation and the
displacement of the outer levels along the wall (rigid near the wall), the
predictor uses discrete velocities, the plane-strain structure is averaged
over its span groups.  FSI-PfS-2a at dt = 1 ms (after the pluck): f = 11.7 Hz
(experiment 11.25, LES 11.73), d_y/D amplitude 0.54 at the paper's probe
(0.648 / 0.629).  PROJECT_RECORD 36.1-36.2.  NOTE: the cloud
workspace of this project was reset on 2026-10-28 and 24b/24c were rebuilt
on top of the update-24 package from the record of the work; the results of
36.2 are those of the rebuilt code.
Update 23 added the on-the-fly turbulence statistics (`statistics.F90`,
`$Statistics` with `SW_statistics`, `Statistics_type1..4`,
`Statistics_first_group / _second_group`, the window / interval keys and
the `Statistics_box / _cyl / _plane` lattice groups): raw moments with
per-point counters (blanked steps skipped, `Coverage`), mean / rms /
Reynolds stresses of every field and the velocity-gradient group (ε_res
from the 12 fluctuating-gradient components, modelled ε and k with a
turbulence model, η, λ_g, Re_λ, Δ/η, vorticity; the gradient from the
face differences of the staggered velocities, accumulated at the cells)
on the original grid
(part frames for moving blocks) and on user-defined boxes, cylinders
(polar components) and planes following the background or a part;
files in `Statistics_volume_solutions/` with unit-carrying `.nam`
files (the volume / surface `.nam` files got units too); accumulators in
the restart file.  Section 35; cases `cases/Statistics_ActuatorDisk`,
`cases/Statistics_TowedCylinder` (`tools/stat_cases.py`), verification
in `verification/update23_statistics` (`check_stats.py`); 29 source files.
Update 22 added the volume-field probe points (`probes.F90`:
Update 22 added the volume-field probe points (`probes.F90`:
`Probe_point<n> = x y z` in the new `$Probe_background` group - points
following the background - or inside a Body / Child_body /
Sub_child_body group - points following the part; one time-history
file per point in `Output_data/`, `NoData` when the point is in a hole
or outside the grid), the `Shear_Stress` variable of the volume files
(after `TTL_Prs`: μ_eff √(2 S_ij S_ij) at the cells, |τ_w| of the wall
faces at the solid-wall nodes) and the solid-wall surface files
(`SW_Prs_Frc_surf`, `Frame_Prs_Frc_surf` in `$Export_Interval`:
`<case>_surf_NNNNNNNN.grd / .p3d`, pressure and friction stress vectors
[Pa] on every wall patch, node values built so that their trapezoidal
integration equals the force files exactly; `surface_output.F90`, the
force integration's face loop factored into the shared
`wall_face_stress`).  Section 34; cases `cases/Probes_TowedCylinder`,
`cases/Probes_SteeredFoil` (written by `tools/probe_cases.py`),
verification in `verification/update22_probes_surface`
(`check_surface.py` integrates a surface file and compares it with the
force file).
Update 21 added the PID controllers (`pid_control.F90`, compiled before
Update 21 added the PID controllers (`pid_control.F90`, compiled before
`sixdof.F90`): `PID_*` blocks inside a Body / Child_body
/ Sub_child_body group or a Force_model group of `$SixDOF` (`PID_use = 1`
opens one), three types - throttle (rev/s), steering (deg), direct force
(a component of a point force) - acting on a grid part (rotated about
`PID_Axis` through its Rot_Center, its Motion_DOF ignored), on a disk
(rps and / or axis angle about `PID_Axis` through `PID_Axis_Center`) or
on a point force; the measured variable is any part's location or
velocity (`PID_Target_body/_child_body/_sub_child_body`, `PID_Target_type`,
`PID_Target_var`), with sequential targets (`PID_Target_value`,
`PID_Target_tol`, `PID_Target_hold_T`), gains `PID_P/I/D`, `PID_Out_init`,
`PID_Out_limit`, `PID_Rate_limit`, `PID_Start_T`; output
`Output_data/<case>_<part>_pid<i>.dat`; the states are in the restart
file.  Section 33; cases `cases/PID_Steer`, `PID_Throttle`, `PID_Direct`,
`PID_Disk`, `PID_DiskSteer` (written by `tools/pid_cases.py`),
verification in `verification/update21_pid` (`check_pid.py` summarises a
run); the variables are annotated as commented-out examples in
`probCFD_template.nml`.  Three fixes came with the tests: the wrenches
of several force models on one mother are added (the last overwrote the
others since update 19), a restart re-evaluates the force models from
the restored state and restores the half-kick flag (restarted free
bodies now reproduce the continuous run exactly).
Update 20 added the moving background (`Moving_BGR`, `Moving_BGR_DoF`,
Update 20 added the moving background (`Moving_BGR`, `Moving_BGR_DoF`,
`BGR_Following_Body`, `Moving_BGR_Filter_T` in `$Case_title`; the
Object_ID-0 blocks follow a Body as their own rigid overset group, earth-
frame far-field values), per-face wave relaxation zones (`Relax_L`) and
the overset re-assembly tolerance (`Ovst_Reassemble_Tol`); section 32,
cases `cases/Cylinder_MovingBGR`, `cases/WaveTank_MovingBGR`,
verification in `verification/update20_moving_background`
(`tools/bgr_cases.py` generates the cases).  Six fixes came with the
validation (32.1 a-f): the restart file carries the background record;
the force integration / wall functions use the ramped Dirichlet wall
velocity of a towing-tank start; a restart applies the boundary values of
the restored time (not of t = 0); the ALE grid velocity is the discrete
displacement velocity of the grid point; Child_body groups may follow
Body_End (the template order now parses completely); the wave probes use
a point-in-cell test (a yawed cell was counted twice).  Open item: a
frame-dependent offset of a few per cent of ½ρU²DL in the pressure force
of a towed body after its start-up, in the moving-grid path since update
12 (verification README, section 1).
Update 19 added the force models (`force_model.F90`): Hough-Ordway and
BEM actuator disks and user-formula point forces as "bodies without a
grid" of the 6-DoF hierarchy (`Force_model` groups of `$SixDOF` with
`Mother_body` / `Mother_child_body` / `Mother_sub_child_body`, `_End = -1`
marks); section 31; cases `cases/ActuatorDisk`, `cases/Foil_Spring`,
verification in `verification/update19_force_models`; 25 source files
(force_model.F90 after sixdof.F90).
Update 18 added the incident-wave model (`$Wave`, `waves.F90`,
`wave_fields.F90`): regular / irregular waves defined in the earth-fixed
frame with the grid-frame shift, the `.bcs` far-field flag 12 (u, v, w,
p, alpha together), relaxation zones, ramp, probes; section 30; case
`cases/RegularWave_Tank`, verification in `verification/update18_waves`
(`tools/wave_cases.py` generates the tank cases).  The build now has 24
source files (waves.F90 before curvi.F90, wave_fields.F90 after
poisson.F90).
Update 17 added the compressible liquid-vapour mixture (`Cav_compressible`,
Wood law, implicit in the pressure equation; verified against the Wood
sound speed and the acoustic impedance) and the thermal cavitation
coupling (`Cav_thermal`: Clausius-Clapeyron p_v(T), rho_v(T), B-factor
delay, implicit latent heat; LN2 vs water); section 29,
`verification/update17_compressible_thermal`.
Update 16 added the heat / energy equation (`$Energy`, `energy.F90`,
temperature at the cells with Dirichlet / flux / convective / adiabatic
BCs through a 6th `.bcs` line, thermal wall function, Boussinesq
buoyancy, latent heat and p_v(T) with cavitation, `Output_data/<case>_heat.dat`;
section 28; case `cases/NaturalConvection_Cavity` = de Vahl Davis Ra 1e5
benchmark, verification in `verification/update16_energy`), the
incremental pressure form `Prs_form = 2` (needed for gravity-driven
flows on wall-clustered grids, section 28.3), and fixed the vapour
ghost cells at Dirichlet-p faces, the turbulence-solver deadlock with
uneven domain counts and the wall distance units (section 28.6).
Update 15 added the cavitation models (`$Cavitation`, `cavitation.F90`,
three-phase water / air / vapour; section 27; cases `cases/NACA66mod_cav`,
`cases/DelftTwist_cav`).  Update 14b added `Mode = 2` (grid-movement check of the impose motion
files into `Check_grid_movement/`) and the missing-file stops (section
26.7).  Update 14 added child / sub-child bodies to the 6-DoF module (hull >
stern drive > propeller hierarchy, relative motion, composite-mass
coupling, `Output_data/` folder, `cases/GPPH` with an annotated name
list): README.md "Update 14" and PROJECT_RECORD.md section 26.  Updates
12-13 (6-DoF rigid bodies, dynamic overset, SSP scheme, runtime name
list, Reason_for_stop.md, restart) are in sections 24-25.  Build with
`make clean && make` (gfortran + OpenMPI, 29 source files).

## 1. What to do first in the new room
1. Unzip, `make clean && make` (gfortran + OpenMPI; the Makefile does
   not track dependencies — always `make clean && make`).
2. Sanity check: `cd regression_BFS && cp ../probCFD . &&
   mpirun -np 2 ./probCFD` and compare `Curvilineargrid_MultiBlockTest_monitor.dat`
   with the reference lines in `regression_BFS/README.md`
   (FGMRES digits are exact for np = 2, otherwise equal to the solver
   tolerance; Jacobi paths are bit-identical for any np).
3. Read `README.md` (what the three 부산 updates added; update 3 =
   overset grids), then `PROJECT_RECORD.md` sections 12-15 for the
   reasoning and the validation numbers; sections 1-11 are the earlier
   rooms.  Overset grid check: set `Ovst_only = 1` in the name list,
   run, read `<case>.ovst` and plot `<case>_iblank.q` (grid-check
   mode needs no field arrays: it works for grids far larger than
   the machine could time-step).

## 2. Package map
| path | content |
|---|---|
| `*.F90`, `Makefile` | solver (15 files, compile order = the SRCS list; `ovs_geom.F90` + `overset.F90` = the overset assembly) |
| `NACA66mod_Ovst_case/` | **update 3**: 3-block OVERSET grid of the foil section (foil O-grid / refinement box / background), README with the assembly numbers |
| `Cylinder_case/` | two-phase surface-piercing cylinder, 4-block O-H grid, fully annotated `probCFD.nml` (use it as the template for every variable) |
| `NACA66mod_case/` | single-phase foil section, O-grid with self-interface cut |
| `PoleCyl_case/` | single-phase pipe, O-type grid with a POLE (axis) face |
| `regression_BFS/` | 10-block backward-facing step, reference monitor lines |
| `reference/` | the user's reference name list and log format, the inherited project record of 2026-09-02, the curvilinear study notes |
| `tools/` | switch-matrix scripts, node single-valuedness check, monitor plot, the force/slice plots of the user's 30 s cylinder run |

## 3. Solver in one paragraph
Fortran 90 + MPI, fully staggered finite differences on multi-block
3-D curvilinear grids (half-index metrics, Park m-matrix, `curvi_form`
1 = unweighted U^q or 2 = J-weighted flux), projection method with
`Iner_Itn` Picard inner iterations per step, momentum by weighted
Jacobi with a true-residual criterion (coefficients assembled once per
solve; implicit convection = linearized 1st-order upwind + frozen
deferred correction of the chosen scheme), convection schemes
`Cnvt_scheme` 1 upwind-1 / 2 upwind-2 / 3 QUICK / 4 ENO-2 — all in
CONSERVATIVE FLUX FORM (face value x face flux), explicit / implicit /
hybrid per direction (`Cnvt_mode`, `Dffs_mode`, `cfl_switch`,
`vn_switch`), Poisson by weighted Jacobi (update-norm criterion) or
FGMRES + 7-point PCG with Jacobi or line-SGS preconditioner (residual
criterion; both norms reported), flux-form operator with exchanged face
fluxes (discrete continuity to solver tolerance on every cell of every
block), two-phase VoF (explicit or implicit CICSAM, sub-cycled,
balanced-force gravity), unified ghost exchange for intra-block halos
and block-to-block interfaces of any orientation incl. partial windows
and self-interfaces (O-cuts), periodic BCs, pole (degenerate) faces,
mass-corrected outflow, ramped Dirichlet BCs, node-averaged
single-valued PLOT3D output, reference-format grouped name list and log.

## 4. Decisions and answers given in this room (so they are not re-asked)
* Momentum/VoF convergence = true algebraic residual, relative to
  ||b||_inf (momentum) or the scaled cell equation (VoF); pressure
  reports update norm AND residual, Jacobi stops on `Trl_Lmt_Prs`,
  FGMRES on `Rsd_Lmt_Prs`.  `Min_Itr_*` are enforced everywhere.
* `SW_Ref_Prs = 1` / `RefPrs`: the all-Neumann Poisson problem is solved
  level-free (compatible RHS) and the field is SHIFTED so the reference
  cell equals `RefPrs`; the value never enters the operator.
* `N_max` and `Initial_Value_fct` were removed on request (update 2);
  `a_cel_initial` = initial alpha of the cells above `FreeSurf_Z`.
* `.bcs` files may have 4 or 5 lines per face for any `Num_of_Phase`
  (a 5th line is auto-detected; ignored with a NOTE for single phase).
* Explicit higher-order convection (upwind-2, QUICK, ENO) is stable
  only below CFL 1/2 per direction: keep `cfl_switch <= 0.3` in hybrid
  mode (the NACA case blew up at the trailing edge with 1.0); QUICK
  with forward-Euler explicit directions is unstable for pure
  convection — use `Cnvt_mode = 2` (+ `Iner_Itn = 2`).  The reader
  prints warnings for these combinations.
* ENO-2 is implemented in flux form with the minmod of the two first
  differences at the donor (conservative, 4-point support); the
  derivative form of the user's formula (minmod of second differences)
  is not conservative on the staggered operator and would need a
  6-point support (PROJECT_RECORD 13.3).
* Pole faces: zero-area boundary (Jh = mh = 0 on the pole plane,
  automatic zero-gradient patch, type-93 record in `.bcs_b2b`); ghost
  values are zero-gradient copies of the same column, not of the
  opposite cell (the PoleCyl grid has 39 azimuthal cells).  An
  interpolated opposite-cell exchange is a possible refinement.
* Junction ghost ambiguity of the O/H grid is resolved at output by
  node averaging over all block copies (verified: zero difference).
* SOR_fct_Vel 0.8 halves the momentum sweeps on the cylinder (91 -> 46)
  but is slower on the NACA case (24 -> 31): case dependent.

* Overset (update 3): overset boundary = uncovered face; wall = .bcs
  patch with Dirichlet u = v = w = 0; bodies closed through surface
  overlap / b2b interfaces / physical boundary patches within
  `Surfovl_tol`; finest grid active, `Ovst_fringe = 2` layers; fringe
  cells implicit in the coupled pressure system (inner BiCGSTAB on the
  coupled operator - the single-grid block-Jacobi PCG stalls at 1e-2);
  iblank appended to the .q files; report `<case>.ovst`.

## 5. Known limitations / natural next steps
* Overset pressure solve costs 3-5x a single-grid solve of the same size
  (coupled inner BiCGSTAB, ~3000-5500 inner iterations per step on the
  NACA66mod_Ovst case): a multigrid / coarse-grid component in the
  preconditioner is the next upgrade.  6-DoF motion (re-assembly every
  step: call ovst_setup again after moving blk(b)%xn/yn/zn and
  rebuilding the metrics) is not implemented.  Two-phase overset (alpha
  interpolated explicitly at fringes) is untested.  Fringe cells that
  interpolate from fringe cells of another block are allowed (counted
  in the report); orphans are reported, not repaired.
* Pressure solve dominates on strongly stretched grids (NACA: ~1000 CG
  iterations, 80 % of the step) -> geometric multigrid preconditioner.
* First-order (backward Euler) time integration; 2nd-order (Crank-
  Nicolson / BDF2) would pair well with `Iner_Itn`.
* Laminar only (Re 1e6 cases are under-resolved); no turbulence model,
  no surface tension, no wall functions, no moving grids.
* Pole ghosts: same-column zero gradient (see above).
* FGMRES results are deterministic per rank count only.
* Run inside scripts with `mpirun --allow-run-as-root --oversubscribe`
  and `< /dev/null`; never `pkill -f probCFD` (kills the shell) — use
  `pkill -x probCFD`.

## 6. Test status at packaging (dev VM, 2 cores)
* BFS: 5 steps x {Cnvt_scheme 1..4} x {Cnvt_mode 1,2,3} all run,
  residuals 1e-9, itp 755-766 (line SGS).
* Two-phase cylinder: 10 steps with the shipped name list (|div|max
  ~2e-7 1/s, liquid volume constant, itp 250-410); QUICK/ENO 1 step
  implicit and hybrid; QUICK implicit with Iner_Itn 2 (itp 443 -> 234).
* NACA66mod: 80 steps stable (cfl_switch 0.3), 8 s/step on 4 ranks.
* PoleCyl: 200 steps, plug core 1.0007 m/s, axis single-valued,
  |div|max 1e-8 1/s, 2.8 s/step.
