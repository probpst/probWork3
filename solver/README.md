# probCFD_부산 — multi-block curvilinear incompressible NS solver, single/two-phase

**New room? Start with `HANDOVER.md`.**

This package is the complete solver source plus three runnable cases:

| folder / file | content |
|---|---|
| `*.F90`, `Makefile` | full solver source (32 files, `fsi_shell.F90`, `fsi.F90`, `fsi_test.F90` added in update 24); build with `make clean && make` (gfortran + OpenMPI); `make fsi_test` builds the stand-alone structural tests |
| `Cylinder_1phase/`, `Cylinder_2phase/` | **update 4 tutorial**: surface-piercing cylinder D = 1 m on the 4-block single grid (O-grid halves + two outer blocks), water only / water-air with the free surface at z = 0 |
| `Cylinder_Overset_1phase/`, `Cylinder_Overset_2phase/` | the same cylinder on the 3-block OVERSET grid (O-grid halves inside a Cartesian background), water only / water-air |
| `NACA66mod_1phase/` | NACA66(mod) section on the 2-block C-grid (C-grid + wake slit), water |
| `NACA66mod_Ogrid_1phase/` | the same section on the 1-block O-grid, water |
| `NACA66mod_Overset_1phase/` | the same section on the 3-block OVERSET grid (foil O-grid + refinement box + background), water |

Run a case from its folder with the solver copied in:

    make                                  # in the package root
    cp probCFD NACA66mod_case/ && cd NACA66mod_case
    mpirun -np 8 ./probCFD                # np <= the CPU total of Input_polyMesh_probCFD.csv

## Update 24 (this delivery, 2026-10-19): fluid-structure interaction of flexible parts

Details: PROJECT_RECORD.md section 36; verification in
`verification/update24_fsi` (cases written by `tools/fsi_cases.py`, plots by
`plot_fsi.py`, structural tests by `fsi_test.F90`); all inputs annotated in
`probCFD_template.nml` (the FSI note before the `$SixDOF` bodies and the key
block in the propeller sub-child body).  New source files `fsi_shell.F90`,
`fsi.F90`, `fsi_test.F90` (32 files).

* **`SW_FSI = 1`** in a Body / Child_body / Sub_child_body group makes the part
  FLEXIBLE.  No structural grid: the part's closed wall surface (its no-slip
  patches) becomes a shell MID-SURFACE whose thickness is measured from the
  surface grid itself (a ray from every upper-side node to the opposite side).
  Inputs: `Youngs_modulus`, `Poisson_ratio`, `Density` (or `Mass`; the other
  follows from the enclosed volume), `Num_NoFSI_Group` + `NoFSI_block_num<n> /
  NoFSI_i_range<n> / NoFSI_j_range<n> / NoFSI_k_range<n>` (clamped wall-node
  ranges), `Num_FSI_Submat_regions` + `Submat_block_num<n> / Submat_i/j/k_range<n>
  / Submat_Density<n> / Submat_Youngs_modulus<n> / Submat_Poisson_ratio<n>`
  (a second material, e.g. a steel mass on a rubber plate); optional
  `FSI_Start_T`, `Damping_ratio`, `FSI_Min_Thickness`, `FSI_Max_Iter`,
  `FSI_Tol`, `FSI_Relax`, `FSI_Mode`, `FSI_Monitor_node<n> = block i j k`;
  24c: `FSI_Plane_Strain` (span axis of a quasi-2-D case: cylindrical
  bending) and `FSI_Test_Motion = A, f [, mode [, shape]]` (an imposed bending
  mode for testing the moving-grid fluid alone; mode 2 = a pluck that releases
  the structure after 2.25 periods - the start-up of the benchmark), and in
  `$Numerical_Scheme` `SW_Quasi_2D = 3` (strictly 2-D across the span of an
  extruded grid).
* **Structure** (`fsi_shell.F90`): geometrically non-linear (co-rotational,
  large rotation / small strain) 4-node shell - membrane + MITC4 plate +
  drilling term - with per-element material and thickness, lumped mass,
  Rayleigh damping, generalized-alpha time integration with Newton iterations,
  banded Cholesky after a reverse Cuthill-McKee numbering.  Verified alone
  (`fsi_test`): linear cantilever (beam theory to 1.2 %), free-vibration
  frequency (0.6 %), Turek-Hron CSM3 (d_y mean -63.5 / amplitude 65.1 mm,
  f 1.096 Hz against -63.6 / 65.2 mm, 1.0995 Hz).
* **Coupling** (`fsi.F90`): wall node forces (the force-file values) to the
  shell nodes and the shell motion back to the wall nodes through the same
  attachment weights (energy-conserving), the part frame's rigid motion and
  inertial loads kept; mesh MORPHING of the part's blocks (boundary layer rigid
  with the wall, the distortion absorbed towards the overset outer boundary),
  rebuilt metrics, discrete grid velocity in the ALE fluxes, re-projected face
  fluxes, geometric-conservation term; STRONG coupling - sub-iterations inside
  every time step with Aitken / interface quasi-Newton relaxation (the added
  mass of water on a rubber plate is handled), restarted from the level-n field.
  A wall patch that is not sheet-like or overlaps a larger one (the KP505
  tip caps) is passive and follows the other patches' mid-surface.
* **Files**: `Output_data/<case>_fsi_<part>_node<n>.dat` (displacement /
  velocity / von Mises history of a monitor node), `<case>_fsi_<part>_shell.*`
  (the mid-surface with thickness, material and clamp flags), the wall surface
  files with `Disp_x/y/z` and `vonMises`; the structural states in the restart
  file.
* **Verification** (`cases/FSI_PfS2a`): the FSI-PfS-2a benchmark of Kalmbach &
  Breuer (2013) as computed by Fabbri et al. (2023) - a fixed cylinder D = 22 mm
  with a para-rubber plate (50 x 2 mm, E 4.1 MPa) and a steel rear mass in
  water at 1.385 m/s (Re 30 470): see section 36.2 of PROJECT_RECORD.md and
  `verification/update24_fsi/README.md` for the frequency, mode and amplitude
  obtained on the quasi-2-D sandbox grid against the experiment (11.25 Hz,
  d_y/D = 0.63 ... 0.65).  `cases/KP505_FSI`: the KP505 propeller with the five
  blades flexible (one child body per blade, clamped inside the hub).

## Update 23 (2026-10-18): on-the-fly turbulence statistics

Details: PROJECT_RECORD.md section 35; verification in
`verification/update23_statistics` (cases written by `tools/stat_cases.py`,
checks by `check_stats.py`); all inputs annotated in `probCFD_template.nml`
(`$Statistics`).  New source file `statistics.F90` (29 files).

* **`$Statistics`** (`SW_statistics = 1`): raw moments accumulated every
  `Statistics_Sample_Intvl` steps inside `Statistics_Start_T ..
  Statistics_End_T` (0 = End_T), the fluctuation statistics formed at
  output (mean, rms = √(⟨x²⟩ − mean²), ⟨x′y′⟩ = ⟨xy⟩ − mean_x mean_y) with
  per-point sample counters that skip the steps at which a point is
  blanked (hole of a moving grid); `Coverage` = sampled fraction of the
  window, unsampled points are written as 0 (numbers only in the PLOT3D
  files).  Accumulators in the restart file.
* **Where** (`Statistics_type1..4`): the original grid blocks (node based,
  on the t = 0 grid; `Statistics_Type1_Frame = 1` turns the vectors of a
  moving block into its part's frame - blade-fixed statistics on a
  propeller block), and user-defined **boxes**, **cylinders** and
  **planes** (`Statistics_box / _cyl / _plane` groups with `Origin` /
  `Facecenter`, `Axis1` / `Axis`, `Axis2` / `Ref`, `Length`, `nPoints`,
  `Radius`, `Theta`, `Polar_Vel`, attached to the background or to a part
  through `Mother_body / _child_body / _sub_child_body`), sampled in the
  moved grid by the probe-point machinery of update 22.
* **What** (`Statistics_first_group`): Coverage, mean and rms of u v w,
  TTL_Prs, Dym_Prs, Shear_Stress, VoF, Vapor, Temperature, ν_t and the
  model variables when present, the Reynolds stresses and k_res;
  (`Statistics_second_group`, velocity only): the 9 mean gradients, the 9
  rms fluctuating gradients, the 3 cross terms, ε_res = ν⟨(∂u′ᵢ/∂xⱼ +
  ∂u′ⱼ/∂xᵢ)²⟩/2 from the **face differences of the staggered
  velocities** (∂u/∂ξ, ∂v/∂η, ∂w/∂ζ compact over the cell's own faces,
  the transverse derivatives the central average of the edge
  differences; accumulated at the cells for the grid output), with a
  turbulence / IDDES model ε_model (the model's
  destruction term: k-ε ε, k-ω/SST β*kω, SA/WMLES ν_t|S|²), ε_total,
  k_model, k_total, the Kolmogorov scale η, the Taylor scale λ_g =
  √(10νk/ε), Re_λ, Δ/η and the mean vorticity.  A cylinder with
  `Statistics_cyl_Polar_Vel = 1` writes axial / radial / tangential
  components.
* **Files** in `Statistics_volume_solutions/`: `<case>_stat_grid_g1/g2_NNNNNNNN.p3d`
  on `<case>_stat_grid.grd`, `<case>_stat_box<n>_g1/g2_NNNNNNNN.p3d` on
  `<case>_stat_box<n>.grd` (cylinders / planes alike), one `.nam` per
  family with the SI units (`u_mean [m/s]`).  The volume and surface
  `.nam` files now carry units too (`u [m/s]`, `TTL_Prs [Pa]`).
* Verification (`cases/Statistics_ActuatorDisk`, `cases/Statistics_TowedCylinder`):
  the grid statistics equal the statistics of the per-step volume files
  to the file precision (8·10⁻⁸ of 1.29 m/s); the lattice points equal the
  probe time histories; ε_res recomputed from the exported components
  to 10⁻¹³ (Cartesian) and to the single precision of the files on the
  O-grid of the cylinder, seam nodes included; a restart continues
  bitwise (grid) / to the last bit of the files (lattices); the
  hole of the towed cylinder appears as Coverage 0; the ring around it
  gives the antisymmetric tangential velocity.  Regression: the probe and
  force files of the towed cylinder are bitwise identical to update 22.

## Update 22 (2026-10-13): probe points, shear-stress field, wall surface files

Details: PROJECT_RECORD.md section 34; verification in
`verification/update22_probes_surface` (cases written by
`tools/probe_cases.py`, the surface integral checked by
`check_surface.py`); all inputs annotated in `probCFD_template.nml`.
New source files `probes.F90` and `surface_output.F90` (28 files).

* **Probe points** (`Probe_point<n> = x y z`, grid frame at t = 0):
  in the new group `$Probe_background` the point follows the background
  (fixed, or moving with the background blocks when `Moving_BGR = 1`);
  inside a `Body` / `Child_body` / `Sub_child_body` group of `$SixDOF`
  it follows that part's 6-DoF motion (a point on a propeller rotates
  with it).  One file per point in `Output_data/`: `<case>_probe<n>.dat`,
  `<case>_body1_probe<n>.dat`, `<case>_body1_child1_probe<n>.dat`,
  `<case>_body1_child1_sub1_probe<n>.dat`, every `Data_Intvl` steps:
  the current position, u v w, TTL_Prs, Dym_Prs, VoF / Vapor /
  Temperature when present, Shear_Stress, the cell (block i j k) and its
  overset status.  The values are the trilinear interpolation of the
  volume-file node values in the cell that holds the point (located in
  the moved grid every export step; the active cell wins over a fringe
  cell); a point in a hole (inside a solid, behind a moving body),
  outside the grid or in a cell no rank holds writes `NoData`.
* **`Shear_Stress`** [Pa], a new variable of every `.p3d` directly after
  `TTL_Prs` (the `.nam` file lists it): μ_eff √(2 S_ij S_ij) at the cells
  from the cell-centred velocity gradient, μ_eff = laminar (phase
  weighted) + eddy viscosity of the turbulence / IDDES model (none
  without a model); at every solid-wall node the friction stress
  magnitude |τ_w| of the wall faces (the wall-function value where one
  is used) - the same node average that integrates to the friction
  force.  No smoothing is applied near the wall (the field shows the
  discrete gradients; the wall value is untouched).
* **Wall surface files** (`SW_Prs_Frc_surf = 1` in `$Export_Interval`):
  `Volume_solutions_grids/<case>_surf_NNNNNNNN.grd` / `.p3d` with every
  volume solution, one PLOT3D block (ni × nj × 1) per solid-wall patch of
  every block, six variables [Pa]: `Prs_x Prs_y Prs_z` = −p n (the
  pressure force per unit area on the wall), `Frc_x Frc_y Frc_z` = the
  friction stress vector on the wall; `<case>_surf.nam`,
  `<case>_surf_patches.txt` (zone → block, patch, body).  The per-face
  values are exactly those of the force integration (its face loop is
  now the shared `wall_face_stress`), and the node values are dual-area
  weighted averages (overset wall weights included) built so that the
  trapezoidal integration over a patch - mean of the four node values ×
  face area, e.g. Tecplot "Integrate" - returns the pressure and friction
  forces of the force files exactly (checked to the single-precision
  storage of the file: 9·10⁻⁸ relative on the towed cylinder, 9·10⁻⁹ on
  the overset foil).  `Frame_Prs_Frc_surf = 0`: the grid axes of the
  volume files (the moved grid when the background moves; the force
  files use the SNAME earth axes X = −x, Y = y, Z = −z); `1`: the
  background frame (the surface at its original grid position).
* Verification (`cases/Probes_TowedCylinder`, `cases/Probes_SteeredFoil`):
  probe points following the towed cylinder and the steered foil move
  exactly with their parts (x = 0.6 − 0.5 t²/0.2 behind the cylinder;
  z = −0.97 − 0.8 sin θ at the foil's trailing edge) and stay in their
  cells, background points follow the moving background, a point inside
  the cylinder writes `NoData`; the surface files of both cases integrate
  to the force files to 10⁻⁸; `Frame_Prs_Frc_surf = 1` keeps the towed
  surface at its original position.
* Regression: Cylinder_MovingBGR (5 steps, surface files on): monitor,
  force and motion files bitwise identical to update 21; the volume file
  gains `Shear_Stress`, every other variable is bitwise identical.

## Update 21 (2026-10-12): PID controllers

Details: PROJECT_RECORD.md section 33; verification in
`verification/update21_pid` (cases written by `tools/pid_cases.py`,
summaries by `check_pid.py`); all variables annotated (commented-out
examples) in `probCFD_template.nml`.  New source file `pid_control.F90`
(26 files; it compiles before `sixdof.F90`).

* **A controller** is a block of `PID_*` lines inside a `Body` /
  `Child_body` / `Sub_child_body` group or inside a `Force_model` group of
  `$SixDOF`, opened by `PID_use = 1` (one per grid part, up to four per
  force model).  `PID_type`: 1 **throttle** (output = rotation rate
  [rev/s]), 2 **steering** (output = angle [deg]), 3 **direct force**
  (output = one component `PID_Component = Force_Fx .. Force_Mz` of a
  `Force_type = 3` point force, added to the formula value).
* **What the output drives**: a grid part (hull, rudder, drive, shaft)
  is rotated about `PID_Axis` (grid frame at t = 0) through its
  `Rot_Center` by the controller - throttle: a rotation at the output
  rate, steering: the output angle - and its `Motion_DOF` is ignored;
  the motion files, the mother / child coupling and the force output
  are those of update 14.  A Hough-Ordway / BEM disk: a throttle
  controller replaces `Prop_RPS` (`Prop_RPS` = initial value), a steering
  controller turns the disk axis about `PID_Axis` through
  `PID_Axis_Center` (default: the disk centre); a disk may carry both.
  A point force: the direct-force output is added to the component.
* **What is measured**: a 6-DoF quantity of any part - `PID_Target_body`
  [, `PID_Target_child_body` [, `PID_Target_sub_child_body`]],
  `PID_Target_type` 1 = location (`PID_Target_var` = x y z [m] roll pitch
  yaw [deg]; a Body in the earth frame, a child relative to its mother)
  or 2 = velocity (u v w [m/s], roll_rate pitch_rate yaw_rate = p q r
  [deg/s] in the body frame).  `PID_Target_value` lists up to 32 targets
  reached in sequence (the next one when |error| <= `PID_Target_tol` for
  `PID_Target_hold_T` seconds); angle errors take the short way round.
* **Law**: `out = PID_Out_init + P e + I ∫e dt − D d(measured)/dt`,
  e = target − measured (`PID_P`, `PID_I`, `PID_D`; the gains carry the
  sign of the actuator), clipped to `PID_Out_limit` (min max, with
  anti-windup) and to `PID_Rate_limit` per second, held at the initial
  output before `PID_Start_T`.  `PID_Out_init` is ignored for a disk
  throttle (`Prop_RPS` is the initial rps).
* **Output** `Output_data/<case>_<part>_pid<i>.dat` (step, time, target
  index, target, measured, error, P / I / D terms, output, angle [deg]);
  the summary lists every controller; the controller states are in the
  restart file (an older file starts them fresh).
* Verification (`verification/update21_pid`, cases `cases/PID_Steer`,
  `PID_Throttle`, `PID_Direct`, `PID_Disk`, `PID_DiskSteer`):
  a foil steered by a controller
  on its own pitch reaches +4 and then −4 deg with rate-limited ramps
  (0.6 deg per step), its motion file reproducing the controller's angle
  and rate exactly; a cylinder spun by a throttle controller holds its
  yaw rate at 60.000 deg/s (output exactly 1/6 rev/s); on a damped
  100-kg mass (the ActuatorDisk box declared body 1) a direct-force
  controller holds +5 / −5 mm, a disk throttle (rps in 0..12) and a disk
  steering (±40 deg, 200 deg/s) push it to +5 and +10 mm; a restart with
  a controller continues identically.  All five actuation paths and the
  sequential targets, tolerances, hold times, output and rate limits
  and the anti-windup are exercised; `check_pid.py` summarises a run.
* Fixes that came with the tests (33.1): several force models on one
  mother now all act on it (the external wrench was overwritten by the
  last model since update 19); a restart re-evaluates the force models
  from the restored field, poses and controller states (before, the
  first restarted step used the sources and wrenches of the initial
  field); the velocity-Verlet half-kick flag is restored (a 0.3 %
  position offset after a restart of a free DOF); a restarted box
  reproduces the continuous run to all printed digits.
* Regression: Cylinder_Overset_2phase (10 steps) and ActuatorDisk
  (5 steps, a disk without a controller) are bitwise identical to
  update 20 (monitor, `.p3d`, force-model file).

## Update 20 (2026-10-06 / 11): moving background

Details: PROJECT_RECORD.md section 32; verification in
`verification/update20_moving_background` (cases written by
`tools/bgr_cases.py`); all variables annotated in `probCFD_template.nml`.

* **`Moving_BGR = 1`** (`$Case_title`, next to `Mode`): the background
  blocks (`Object_ID = 0`) follow the Body `BGR_Following_Body` as one
  rigid group in the DOFs permitted by `Moving_BGR_DoF` (six 0/1 flags,
  e.g. `1 1 0 0 0 1` = planar motion for manoeuvring in waves, `1 1 1 1 1
  1` = everything), optionally low-pass filtered (`Moving_BGR_Filter_T`
  [s]) - CFDShip-Iowa's `mesh_moving = 2` / `is#_o` mechanism.  The hull
  blocks move relative to the background only in the other DOFs; the far
  field (still water, incident wave, bottom) stays earth-fixed; the
  background is its own rigid motion group of the overset with its own
  donor cache.  Dirichlet velocities on the background are earth-frame
  values (a free-running ship has `U_inf = 0`).  Output
  `Output_data/<case>_bgr_motion.dat`; restart carries the background
  state (same `Moving_BGR` as the file).
* **Waves on a moving background**: the moved grid coordinates are the
  earth positions, so the flag-12 values need nothing; the relaxation
  zones stay attached to the grid faces and are now given per face,
  `Relax_L = xmin xmax ymin ymax zmin zmax` (the update-18 keys remain as
  aliases); the probes are earth-fixed.
* **`Ovst_Reassemble_Tol`** (`$Overset`, [m], default 0): the dynamic
  overset is re-assembled only when a moving part has moved more than the
  tolerance relative to the background since the last assembly (the
  rigid-group donor caches and the warm-started donor search of update 12
  stay as they are).
* Validation (`verification/update20_moving_background`, cases
  `cases/Cylinder_MovingBGR`, `cases/WaveTank_MovingBGR`): the towed
  cylinder with the background following it reproduces the inertial
  forces of the ramp in both frames to 0.3 % (added mass C_a = 1.00,
  Froude-Krylov of the accelerating stream), agrees with the body moving
  through the static background to 0.1 %, re-uses the assembly in
  59 of 60 steps with `Ovst_Reassemble_Tol` (3·10⁻⁶), restarts with the
  background state; the wave tank surging / yawing through the
  earth-fixed Airy wave keeps the wave at the earth-fixed probes
  (heights within 2-6 % of the Airy wave, the static tank: +2-5 %).  Open item: a frame-dependent offset of a few
  per cent of ½ ρ U² D L in the pressure force of a towed body after the
  start-up transient, present in the moving-body path since update 12
  (section 1 of the verification README).
* Fixes that came with the validation (32.1): restart record of the
  background, ramped wall velocity in the force integration and the
  wall functions, boundary values of the restored time on a restart,
  discrete (displacement) grid velocity of the ALE fluxes.
* Regression: with `Moving_BGR = 0` Cylinder_Overset_2phase (10 steps) is
  bitwise identical to update 19.

## Update 19 (2026-10-06): force models - Hough-Ordway disk, BEM disk, point force

Details: PROJECT_RECORD.md section 31; cases `cases/ActuatorDisk`,
`cases/Foil_Spring`; verification in `verification/update19_force_models`;
every variable is listed and annotated in the `$SixDOF` group of
`probCFD_template.nml`.

* **Force models** (`force_model.F90`): "bodies without a grid" in the 6-DoF
  hierarchy - `Force_model = n ... Force_model_End = -1` groups of `$SixDOF`
  that hang on a Body / Child_body / Sub_child_body (`Mother_body`,
  `Mother_child_body`, `Mother_sub_child_body`) or on nothing
  (`Mother_body = 0`, fixed in space); up to 64, any mix, no children.
  `Body_End`, `Child_body_End`, `Sub_child_body_End = -1` close the parts.
  The geometry points move with the mother; the wrench on the mother
  (about its Rot_Center, thrust and torque together) is transferred up the
  hierarchy like the fluid forces; the disk forces are momentum sources of
  the fluid, weighted by the liquid fraction.
* **`Force_type = 1`** Hough-Ordway actuator disk: `Prop_P1/P2, Prop_D,
  Prop_Hub, Prop_RPS, Prop_KT, Prop_KQ` (K_T, K_Q quadratic in J; J from
  the disk-averaged inflow relative to the moving disk); the discrete disk
  cells carry exactly T and Q.
* **`Force_type = 2`** pure blade-element momentum disk: `Prop_Blades,
  Prop_Sections, Prop_rR, Prop_cD, Prop_PD, Prop_fc, Prop_CD0` (+ optional
  `Prop_CLa, Prop_CDk`); local inflow angle per cell, thin-airfoil lift
  with camber, parabolic drag, Prandtl tip / hub loss; K_T, K_Q reported.
* **`Force_type = 3`** point force: `Force_Point`, `Force_Fx/Fy/Fz` and
  optional `Force_Mx/My/Mz` as quoted formulas of `x y z roll pitch yaw u v
  w p q r t` of the mother (SNAME body frame; springs, moorings, control
  forces, empirical thrusters ...); without a mother the force goes into
  the fluid at the point.
* Output `Output_data/<case>_force_model<n>.dat` (mother's body frame,
  moments about the mother's Rot_Center and about the model's own centre;
  disks: n, J, K_T, K_Q, T, Q, liquid fraction, cells), `[ Force model ]`
  monitor line.  `tools/force_cases.py` writes the test cases.
* Validation: a Hough-Ordway disk fixed in a uniform stream converges to
  the momentum-theory solution (T = 43.0 N applied vs 43.7 N from
  2 ρ A U_d (U_d − U) at the computed disk inflow); the pure BEM disk
  reports its K_T, K_Q (0.33 at J = 0.79 with C_La = 2π, 0.28 at J = 0.75
  with 3.8 — about twice the B4-70 series: the section data of a pure BEM
  must be set against open-water data); a stiff point-force spring on the
  overset foil (m = 10⁴ kg, period 0.1 s) oscillates with the analytic
  amplitude (9.97 / 10.0 mm) and period (+0.7 %, the foil's added mass).
  The velocity-Verlet half-kick added to the body integrator for this
  test changes moving-body results of earlier updates at the O(Δt) level.
* Regression: without a `Force_model` group Cylinder_Overset_2phase (10
  steps) is bitwise identical to update 18.

## Update 18 (2026-10-04): incident waves (regular and irregular)

Details: PROJECT_RECORD.md section 30; case `cases/RegularWave_Tank`;
verification in `verification/update18_waves`.

* **Incident waves** (`$Wave`, `waves.F90`, `wave_fields.F90`): regular
  (Airy or Stokes 2nd order, deep or finite depth) and irregular seas
  (Bretschneider / JONSWAP / Pierson-Moskowitz, long- or short-crested
  with cos^2 spreading, or a component file) for two-phase VoF runs,
  following the CFDShip-Iowa treatment of the exact potential wave on the
  far-field boundaries with the three changes agreed in the discussion:
  the wave is defined in the earth-fixed frame and the grid frame is
  converted with one shift X = x - integral(U_mean) (encounter frequency
  for a body held in a stream automatically; the later moving-background
  mode only adds its displacement to `wave_shift_extra`); one far-field
  flag for all variables; relaxation zones towards the incident wave.
* **`.bcs` flag 12** on the u, v, w, p AND alpha lines of a surface (any
  face: inlet, outlet, sides, bottom; the u/v/w values are the mean
  velocity): velocity = ramped mean + orbital velocity (frozen at the
  surface value in the air), alpha ghost = fill below the exact elevation,
  pressure ghost = interior + exact wave-pressure difference (the face
  stays a fixed-flux face of the projection).
* **Relaxation zones** `Relax_In_L / Relax_Out_L / Relax_Side_L` [m]
  (Jacobsen et al. 2012, waves2Foam weight) from the x-min / x-max / y
  faces of the bounding box of the wave blocks: velocity (water only)
  and alpha are blended towards the incident wave after every step; they
  absorb body-generated waves at the inlet and let the incident wave out.
* `Wave_Ramp_T` (half-cosine ramp), `Wave_Init = 1` (exact field at
  t = 0), wave probes (`Output_data/<case>_wave.dat`: computed and exact
  elevation), the component list in `Output_data/<case>_components.wav`,
  a `[ Wave ]` monitor line; `Wave_Ramp_T` and the zone lengths are
  runtime-changeable.  `tools/wave_cases.py` generates the tank cases.
* Validation (2-D tank, λ = 1 m, a = 1 cm, 32 cells / λ, 200 steps /
  period): flat start with a two-period ramp - the wave develops at the
  group velocity, the amplitude is kept within 2.5 % over three
  wavelengths of free interior, the numerical phase speed is 2-3 % low
  on 32 cells per wavelength (a spatial error: the phase lag after two
  periods falls 20° / 13° / 10° on 16 / 32 / 64 cells per wavelength and
  is independent of Δt; the steady lag behind the inlet zone falls from
  16.5° to 9.9° when the grid is doubled); exact initial field (`Wave_Init`)
  kept for six periods; mean current 0.5 m/s - the probes oscillate at
  the encounter period; JONSWAP sea (Hs 2 cm, Tp 0.8 s, 60 components) -
  see `verification/update18_waves/README.md`.
* Not for waves: `Time_scheme = 2` with the geometric VoF (the stage
  average smears the interface: 12 % amplitude loss in two periods).
* Regression: with `SW_Wave = 0` Cylinder_Overset_2phase (10 steps) is
  bitwise identical to update 17.

## Update 17 (2026-10-04): compressible mixture, thermal cavitation

Details: PROJECT_RECORD.md section 29; verification in
`verification/update17_compressible_thermal`.

* **`Cav_compressible = 1`** (`$Cavitation`): the liquid-vapour mixture is
  compressible with the Wood / Wallis sound speed built from `Snd_liquid`
  (1480) and the vapour as an isothermal ideal gas (c_v^2 = p/rho_v, about 370 m/s at p_v): 1/(rho c^2) = sum alpha_k/(rho_k c_k^2),
  so c drops from 1480 m/s in water to a few m/s in the 10-90 % mixture.
  The pressure equation becomes (1/(rho c^2)) dp/dt + div u = S_cav, the
  compressibility term implicit on the Poisson diagonal (next to the
  pressure-implicit mass-transfer term of update 15); the vapour volume
  shrinks with the pressure rise (dV_v/V_v = -dp/(rho_v c_v^2)); an
  all-Neumann pressure problem no longer needs the level pinned.  The air
  of a two-phase cell stays incompressible.  Isothermal - no energy
  equation needed.  Monitor: max acoustic CFL, min mixture sound speed.
  Verified: a closed tube of a 30 % mixture compressed by a piston rises
  at 93 Pa/s (Wood law 92.3), the piston launches a pressure step of
  rho c u = 25.43 Pa (25.41) that travels at 3.63 m/s (3.63).
* **`Cav_thermal = 1`** (needs `SW_Energy = 1`, any fluid): p_v(T) and
  rho_v(T) by Clausius-Clapeyron through the user's (`Prs_Vapor`, `Lou_v`,
  `Latent_Heat`) at `T_init`, the thermodynamic delay of the evaporation
  (B-factor: the bubble of radius R_B is cooler than the liquid by
  dT* = 1/2 (rho_v L/(rho_l c_p)) R_B/sqrt(alpha_l t), t = `Cav_thermal_t`
  (default 1 ms); evaporation driven by p_v(T - dT*)), and the latent-heat
  sink of the energy equation with its derivative dp_v/dT taken
  implicitly.  `Pv_of_T = 1` alone keeps the water (Buck) law of update
  16.  `vap_cel_initial`: initial vapour fraction (tests).  LN2 at 20 m/s:
  the cavity cools by 0.3-0.9 K and loses 9 % of its volume; cold water:
  0.01 K, 0.7 %.
* Use the compressible mixture after the start-up transient (restart, or
  `Cav_Start_T` after it): from an impulsive start the far field dips
  below p_v and the vapour fog expands with the pressure.  A cavitating
  cylinder continued from a restart follows the incompressible solution
  within 1-5 %.
* Regression: with the new switches off the cavitating cylinder is
  bitwise identical to update 16b.

## Update 16 (2026-10-03): heat / energy equation

Details: PROJECT_RECORD.md section 28.  Case: `cases/NaturalConvection_Cavity`
(de Vahl Davis square cavity, Ra = 1e5, Pr = 0.71, with the benchmark
comparison in its README); verification runs in
`verification/update16_energy` (1-D conduction against the series
solution, heat-flux / convective boundary conditions against the exact
linear profile, step advection with the four schemes).

* **Temperature equation** (`$Energy`, new module `energy.F90`):
  `SW_Energy = 1` solves dT/dt + u.grad T = div(k_eff grad T)/(rho c_p)
  for the temperature T [K] at the cell centres (field id 20,
  'Temperature' in the .p3d files, min / max in the monitor and the
  monitor file).  One-fluid mixture: rho c_p and k volume-weighted over
  liquid / air / vapour (`Cp_l`, `k_l`; `Cp_a`, `k_a` read with two
  phases; `Cp_v`, `k_v` with cavitation), k_eff = k + mu_t c_p/`Pr_t`
  (default 0.85).  `T_init` is the initial and the Boussinesq reference
  temperature.
* **Numerics** with their own controls: `Max_Itr_Engy`, `Min_Itr_Engy`,
  `Rsd_Lmt_Engy`, `SOR_fct_Engy` (SOR sweeps), `Cnvt_scheme_Engy` 1
  upwind-1 / 2 upwind-2 / 3 QUICK / 4 ENO-2 (minmod; bounded),
  `Cnvt_mode_Engy` and `Dffs_mode_Engy` 1 explicit (sub-cycled) / 2
  implicit / 3 hybrid per direction by the local CFL / Von Neumann
  number (`cfl_switch_Engy`, `vn_switch_Engy`), all runtime-changeable.
  Non-conservative convective form (bounded with the cavitation
  dilatation and at overset fringes), orthogonal diffusion implicit /
  explicit, cross-derivative diffusion explicit, implicit Euler per inner
  iteration, SSP stages of `Time_scheme = 2`.
* **Boundary conditions**: one more line per surface of the `.bcs` file
  (after the alpha line, which is then required in a two-phase run):
  `0` zero gradient (adiabatic), `1 T_w` Dirichlet [K] (isothermal wall or
  inlet; on a flag-11 wall the thermal wall function of Kader gives the
  wall heat flux), `3 q_w` heat flux [W/m2] into the fluid, `4 h T_inf`
  convective (Robin).  The heat flow [W] and the mean wall temperature of
  every such surface are written to `Output_data/<case>_heat.dat` every
  `Data_Intvl` steps; the total wall heat is in the energy monitor block.
* **Buoyancy**: `Beta_l`, `Beta_a` [1/K] (default 0 = passive scalar) add
  the Boussinesq acceleration -beta (T - T_init) g to the momentum
  equation.
* **Cavitation coupling** (`$Cavitation`): `Latent_Heat` (sink -mdot L of
  the energy equation, default water) and `Pv_of_T = 1` (vapour pressure
  from the saturation law p_v(T) of water at the cell temperature instead
  of the constant `Prs_Vapor`).
* **Restart** file v3 (older files are read; a file without T starts the
  energy equation from `T_init`, so a flow-only run can be continued with
  `SW_Energy = 1`).
* **`Prs_form = 2`** (`$Numerical_Scheme`, default 1 = unchanged):
  incremental pressure form of the fractional step (the gradient of the
  latest pressure in the predictor, the increment projected).  Found
  while validating the cavity: with the non-incremental form a fluid at
  rest under gravity develops spurious velocities at no-slip walls
  wherever the wall-normal diffusion is implicit (VN >= `vn_switch`, i.e.
  wall-clustered grids), because the predictor diffuses the body-force
  velocity dt g; form 2 keeps the fluid exactly at rest and is the right
  choice for gravity / buoyancy-driven flows.  With form 2 the initial
  pressure of a single-phase run with gravity is hydrostatic.
* **Blank variable** (16b): the `Blank` of the `.p3d` solution files is now
  the exact overset cell status without any combination of neighbouring
  cells - node (i,j,k) carries the status of cell (i,j,k) (the max-face
  nodes repeat the last cell): 1 active, 0 fringe, -1 hole.  The separate
  `<case>_NNNNNNNN_cblank.p3d` file is no longer written (`tools/blank_slice.py`
  works without it).
* Fixes found on the way: the vapour-fraction ghost cells at Dirichlet-p
  faces read the metric out of bounds (field 19 fell through to the
  pressure flags of `apply_bc`; now zero gradient like the other
  scalars, together with the temperature); the SOR sweeps of the
  turbulence transport ran per domain with a collective inside - a
  deadlock when the ranks own different numbers of domains (fewer ranks
  than csv sub-domains) - now one collective per sweep over all owned
  domains (results unchanged); the wall distance of the turbulence
  module was in metres while y+ assumed Lch units (exact only for
  Char_L = 1; now divided by Lch).

## Update 15 (2026-10-02): cavitation models

Details: PROJECT_RECORD.md section 27.  Cases: `cases/NACA66mod_cav`
(NACA66(mod) section, sigma = 0.91) and `cases/DelftTwist_cav` (Delft
twist-11 hydrofoil, sigma = 1.07, SST + Reboud).

* **Three transport-equation models** (`$Cavitation`, new module
  `cavitation.F90`): `Cav_model` 1 Schnerr-Sauer, 2 Zwart-Gerber-Belamri,
  3 Kunz, with their literature constants as name-list defaults
  (`SS_Nuclei`, `SS_Radius`; `ZGB_Fvap`, `ZGB_Fcond`, `ZGB_Rnuc`, `ZGB_Rb`;
  `Kunz_Cdest`, `Kunz_Cprod`, t_inf = Char_L/Char_U).  `SW_Cav` switches the
  model on, `Prs_Vapor`, `Lou_v`, `mu_v` give the vapour, `Cav_Start_T` the
  time the mass transfer starts (50-step ramp, runtime-changeable).
* **Three-phase formulation**: the vapour is the volume fraction of the
  liquid-vapour part of a cell; the water/air VoF `alpha` (CICSAM or PLIC)
  is untouched, rho and mu are the three-phase mixture values at every
  half-index location.  Mass transfer only in water (none where the air
  fraction exceeds 0.5); a cavity reaching the free surface condenses by
  the pressure itself.
* **Vapour transport** (`Vapor_scheme`): 1 explicit CICSAM, 2 implicit
  CICSAM (default), 3 second-order ENO (minmod) face values with the time
  integration of `Cnvt_mode` (explicit, implicit or hybrid by the Courant
  number) - all through the VoF transport of `vof.F90`, with the
  evaporation explicit and the condensation implicit so that the fraction
  stays in [0,1] (the residual clipping is reported).
* **Continuity**: div u = mdot (1/rho_v - 1/rho_l) enters the Poisson
  equation with its pressure derivative taken implicitly (a diagonal
  term that pins the cavity pressure at p_v - without it the cavity
  flashed on and off from step to step), the vapour fraction is made
  consistent with the implicit part after the pressure solve, the
  explicit production is capped at 0.2 per step; the divergence monitor
  and the outflow mass-flux correction include the source.  Use
  `Iner_Itn = 2` with cavitation (Picard iteration of the coupling).
* **Absolute pressure**: with `SW_Cav = 1` the solver pressure is the
  absolute pressure - the level is the Dirichlet p of the `.bcs` file (or
  `RefPrs`) and `p_cel_initial`; sigma is printed at start-up, `Dym_Prs` of
  two-phase files has the level removed.
* `Reboud = 0/1` (turbulence variables): Reboud's eddy-viscosity
  correction (n = 10) in cavitating cells.
* **Validation** (PROJECT_RECORD.md 27.5, `verification/update15_cavitation`):
  cylinder at sigma 0.7 for the three models, three vapour schemes, SSP
  and restart (smooth cavity growth, div u = S to 1e-3); NACA66mod at
  sigma 0.3 (thin sheet along the suction side where Cp < -sigma, the
  pressure pinned at p_v); DelftTwist start-up.
* **Output**: `Vapor` after `VoF` in the `.p3d` files (name file updated),
  vapour volume and max fraction in the monitor line and the log block
  (with the evaporating / condensing cell counts), vapour fraction in the
  restart file (`.rsto` v2; v1 files still read).  Overset: the vapour
  fraction (and now also k, omega, mu_t) is interpolated into the fringe
  cells.

## Update 14b (2026-10-02): Mode = 2 grid-movement check, missing-file stops

Details: PROJECT_RECORD.md section 26.7.  Case: `cases/GPPH` delivered with
`Mode = 2` (stern drive steering to +45 deg about its raked axis, propeller
turning on its shaft, both from impose motion files).

* **`Mode = 2`** (`$Case_title`): grid-movement check of the imposed
  motions before the CFD run.  Right after the name list, the grid and the
  csv are read (no overset assembly, no field arrays: any grid size), every
  part with imposed DOFs follows its impose motion file (free DOFs held,
  fixed children ride with their mothers), and the grid is written every
  `Output_Intvl` steps of `Del_T` up to `End_T` (and at the last step) to the
  sub-folder **`Check_grid_movement/`** as `<case>_NNNNNNNN.grd` (PLOT3D,
  binary single precision) **plus a companion PLOT3D solution file
  `<case>_NNNNNNNN.q`** whose block headers carry `TIME = t` (and `RE` =
  step): Tecplot's PLOT3D loader takes the solution time of a grid /
  solution pair from that header, so loading all `.grd` + `.q` pairs
  together activates the solution-time bar and animates the movement
  (without the `.q` files all grids are shown at once).  Its five
  variables (`<case>.nam`: Part_ID, dX, dY, dZ, Disp) are the 6-DoF part
  of every block and the node displacement from the initial grid, so the
  moved parts can be coloured.  `<case>_movement.dat` lists the relative
  and absolute attitude and the rotation-centre position of every imposed
  part at the exported steps; the log prints one line per export.  The run
  then stops (`Reason_for_stop.md`: "Mode = 2 ... finished").
* **Stops with a message** (log + `Reason_for_stop.md`, clean MPI exit):
  "There is no name list file" (no `probCFD.nml`), "There is no grid file"
  (no `<case>.grd`), and for Mode = 2 "There is no impose motion" (no part
  with `Motion_DOF = 1`) or "There is no impose motion file: <name>".
* `tools/make_axis_rotation.py`: `--stop_deg A` turns to the angle A and
  holds (rudder / stern-drive steering).

## Update 14 (2026-10-01): child and sub-child bodies of the 6-DoF motion

Details: PROJECT_RECORD.md section 26.  Case: `cases/GPPH` (hull = body,
stern drive = child body, propeller + hub = sub-child body, annotated
`probCFD.nml`).

* **Three-level hierarchy** body > child body > sub-child body (several
  bodies, several children per body, several sub-children per child).  A
  child / sub-child has the same properties as a body (`Mass`, `CG`,
  `Rot_Center`, `Motion_DOF`, `Inertia`, `Prod_Inertia`, `Added_Mass`,
  `Motion_Damp`, `Motion_File`) in its OWN axes, plus `Mother_body` (and
  `Mother_child_body` for a sub-child).
* **csv**: rows `Child_Object_ID` and `Sub_Child_Object_ID` after
  `Object_ID` (both optional) give the child / sub-child a block belongs
  to (0 = the body / child itself).  Every part is its own rigid overset
  motion group.
* **Name list** (`$SixDOF`): `Body = k`, then its children `Child_body = c`
  each followed by its sub-children `Sub_child_body = s`, then the next
  body (`[body1] [Child_body1] [Sub_child_body1] ... [Child_body2] ...
  [body2]`); `Mother_body` / `Mother_child_body` are checked against the
  enclosing lines.
* **Motion of a child / sub-child is RELATIVE to its mother** (translation
  in the mother frame, Euler angles relative to the mother frame, DOF 0 =
  rigid with the mother); the absolute pose, velocity and acceleration
  follow by composition (`update_transform`).  Imposed motion from
  `<case>_body<k>_child<c>_imposed.dat` / `<case>_body<k>_child<c>_sub<s>_
  imposed.dat` (same format as the body's file): any rotation /
  translation through the Euler angles - `tools/make_axis_rotation.py`
  writes the table of a constant-speed rotation about an arbitrary
  (inclined) axis, e.g. a propeller on a 2.27-degree shaft.
* **Coupling**: a part moves with the composite mass / CG / inertia of its
  subtree at the current relative poses; the fluid wrench of a child (and
  its sub-children) is transferred to the mother except the components
  along the child's FREE DOFs (a free-wheeling propeller passes its thrust,
  not its torque; an imposed rotation is a motor-driven shaft and passes
  everything); the relative equations of a child include the apparent
  forces of the mother's motion (centripetal, Coriolis, Euler terms).
  Neglected: the inertial reaction of the relative motion of a light
  appendage on its mother.
* **Output** in the new folder **`Output_data/`**: the body files as before
  (`<case>_{force,moment}_{body,earth}<k>.dat`, `<case>_motion<k>.dat`,
  own walls) and `<case>_body<k>_child<c>[_sub<s>]_{force,moment,motion}.dat`
  with the force / moment in the MOTHER's frame (moment about the child's
  own O) and the motion relative to the mother (plus the absolute attitude
  and angular velocity).  The log lists the children indented under their
  body with their relative motion.  Restart files hold every part.
* **Validation** (PROJECT_RECORD.md 26.5): Box_6DoF with one block of body
  1 declared a fixed child (mass 60 + 40 kg, offset rotation centre) =
  the single-body run (identical motion, forces to 2e-9 when the rounding
  is identical; the remaining 4e-5 is the grid's own round-off sensitivity,
  shown by a 1e-10 m perturbation of the rotation centre); a child with
  imposed yaw + translation carrying a sub-child with imposed roll +
  translation = a single body with the composed Euler angles (grid to the
  file precision, absolute angular velocity exact); restart of the
  hierarchy exact; GPPH: start-up with the three parts (the moving steps
  need more memory than the sandbox has - run on your machine).

## Update 13 (13b of 2026-10-01): SSP time scheme, dynamic pressure, runtime name list, stop report, log layout, output folder, Tecplot variable names, restart

Details: PROJECT_RECORD.md section 25.

* **`Time_scheme`** (`$Solver_Iteration`, default 1): 1 = the existing scheme
  (first order in time), 2 = **2nd-order SSP (Heun)**: the fractional step
  (VoF transport, momentum predictor, pressure projection, velocity
  correction, turbulence) is applied twice, U(1) = S(U^n), U(2) = S(U(1)),
  and U^(n+1) = (U^n + U(2))/2 for every transported field; the pressure is
  that of the second stage.  Twice the cost per step.
* **Dynamic pressure**: two-phase solution files carry `p_dyn = p - p_hyd`
  as the variable after `alpha` (`p_hyd = rho_ref(z) |g| (z_fs - z)`, water
  below the still free surface, air above); variables: u, v, w, p, alpha,
  p_dyn (+ iblank for overset cases).
* **Runtime name-list updates**: `probCFD.nml` is re-read every step (rank 0,
  content hash); when it changed, the runtime-changeable values are applied
  to the running calculation and listed in the log: `Del_T`, `End_T`,
  `Iner_Itn`, `Time_scheme`, the momentum / pressure / VoF / turbulence
  iteration controls, `Output_Intvl`, `Log_Intvl`, `Data_Intvl`,
  `Cnvt_scheme`, `Cnvt_mode`, `Dffs_mode`, `cfl_switch`, `vn_switch`,
  `vof_cfl_max`, `t_ramp`, `Stop_Rsd_Max`, `Stop_CFL_Max`, `poi_verbose`,
  `Motion_Start_T`.  The step count follows a new `Del_T` / `End_T` from the
  current time on (the time is accumulated, not `step x Del_T`).  Static
  inputs (grid, phases, fluids, overset, turbulence model, VoF scheme, 6-DoF
  body data) are ignored until a restart.
* **`Reason_for_stop.md`**: written in the working folder whenever the run
  ends - the three stop criteria (NaN, residual / update norm, CFL), input
  errors and the normal end.  It holds the reason, the cells / values the
  criterion reported (all ranks), a table of the last 40 steps (CFL,
  residuals, iterations, liquid volume, wall time), a diagnosis with the
  settings to change, and the settings of the run.
* **Log**: the 6-DoF block lists force, moment, position, angles, linear
  and angular velocity one per line with aligned labels, an empty line
  between bodies, then the overset re-assembly counts; 21 empty lines
  between time-step groups (3 x before).
* **Variable names (Tecplot)**: the solution variables are written in the
  order u, v, w, TTL_Prs, Dym_Prs, VoF, Blank (single phase: u, v, w,
  TTL_Prs, Blank; Blank only for overset cases) and the PLOT3D name file
  `Volume_solutions_grids/<case>.nam` lists them: load it as the "Name
  file" of Tecplot's PLOT3D loader together with the grid and the solution
  and the variables appear with these names instead of F1V1, F1V2, ...
* **Restart**: `Mode` (first line of `$Case_title`): 0 = standard run from
  step 0, 1 = restart from the latest `<case>_NNNNNNNN.rsto` (the latest is
  the last entry of `<case>_restart.log` whose file exists).  `Restart_Intvl`
  (`$Export_Interval`, default 0 = none) writes a restart file every N steps
  and at the last step.  The file (binary, rank 0) holds the step, the time,
  every sub-domain's fluxes, face velocities, pressure, VoF and turbulence
  fields with their ghost layers, and the 6-DoF state (pose, velocities,
  accelerations, last forces) of every body; it is independent of the number
  of MPI ranks.  A restarted run moves the grids to the restored poses,
  re-assembles the overset, continues the monitor and the force / motion
  files (append) and reproduces the uninterrupted run (cylinder: to all
  monitor digits; Box_6DoF: forces to 1e-8 relative).  `Restart_Intvl` can be
  changed at run time; the restart files stay in the working folder.
* **Volume files**: `Volume_solutions_grids/<case>_NNNNNNNN.grd` (grid, also
  the initial one as `_00000000`) and `<case>_NNNNNNNN.p3d` (solution),
  (update 16b: the exact cell status is the `Blank` variable itself; no separate `_cblank.p3d` file any more);
  the grid check (`Ovst_only = 1`) writes `<case>_00000000.grd`,
  `_00000000_iblank.p3d` there (the `_cblank.p3d` of updates 12c-16 is gone since 16b).
  `tools/blank_slice.py` reads the new layout (and the old one).

## Update 12 (final 12c state of 2026-09-30, incl. the exact cell blank, self-overlap donors): 6-DoF rigid-body motion, dynamic overset

New module `sixdof.F90` (details: PROJECT_RECORD.md section 24).

* **Objects**: row 6 of `Input_polyMesh_probCFD.csv` (`Object_ID, o1, ..., oN`)
  gives the 6-DoF object of every block: 0 background, 1 first body and its
  refinement grids, 2 second body, ...; any number of bodies (up to 64).
* **`$SixDOF` group** at the end of `probCFD.nml`, per body after `Body = k`:
  `Mass` [kg], `CG` and `Rot_Center` (x y z, grid coordinates [m]),
  `Motion_DOF` (surge sway heave roll pitch yaw: 0 fixed, 1 imposed, 2 free),
  `Inertia` (Ixx Iyy Izz about the CG) and `Prod_Inertia` (Ixy Ixz Iyz)
  [kg m2], optional `Added_Mass` (6, stabilizes the explicit fluid-body
  coupling), `Motion_Damp` (6, artificial damping [1/s]), `Motion_File`;
  global `Motion_Start_T` (free DOFs released from this time on).
* **Convention: SNAME / Fossen.**  Earth frame x forward (= -x grid),
  y starboard (= y grid), z down (= -z grid); body frame aligned with it at
  t = 0, origin at the centre of rotation; Euler angles Z-Y-X (yaw, pitch,
  roll); surge/sway/heave = displacement of the centre of rotation in the
  earth frame (heave positive DOWN), roll positive starboard down, pitch
  positive bow up, yaw positive bow to starboard.  The rigid-body equations
  are Fossen's general form about an arbitrary origin (Eq. 2.89 / 2.91, with
  the CG offset and the full inertia tensor), solved directly (6x6) in the
  earth-frame representation so that fixed / imposed DOFs are simple rows.
* **Imposed motion** file `<case>_body<k>_imposed.dat`: rows `t x y z phi
  theta psi` ([s], [m], [deg]); cubic Hermite interpolation in time.
* **Moving grids**: the blocks of a moving object are moved rigidly from the
  original coordinates every step (metric vectors rotated, J unchanged),
  the advecting fluxes become relative fluxes U - Ug (ALE), the no-slip wall
  velocity is the grid velocity, the overset assembly is redone from the
  new coordinates (re-entrant `overset.F90` with caches of the motion-
  invariant parts: self hole-cut, surface trimming, same-object donor
  search), uncovered cells are filled before the fringe interpolation.
  The moved grid is written as `<case>_NNNNNNNN.x` with every `.q`.
* **Forces**: pressure and friction on the wall patches of every body
  (weighted by the surface-overset trimming so overlapping patches count
  once, hole cells skipped), moments about the centre of rotation.
  Files every `Data_Intvl` steps (`$Export_Interval`):
  `<case>_force_body<k>.dat`, `<case>_moment_body<k>.dat`, `<case>_force_earth<k>.dat`,
  `<case>_moment_earth<k>.dat` (friction / pressure / total [N], [N m]) and
  `<case>_motion<k>.dat` (position, velocity, acceleration), k = body number.
* `Grid_scale` (`$Numerical_Scheme`): multiplies the `.grd` coordinates
  (used for `Box_6DoF`: file units read as decimetres).
* **Cases**: `NACA66mod_Overset/` (imposed pitch 2 deg + heave 0.02 m,
  period 4 s), `cases/KP505/` (imposed rotation, J = 0.7), `cases/JBC/`
  (free heave + pitch = sinkage and trim), `Box_6DoF/` (two free-floating
  boxes, all six DOFs, two-phase).
* A body at rest reproduces the static solver to round-off; a cylinder grid
  rotated at 0.001 rad/s reproduces the static drag to 1e-6.
* **Solver robustness** (found on `Box_6DoF`, PROJECT_RECORD.md 24.3): inner
  GMRES instead of BiCGSTAB for the coupled overset pressure operator
  (two-phase overset systems no longer stall; the cylinder case needs
  84 instead of 237 iterations per step), hydrostatic-consistent pressure
  interpolation across the free surface, orphan fringe cells / faces closed
  from a same-block neighbour, prescribed Cartesian velocities on Dirichlet
  faces (no reconstruction loop at wall edges), donor promotion and a
  neighbour-average term in the fringe rows (no singular interpolation
  cycles).  `<case>.ovst` reports the orphans closed and the donor
  promotions.
* **Update 12c - no empty spaces in the assembly** (PROJECT_RECORD.md 24.4):
  the "active cells only" slices of Box_6DoF and KP505 showed fluid regions
  solved by no grid (tube corners under the refinement boxes, the layers
  between the blade grids and the refinement, the cut layer where two body
  grids meet).  `classify` now (a) keeps an active cell next to a hole
  active when its donor stencil holds no solved cell (the dominated holes
  around it become fringes; if only solid holes remain, they act as the
  wall - "solid-bounded"), (b) re-activates dominated holes that lie in no
  solved cell, (c) runs a coverage pass: every non-active fluid cell must
  have its centre, nodes, edge midpoints and face centres (27 points)
  inside solved cells of other grids, else it is made active; a FORCED
  cell (behind an overset face, never solvable)
  is covered instead by activating the finest solvable cell of another grid
  around each uncovered node.  Solid holes are never reopened: no cell
  inside a body is ever solved.  (d) The winding-number hole cut is now a
  fast winding number after Barill et al. (per-body cluster trees, dipole +
  second-moment far field, exact near field): the DelftTwist foil O-grid
  lost 1865 spurious "inside" cells, and the hole cut of GPPH is 8 times
  faster (assembly 243 s instead of 704 s).  (e) A body thinner than a
  cell (KP505 blade tips in background cells) now cuts the cells it passes
  through, and the coverage pass never activates a cell with a test point
  inside a body.  (f) A permeation check (grid-check mode, `Ovst_only = 1`)
  counts the active cells with any of 27 test points inside a body
  ("permeation check" line in the log / `.ovst`): no case has an active
  cell with its centre inside a body; the few cells with an edge / node
  point inside (hull cases, 0.01 - 0.03 %) sit at imperfect joints of two
  wall patches (centre plane of the ONRT / JBC hull halves, appendage
  roots) - see PROJECT_RECORD.md 24.4.
  New report lines "coverage :" and "solid-bounded :" (log and `.ovst`).
  Checked on all 18 overset case folders with `tools/blank_slice.py`
  (`blank_check/*.png` in every folder: y = 0 or z = 0 slices, blank status
  and active-only): no gap anywhere except the body interiors.  Per-case
  counts (coverage, solid-bounded cells, orphans, assembly time) in
  PROJECT_RECORD.md 24.4; `<case>.ovst` and `ovst_check.log` of every
  folder hold the full report.  Usage of the plot tool:
  `python3 tools/blank_slice.py <case_dir> y 0 [--zoom x0,x1,y0,y1]`
  (needs the `<case>.x` / `<case>_iblank.q` of a run with `Ovst_only = 1`;
  with the cell-centred `<case>_cblank.q` beside them - written by the
  solver with every `_iblank.q` and every solution `.q` as
  `<case>_NNNNNNNN_cblank.q` - the plot uses the exact cell status).
* **Blank output**: the node `[Blank_var]` of the solution `.q` marks a node
  −1 when any adjacent cell is a solid hole (thin bodies show); the exact
  cell status is in `<case>_NNNNNNNN_cblank.q` (cell-centred Plot3D q,
  dimensions nx ny nz).

## Update 11: surface tension

* `Surf_tension = 1` / `Sigma_ST` in `$VoF`: CSF surface tension (curvature from
  the smoothed volume fraction, balanced-force face gradient) for two-phase runs
  with any of the three VoF schemes; monitor with max curvature and the capillary
  time-step limit.  Also three latent defects fixed (output array bounds, cell-
  field classification of ids >= 14 in the block-to-block exchange, lower bounds
  of the gradient stencils).  PROJECT_RECORD.md section 23.

## Update 10: calculation-stop criteria

* `stop_criteria.F90`: the run stops with the grid block, index and position
  written to the log whenever a NaN appears in any field, whenever a residual
  |b-Ap|/|b| or the pressure update norm max|dp| exceeds `Stop_Rsd_Max` (1e9), or
  whenever the maximum local CFL number exceeds `Stop_CFL_Max` (1e4).
  PROJECT_RECORD.md section 22.

## Update 9: turbulence models, IDDES, wall function

* `TURB_model` (1 SA, 2 k-epsilon, 3 k-omega, 4 k-omega SST), `IDDES_model`
  (hybrid RANS/LES length scale on top of the chosen model, WMLES alone with
  model 0) and the nine turbulence-solver controls (`Max/Min_Itr_Turb`,
  `Rsd_Lmt_Turb`, `SOR_fct_Turb`, `Cnvt_scheme/mode_Turb`, `Dffs_mode_Turb`,
  `cfl/vn_switch_Turb`) in the name list before `$Overset`; module
  `turbulence.F90`; the eddy viscosity enters through `mu_h` (convection scheme
  unchanged).
* bcs flag 11 = no-slip wall + wall function (u, v, w all 11, otherwise the run
  stops); wall-function type per model, SST wall function by default.
* `NACA66mod_3D` (10 point-matched blocks, two-phase, wall function) and
  `DelftTwist` (3 overset blocks, single phase, wall function) case folders.
* `FlatPlate_WF`: verification case - SST + wall function at y+ 29-50, Cf within
  3 % of the Schultz-Grunow correlation over Re_x 2.5e5-2e6 (see its README).
  PROJECT_RECORD.md section 21.

## Update 8: divergence guards, active overlap, fringe reconstruction

* Pressure solver: breakdown / NaN guards in the inner BiCGSTAB and FGMRES
  (the sudden NaN at t ~ 1 s of the cylinder runs was a Krylov breakdown).
* Overset: forced fringe cells no longer dominate the other grid, and the
  active regions of overlapping grids overlap by one cell
  (`Ovst_active_overlap = 1`) - no blanked gaps in Tecplot plots.
* Geometric VoF: fringe reconstruction split per donor cell (no plane
  extension) - no spurious bubbles at start-up.  PROJECT_RECORD.md section 20.

## Update 7: geometric VoF

* `vof_scheme = 3` (name list `$VoF`): geometric PLIC volume-of-fluid
  transport - Youngs normal, Scardovelli-Zaleski plane constant,
  direction-split flux slabs with the Weymouth-Yue divergence correction
  (conservative, bounded for face Courant <= 0.5, automatic sub-cycling),
  ghost / overset refresh after every sweep, and GEOMETRIC reconstruction of
  alpha in the overset fringe cells (liquid fraction of the receiver cell on
  the donor's PLIC plane).  Schemes 1 and 2 (CICSAM) are unchanged and remain
  the default.  See `PROJECT_RECORD.md` section 19.
* Cylinder_Overset_1phase / _2phase grids: z distribution clustered at the
  free surface (dz 0.04 m at z = 0, node exactly at z = 0, same node count).

## Update 6: solid-body detection for overlapping surfaces

* Hole cutting rebuilt: overlap trimming (surface overset) -> generalized
  winding number of the trimmed, symmetry-mirrored body surface for the
  cut layer -> flood fill of every block from its non-wall faces (body
  interiors, appendage roots, duct penetrations).  `Surfovl_tol` no
  longer shrinks the solid.  A cell of a wall block inside the grid of
  another wall block whose wall is a different surface there (duct
  mouth, region above an appendage) is never solved by that block,
  regardless of the cell sizes; cells next to their own wall are never
  dominated.  See `PROJECT_RECORD.md` section 18 and `OVERSET_CASES.md`.
* Nodes inside a body are always holes in `<case>_iblank.q`; boundary-layer
  cells (`Ovst_wall_layers`, default 8, name list `$Overset`) stay active
  unless cut by another body; KCS deck-region refinement.
* `<case>_iblank.q`: a fringe node directly adjacent to an active node is
  written as 0.5, a fringe node adjacent to a 0.5 node as 0.25, one adjacent
  to a 0.25 node as 0.125 (1 active, 0 other fringe, -1 hole); see
  `tools/iblank_half_fringe.py` for the equivalent post-processor.
* Orphan accounting: active cells bounded by solid holes only (no other grid
  there) and fringe faces on covered block faces are not orphans; donor-less
  fringe faces next to solid holes become wall faces; `poi_verbose = 1` lists
  the first orphan faces with their positions.
* Cells inside the near-wall band (`Ovst_wall_layers`) of another wall block
  are holes, never fringe cells (always on): no background / refinement fringe
  cell touches a body's boundary layer.
* NSWC15E blocks 1-4 rebuilt in one step from the original grid so that the
  L-shaped block 2 / block 3 interface (i-min + k-min windows) is point-matched
  again and its corner region is active, not fringe; the interface detector now
  warns when two faces nearly coincide (within 1000 x `b2b_tol`) but are not
  connected.
* Block-to-block interfaces of the 12 cases repaired (common wall-normal
  point count per case), a start-up warning flags partly-b2b faces;
  block order of the case folders restored (background last); every
  case folder carries `<case>.x` and `<case>_iblank.q` from the grid
  check of this code.

## Update 5: overset orphan reduction, 12 overset cases

* **Solver change (overset.F90)**: nearest-valid-donor fallback.  A cell
  next to a hole whose containing cell in every other block is a solid
  hole (the coarse background cell straddling another body's surface)
  had no donor at all and stayed an orphan; it now takes the nearest
  non-solid cell of the finest other block within two index layers as a
  nearest-neighbour donor (never used for the domination test), and
  fringe faces get the same fallback on the nearest non-hole cell.
  Across the 12 overset cases the orphan counts drop by 96 % (cells) /
  98 % (faces); see `OVERSET_CASES.md` for the per-case numbers.
* **Solver fix (curvi.F90, `build_metrics`)**: degenerate metrics at
  corner-ghost, wedge-edge and boundary node/face positions no longer
  abort the run (they are filled from the nearest interior cell); only a
  degenerate interior cell centre is treated as a grid defect.  All ship
  grids (two block faces on the symmetry plane) needed this to start.
* 12 overset case folders added (`cases/`): KCS, DTMB, GPPH, DCWJ, ONRT,
  JBC, KP505, NUMRT, LCS, NSWC15E, JUSKFoil, SPSCFoil - each with grid,
  bcs, name list, csv and the assembly report of the new code.

## Update 4: tutorial cases, PLOT3D-ASCII grids only

* The grid file is read as PLOT3D-ASCII (`<case>.grd`) only; the binary
  `.xyz` input of update 3 was removed (`tools/cfdship2probcfd.py`
  converts binary CFDShip-Iowa grids to ASCII when needed).
* Seven tutorial case folders replace the earlier case folders (which
  are in the previous deliveries): cylinder single-grid / overset, each
  single-phase and two-phase, and the NACA66(mod) section on a C-grid,
  an O-grid and an overset grid.  Every case: free stream 1 m/s along
  +x, 2nd-order upwind convection, towing-tank start.  Each folder holds
  the grid, the final `.bcs`, the csv, the annotated `probCFD.nml`, a
  README and `monitor_reference.dat` (the first 10 steps).
* `tools/make_cases.py` (writes the seven `.bcs` / csv files from the
  supplied grids), `tools/compare_naca.py`, `tools/compare_cyl.py` (the
  cross-grid comparisons of section "Verification" in
  `TUTORIAL_CASES.md`).
* **Solver fix (boundary.F90, `apply_bc`)**: at a zero-gradient /
  outflow boundary of the face-normal velocity family the STAR
  (predictor) field had its boundary face overwritten by the first
  interior face.  The momentum equation is solved at that face, so the
  projection had to remove the difference between the extrapolated and
  the momentum-consistent flux every step - at the cost of
  rho du dx / (2 dt) in the outlet cell.  On grids with large outlet
  cells this left a uniform pressure offset of kPa size and a radial
  pressure oscillation in the last cells (the C-grid: -80 kPa level with
  4.4 m outlet cells; the overset background: -6.5 kPa; the O-grid:
  -1 kPa), i.e. the "pressure level" and "wiggle" problems of the NACA
  tutorial cases.  The boundary face is now left to the momentum
  solution (+ projection / mass-flux scaling) and both ghost faces are
  copied from it.  The star outlet flux now matches the inflow to
  3e-4 instead of 3-5 %.  The BFS regression (previous deliveries)
  changes only in the last digits of |div|max and by a few percent in
  the iteration counts (its outlet cells are tiny).
* **Solver fix 2 (boundary.F90)**: the Dirichlet ghost values of the
  tangential velocity families were mirrored in TRANSFORMED space with
  the extrapolated ghost metrics; on curved / stretched boundaries (the
  O-grid outer circle, the C-grid arc) the Cartesian ghost velocity was
  off by percents (0.975 instead of 0.9998 on the O-grid inlet), which
  is convected into the first cell and answered by the projection with
  a pressure wiggle.  The mirror now acts on the Cartesian velocity
  vector (u_g = 2 u_bc - u(b), assembled at the interior face as in
  cart_to_contra) and is projected with the ghost metrics.
* **Solver fix 3 (convection.F90, `face_setup`)**: the advecting volume
  fluxes at the half positions were J_cell x mean(U^q) of the two
  faces; on stretched grids that is not the mean of the two face
  fluxes continuity is built from (6 % off in the first inlet cell,
  where the ghost spacing is a copy of the first cell), so a uniform
  inflow got a spurious convective source in the first cell (the
  "too low inlet pressure" of the overset case: -111 Pa in the first
  background cell).  Now the mean of the two face fluxes J_f U_f/fvw_f
  is used: a uniform field convects nothing to machine precision.
  The three fixes together move the BFS regression in the 3rd digit of
  the CFL and |div| columns; its reference lines with this solver
  (np = 2, FGMRES + line SGS, Cnvt_scheme 2) are now:

          1   2.0000E-06  0.0218  0.0199  0.0000  0.0092  5.182E-12  2.905E-04  1.122E-05    5 4.09E-10  1119 1.84E-08 5.92E+03
          2   4.0000E-06  0.0202  0.0171  0.0000  0.0074  7.233E-13  9.064E-04  3.587E-05    5 9.12E-10  1204 3.71E-08 5.89E+03
          3   6.0000E-06  0.0217  0.0166  0.0000  0.0077  2.699E-14  3.809E-04  1.283E-05    5 2.87E-10   878 2.09E-08 7.53E+01
          4   8.0000E-06  0.0225  0.0162  0.0000  0.0083  3.730E-15  3.771E-04  1.444E-05    5 4.07E-10   768 2.75E-08 1.59E+01
          5   1.0000E-05  0.0229  0.0159  0.0000  0.0086  1.805E-15  1.559E-04  5.088E-06    5 4.43E-10   750 1.47E-08 8.35E+00

## Update 3 (previous delivery): OVERSET (Chimera) grids

Sources: `ovs_geom.F90` (geometry kernels), `overset.F90` (assembly and
run-time fringe exchange) + hooks in `mod_global`, `read_input`,
`bcs_b2b`, `b2b_exchange`, `curvi`, `convection`, `u_momentum`,
`poisson`, `vof`, `output`, `Main` (15 files now; compile order = the
SRCS list).  New cases `NACA66mod_Ovst_case/` and `SPSCFoil_case/`;
the former is (3-block overset grid of the foil section:
foil O-grid + refinement box + background, all overset-connected, water,
1 m/s along +x).  Every earlier case is unchanged (the BFS regression is
bit-identical; the overset code is inactive when all block faces are
covered).

1. **Overset boundary faces.**  A block face (or a window of it) that
   is neither in `<case>.bcs` nor a block-to-block / periodic / pole face
   of `<case>.bcs_b2b` is an overset boundary: the block ends inside
   another block.  Detected automatically (`SW_Ovst = -1`); `1` forces
   the assembly, `0` switches it off (with a warning if such faces
   exist).  Automatic zero-gradient patches keep the ghost values behind
   such faces finite.
2. **Solid bodies from the no-slip walls.**  A `.bcs` patch with
   Dirichlet u = v = w = 0 is a wall.  Walls are triangulated (normals
   into the fluid) and grouped into bodies through their boundary edge
   nodes: an edge node is connected when it lies within `Surfovl_tol`
   of another wall patch (closed surface WITH SURFACE OVERLAP, hydrofoil
   / KCS type, task 4.1), on a block-to-block interface of its block
   (the wall continues in the neighbour block, e.g. the O-cut), or
   within `Surfovl_tol` of a physical boundary patch of any block
   (symmetry plane / free surface / end of the background grid, NSWC15E
   and NACA66mod type, task 4.2).  A body without unconnected edge
   nodes is CLOSED.  `<case>.ovst` reports the counts per body.
3. **Hole cutting.**  All cells of all blocks whose centre or any node
   lies inside a closed body are HOLES (closest-point test on the
   triangulated walls with angle-weighted pseudo normals: robust for
   overlapping patches and sharp trailing edges).  Open surfaces cut no
   holes (reported).
4. **Active / fringe / hole classification** ("implicit hole cutting"):
   the centre of every cell is located in every other block (bucket
   search + Newton inverse trilinear map); the finest containing
   non-hole cell is its donor; a cell is dominated when its donor is
   finer (relative volume tolerance `Ovst_vol_tol`, ties keep both
   active -> a symmetric grid gives a symmetric blanking, task 5).
   ACTIVE = not dominated, not within `Ovst_fringe` (= 2) layers of an
   overset face, not hole.  FRINGE = non-active non-hole cells within a
   Chebyshev index distance 2 of an active cell (the 2nd-order upwind /
   QUICK / ENO stencils) plus every donor-stencil member of a fringe;
   a dominated fringe whose donors are boundary fringes of the other
   block becomes active (no fringe-from-boundary-fringe cycles); the
   rest are HOLES.  Result on the NACA66mod_Ovst grid: the SUGGAR-like
   pattern of the CFDShip-Iowa reference (foil grid blanked in its outer
   layers with a 2-layer fringe band, foil-shaped hole in the
   refinement block, rectangular hole in the background), no orphans.
5. **Coupling in the solver.**  Faces with an active neighbour cell are
   solved; faces with only fringe neighbours are FRINGE FACES whose
   Cartesian velocity vector is interpolated from the donor block's
   staggered face lattices (trilinear in the donor's computational
   space; u from i-faces, v from j-faces, w from k-faces) and projected
   onto the transformed component with the local metrics - in every
   momentum sweep, after the correction and at start-up.  Fringe cells
   get p (implicitly) and alpha (explicitly) from the donor cell-centre
   lattice.  The pressure equation is solved as ONE coupled system: the
   fringe rows are `p_f - I p = 0` inside the FGMRES operator, and the
   FGMRES preconditioner is an inner BiCGSTAB on the coupled 7-point +
   interpolation operator (`Ovst_inner_max`, `Ovst_inner_tol`), block-
   Jacobi preconditioned (line SGS / Jacobi per sub-domain, Dirichlet at
   the fringes).  Hole cells and faces hold zero and are skipped
   everywhere (momentum, projection, VoF, monitors).
6. **Output.**  The node iblank (1 active, 0 fringe, -1 hole; the
   CFDShip-Iowa convention) is appended as the LAST variable of every
   `.q` file when the overset assembly is active (nvar 5 or 6);
   `<case>.ovst` is the assembly report; `Ovst_only = 1` stops after
   the assembly (grid check mode).  Assembly time on the NACA66mod_Ovst
   grid: 1.3 s (replicated on every rank, deterministic).
7. **Name list** group `$Overset`: `SW_Ovst` (-1/0/1), `Surfovl_tol`
   [m], `Ovst_fringe` (2..3), `Ovst_vol_tol`, `Ovst_only`,
   `Ovst_inner_max` (150), `Ovst_inner_tol` (1e-2), `poi_verbose`
   (FGMRES / BiCGSTAB history).  6-DoF grid motion (re-assembly per
   step) is future work: the assembly is a single `ovst_setup` call
   whose inputs are the node coordinates.
8. **Grid-check mode**: `Ovst_only = 1` runs the assembly BEFORE the
   field arrays are allocated (an 11 M cell grid fits in ~1.5 GB) and
   writes `<case>.ovst`, `<case>.x` and `<case>_iblank.q` (one
   variable: 1 active, 0 fringe, -1 hole).  Reading the reference
   NACA66mod grid (converted with the tool below), the blanking agrees
   with the SUGGAR iblank of the reference solution on 99.6 / 97.6 /
   97.0 % of the nodes of its three blocks.  The grid file is read as
   PLOT3D-ASCII (`<case>.grd`) only.
9. **Tool** `tools/cfdship2probcfd.py`: converts a CFDShip-Iowa binary
   PLOT3D grid + `.bcs` (type 20 = no-slip) into the probCFD `.grd` /
   `.bcs` format, optionally coarsened; `tools/hydrofoil_assembly_test/`
   = the assembly-only test on the coarsened hydrofoil example (7
   blocks, 6 OVERLAPPING wall patches -> one closed body, task 4.1);
   `tools/readp3d.py`, `tools/readq.py` = Python readers of the example
   files and of the solver's `.x` / `.q` output (iblank = last variable).

## Update 2 (previous delivery): cleanup, log spacing, QUICK / ENO, pole boundary, PoleCyl case

1. `N_max` and `Initial_Value_fct` are removed: the run length is
   `End_T / Del_T` steps, the pressure solve always starts from the
   previous pressure.
2. Seven empty lines separate the time-step blocks of the screen log.
3. Convection schemes, `Cnvt_scheme` (renamed from `conv_scheme`):
   1 = 1st-order upwind, 2 = 2nd-order upwind, 3 = QUICK, 4 = 2nd-order
   ENO.  All four work with `Cnvt_mode` 1 (explicit), 2 (implicit) and
   3 (hybrid per direction): in the explicit directions the scheme's
   face value enters the explicit convection term, in the implicit
   directions the 1st-order upwind part is in the matrix and the
   (scheme - 1st order) difference is the deferred correction from the
   latest iterate (`Iner_Itn = 2` converges it).  Verified on the BFS
   case (all 12 scheme x mode combinations) and the two-phase cylinder.
   Stability notes printed by the reader: explicit 2nd-order upwind /
   QUICK / ENO need `cfl_switch <= 0.3` in hybrid mode; QUICK with
   forward-Euler explicit directions is unstable for pure convection
   (use `Cnvt_mode = 2`).  The ENO scheme is the conservative flux form
   `phi_f = phi_D + 1/2 minmod(phi_A - phi_D, phi_D - phi_U)` (the
   smoother of the central and upwind slopes; the derivative form of
   the attached formula is not conservative and needs a 6-point stencil,
   see PROJECT_RECORD 13.3).
4. Pole boundary condition: the interface search tests every block
   face for degeneracy (all nodes along one index coincide within
   `b2b_tol`); such a face is a pole (axis line or point), excluded
   from the interface matching and written to `<case>.bcs_b2b` as a
   type-93 record next to the type-92 interfaces.  The solver zeroes
   the face Jacobian / metrics on the pole plane (zero area: no mass,
   momentum or pressure flux), adds an automatic zero-gradient patch
   for all variables, and averages the azimuthal node copies in the
   output.  Nothing is needed in the `.bcs` file.
5. `PoleCyl_case/`: single-block pipe, 1 m/s along +i, water, no
   gravity — see its README.

## Update 1 (earlier delivery): the "probCFD_부산 Update Task List"

1. **Inner iterations** `Iner_Itn` (`$Solver_Iteration`): the momentum
   predictor, the Poisson solve and the correction (and, two-phase, the
   VoF transport from alpha^n) are repeated `Iner_Itn` times per time
   step, with the advecting fluxes, the explicit terms and the 2nd-order
   deferred correction taken from the latest iterate.  `Iner_Itn = 1`
   is the classic fractional step; 2 gives a Picard-converged, fully
   implicit (backward Euler) step (the pressure of iteration 2 starts from a converged field and the
   divergence after the step drops 2-3x on the BFS case, ~3x on the
   two-phase cylinder; see `regression_BFS/README.md`).
2. **Convergence criteria**:
   * momentum predictor: true algebraic residual
     `||b - A u*||_inf / ||b||_inf < Rsd_Lmt_Vel`, `Min/Max_Itr_Vel`,
     weighted Jacobi `SOR_fct_Vel`; max and average residuals reported
     per component u, v, w;
   * VoF (`vof_scheme = 2`, implicit CICSAM): true residual of the
     scaled cell equation `max |alpha^n - alpha - dt/V sum F alpha_f| <
     Rsd_Lmt_VoF`, `Min/Max_Itr_VoF`, `SOR_fct_VoF`;
   * pressure: BOTH the update norm `max|p_new - p_old|` and the
     residual `||b - A p||_inf / ||b||_inf` are computed and printed
     every step; `poi_solver = 1` (Jacobi) stops on the update norm
     (`Trl_Lmt_Prs`), `poi_solver = 2` (FGMRES + PCG) on the residual
     (`Rsd_Lmt_Prs`); `Min/Max_Itr_Prs`, `SOR_fct_Prs` (Jacobi),
     (`Initial_Value_fct` and `N_max` were removed again in update 2).
3. **Renewal upgrades**:
   * junction ghost-cell ambiguity: the solution writer now averages
     every field over all block copies of each coincident block-boundary
     node (coincidence groups from the node coordinates, tolerance
     `b2b_tol`), so u, v, w, p (and alpha) are single-valued on the O/H
     interface INCLUDING the three-block junction lines (verified: zero
     difference between block copies on all 6580 shared nodes of the
     cylinder grid);
   * line-relaxation preconditioner `poi_precond = 2` (line symmetric
     Gauss–Seidel with tridiagonal solves along the strongest-coupled
     direction of each block) for the PCG inside `poi_solver = 2`:
     2670 -> 411 CG iterations on the first two-phase cylinder step,
     2000+ (cap) -> 766-1122 on the BFS case (details in PROJECT_RECORD 12.3).
4. **Name list in the reference format** (`$Case_title`, `$Control`,
   `$Fluid_Properties`, `$Reference_Pressure`, `$Cell_Initial_Condition`,
   `$Solver_Iteration`, `$Export_Interval`, each closed by `/end`) with
   the reference variable names (`Case_name`, `Char_L`, `Char_U`, `Del_T`,
   `End_T`, `FreeSurf_Z`, `Grvty`, `Num_of_Phase`, `Lou_w`, `mu_w`, `Lou_a`,
   `mu_a`, `SW_Ref_Prs`, `Loc_X/Y/Z`, `RefPrs`, `u/v/w/a/p_cel_initial`,
   `Iner_Itn`, `Max/Min_Itr_Vel`, `Rsd_Lmt_Vel`, `SOR_fct_Vel`,
   `Max/Min_Itr_Prs`, `Trl_Lmt_Prs`, `SOR_fct_Prs`,
   `Max/Min_Itr_VoF`, `Rsd_Lmt_VoF`, `SOR_fct_VoF`,
   `Output_Intvl`) plus the probCFD-specific `Rsd_Lmt_Prs`, `Log_Intvl`
   and the `$Numerical_Scheme` group (discretization / solver switches).
   Unknown or misspelled names stop the run with the line number.
5. **Log output in the reference layout**: banner + code name
   `probCFD_부산` first, the `PC_# / BLK_# / i_th j_th k_th / node ranges`
   splitting table, and per time step: current time, `Max Local Total /
   <i> / <j> / <k> CFL#`, `Cell abs(DIV) {min/max/avg}`, the `[ I ] [ J ]
   [ K ] Face Monitor` blocks (sweeps, max / avg residual), the
   `[ Cell ] Center Monitor : Pressure - Inner Itr#` blocks (iterations,
   residual, update norm), the `CICSAM-VoF` monitor (two-phase) and the
   wall-clock time of the step with its VoF / momentum / Poisson split.
   A one-line-per-step history is also written to `<Case_name>_monitor.dat`.
6. **NACA66mod case** (`NACA66mod_case/`), see its `README.md`.

Also in this update: the momentum sweep assembles its coefficients once
per solve (3x faster predictor); the 2nd-order deferred correction is
frozen during the sweeps (the previous in-sweep update diverged for
CFL > 1/2 — `Cnvt_scheme = 2` now runs on the two-phase cylinder);
`.bcs` files may have 4 or 5 lines per face for any `Num_of_Phase`.

## Interpretation of one reference variable
* `a_cel_initial` is the initial volume fraction of the cells ABOVE
  `FreeSurf_Z` (the cells below are liquid, alpha = 1); 0 = air.

## Log output example (two-phase, one step)

     Step        2   Current Time(real) :             0.0100000000 s
     Max Local Total CFL#        :                                7.0201
     Max Local <i> CFL #         :                                0.0457
     Max Local <j> CFL #         :                                7.0128
     Max Local <k> CFL #         :                                0.1098
     Cell abs(DIV) {min/max/avg} :    4.52E-15   9.19E-07   1.28E-08

     [ I ] Face Monitor :          < u Vel >
      --> TTL # of Itr  :              46
      --> Max Residual  :          7.1400E-10
      --> Avg Residual  :          3.2100E-12
     ... [ J ], [ K ] ...
     [ Cell ] Center Monitor : Pressure - Inner Itr#:            1
      ---> Total Number of Iteration                             363
      ---> Pressure Residual |b-Ap|/|b|                 1.810000E-09
      ---> Pressure Update norm max|dp|                 2.440000E+00

     [ Cell ] Center Monitor :    < CICSAM-VoF >
      --> TTL # of Itr  :             180
      --> Max Residual  :          4.5800E-10
      --> Avg Residual  :          1.2000E-12
      --> sub-cycles    :              15

     Cell abs(DIV) {min/max/avg} :    0.00E+00   3.21E-08   4.94E-10  (after correction)
     Liquid volume [m3] :   1.90436E+02   alpha min/max :  0.00000  1.00000   clipped [m3] :  4.175E-09 (cumulative)

     ClockTime[s] for ONE time step approach:           23.394 s
       (VoF / momentum / Poisson+correction :      0.460    12.383    10.505 s)

The CFL numbers are the computational ones, `dt |U^q|` per grid
direction (unit index spacing), evaluated on the field the step starts
from; `Cell abs(DIV)` is printed in 1/s, before and after the step.
