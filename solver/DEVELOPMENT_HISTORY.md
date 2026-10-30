# probCFD_부산 — Development History

*A summary of the code development up to update 13 (2026-09-30). The full technical record, including every formula, test number and fix, is in `PROJECT_RECORD.md` (sections 1–25). Usage is in `README.md`.*

---

## 1. What the code is today

**probCFD_부산** is a 3-D incompressible Navier–Stokes solver written in Fortran 90 with MPI. It is built for ship and marine hydrodynamics and has these parts:

| Area | Capability |
|---|---|
| Discretisation | Finite differences on **fully staggered, structured, multi-block grids** with a general 3-D curvilinear transformation (half-index metrics). |
| Grids | Multi-block PLOT3D-ASCII input. Block-to-block interfaces are **detected automatically** (all 24 orientations, partial faces, several windows per face, pole faces). Blocks are split across MPI ranks through a csv file. |
| Overset (Chimera) | Hole cutting of closed solid bodies, fringe/donor search, orphan handling, and gap-free coverage with no solved cell inside a body. |
| Flow | Projection (fractional-step) method; Jacobi momentum predictor; FGMRES + PCG pressure solver with an inner GMRES for overset systems. |
| Convection | 1st/2nd-order upwind, QUICK, 2nd-order ENO. Explicit, implicit or hybrid per term and per direction. |
| Time | Inner iterations; 1st-order scheme or **2nd-order SSP (Heun)**. |
| Two-phase | Air–water VoF: explicit CICSAM, implicit CICSAM or **geometric PLIC**. Gravity, surface tension (CSF), and dynamic pressure output. |
| Turbulence | Spalart–Allmaras, k-ε, k-ω, k-ω SST; IDDES; wall functions. |
| Motion | **6-DoF rigid-body motion** (SNAME/Fossen), imposed or free per DOF, for any number of bodies, with dynamic overset re-assembly each step. |
| Robustness | Stop criteria (NaN, residual, CFL), a `Reason_for_stop.md` report, and name-list changes applied at run time. |

Quantities are non-dimensionalised internally by `Char_L`, `Char_U` and the water density. All inputs and outputs are dimensional.

---

## 2. Source structure (32 modules)

| File | Role |
|---|---|
| `mod_global.F90`, `mod_mpi.F90` | Global data (blocks, sub-domains, fields), MPI setup, domain assignment |
| `stop_criteria.F90` | Stop criteria and the `Reason_for_stop.md` report |
| `bcs_b2b.F90`, `b2b_exchange.F90` | Detection of block-to-block interfaces (`<case>.bcs_b2b`) and the unified ghost-exchange engine |
| `ovs_geom.F90`, `overset.F90` | Geometry kernels (buckets, point location, surfaces, fast winding number) and the overset assembly |
| `curvi.F90` | Curvilinear metrics and conversions between Cartesian and contravariant velocities |
| `sixdof.F90` | 6-DoF bodies, child and sub-child bodies (update 14): forces, equations of motion (relative, composite mass), grid motion, output files in `Output_data/` |
| `statistics.F90` | on-the-fly turbulence statistics on the grid and on user-defined lattices (update 23) |
| `fsi_shell.F90`, `fsi.F90`, `fsi_test.F90` | fluid-structure interaction (update 24): co-rotational shell solver (stand-alone, tested by `fsi_test`), mid-surface from the wall surface grid, load / motion transfer, mesh morphing, strongly coupled sub-iterations |
| `surface_output.F90`, `probes.F90` | solid-wall surface stresses and files, cell shear stress (update 22); volume-field probe points (update 22) |
| `pid_control.F90`, `force_model.F90` | PID controllers (update 21: throttle / steering / direct force on grid parts and force models) and force models (update 19: Hough-Ordway / BEM disks, point forces) |
| `waves.F90`, `wave_fields.F90`, `energy.F90` | incident waves (update 18), energy equation (update 16) |
| `restart.F90` | restart files `.rsto` (update 13b) |
| `cavitation.F90` | cavitation models: vapour fraction sources, continuity source, monitor (update 15) |
| `read_input.F90` | Name list (start-up and run-time re-read), csv, grid, `.bcs` |
| `vof.F90`, `surface_tension.F90` | Volume-fraction transport and CSF surface tension |
| `turbulence.F90` | RANS/IDDES models, wall distance, wall functions |
| `boundary.F90`, `convection.F90`, `u_momentum.F90`, `poisson.F90` | Boundary conditions, convection and diffusion, momentum predictor, pressure solver |
| `output.F90` | PLOT3D grid and solution export, node and cell blank files |
| `Main.F90` | Start-up, the time-marching loop, logging |

---

## 3. Timeline of the development

| When | Update | Main content |
|---|---|---|
| ≤ 2026-07-23 | Earlier rooms (M1–M7) | b2b interface detection, multi-block solver, csv MPI splitting, periodic BCs, 3-D curvilinear transformation, backward-facing-step case, `curvi_form` |
| 2026-09-02 | Two-phase room (M8–M12) | `Cnvt_mode`/`Dffs_mode` time integration, CICSAM VoF, gravity, consistent projection at boundaries and interfaces, FGMRES+PCG, BC ramp `t_ramp`, surface-piercing cylinder case |
| 2026-09-03 | Update 1 | Inner iterations, residual-based convergence criteria, dual pressure norms, single-valued node output, line-SGS preconditioner, grouped name list, reference log format, implicit CICSAM, NACA66mod case |
| 2026-09-04 | Update 2 | QUICK and ENO convection, pole (degenerate-face) boundary, PoleCyl case |
| 2026-09-05/06 | Update 3 | **Overset (Chimera) grids**, SPSCFoil case, the `Ovst_only` grid-check mode |
| 2026-09-07/08 | — | ASCII-only grid input (a converter tool for CFDShip grids); **outlet pressure offset fix** (the boundary face is no longer overwritten) |
| 2026-09-13 | Update 5 | Overset orphan reduction (nearest-valid-donor fallback); start-up checks of 12 ship and foil cases |
| 2026-09-15…20 | Update 6 (6b–6h) | **Solid-body detection rebuilt**: surface trimming, winding number, flood fill; symmetric blanking; boundary-layer holes; fringe layers in the blank file |
| 2026-09-20/21 | Update 7 | **Geometric PLIC VoF** (`vof_scheme = 3`), geometric fringe reconstruction of alpha |
| 2026-09-21 | Update 8 | Krylov divergence guards, active overlap, gap fixes in the blanking |
| 2026-09-22/23 | Update 9 (9b, 9c) | **Turbulence models, IDDES, wall functions**, flat-plate verification, NACA66mod_3D and DelftTwist cases |
| 2026-09-23 | Update 10 | **Stop criteria** (NaN, residual/update norm, CFL) with the location of the offending cell |
| 2026-09-23 | Update 11 (11b) | **Surface tension (CSF)**, initial drop shapes, drop test |
| 2026-09-28 | Update 12 (12b) | **6-DoF motion, dynamic overset**, force/moment/motion files, Box_6DoF, KP505 and JBC setups, solver robustness fixes |
| 2026-09-29/30 | Update 12c | **Gap-free overset assembly with no solved cell inside a body**, fast winding number, permeation check, exact cell blank, self-overlap donors |
| 2026-09-30 | Update 13 | **2nd-order SSP time scheme**, dynamic pressure, run-time name-list updates, `Reason_for_stop.md`, new log layout, `Volume_solutions_grids/` with `.grd`/`.p3d` |
| 2026-10-01 | Update 13b | Tecplot variable names (`.nam` file), **restart** (`Mode`, `Restart_Intvl`, `.rsto`) |
| 2026-10-01 | Update 14 | **Child and sub-child bodies of the 6-DoF motion** (hull > stern drive > propeller), relative motion, composite-mass coupling, `Output_data/`, GPPH hierarchy case, `tools/make_axis_rotation.py` |
| 2026-10-02 | Update 14b | **`Mode = 2` grid-movement check** (`Check_grid_movement/`), stops for a missing name list / grid / impose motion, GPPH with steering drive + turning propeller |
| 2026-10-02 | Update 15 | **Cavitation models** (Schnerr-Sauer, Zwart-Gerber-Belamri, Kunz; `cavitation.F90`), three-phase water/air/vapour, `Vapor_scheme`, Reboud correction, NACA66mod_cav and DelftTwist_cav cases |
| 2026-10-03 | Update 16 | **Energy equation** (`energy.F90`), thermal wall function, buoyancy, latent heat |
| 2026-10-04 | Update 17 | Compressible liquid-vapour mixture, thermal cavitation coupling |
| 2026-10-05 | Update 18 | **Incident waves** (`waves.F90`, `wave_fields.F90`): regular / irregular, flag 12, relaxation zones, probes, RegularWave_Tank |
| 2026-10-06 | Update 19 | **Force models** (`force_model.F90`): Hough-Ordway / BEM disks, point forces, velocity-Verlet half-kick |
| 2026-10-06/11 | Update 20 | **Moving background** (`Moving_BGR`, background follows a Body), `Relax_L` per face, `Ovst_Reassemble_Tol`, discrete ALE grid velocity, restart / probe / name-list fixes, Cylinder_MovingBGR and WaveTank_MovingBGR cases |
| 2026-10-12 | Update 21 | **PID controllers** (`pid_control.F90`): throttle / steering / direct-force controllers on grid parts and force models, sequential targets, output and restart state; PID_Steer, PID_Throttle, PID_Direct, PID_Disk, PID_DiskSteer cases |
| 2026-10-19 | Update 24 | **Fluid-structure interaction** (`fsi_shell.F90`, `fsi.F90`): flexible Body / Child / Sub-child parts (`SW_FSI = 1`) with a geometrically non-linear shell built from the part's own wall surface grid (thickness by ray pairing, clamped NoFSI groups, Submat regions), mesh morphing of the part's blocks, strongly coupled sub-iterations; CSM3 structural test, FSI-PfS-2a benchmark (cylinder + rubber plate + rear mass), KP505_FSI run-through case |
| 2026-10-28 | Update 24b/24c | **FSI verification work** on FSI-PfS-2a: level-n restart of the coupling sub-iterations (fluxes + Cartesian field, old metrics kept through the step, full pressure kept as the Poisson guess), scaled IQN, rigid morphing columns, structural sub-step fallback; fix of the re-projection's old metrics (move_grids overwrote them with the undeformed ones); `SW_Quasi_2D`, `FSI_Plane_Strain`, `FSI_Test_Motion` (incl. the pluck start-up, `fsi_shell.shell_mode`); slab end caps in the hole cutting; one-body elliptic grid of the benchmark (convex outer boundary, seam at the cylinder front, no sliver column); rebuilt after the workspace reset of 2026-10-28 |
| 2026-10-13 / 10-18 | Update 23 | **On-the-fly turbulence statistics** (`statistics.F90`): mean / rms / Reynolds stresses and the gradient group (ε, η, λ, Re_λ) on the original grid and on boxes / cylinders / planes following a frame, blanked steps skipped, restart; 23b: the gradient from the face differences of the staggered velocities (accumulated at the cells); Statistics_ActuatorDisk, Statistics_TowedCylinder cases |
| 2026-10-13 | Update 22 | **Probe points** (`probes.F90`: background- and part-following points, one time history per point, NoData in holes), **Shear_Stress** variable of the volume files (|τ_w| at the wall nodes), **wall surface files** (`surface_output.F90`: pressure / friction stress vectors whose trapezoidal integral equals the force files); Probes_TowedCylinder, Probes_SteeredFoil cases |

---

## 4. The development in detail

### 4.1 Foundations (earlier rooms)

- **Block-to-block interfaces** are detected automatically and written to `<case>.bcs_b2b` as type-92 records with affine index maps. Pole faces use type-93 records.
- A **single ghost-exchange engine** handles intra-block halos, b2b patches and periodic faces. Results do not depend on the number of MPI ranks.
- **Curvilinear formulation**: contravariant face velocities U¹, U², U³ are the exchanged variables. Half-index metrics are used, and `curvi_form` selects an unweighted or a J-weighted flux.
- The **projection was made exactly consistent** at non-orthogonal boundaries, pressure outlets and multi-block interfaces (a flux-form Poisson operator with synchronised interface Jacobians). This fixed a divergence that grew at the boundaries.

### 4.2 Solver controls and schemes (updates 1–2)

- `Iner_Itn` inner iterations; residual-based stop criteria for the momentum, pressure and VoF solvers.
- `Cnvt_scheme` 1–4 (upwind 1/2, QUICK, ENO) and `Cnvt_mode`/`Dffs_mode` (explicit, implicit, hybrid).
- A grouped name list (`$Control`, `$Fluid_Properties`, …) and a fixed log format.
- Node output made single-valued at block interfaces.

### 4.3 Overset grids (updates 3, 5, 6, 8)

- **Overset boundaries**: any block face that is not in `.bcs` or `.bcs_b2b` is an overset boundary. The cells behind it are forced fringes.
- **Solid bodies**: the no-slip patches are grouped into closed bodies. Closure uses the tolerance `Surfovl_tol`, which covers overlapping patches (hydrofoil, KCS) and wall edges that end on the domain boundary (NSWC15E, NACA66mod).
- **Hole cutting (update 6)**: triangles hidden by a finer overlapping patch are trimmed, the cut layer is decided by a winding-number test, and a flood fill from the non-wall faces finds body interiors, appendage roots and duct penetrations. Symmetric grids give symmetric blanking.
- **Boundary-layer holes**: cells inside another wall block's near-wall band are holes, never fringes.
- **Orphans** are reduced by a nearest-valid-donor fallback and closed from a neighbour in the same block.
- **`Ovst_only = 1`** runs the grid check only: the `.ovst` report plus the grid and blank files, without allocating the flow fields.
- Assembly of SPSCFoil was compared against SUGGAR on NACA66mod: the node blank is identical on 97–99.6 % of the nodes.

### 4.4 Two-phase flow (two-phase room, updates 7, 11)

- One-fluid formulation with a balanced-force pressure gradient; alpha boundary conditions; the liquid volume is monitored.
- VoF schemes: explicit CICSAM (sub-cycled), implicit CICSAM (iterative), and **geometric PLIC** (direction-split, Weymouth–Yue) with geometric reconstruction of alpha in overset fringes.
- **Surface tension (CSF)**: curvature from a smoothed alpha field, balanced-force face gradient, and a capillary time-step monitor.

### 4.5 Turbulence (update 9)

- SA, k-ε, k-ω (Wilcox) and k-ω SST, each with an optional IDDES length scale (WMLES alone with model 0).
- Wall function via `.bcs` flag 11; exact wall distance; verified on a flat plate (Re 10⁵–10⁸).
- New cases: NACA66mod_3D (10 point-matched blocks, two-phase) and DelftTwist (overset, wall function).

### 4.6 Robustness (updates 8, 10, 12)

- **Krylov guards**: breakdown detection in BiCGSTAB and FGMRES, and NaN-safe back substitution.
- **Stop criteria** (update 10): NaN, residual or update norm above `Stop_Rsd_Max`, CFL above `Stop_CFL_Max`. Each reports the block, index and position of the offending cell.
- **Box_6DoF fixes** (update 12):
  - an inner GMRES replaced BiCGSTAB for the two-phase overset pressure operator;
  - the pressure interpolation across the free surface is hydrostatic-consistent;
  - orphan cells and faces are filled from a same-block neighbour;
  - Cartesian velocities of Dirichlet faces are prescribed directly, which removed a ×20 amplification per sweep at wall edges;
  - donor promotion and a neighbour-average term remove singular fringe cycles.

### 4.7 6-DoF motion (update 12)

- **Convention: SNAME/Fossen.** The earth frame has x forward, y to starboard and z down (it is the grid frame with x and z reversed). The body frame is centred at the centre of rotation. Euler angles are Z-Y-X. Surge, sway and heave are earth-frame displacements; u v w and p q r are body-frame velocities.
- The equations are Fossen's rigid-body equations about an arbitrary origin, with the CG offset and the full inertia tensor, solved as a 6×6 system.
- **Inputs**:
  - row 6 of the csv, `Object_ID`, assigns every block to a body;
  - a `$SixDOF` group per body: `Mass`, `CG`, `Rot_Center`, `Motion_DOF` (0 fixed / 1 imposed / 2 free), `Inertia`, `Prod_Inertia`, `Added_Mass`, `Motion_Damp`;
  - imposed-motion files `<case>_body<k>_imposed.dat`.
- **Moving grids**: the grid moves rigidly; the metrics are rotated; ALE relative fluxes are used; the wall velocity is the grid velocity; the overset assembly is redone every step, with caches for the parts that do not change with motion.
- **Forces**: pressure and friction are counted once where patches overlap; moments are taken about the centre of rotation. Output files per body `k`: `_force_body<k>`, `_moment_body<k>`, `_force_earth<k>`, `_moment_earth<k>`, `_motion<k>.dat`.
- **Cases**: NACA66mod_Overset (imposed pitch and heave), KP505 (imposed rotation, J = 0.7), JBC (free sinkage and trim), Box_6DoF (two free boxes, two-phase).

### 4.8 Gap-free and permeation-free overset (update 12c)

This update came from your slice plots of Box_6DoF and KP505, which showed empty regions, and from your condition that **no grid point may ever lie inside a solid body**.

- **Safety rule**: an active cell next to a hole becomes a fringe only if its donor stencil contains a solved cell. Otherwise it stays active, and solid holes act as its wall.
- **Coverage pass**: every non-active fluid cell must have all 27 test points (nodes, edge midpoints, face centres, centre) inside solved cells of other grids. If not, the cell is activated. If it is a forced cell, the coarser cell underneath is activated instead.
- **Thin bodies**: a cell that a surface passes through is cut, even if its centre and nodes are outside (KP505 blade tips). A cell with a test point inside a body is never activated.
- **Fast winding number** (Barill et al.): per-body cluster trees with dipole and second-moment terms. It removed 1865 false "inside" cells on DelftTwist and made the GPPH hole cut 8× faster.
- **Permeation check** in grid-check mode: 0 active cells with a centre inside a body across all 18 overset cases.
- **Exact cell blank files** (`_cblank`) and a node blank that respects solid holes, for correct plots.
- **Self-overlap donors**: a block can take donors from the far part of itself (KP505 tip caps).
- `tools/blank_slice.py` draws sliced blank-status plots. It was used on all 18 overset cases.

### 4.9 Usability (update 13)

- **`Time_scheme = 2`**: 2nd-order SSP (Heun). Each step runs two stages and averages them.
- **Dynamic pressure** `p_dyn = p − p_hyd` is added to two-phase solution files.
- **Run-time name-list updates**: `probCFD.nml` is re-read every step. `Del_T`, `End_T`, iteration controls, intervals, schemes and stop limits apply to the running calculation, and the changes are logged.
- **`Reason_for_stop.md`** is written on every end (divergence, input error or normal end). It contains the reason, the offending cells, the last 40 steps, a diagnosis with the settings to change, and the run settings.
- **Log**: the 6-DoF block is laid out one quantity per line, and step groups are separated by three times the previous spacing.
- **Output**: files go to `Volume_solutions_grids/` as `<case>_NNNNNNNN.grd` (including `_00000000` for the initial grid) and `<case>_NNNNNNNN.p3d`.

### 4.10 Cavitation (update 15) and the energy equation (update 16)

- **Cavitation** (`$Cavitation`, `cavitation.F90`): Schnerr–Sauer, Zwart–Gerber–Belamri and Kunz as transport-equation models for a vapour fraction carried next to the water/air VoF; continuity source with a pressure-implicit part on the Poisson diagonal, `Vapor_scheme`, `Reboud`, absolute pressure level from the BCs; cases NACA66mod_cav, DelftTwist_cav.
- **Energy equation** (`$Energy`, `energy.F90`): temperature at the cells with its own solver controls (`*_Engy`), non-conservative bounded convection with four schemes, implicit/explicit/hybrid diffusion, Dirichlet / heat-flux / convective / adiabatic BCs through a 6th `.bcs` line, Kader thermal wall function on flag-11 walls, Boussinesq buoyancy (`Beta_l`, `Beta_a`), latent heat and p_v(T) with cavitation, `Output_data/<case>_heat.dat`, restart v3.
- **`Prs_form = 2`**: incremental pressure form (gradient of the latest pressure in the predictor). Needed for gravity-driven flows on wall-clustered grids: the non-incremental predictor diffuses the body-force velocity dt g at no-slip walls where the diffusion is implicit, and a fluid at rest starts to move.
- **Fixes**: vapour ghosts at Dirichlet-p faces, turbulence-solver deadlock with fewer ranks than sub-domains (collectives inside a per-domain loop), wall distance in Lch units.
- **Update 17**: `Cav_compressible = 1` (Wood-law mixture compressibility implicit in the pressure equation, vapour as an isothermal gas; verified on a piston-compressed tube and a travelling pressure step) and `Cav_thermal = 1` (Clausius-Clapeyron saturation properties, B-factor thermodynamic delay, implicit latent heat; LN2 cavity -9 %, water 0.7 %).

### 4.11 Incident waves (update 18)

- **Wave model** (`$Wave`, `waves.F90`): regular (Airy / Stokes 2nd order, deep or finite depth) and irregular seas (Bretschneider, JONSWAP, Pierson-Moskowitz; cos² spreading; component file), defined in the earth-fixed frame; the grid frame is converted with the shift X = x − ∫U_mean dt (encounter frequency automatically; `wave_shift_extra` is the hook of the future moving background).
- **Boundary condition**: `.bcs` flag 12 on u, v, w, p, alpha together — velocity = mean + orbital velocity (frozen at the surface value in the air), alpha = fill below the exact elevation, pressure ghost from the exact wave pressure; usable on every far-field face (the CFDShip-Iowa type-17 condition without the frame switch).
- **Relaxation zones** (`wave_fields.F90`): waves2Foam blend of velocity (water only) and alpha towards the incident wave over `Relax_In_L / Relax_Out_L / Relax_Side_L`; ramp, `Wave_Init`, probes (`Output_data/<case>_wave.dat`).
- **Validation**: 2-D tank, λ = 1 m, a = 1 cm: amplitude within 2.5 %, phase speed 2–3 % low on 32 cells per wavelength (first order in Δx, independent of Δt); encounter period with a mean current; JONSWAP sea.

### 4.12 Force models (update 19)

- **Three force models** (`force_model.F90`, validated: Hough–Ordway disk vs momentum theory within 1.5 %, point-force spring vs the analytic oscillation; the pure BEM needs its section data calibrated against open-water K_T/K_Q) as "bodies without a grid" in the 6-DoF hierarchy: `Force_model = n ... Force_model_End = -1` groups of `$SixDOF` hanging on a body, child or sub-child (or fixed in space). Hough–Ordway actuator disk (K_T/K_Q curves, discrete-exact thrust and torque, J from the disk inflow relative to the moving disk), pure blade-element momentum disk (station geometry, thin-airfoil lift with camber, parabolic drag, Prandtl tip loss), and a point force / couple given by user formulas of the mother's state (compiled expression evaluator). Disk forces enter the momentum equation weighted by the liquid fraction; the reaction (thrust + torque) acts on the mother and is transferred up the hierarchy; output files per model in the mother's body frame.

### 4.13 Moving background (update 20)

- **`Moving_BGR`**: the background blocks follow a Body in the permitted DOFs (`Moving_BGR_DoF`), optionally filtered, as their own rigid overset group; far-field Dirichlet values are earth-frame; the wave BCs need nothing (moved coordinates = earth positions), the relaxation zones stay attached to the grid faces (`Relax_L` per face); `Ovst_Reassemble_Tol` skips the re-assembly for sub-tolerance relative motion.
- Fixes: the restart file carries the background record; the force integration and the wall functions use the ramped Dirichlet wall velocity of a towing-tank start (`u_init (1 - ramp)`), as the boundary conditions do; a restart of a case with `t_ramp > 0` applies the boundary values of the restored time (before: those of t = 0 for the first step); the ALE grid velocity and the moving-wall velocity are the discrete displacement velocity of the grid point (exactly zero relative flux through a moving wall; moving-body results change at O(a dt) during accelerations); `Child_body` / `Sub_child_body` groups may follow the `_End` mark of their mother (the template order); the wave probes use a point-in-cell footprint test (a yawed cell was counted twice).
- Open item: a frame-dependent offset of a few per cent of ½ρU²DL in the pressure force of a towed body after its start-up transient (moving-grid path since update 12; the moving background reproduces the moving-body path to 0.1 %).
- Validation: towed cylinder with the background following it against the fixed cylinder in the stream (towing-tank start) and against the body moving through the static background; restart; wave tank surging / yawing through the earth-fixed wave (`verification/update20_moving_background`).

### 4.14 PID controllers (update 21)

- **`PID_*` blocks** (`pid_control.F90`, before `sixdof.F90`) inside a Body / Child_body / Sub_child_body group or a Force_model group, opened by `PID_use = 1`. Three types: throttle (output = rev/s), steering (output = deg), direct force (output = a component of a point force). A controlled **grid part** is rotated about `PID_Axis` through its `Rot_Center` (Rodrigues rotation in place of the Euler-angle pose; the equivalent Euler angles / rates are written into the 6-DoF state so the motion files and the mother / child coupling are unchanged; its `Motion_DOF` is ignored). A controlled **disk** gets its rps from a throttle controller and / or its axis turned about `PID_Axis` through `PID_Axis_Center` by a steering controller; a **point force** gets the direct-force output added to `PID_Component`.
- **Measurement**: any part's location (eta, deg for angles) or velocity (nu, deg/s for rates) by `PID_Target_body / _child_body / _sub_child_body`, `PID_Target_type`, `PID_Target_var`; up to 32 sequential targets (`PID_Target_value`, advanced when the error stays within `PID_Target_tol` for `PID_Target_hold_T`); angle errors wrapped.
- **Law**: `out = out0 + P e + I ∫e dt − D d(meas)/dt`, output limits with anti-windup, rate limit, start time; the state is replicated on all ranks and stored in the restart file; `Output_data/<case>_<part>_pid<i>.dat`.
- Fixes: the external wrenches of several force models on one mother are added (the last one overwrote the others since update 19); a restart re-evaluates the force models from the restored state (the first restarted step used the initial field's sources and wrenches) and restores the velocity-Verlet half-kick flag (0.3 % position offset of a restarted free DOF); the PID states are read before the file is closed.
- Verification (`verification/update21_pid`, cases written by `tools/pid_cases.py`): a steered foil (grid part, sequential pitch targets, rate limit), a spun cylinder (grid-part throttle on its own yaw rate), a direct-force heave controller, a disk throttle and a disk steering controller on a damped 1-DOF mass (the ActuatorDisk box declared body 1), a restart with controllers, and the regression of cases without controllers.

### 4.15 Probe points, shear stress, wall surface files (update 22)

- **Probe points**: `Probe_point<n> = x y z` in `$Probe_background` (following the background: fixed or moving with the background blocks) or inside a Body / Child_body / Sub_child_body group (following the part). Located every `Data_Intvl` step in the moved grid (last-cell neighbourhood, then the bucket search of every block; the active cell preferred); trilinear interpolation of the volume-file node values; one file per point (`<case>_[body1_[child1_[sub1_]]]probe<n>.dat`); `NoData` in holes / outside the grid.
- **`Shear_Stress`** in the volume files after `TTL_Prs`: μ_eff √(2 S_ij S_ij) at the cells (laminar + eddy viscosity), |τ_w| of the wall faces at the solid-wall nodes (wall function included); no near-wall smoothing.
- **Wall surface files** (`SW_Prs_Frc_surf`, `Frame_Prs_Frc_surf`): `<case>_surf_NNNNNNNN.grd / .p3d`, one block per wall patch, pressure stress −p n and friction stress vectors [Pa]; the face loop of the force integration became the shared `wall_face_stress`, and the node values are dual-area averages whose trapezoidal integral reproduces the force files exactly (`check_surface.py`).

### 4.17 Fluid-structure interaction (update 24)

- **Input**: `SW_FSI = 1` inside a part group; `Youngs_modulus`, `Poisson_ratio`, `Density` or `Mass`; `Num_NoFSI_Group` + `NoFSI_block_num<n> / NoFSI_i/j/k_range<n>` (clamped wall nodes); `Num_FSI_Submat_regions` + `Submat_*<n>` (a second material); optional `FSI_Start_T`, `Damping_ratio`, `FSI_Min_Thickness`, `FSI_Max_Iter`, `FSI_Tol`, `FSI_Relax`, `FSI_Mode`, `FSI_Monitor_node<n>`.
- **Structure from the surface grid**: every no-slip wall patch of the part is a sheet (around / span directions, upper / lower / cap nodes by the section's principal axis); a ray from each upper node to the opposite side gives the mid-surface node and the local thickness; sheets of one part are merged at coincident nodes; patches that are not sheet-like or overlap a larger one are passive (their nodes follow the mid-surface). Every wall node is attached to the mid-surface with bilinear weights and a rigid offset: loads go in, motion comes out through the same weights.
- **Shell** (`fsi_shell.F90`): co-rotational flat 4-node shell (membrane + MITC4 plate + drilling), per-element material / thickness, lumped mass, Rayleigh damping, generalized-alpha with Newton iterations, RCM-banded Cholesky. Stand-alone tests (`make fsi_test`): cantilever statics, free vibration, large deflection, Turek-Hron CSM3.
- **Coupling**: predictor, mesh morphing of the part's blocks (rigid boundary layer, cosine decay of the rotation to the outer boundary, rigid pose on top), rebuilt metrics, discrete grid velocity in the ALE fluxes, re-projected face fluxes, GCL term; the inner iterations of the fractional step are the FSI sub-iterations (wall forces -> structure -> IQN-ILS relaxation -> morph), restarted from the level-n field; `FSI_Start_T` holds the structure while the flow develops.
- **Output**: monitor-node histories, the mid-surface mesh files, `Disp_x/y/z` and `vonMises` on the wall surface files; structural states in the restart file.
- **Verification**: CSM3 (f 1.096 Hz, amplitude 65.1 mm against 1.0995 Hz, 65.2 mm), FSI-PfS-2a on a quasi-2-D overset grid (section 36.2 of PROJECT_RECORD.md), KP505 with flexible blades (run-through).

### 4.16 On-the-fly turbulence statistics (update 23)

- **`$Statistics`**: raw moments (sums of x, x², xy and of the 9 gradients, their squares and 3 cross products, the modelled ε and k, ν, Δ) accumulated per point with counters that skip blanked steps; the fluctuation statistics are formed at output. Group 1: Coverage, mean / rms of every field, Reynolds stresses, k_res. Group 2: the gradient statistics, ε_res from the 12 components, ε_model (model destruction term), totals, η, λ_g, Re_λ, Δ/η, vorticity.
- **Where**: the original grid (group 1 at the nodes, group 2 at the cells with the face-difference gradient, t = 0 grid, part frames for moving blocks) and boxes / cylinders / planes attached to the background or a part, sampled through the probe machinery (`pt_locate` / `pt_owner_key` / `pt_sample` of `probes.F90`, now a fixed 24-value layout shared by the probe files and the statistics). Cylinders can write polar components.
- Restart: the accumulators follow the PID states in the restart file; the lattice sums restored on rank 0 combine with the new partial sums.
- Verification: grid statistics against the per-step volume files (file precision), lattice points against probe histories, ε_res from the exported components, restart, SST path, coverage of a hole.

---

## 5. Test and validation cases

| Group | Cases |
|---|---|
| Tutorial / verification | Backward-facing step, Cylinder 1/2-phase, FlatPlate Re 10⁵–10⁸ and WF, NACA66mod (1-phase, O-grid, 3D), PoleCyl, Drop2D_SurfTension, Cylinder sigma 0.7 (cavitation), 1-D conduction slabs and step advection (energy), NaturalConvection_Cavity (de Vahl Davis Ra 1e5) |
| Overset, small | Cylinder_Overset_1/2phase, NACA66mod_Overset (+ 1phase), JUSKFoil, DelftTwist |
| Overset, ships and propulsors | KCS, JBC, ONRT, DTMB, LCS, DCWJ, GPPH, NUMRT, NSWC15E, SPSCFoil, KP505 |
| 6-DoF | Box_6DoF (free), NACA66mod_Overset (imposed), KP505 (imposed rotation), JBC (free heave and pitch), GPPH (hull body + stern-drive child + propeller sub-child, update 14) |
| PID controllers (update 21) | PID_Steer (foil steered to +-4 deg), PID_Throttle (cylinder spun to a yaw rate), PID_Direct / PID_Disk / PID_DiskSteer (point force, disk throttle, disk steering holding the heave of a damped mass) |
| FSI (update 24) | FSI_PfS2a (cylinder with a flexible rubber plate and a steel rear mass, Kalmbach & Breuer 2013 / Fabbri et al. 2023), KP505_FSI (five flexible blades as child bodies of the hub, run-through) |
| Statistics (update 23) | Statistics_ActuatorDisk (grid statistics = per-step files, box / cylinder / plane, restart, SST), Statistics_TowedCylinder (body-following lattices = probe histories, hole coverage, polar ring) |
| Probes / surface files (update 22) | Probes_TowedCylinder (body- and background-following points, NoData inside the cylinder, surface integral = force file), Probes_SteeredFoil (points turning with the foil, overset wall weights) |

Reference checks:

- A body at rest reproduces the static solver to round-off.
- A slowly rotating cylinder grid reproduces the static drag to 10⁻⁶.
- Box_6DoF gives a hydrostatic force of 1034 N per box (ρgV/2 = 1030 N).
- The Cylinder_Overset_2phase monitor is unchanged across updates 12b → 13 (with `Time_scheme = 1`).
- Update 14: a block of a body declared a fixed child reproduces the single-body run; a child + sub-child with imposed motions reproduce the single body with the composed Euler angles (PROJECT_RECORD.md 26.5).

---

## 6. Conventions

- **Input** (in the case folder): `probCFD.nml`, `Input_polyMesh_probCFD.csv` (rows: i/j/k split, number of blocks, Object_ID), `<case>.grd` (PLOT3D ASCII), `<case>.bcs`; `<case>.bcs_b2b` is generated.
- **Output**:
  - `Volume_solutions_grids/<case>_NNNNNNNN.grd` and `.p3d`, with variables u, v, w, p (+ alpha, p_dyn) (+ Vapor) (+ Temperature) (+ Blank = exact cell status 1 / 0 / -1);
  - `<case>_monitor.dat`, `<case>.ovst`, the 6-DoF data files, and `Reason_for_stop.md`.
- **Build**: `make clean && make` (gfortran + OpenMPI). The Makefile does not track dependencies.

---

## 7. Open items

- **KP505 tip caps**: the two ends of each C-type cap overlap by only 2 cells, so a thin band at the seam is interpolated but never solved. The fix is in the grid: point-match the ends, or extend one end by at least 5 layers.
- **Hull cases**: many "solid-bounded" cells remain where a wall block ends with only forced layers. Extending those wall blocks by the fringe width would remove them.
- **Coarse Box_6DoF grid**: it fails after a few steps because the box edges are almost entirely forced fringes, which is a limitation of the grid.
- **Impulsive starts** of imposed motion: use `t_ramp` or a smooth motion start.
- **Fluid–body coupling** of the rigid 6-DoF is explicit; `Added_Mass` stabilises light bodies. The flexible parts (update 24) are coupled strongly (sub-iterations inside the step).
- **Towing-tank start with tiny wall cells**: during `t_ramp` the ramped wall velocity has a normal component on walls not aligned with the flow, so fluid passes through them; the implicit upwind of that inflow (whose ghost value is the mirror of the interior value) makes the momentum sweeps diverge when the wall-normal cell is much smaller than U dt (the 0.1-mm root cells of KP505: the delivered KP505 cases of updates 12-23 could not take a first step in the sandbox for this reason). `cases/KP505_FSI` uses `t_ramp = 0`. A clean fix would keep the wall impermeable during the ramp (ramp the tangential part only) or treat the mirrored ghost implicitly.
- **Large cases** (KCS, ONRT, JBC, …) were checked for assembly and start-up in the sandbox only. Full runs need a larger machine.
- **Explicit diffusion in several directions**: the hybrid switch tests each direction against `vn_switch` separately; two explicit directions of VN < 0.5 each can exceed the 2-D explicit limit (0.5 in total) and produce an odd-even instability (seen in the air-filled cavity on its wall-clustered grid). Use `Dffs_mode = 2` or a smaller `vn_switch` where the Von Neumann numbers are moderate in more than one direction.
- **Pure BEM loading**: with the face pitch and a thin-airfoil polar the BEM disk over-predicts the angle of attack of a B-series blade (K_T about twice the series value at the same J, even with a finite-span lift slope). Set `Prop_CLa` / `Prop_fc` / the pitch against known open-water K_T, K_Q (reported every step), or use the Hough–Ordway model when only the curves are known.
- **Body integrator**: the velocity-Verlet half-kick of update 19 changes moving-body results of earlier updates at O(Δt); the stored restart state has no half-kick flag, so the first step after a restart skips it (negligible).
- **Wave dispersion**: the numerical phase speed of a wave on 32 cells per wavelength is 2–3 % low and the error falls only at first order with Δx (Δt has almost no influence) — the spatial accuracy of the free-surface / pressure coupling at the interface, not the time scheme. Use 48–64 cells per wavelength when the phase matters (wave-force phases, added resistance).
- **`Time_scheme = 2` and `Iner_Itn = 2` with the geometric VoF** damp propagating waves by 12 % in two periods (the SSP stage average of two PLIC fields is not a PLIC field; the Picard pass re-advects alpha with the corrected fluxes). Use `Time_scheme = 1`, `Iner_Itn = 1` for wave runs.
- **Relaxation zones** blend alpha, so the liquid volume is not conserved inside the zones (monitor `Vliq`); the velocity blend leaves a divergence of order 0.1 (non-dimensional) in the zone cells that the next projection removes.
- **Waves on a mean current**: the amplitude grew 5–13 % along two wavelengths at U = 0.5 m/s (U/c = 0.4) with the 2nd-order upwind scheme — to be checked with grid refinement before relying on wave-current cases.
- **Statistics on moving grids**: a fixed grid index of a moving block averages in that block's part frame (`Statistics_Type1_Frame = 1`); earth-fixed statistics in a region swept by a moving part come from a background-attached lattice, whose points are skipped (Coverage < 1) while they lie in the moving hole. The gradient of ε_res uses the face differences of the staggered velocities for the three own-direction derivatives (compact, accumulated at the cells) and the central average of the edge differences for the transverse ones; ε_res is a resolved-scale quantity and only physical where Δ/η is small. The SA / WMLES modelled ε is the equilibrium estimate ν_t|S|²; the SST-IDDES one uses the RANS destruction β*kω also in the LES region.
- **PID controllers on grid parts** are kinematic: the controller sets the rotation (rate or angle) of the part directly, without the part's inertia, and a controlled part has no free DOF of its own (a steered rudder is a child whose hull carries the free DOFs). Because the loop is algebraic for a steering controller (the angle IS the output), the proportional gain must stay below ~0.7 (deg per deg) or the discrete loop oscillates — put the authority in the integral gain. A controller acting through a force model (disk rps / angle, point force) drives the plant through the 6-DoF dynamics as usual.
- **Controller gains are case-specific**: the delivered PID cases are tuned for their plants (a damped 100 kg mass, the kinematic foil / cylinder); start a new application from the output of the first steps (`_pid<i>.dat`) and the sign rule (positive gains when a positive output drives the measured variable up).
- **`Prs_form = 1`** (default, for backward compatibility) carries the spurious-velocity defect at no-slip walls under gravity described in 4.10; consider `Prs_form = 2` for new gravity-driven cases.
