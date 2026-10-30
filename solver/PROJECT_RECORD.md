# probCFD PROJECT RECORD — complete work history and handover document
Date of record: 2026-09-04 (probCFD_부산 update 2, section 13).  Previous
records: 2026-09-03 (section 12), 2026-09-02 (sections 6-11),
2026-07-23 (sections 1-5).
Purpose: full record of the project so it can continue seamlessly in a
new room.  This folder is the CURRENT, VALIDATED source code plus the
two-phase cylinder example case; `regression_BFS/` holds the previous
single-phase case with its reference values.  Sections 1-5 are the
history of the previous rooms (condensed); sections 6-11 are this room.

===========================================================================
1. WHAT probCFD IS (final state)
===========================================================================
A 3-D incompressible Navier–Stokes solver, Fortran 90 + MPI (OpenMPI),
finite differences on FULLY STAGGERED structured multi-block grids with
a GENERAL 3-D CURVILINEAR COORDINATE TRANSFORMATION; since this room
also TWO-PHASE (air–water, algebraic VoF, one-fluid formulation) with
gravity, selectable explicit/implicit/hybrid time integration of the
convection and diffusion terms, and a Krylov pressure solver.
Projection (fractional-step) method; Jacobi iterations for the
momentum predictor; 1st/2nd-order upwind convection (`conv_scheme`).
Non-dimensionalized internally by Lch, Uch and rho1 (fluid 1);
inputs/outputs are dimensional.

Capabilities, in order of development:
  M1  b2b interface auto-detection (`bcs_b2b.F90` -> `<case>.bcs_b2b`,
      type-92 records, affine index maps, all 24 orientations, partial
      faces; since this room ALSO several disjoint windows per face pair).
  M2  Multi-block solver with the unified ghost-exchange engine
      (intra-block halos + b2b patches, record-ordered unpacking ->
      results independent of the rank count for the Jacobi paths).
  M3  Per-block MPI splitting via csv.
  M4  Periodic BCs (velocity flag 3).
  M5  3-D curvilinear transformation (half-index metrics, Park form).
  M6  Backward-facing-step example (Pitz & Daily grid) + `make_case.py`.
  M7  `curvi_form` switch (1 unweighted U^q, 2 J-weighted flux).
  M8  (this room) `Cnvt_mode` / `Dffs_mode` explicit-implicit-hybrid
      time integration per term and per direction.
  M9  (this room) two-phase VoF (CICSAM), gravity, variable properties,
      balanced-force pressure gradient, alpha in the .q files, alpha BCs.
  M10 (this room) projection made exactly consistent at non-orthogonal
      physical boundaries, pressure outlets and multi-block interfaces
      (flux-form Poisson operator with exchanged face fluxes, synchronized
      interface Jacobians); FGMRES + PCG pressure solver (`poi_solver=2`).
  M11 (this room) start-up ramp of Dirichlet velocity BCs (`t_ramp`);
      effective per-family BC flags on reoriented blocks.
  M12 (this room) surface-piercing cylinder example (4-block O–H grid)
      + `make_case_cylinder.py`.
  M13 (probCFD_부산 update, section 12) inner iterations, residual-
      based criteria, dual pressure norms, single-valued node output,
      line-SGS preconditioner, grouped name list, reference log format,
      implicit CICSAM option, NACA66mod case.
  M14 (update 2, section 13) QUICK and 2nd-order ENO convection
      (Cnvt_scheme 3/4), pole (degenerate face) boundary condition with
      type-93 records in .bcs_b2b, PoleCyl pipe case; N_max and
      Initial_Value_fct removed.

===========================================================================
2. FILE FORMATS AND CONVENTIONS
===========================================================================
* Grid `<case>.grd`: multi-block ASCII PLOT3D; nblk; ni nj nk per
  block; then x,y,z arrays per block, i fastest. Meters.
* BCs `<case>.bcs`: nblk; per block: nfaces; per exterior face:
  direction d (+-1,+-2,+-3), i-range, j-range, k-range (NODE indices),
  then per variable u,v,w,p "flag value" lines and, ONLY when
  Num_of_Phase = 2, a 5th line for alpha.  Comment lines starting with
  '!' and trailing text after the numbers are allowed (the shipped
  Cylinder.bcs is fully annotated).
    velocity/p flags: 1 = Dirichlet, 0 = zero gradient, 3 (u line) =
    periodic; alpha flags: 0 = zero gradient, 1 = Dirichlet value in
    [0,1], 2 = still-water profile (liquid below z_free_surface; use at
    inlets).
  Dirichlet velocity values are GLOBAL CARTESIAN components [m/s].
  The condition applied to a TRANSFORMED face family U^q is the flag
  of the Cartesian axis the grid direction q points along on that
  patch (`bc_eflag`, computed from the patch geometry in read_bcs) -
  e.g. on the cylinder O-grid, whose i-axis is vertical, the bottom
  condition of U1 is the w condition.  Interior (b2b) faces get NO
  record.
* `<case>.bcs_b2b`: auto-generated at first start; do not ship.
* csv: rows i_split/j_split/k_split/n_of_BLK per block column.
* Namelist `probCFD.nml`: the shipped one is annotated variable by
  variable.  Pitfall: a value placed after a '!' is silently ignored.
  New this room: Cnvt_mode, Dffs_mode, cfl_switch, vn_switch,
  Num_of_Phase, rho1, mu1, rho2, mu2 (rho, mu accepted as aliases of
  rho1, mu1), grav(3), z_free_surface, vof_cfl_max, t_ramp, poi_solver.
* Outputs: binary PLOT3D `<case>.x` (once) and `<case>_NNNNNNNN.q`
  (float32; u,v,w [m/s], p [Pa], + alpha when Num_of_Phase = 2):
  node-interpolated, GLOBAL-frame components for all blocks.
* Monitor line per step: max|div| (nondimensional, after correction),
  itm (momentum Jacobi sweeps), itp (Poisson iterations; solver 2:
  total inner CG iterations), res(p) (solver 1: max change per sweep;
  solver 2: ||b - A p||_inf/||b||_inf of the full operator); two-phase
  adds nsub (VoF sub-cycles), Vliq [m3], amin, amax, clipped [m3]
  (cumulative volume removed by the [0,1] clipping).

===========================================================================
3. THE CURVILINEAR FORMULATION (unchanged core; see the previous record
   and `curvilinear_study_notes.md` for the derivation)
===========================================================================
* Primary face unknowns: transformed components U^q = u.grad(xi_q)
  (curvi_form 1) or volume fluxes J U^q (form 2) at the q-faces,
  exchanged across halos/b2b/periodic patches; p (and alpha) at cells.
* Cartesian face velocities u,v,w (x/y/z components at the i/j/k face
  families) are local derived fields; the momentum predictor solves
  them (no Christoffel terms): cart_to_contra -> refresh (exchange +
  BCs in the stored-variable space) -> contra_to_cart(ghosts only).
* Metrics on the half-index lattice from ghost-extended node
  coordinates (Xh, Jh, mh = d(xi_q)/d(x_c)); accessors mq, gql, Gcap,
  tq, fvw; NEW: alph_h, rho_h, mu_h (fluid properties at any half point,
  arithmetic averages of the touching cells' alpha).
* Operators: conservative flux-form convection with upwinded values;
  diffusion (1/(J rho Re)) sum_q D_q[mu_h sum_l G_ql d(phi)/d(xi_l)],
  compact parts implicit, cross parts lagged; Poisson and corrector in
  flux form (section 7).

===========================================================================
4. KEY DECISIONS OF THE PREVIOUS ROOMS (kept)
===========================================================================
D1 face unknowns = transformed components (both forms behind fvw);
D2 momentum for global-Cartesian components; D3 conservative operators
(rectilinear limit reproduces the old Cartesian code); D4 BC values
global Cartesian; D5 pointwise half-index metrics (GCL to O(dx^2));
D6 reoriented blocks are legitimately different discretizations.

===========================================================================
5. VALIDATION OF THE PREVIOUS ROOMS (still true for the paths that were
   not changed; the digits of R5 are superseded, see R7)
===========================================================================
R1 rectilinear 27-block cavity reproduces the pre-curvilinear code;
R2 rank-count invariance (Jacobi paths); R3 orientation tests;
R4 curved-grid GCL residual ~2e-4; R5 BFS first-five-step digits
3.403E+01 ... 1.026E+01 (form 1) — reproduced bit for bit by the
refactored operators of this room BEFORE the projection fix (section
7.3), i.e. the refactoring itself changed nothing; R6 stability bound.

===========================================================================
6. THIS ROOM: REQUIREMENTS AND WHAT WAS BUILT
===========================================================================
Requirement specification (user): (a) Cnvt_mode / Dffs_mode switches
(1 explicit, 2 implicit, 3 hybrid per i/j/k direction by the local
CFL / Von Neumann number); (b) two-phase capability with an algebraic
VoF method: new variable alpha in the solution file, BC and IC support,
Num_of_Phase switch, heavy liquid filling z < 0 initially, both fluids'
density/viscosity and the gravity vector from the namelist; (c) an
example case on the supplied 4-block cylinder grid with inlet 1 m/s at
x-min, outlet at x-max, zero-gradient side walls, p = 0 at the top,
interface at z = 0; (d) a case zip with sources, Makefile, grid,
annotated namelist, per-face documented bcs.
Decisions taken with the user: BiCGSTAB-class Krylov pressure solver
allowed (became FGMRES+PCG, see 7.5); CICSAM as the VoF scheme; free-slip
bottom; surface tension omitted (We ~ 1e4).

6.1 Time-integration switches (convection.F90 `mom_ops`, u_momentum.F90)
  * face_setup(c, pt, ...) returns for every face unknown the six
    advecting fluxes J U^q (time level n), and the masks cm(q), dm(q)
    (1 = implicit, 0 = explicit).  Hybrid: cm(q)=1 iff
    CFL_q = dt |U^q| >= cfl_switch (computational-space CFL, evaluated
    from the advecting fluxes divided by the local J), dm(q)=1 iff
    VN_q = dt (mu_h/(rho_f Re)) g^{qq} >= vn_switch.  Defaults: 1.0 and
    0.5 (classic 1-D criteria; in 3-D the explicit directions must
    together satisfy sum CFL <= 1, sum VN <= 1/2 - tighten the
    thresholds if needed).
  * compute_rhs_momentum: rhs = u^n/dt - CONV_exp(u^n) + DIFF_exp(u^n)
    + g*, every term multiplied by its mask; the default path
    (Cnvt 1, Dffs 2, single phase) is bit-identical to the old code.
  * solve_momentum: implicit compact diffusion in the Jacobi diagonal,
    cross terms lagged; implicit convection = linearized FIRST-ORDER
    upwind in the matrix (diagonal max(a+,0)-min(a-,0) => diagonal
    dominance), and for conv_scheme = 2 the difference between the 2nd-
    and 1st-order upwind fluxes as a deferred correction on the current
    iterate (converged result = fully implicit 2nd-order upwind).  With
    everything explicit the loop converges in 2 sweeps.
  * All three families share one generic code path (face_range,
    face_active, cross_sum, diag_diff, conv_sum, xder); the loop order
    of the old code is preserved so the default path stays bit-exact.

6.2 Two-phase VoF (vof.F90, mod_global.F90, boundary.F90, output.F90)
  * alpha at cell centers (field id 8), exchanged like p; ghost/BC
    through apply_bc (flags above; flag 2 = fill_fraction of the ghost
    cell from its node z-range vs zfs = z_free_surface/Lch).
  * Initial condition: alpha = fill_fraction(cell) (partial fill for
    cells cut by the plane), then rho, mu everywhere from alpha.
  * Transport once per step BEFORE the momentum step with the
    divergence-free fluxes F^q = J U^q of time level n: explicit,
    conservative cell update; face values by CICSAM (HYPER-C / ULTIMATE-
    QUICKEST blended by gamma = cos^2 of the angle between grad(alpha)
    at the donor cell and the face normal grad(xi_q)); normalized-
    variable band check; donor-cell Courant number; result clamped to
    the donor/acceptor band; automatic sub-cycling nsub =
    ceil(max outflow Courant / vof_cfl_max); final [0,1] clipping with
    the clipped volume reported.  Work arrays rhsu/v/w hold the face
    values (free at that point of the step).
  * Properties: rho* = alpha + (rho2/rho1)(1-alpha), mu* likewise, at
    any half point via alph_h (average of the touching cells).  The
    momentum equation is the one-fluid equation divided by the face
    density: du/dt + CONV = -(1/rho) grad p + (1/(rho Re)) div(mu grad u)
    + g*, g* = g Lch/Uch^2.  (The transposed-gradient part of the
    viscous term, nonzero only where mu varies, is omitted - a common
    simplification; document if it is ever added.)
  * Output: nvar = 5 in the PLOT3D function files when Num_of_Phase=2.

6.3 Example case: see README.md in this folder.

===========================================================================
7. THIS ROOM: PROBLEMS FOUND IN THE INHERITED CODE AND THEIR FIXES
   (each was found by a diagnostic, all fixes are in the shipped code)
===========================================================================
7.1 Dead code / stale comments: the 1-D arclength metrics (xf/xc/dxc/
    dxu, build_block_metrics, fill_ghost_metric) were no longer read by
    any kernel - removed; stale "local frame" comments corrected.
7.2 Step-1 ghost inconsistency: the initial Cartesian ghost values were
    never reconstructed from the BC-applied U -> contra_to_cart added
    to initialize_fields (affects only the very first RHS).
7.3 INCONSISTENT PROJECTION AT BOUNDARIES (the reason max|div| of the
    BFS case never fell below O(1)): on non-orthogonal walls/inlets the
    Poisson operator carried a cross-term pressure flux g^{ql} p_l
    through faces whose velocity is fixed (the corrector never touched
    them), and at the outlet (Neumann velocity + Dirichlet p) it
    carried the Dirichlet flux through a face the corrector left at its
    star value.  Fix: block-boundary face TYPES (tu_lo.. arrays,
    face_type()): 0 interface, 1 fixed-flux physical boundary (no
    pressure flux, corrector keeps the star value), 2 pressure outlet
    (Dirichlet p + zero-gradient normal velocity: corrected by the
    projection like an interior face; the zero-gradient copy is applied
    to the STAR field only, the corrected field keeps the projected
    face value and copies it outward).  BFS max|div| went from 34 to
    ~1e-5 (with a converged Poisson solve).
7.4 MULTI-WINDOW INTERFACES (found on the cylinder grid):
    a) bcs_b2b: the matched node set of a face pair is now split into
       4-connected components, each a rectangle -> a face that meets
       the same neighbour in several disjoint windows (blocks 3<->4 on
       the y=0 plane, separated by the O-grid) is exported as several
       records (12 records on the cylinder grid instead of 8).
    b) Exchange windows widened by ng tangentially were writing ghost
       cells beyond a window junction with data that is the geometric
       continuation of the WRONG neighbour (block 4's face bends onto
       the O-grid there while block 3's ghost cells belong to block 1).
       Fix: every exchange pass runs twice - widened windows first (to
       reach edge/corner ghosts no window covers), then the un-widened
       CORE windows overwrite (authoritative).  The ghost coordinates
       (fill_ghost_coords) use the same rule.  The widened sub-pass never
       writes the shared face plane itself (owned dofs of another
       window were being clobbered).
    c) At a window junction the ghost-node column is claimed by two
       windows, so the receiver's edge ghost cell is distorted and the
       two blocks computed DIFFERENT Jacobians for the same shared face
       (2.35e-3 vs 3.70e-3): the volume flux J U^q through the face
       differed between the sides (continuity and VoF mass conservation
       broken there).  Fix: sync_face_jacobian - non-owned interface
       faces take the owner's J_f (exchanged once after build_metrics).
    d) Even so, the owner and the receiver discretize the cross terms of
       a shared face with different tangential neighbours at a window
       edge.  Fix (general, also covers 7.3): the Poisson operator is
       assembled in FLUX FORM - face_grad computes the transformed
       pressure-gradient flux G^q at every face this domain computes
       (ids 11..13 in rhsu/v/w), exchange_faces delivers the OWNER's
       G^q to non-owned interface faces and halo faces, cell_div_G
       forms (A p)_P = (1/J) sum(G+ - G-); the corrector uses the same
       exchanged G^q (U^q -= dt fvw G^q/J).  Discrete continuity of the
       corrected velocity now equals the Poisson residual on every cell
       (verified: expected vs actual divergence identical to all digits).
7.5 PRESSURE SOLVER: Jacobi-preconditioned BiCGSTAB (first choice) is
    fine single-phase but wanders and breaks down on the 1000:1
    coefficient jump; a deferred-correction PCG converged only at ~0.3
    per outer sweep on the skewed BFS grid.  Final: FGMRES(20) on the
    full 19-point operator (true residual at every restart), right-
    preconditioned by an inexact Jacobi-PCG solve (1e-2 reduction, <=200
    iterations, one-hop ghost refresh) of the SYMMETRIC 7-point part,
    whose face coefficients are exchanged (owner values) so the matrix
    is exactly symmetric across interfaces.  Krylov vectors need an
    exactly linear operator: refresh_all(iv, reps=3) repeats the
    exchange+BC cycle so that relayed edge/corner ghosts are not one
    cycle stale (harmless for fixed-point iterations, fatal for Krylov).
    All-Neumann problems: J-weighted mean of b removed, level pinned.
7.6 SPURIOUS CURRENTS (two-phase): the Cartesian face velocities were
    corrected with a corner-averaged dp/dxi divided by the local face
    density -> in the air next to the interface the water's pressure
    gradient produced ~500 g dt of spurious vertical velocity per step
    (60,000 m/s within 3 steps).  Fix: balanced-force form - every
    directional pressure difference is divided by the density of ITS
    OWN face first (Ql), then corner-averaged (Qavg); used both in the
    cross terms of G^q and in the Cartesian correction.  On this
    extruded grid with a flat interface the hydrostatic state is then
    exactly at rest.
7.7 BC FLAG MAPPING ON REORIENTED BLOCKS: apply_bc used the flag of the
    Cartesian component with the SAME INDEX as the transformed family
    (U1 <-> u flag).  On the O-grid (i vertical) the free-slip bottom
    (w = 0 Dirichlet, u,v zero-gradient) therefore became zero-gradient
    for the bottom-normal flux: water fell through the floor and no
    hydrostatic pressure formed.  Fix: bc_eflag(q) = flag of the
    Cartesian axis along which grid direction q points on the patch
    (face normal / grid-line tangents at the patch centre from the node
    coordinates); used by apply_bc and by the face-type classification.
7.8 START-UP: an impulsive start (u = 1 m/s everywhere at t = 0 next to
    the cylinder wall) needs an impulsive pressure ~rho_w U L/dt which
    is transmitted into the 0.7 mm air cells at the contact line ->
    227 m/s in the first step.  Fix: `t_ramp` (linear ramp of all
    Dirichlet velocity values, bc_ramp = 0 at t = 0) and u_init = 0;
    the example starts from rest.  Also the VoF step must see a
    divergence-free field: the BC ramp factor is set BEFORE
    initialize_fields so the initial inlet flux matches the interior.

===========================================================================
8. VALIDATION CATALOG OF THIS ROOM (reference numbers)
===========================================================================
R7 BFS regression (regression_BFS/, np = 2, 5 steps), see its README:
   Jacobi: 4.556E-01 2.351E-01 2.712E-01 2.116E-01 1.770E-01 (bit-
   identical np=1/np=2); FGMRES+PCG: 1.513E-03 2.651E-03 3.631E-06
   2.080E-06 2.154E-07 with res(p) 3.4E-05 ... 7.2E-09.  Form 2 with
   FGMRES: 9.174E-03 ... 3.632E-07.
R8 Switch matrix on BFS (3 steps, FGMRES): Cnvt 2 (itm 7-8), Cnvt 3
   (= explicit at this dt), Dffs 1 and 3 (itm 2), conv_scheme 2 with
   Cnvt 2 (itm 8): all run, max|div| 1e-5 ... 5e-5 (solver tolerance).
R9 Single-phase cylinder (water only, gravity, Cnvt/Dffs 3, FGMRES):
   max|div| 1.7E-04 -> 2E-05 over 6 steps, itm 24 then 2, itp ~1500-2000.
R10 Two-phase cylinder from rest (shipped namelist): step 1-6 max|div|
   1.7E-10 ... 1.6E-08, res(p) 3.5E-10 ... 6E-09, nsub = 1, liquid
   volume 190.436 -> 190.441 m3 (= the ramped inflow), clipped 6E-11 m3,
   p_max ~2.0E4 Pa (hydrostatic), velocities O(0.01-0.05 m/s) while the
   inlet ramps.  Bounds-checked build (-fcheck=all -ffpe-trap) runs the
   same steps with identical digits.  Longer run: see section 10.
R11 Projection consistency check: at the worst cell of the single-phase
   cylinder, expected divergence div(U*) - dt A p and the actual
   divergence agree to all printed digits; max|div| then equals the
   Poisson tolerance times ||b||_inf dt (impulsive-start wall cells).

===========================================================================
9. KNOWN PROPERTIES / LIMITATIONS / GOTCHAS
===========================================================================
* Cost: the two-phase Poisson solve dominates (~2000-2700 inner CG
  iterations per step with the 1e-6 relative tolerance; the Jacobi-PCG
  preconditioner is weak on the wall-clustered anisotropic cells).  A
  line-relaxation or multigrid preconditioner is the natural upgrade.
* poi_tol is relative to ||b||_inf; after an impulsive start ||b|| is
  dominated by the wall cells and the divergence elsewhere is
  correspondingly looser - start from rest with t_ramp.
* FGMRES results are deterministic for a given rank count only (global
  reduction order); the Jacobi paths remain bit-identical for any np.
* Hybrid thresholds are per direction (as specified); the explicit
  directions must jointly satisfy the multi-D stability limits.
* Implicit convection (Cnvt 2/3) uses Jacobi sweeps: at large CFL the
  contraction per sweep is ~CFL/(1+CFL) - raise mom_maxiter.
* CICSAM is explicit: sub-cycling keeps the donor Courant < vof_cfl_max;
  the clipped volume monitor shows any boundedness loss.
* Surface tension not included; viscous transposed-gradient term
  omitted; laminar (no turbulence model) at Re_D = 1e6.
* Mixed (slip-type) BC flags on CURVED patches use one patch-averaged
  axis mapping (bc_eflag); exact on planar / axis-aligned patches.
* The Makefile does not track dependencies: `make clean && make`;
  compile order = the SRCS list (vof.F90 before boundary.F90).
* Run inside scripts with `mpirun --allow-run-as-root --oversubscribe`
  and `< /dev/null`; do not `pkill -f probCFD` from a shell whose command
  line contains the word (it kills the shell) - use `pkill -x probCFD`.

===========================================================================
10. STATE OF THE LONG TWO-PHASE RUN AT PACKAGING TIME
===========================================================================
See the end of README.md / the delivery message for the number of steps
completed on the development machine (2 cores, ~25 s per step) and the
observed behaviour (liquid volume tracking the ramped inflow, bounded
alpha, divergence at the solver tolerance).  The shipped namelist runs
2000 steps (10 s); nothing in the code limits longer runs.

===========================================================================
11. HOW TO RESUME IN A NEW ROOM
===========================================================================
1. Unzip; `make clean && make` (gfortran + OpenMPI).
2. Sanity: `cd regression_BFS && cp ../probCFD . && mpirun -np 2 ./probCFD`
   and compare with R7; then the cylinder case for a few steps (R10).
3. Natural next steps: stronger Poisson preconditioner (line-Jacobi /
   multigrid); implicit VoF or geometric (PLIC) interface reconstruction;
   surface tension (CSF); turbulence model; second-order time stepping;
   moving grids.

===========================================================================
12. probCFD_부산 UPDATE (2026-09-03): THE UPDATE TASK LIST
===========================================================================
Package layout: solver sources at the top level; `Cylinder_case/`,
`NACA66mod_case/`, `regression_BFS/` each with grid, bcs, csv and the
grouped name list.  Everything below is implemented and tested (12.6).

12.1 Inner iterations (`Iner_Itn`, Main.F90)
  Per time step the sequence [VoF transport from alpha^n] -> momentum
  predictor -> outflow correction -> Poisson -> correction is repeated
  Iner_Itn times.  New per-domain arrays un/vn/wn (Cartesian face
  velocities of level n) and alpn (alpha^n).  compute_rhs_momentum uses
  un/dt for the time derivative and the LATEST iterate (u, U1..U3 after
  the previous correction) for the explicit convection/diffusion terms,
  the advecting fluxes ap/am of face_setup and the 2nd-order deferred
  correction; the predictor starts from the latest iterate.  Iner_Itn = 1
  reproduces the previous single-pass scheme.  Two-phase: iteration m > 1
  re-advects alpha from alpn with the corrected fluxes of iteration m-1.
  Cost: ~2x per step; benefit: the pressure of iteration 2 starts from a
  converged field (itp 411 -> 190 on the cylinder) and the divergence
  after the step drops 2-3x (BFS, cylinder).

12.2 Convergence criteria (u_momentum.F90, vof.F90, poisson.F90)
  * Momentum: coefficients assembled once per solve into dom%mco(26,:,:,:,3)
    (diag, 6 compact neighbours, 12 cross-derivative coefficients, 6
    masked advecting fluxes, frozen deferred correction); the sweep is
    a pure stencil (38 s -> 12 s per cylinder step on the dev VM).  The
    residual r = b + off - diag*u* of every face equation is evaluated in
    the sweep; criterion max|r| / ||b||_inf < Rsd_Lmt_Vel (b of all
    three components), Min/Max_Itr_Vel, weighted Jacobi SOR_fct_Vel
    (0.8 halves the sweeps on the cylinder: 91 -> 46; on the NACA case
    1.0 is faster: 24 vs 31).  Max and average residual per component.
  * 2nd-order deferred correction (conv_scheme 2): previously updated
    inside every Jacobi sweep - for CFL > 1/2 that iteration is not a
    contraction (1-D analysis: off-diagonal sum 3a vs diagonal 1/dt + a)
    and the cylinder switch matrix blew up (3.6E+29 at step 1).  Now
    evaluated once per solve from the latest iterate (mc(26)); the
    fully implicit 2nd-order solution is reached through Iner_Itn.
  * VoF: `vof_scheme = 2` (default) implicit CICSAM: face value
    (1-beta) alpha_D + beta alpha_A with beta = (a~_f - a~_D)/(1 - a~_D)
    from the CICSAM normalized face value; beta re-evaluated during the
    first 8 sweeps then frozen (with beta updated every sweep the
    residual stalls at ~1e-7); weighted Jacobi to the true residual of
    the scaled equation alpha + dt/V sum F alpha_f = alpha^n
    (Rsd_Lmt_VoF, Min/Max_Itr_VoF, SOR_fct_VoF); sub-cycling kept
    (vof_cfl_max) so the diagonal stays >= 1 - Co > 0.  `vof_scheme = 1`
    is the previous explicit CICSAM.  Both clip to [0,1] with the clipped
    volume monitored.  Cylinder: 12 sweeps per sub-cycle to 1e-9.
  * Pressure: both solvers return the residual ||b - A p||_inf/||b||_inf
    and the update norm max|dp| (non-dim.).  Jacobi: weighted update
    (SOR_fct_Prs; sign: the Laplacian diagonal is -1/wrk), criterion
    update norm < Trl_Lmt_Prs after Min_Itr_Prs.  FGMRES: criterion
    residual < Rsd_Lmt_Prs checked at every restart (true residual),
    Arnoldi loop left as soon as the recursive residual 2-norm is below
    0.5 Rsd_Lmt_Prs ||b||; update norm of the last cycle reported.
    Initial_Value_fct scales the starting pressure.

12.3 Renewal upgrades
  * Single-valued node output (output.F90, module node_groups): rank 0
    builds once the coincidence groups of all block-boundary nodes
    (heap sort on x + union-find with tolerance b2b_tol) and averages
    every output field over the block copies of each group before
    writing.  Cylinder grid: 6580 groups / 13300 copies (3-member groups
    on the O/H junction lines); verified max difference between block
    copies = 0 for u, v, w, p, alpha.  NACA O-grid cut: 405 groups.
    The solver-side ghosts at junctions are unchanged (wide-then-core
    exchange of the previous room); the ambiguity is resolved where it
    was visible, in the node interpolation.
  * Line symmetric Gauss-Seidel preconditioner (`poi_precond = 2`,
    poisson.F90 line_sgs/line_solve): M = (D_L + L) D_L^-1 (D_L + U)
    with D_L the tridiagonal part along the line direction chosen per
    block (largest sum of face coefficients), Thomas solves per line,
    forward/backward line ordering, block-Jacobi across sub-domains and
    interfaces (neighbours outside the owned range dropped) -> M is SPD,
    valid inside the CG.  Work arrays kpy, kpt.  Effect: cylinder step 1
    itp 2670 -> 411 (Jacobi -> line SGS), BFS 2000+ -> 1116, NACA ~1000
    (the NACA grid, stretched from 6e-4 m to 1.2 m, still needs ~1000
    CG iterations: multigrid is the next step for it).

12.4 Name list (read_input.F90 read_namelist)
  Hand-written parser of the grouped `$Group ... /end` format: comments
  after "!", case-insensitive names, comma/blank separated values,
  quoted strings, unknown names -> error with the line number.  Mapping
  to the internal variables: Case_name->casename, Char_L/U->Lch/Uch,
  Del_T->dt_dim, End_T (+N_max) -> nstep = min(N_max, nint(End_T/Del_T)),
  FreeSurf_Z->z_free_surface, Grvty->grav, Lou_w/mu_w/Lou_a/mu_a ->
  rho1/mu1/rho2/mu2, SW_Ref_Prs (0/1)->p_ref_switch, Loc_X/Y/Z, RefPrs,
  u/v/w_cel_initial, a_cel_initial->alpha_init (alpha of the cells above
  FreeSurf_Z), p_cel_initial, the $Solver_Iteration set (mom_*, poi_*,
  vof_*), Output_Intvl, Log_Intvl, and $Numerical_Scheme: conv_scheme,
  curvi_form, Cnvt_mode, Dffs_mode, cfl_switch, vn_switch, poi_solver,
  poi_precond, vof_scheme, vof_cfl_max, t_ramp, b2b_tol.  The old
  names (casename, Lch, rho1, ...) are accepted as aliases.  Defaults:
  Cnvt/Dffs 3, poi_solver 2, poi_precond 2, vof_scheme 2, Iner_Itn 1.

12.5 Log (Main.F90, mod_mpi.F90)
  Banner with the code name probCFD_부산, input summary, block splitting
  table (`PC_#: rank  BLK_#: b  i_th j_th k_th  i-th: n1~n2 ...`, node
  ranges), then per step (Log_Intvl): time, Max Local Total/<i>/<j>/<k>
  CFL# (dt |U^q| per direction on the field the step starts from; total
  = max over cells of the sum), Cell abs(DIV) {min/max/avg} in 1/s, the
  face monitors [I]/[J]/[K] (sweeps, max/avg residual) per inner
  iteration, the pressure monitor per inner iteration (iterations,
  residual, update norm), the CICSAM-VoF monitor, |div| after the
  correction, liquid volume / alpha range / clipped volume, and the
  step wall-clock time with the VoF / momentum / Poisson split.  The
  history file <Case_name>_monitor.dat holds one line per step.

12.6 Validation of this update (dev VM, 2 cores)
  V1 BFS regression (regression_BFS/README.md): precond 1 reproduces
     the previous divergence levels (6.532E-01 1/s = 1.82E-03 non-dim.
     vs 1.513E-03 before, momentum now converged to 1e-9); precond 2:
     itp 766-1122 with res(p) < 4e-8 where precond 1 hit the 2000 cap.
     Jacobi Poisson path runs (SOR 0.8).  Iner_Itn 2: itp 586+432 and
     |div|max 5.3E-04 ... 1.3E-04 1/s vs 2.9E-04 ... 1.6E-04 (one pass).
  V2 Two-phase cylinder, shipped namelist, 3 steps: |div|max 9.2E-07,
     3.2E-08, 4.2E-08 1/s; itm 91 (46 with SOR 0.8); itp 411, 363, 295;
     VoF 15 sub-cycles, 180 sweeps, residual 4.6E-10; liquid volume
     190.436 m3 constant; clipped 4E-9 m3.  Iner_Itn 2: itp 411+190,
     |div|max 3.6E-07, 1.7E-07, 3.5E-07.  conv_scheme 2: stable
     (|div|max 9.9E-07, 1.4E-07).
  V3 Switch matrix on the cylinder (1 step each, all run): phase 1 with
     the 5-line bcs (NOTE printed, alpha lines ignored), phase 1/2 with
     Jacobi Poisson, precond 1 (itp 2670), explicit VoF, Cnvt 1 + Dffs 2,
     Cnvt 2 + Dffs 1, conv_scheme 2, curvi_form 2, impulsive start,
     start from rest, SOR 0.8, Initial_Value_fct 0, SW_Ref_Prs 1 with a
     Dirichlet p (warning), bad switch value and misspelled variable
     (both abort with a message).
  V4 NACA66mod case, 20 steps: |div|max 1E-06 ... 1E-04 1/s, itm 31-33,
     itp 600-1100, CFL_j 1.8 (wall cells during the ramp), no
     divergence growth; node values single-valued across the O-cut.

12.7 Known properties added by this update
  * The pressure solve dominates on stretched grids (NACA: ~1000 CG
    iterations, 22 of 27 s per step on the dev VM).  Multigrid (or a
    2-D-aware direct line solver for k-thin grids) is the next upgrade.
  * The reported momentum residual with conv_scheme 2 is that of the
    linear system with the frozen deferred correction; the non-linear
    residual converges through Iner_Itn.
  * The implicit VoF residual is the one of the linear system with the
    frozen CICSAM weights (after 8 Picard sweeps).
  * N_max and a_cel_initial are interpreted as stated in README.md.
  * FGMRES results are deterministic per rank count; the Jacobi paths
    and the node averaging are rank-count independent.
  * Iner_Itn > 1 doubles the cost per step; use 2 with conv_scheme 2.

===========================================================================
13. probCFD_부산 UPDATE 2 (2026-09-04)
===========================================================================
13.1 Cleanup: N_max and Initial_Value_fct removed from the reader, the
  globals, poisson.F90 and all name lists; nstep = nint(End_T/Del_T);
  the Poisson solve starts from the previous pressure.  Seven empty
  lines precede every time-step block of the screen log.

13.2 Cnvt_scheme (renamed from conv_scheme; the old name is still
  accepted): 1 upwind-1, 2 upwind-2, 3 QUICK, 4 ENO-2.  All schemes are
  face interpolations in `upw` (convection.F90) on the 4-point support
  (U, D, A, and the next point) that the explicit term (conv_sum) and
  the frozen deferred correction mc(26) of the implicit directions
  share, so explicit / implicit / hybrid modes need no further code.
  QUICK: phi_f = 3/8 phi_A + 3/4 phi_D - 1/8 phi_U.  ENO-2: phi_f =
  phi_D + 1/2 minmod(phi_A - phi_D, phi_D - phi_U).

13.3 ENO form: the attached 2-D formula is the derivative (Hamilton-
  Jacobi / level-set) form  u dphi/dx ~ u [ (phi_i - phi_{i-1})/dx +
  dx/2 m(D--, D+-) ]  with the minmod of the two SECOND differences
  choosing between the central and the one-sided second-order stencil.
  In the code's conservative flux form the same stencil choice would
  need phi_{i-2} at the plus face and phi_{i-3} at the minus face of a
  face unknown (a 6-point support, beyond the 3 ghost layers exchanged
  across interfaces), and the derivative form itself is not
  conservative on the staggered curvilinear operator (the advecting
  fluxes ap/am are face fluxes, not a cell velocity).  The delivered
  ENO-2 therefore selects the smoother of the two FIRST differences at
  the donor (minmod slope): same order, essentially non-oscillatory,
  TVD, conservative, 4-point support.  On smooth data both choices
  coincide (they select the same stencil whenever the second
  differences have equal sign and ordering).

13.4 Pole boundary (bcs_b2b.F90, read_input.F90, curvi.F90, guards):
  * detect_pole: a face is degenerate when all its nodes coincide
    along one tangential index (within b2b_tol) - line pole; along
    both - point.  Pole faces are skipped in the pair matching and
    exported as type-93 records (12-line layout, donor = itself).
  * b2b_read_file separates type 92 (b2b) and 93 (pole: n_pole,
    pole_blk, pole_dir, pole_rng).  read_bcs reads the .bcs_b2b file
    first, sets blk(b)%pole(side,axis) and appends one automatic patch
    per pole face with all flags 0 (zero gradient for u,v,w,p,alpha);
    a .bcs record on a pole face is allowed but overridden (warning).
  * build_metrics: half positions on a pole plane get Jh = 0 and mh = 0
    (instead of the zero-Jacobian abort).  Consequences: the pole face
    unknowns are boundary faces (never swept), the face fluxes
    Jh/fvw U vanish, Gm/Gcap coefficients vanish, the face is type 1
    (fixed zero flux) for the projection, the outflow correction sees
    zero flux.  Guards for the zero Jacobian: fvw returns 1 where
    Jh = 0 (curvi_form 2), gv returns 0, the hybrid CFL test divides by
    max(|Jh|, tiny).
  * Ghost cells beyond the pole come from the linear extrapolation of
    fill_ghost_coords, i.e. the reflection through the axis (the
    opposite side of the pipe), with a negative Jacobian that is only
    used in property averages.  The ghost VALUES are zero-gradient
    copies of the first cell of the same column (not of the opposite
    column): exact for uniform fields and for axial flow, second-order
    otherwise, and independent of the azimuthal cell count (the PoleCyl
    grid has 39 azimuthal cells, so no cell has an exact opposite).
    An opposite-cell exchange would need interpolation across the
    azimuthal cut and is left as a possible refinement.
  * The node output averages all azimuthal copies of each axis node
    (1200 groups on PoleCyl: 40 axis points x ... plus the O-cut).

13.5 PoleCyl case: see PoleCyl_case/README.md.  200-step test: no
  NaN, divergence 1e-9 .. 1e-7 1/s, plug core at 1.0007 m/s, axis
  values single-valued, pressure uniform over the cross section,
  entrance-region cross velocities of 0.02 m/s, acceleration head
  rho L dU/dt during the inlet ramp.

13.6 Tests of this update: BFS 5 steps x {scheme 1..4} x {mode 1,2,3}
  all run (residuals 1e-9, itp 755-766, |div|max 1.6e-4 1/s); two-phase
  cylinder 1 step with QUICK / ENO in implicit and hybrid mode and
  QUICK implicit with Iner_Itn 2 (see the delivery message for the
  digits); PoleCyl 200 steps.

===========================================================================
14. END OF THE 부산 ROOM (2026-09-04): HANDOVER
===========================================================================
Package = this folder: sources, Cylinder / NACA66mod / PoleCyl /
regression_BFS cases, reference/ (user's name-list and log references,
inherited record, curvilinear study notes), tools/ (switch-matrix
scripts, node single-valuedness check, monitor plot).  HANDOVER.md at
the top level summarizes the state, the decisions taken in this room
(residual criteria, SW_Ref_Prs shift, removed variables, CFL limits of
the explicit higher-order schemes, ENO flux form, pole treatment,
output node averaging) and the next steps (multigrid preconditioner,
2nd-order time integration, turbulence model, opposite-cell pole
exchange).  Answer of record: ALL convection schemes are evaluated in
the conservative flux form (face value x face flux, telescoping),
explicit, implicit (1st-order upwind matrix + frozen deferred
correction) and hybrid alike.

===========================================================================
15. probCFD_부산 UPDATE 3 (2026-09-05/06): OVERSET (CHIMERA) GRIDS
===========================================================================
Task (user): implement an overset mesh algorithm. (1) review the example
overset solutions (hydrofoil, KCS, NACA66mod, NSWC15E; CFDShip-Iowa grid
+ solution + .bcs + Info_README), (2) block faces defined neither in
.bcs nor .bcs_b2b are overset boundaries, (3) blanking / fringe / active
cells without discontinuity between blocks (6-DoF relative motion:
future work), (4) solid body = closed set of no-slip surfaces, blank all
grid points inside; 4.1 closed surfaces with surface overlap evaluated
with a tolerance Surfovl_tol (hydrofoil, KCS), 4.2 closed boundary edge
of the wall touching the end of the background grid, same tolerance
(NSWC15E, NACA66mod), (5) symmetric grid -> symmetric blanking, (6) case
folder for NACA66mod_Ovst.grd (foil / refinement / background, all
overset-connected, 1 m/s along +x, water), (7) zip deliverable.

15.1 What the example files contain (examples/readp3d.py of the room)
  Fortran-unformatted little-endian multi-block PLOT3D.  hydrofoil /
  NACA66mod: <case>p3d.xyz (x,y,z double) + <case>p3dpfi.f (3 doubles:
  p, f, iblank = "F1V3").  KCS / NSWC15E: p3dxyzi (x,y,z float32 + int
  iblank) + p3dpuvwf (32-byte header record + 5 float32 variables).
  iblank is NODE based: 1 active, 0 fringe (interpolated), -1 hole.
  CFDShip .bcs types: 20 no-slip, 92 patched b2b, 61 symmetry, 40 zero
  gradient, 9/10 inlet, 11/12 exit, 13/28 far field.
  NACA66mod (3 blocks 344x123x5 / 120x5x120 / 241x5x80): SUGGAR overlap
  minimisation: the outer layers of the foil grid are holes with a
  2-layer fringe band inside them, the refinement block has 2 fringe
  layers on its outer boundary and a foil-shaped hole where the foil
  grid is finer, the background has a rectangular hole under the
  refinement block with 2 fringe layers.  hydrofoil: 7 blocks, 6 wall
  patches whose union is the body - the patches OVERLAP each other
  (task 4.1).  NSWC15E / NACA66mod: the wall edges lie on the symmetry
  / background planes (task 4.2).  KCS: 73 blocks, 58 wall patches.
  Target grid NACA66mod_Ovst.grd: 344x44x5 foil O-grid (O-cut i = 1 /
  344, outer boundary ~0.4 c away), 120x5x120 refinement (x -0.5..1.3,
  z -1.55..-0.15, cells 15 x 12 mm), 241x5x80 background (x -3..10.5,
  z -2.73..2.05), all with the same 5 y planes 0..0.2 m; the foil chord
  0..1 along x at z ~ -1 m.

15.2 Design decisions
  * Definitions: CELL status 1 active / 0 fringe / -1 hole (dominated by
    a finer grid) / -2 hole inside a solid body; FACE status 1 solved
    (>= 1 active neighbour cell, no hole neighbour) / 0 fringe
    (interpolated) / -1 hole.  Node iblank for output = max over the
    adjacent cells (clipped to -1).
  * Coverage per face cell (bcs patches, xp windows = b2b + periodic,
    poles); uncovered => overset boundary; Ovst_fringe (2) layers behind
    it (tangential radius 1) are "forced" (never solved).  Prepended
    automatic zero-gradient patches keep those ghosts finite.
  * Wall = bcs patch with u,v,w flags 1 and values 0 (no explicit type
    code exists in the probCFD .bcs format; inlets have non-zero values).
  * Bodies: union-find over wall patches connected by edge-node
    contacts (a) within Surfovl_tol of another wall patch (closest point
    on the global wall triangle soup, own patch excluded), (c) on a
    b2b/periodic window of the own block face (O-cut), (b) within
    Surfovl_tol of the triangle soup of all non-wall physical patches.
    Closed iff no unconnected edge node.  Only closed bodies cut holes.
  * Inside test: closest point on the body's triangles (buckets, ring
    search) + angle-weighted pseudo normal (Baerentzen-Aanaes); signed
    distance < -1e-9 L => inside; cell is a hole if its centre or any
    node is inside.  Handles overlapping double surfaces (the closest
    point is unique enough) and the sharp trailing edge (vertex pseudo
    normal points downstream).
  * Classification = implicit hole cutting (Lee & Baeder): every cell
    centre located in every other block (per-block Cartesian buckets of
    cell bounding boxes, ~1 cell per bucket, Newton inverse trilinear
    map, the most interior candidate wins - mirror symmetric); finest
    containing non-hole cell = donor; dominated iff vol(donor) <
    vol(cell)(1 - Ovst_vol_tol) (ties: both active -> symmetric).
    Fixed point of: fringe = non-active non-hole within Chebyshev nf of
    an active cell OR needed as donor-stencil member; a dominated fringe
    whose donor stencil contains a FORCED cell of the donor block becomes
    active (breaks mutual boundary-fringe cycles); safety pass: active
    with a hole within nf -> fringe (or orphan, reported).
  * Interpolation: trilinear in the donor's node-index coordinates
    (xi = i_cell + r) on the staggered lattices: offset 0.5 on cell-type
    axes, 1.0 on the face axis (u lattice at (f+1, j+.5, k+.5)); one
    routine for cells, u, v, w.  Invalid stencil members (holes, ghosts
    across uncovered faces) get weight 0 and the rest renormalised;
    ghosts across covered faces (BC mirrors, b2b data) are allowed.
    Fringe faces receive the full Cartesian vector (3 lattices) and
    U^q = fvw * m . u with the RECEIVER's face metrics; the Cartesian
    component of the face family is stored too.
  * Exchange: entries sorted by (sender rank, receiver rank), one message
    per rank pair (tags 100000+), donor domain = owner of the containing
    cell (stencil offsets -1..+1 lie inside its halo), all ranks build
    the identical table.  Face update after contra_to_cart (halo copies
    of u,v,w are then current); in correct_velocity a second refresh
    propagates the interpolated fringe faces to the halos.
  * Momentum: face_active = owned AND face status 1; cart_to_contra skips
    non-solved faces (keeps the interpolated U^q); fringe faces updated
    in every Jacobi sweep from the donor's latest star iterate.
  * Pressure: rows of fringe cells are -diag_P (p_f - I p) (diag_P =
    Jacobi diagonal, same scale as the Poisson rows), hole rows -p_h,
    rhsp = 0 at non-active cells, G = 0 on non-solved faces.  The
    single-grid preconditioner (PCG on the masked 7-point operator,
    Dirichlet at the fringes = block Jacobi across grids) makes FGMRES
    stall at |b-Ap|/|b| ~ 1e-2 (restart 20) resp. converge ~2x per 60
    Arnoldi steps: a minimal-overlap Schwarz iteration.  Fix: inner
    BiCGSTAB on the COUPLED 7-point + interpolation operator (fringe
    rows J diag_P (x_f - I x), hole rows J x, unmasked coefficients in
    dom%kc(:,:,:,4:6)) preconditioned by the block-Jacobi line SGS /
    Jacobi (masked coefficients dom%kc(:,:,:,1:3), diagonal kept exact
    through dom%ovd).  NACA66mod_Ovst step 1: 16-20 Arnoldi steps ->
    3.4e-8, |div|max after correction 8.7e-7 1/s.  Overlapping-Schwarz
    variants (dominated cells as overlap, restricted or not) were
    slower and are not used.  All-Neumann RHS mean removal is skipped
    with overset (the fringe rows are not conservative; pinning stays).
  * Output: iblank as the last variable of the .q (no averaging), the
    .ovst report; Ovst_only = 1 stops after write_solution(0).

15.3 Files
  new: ovs_geom.F90 (bkt_t buckets + locate_point + inv_trilinear +
  cell_volumes; surf_t + surf_add_patch/finish/closest/inside,
  tri_closest), overset.F90 (module overset: ovst_setup = face_coverage,
  add_auto_patches, wall_bodies, cut_holes, mark_forced,
  find_cell_donors, classify, ghost_status, build_entries,
  domain_status, write_report; ovs_finish_status, ovs_update_cells,
  ovs_update_faces, exchange, ovs_node_iblank, ovs_dump_cell).
  changed: mod_global (params, dom%cs/fs1-3/ovd/kc/krr, pointers,
  fstat/cstat), read_input ($Overset), convection (face_active /
  face_owned), curvi (cart_to_contra skip), b2b_exchange (refresh_all ->
  ovs_update_cells), u_momentum (ovs_update_faces(5)), poisson
  (cell_div_G(y, iv) rows, face_coefs / mask_coefs / apply_Mfull /
  bicgstab_c / vptr codes, correct_velocity face update, ovs_diag),
  Main (ovst_setup, Ovst_only stop, initialize_fields, active-only
  monitors), vof (active cells only), output (iblank), Makefile.
  Tool: tools/cfdship2probcfd.py.

15.4 Validation
  * BFS regression bit-identical (overset inactive).
  * NACA66mod_Ovst assembly: 5 uncovered face windows, 1 wall patch ->
    1 body CLOSED (edge nodes: 10 on the O-cut interface, 684 on the
    y = 0 / 0.2 symmetry patches within 1e-4 m, 0 unconnected), 2428
    solid holes; active/fringe/hole = 48916/2824/7256 (foil),
    40092/5940/10612 (refinement), 57144/2512/16184 (background);
    11276 fringe cells + 36611 fringe faces, 0 orphans; identical on
    all 5 spanwise planes; pattern = CFDShip-Iowa reference (room plots
    examples/naca_blank.png vs ovst_blank.png).  1.3 s.
  * Solver: momentum 26 sweeps (5e-10), pressure 3.4e-8 in 5550 inner
    iterations (2 FGMRES restarts), |div|max 8.7e-7 1/s after
    correction, 450 s/step on the 1-core sandbox (np = 1, 6 sub-
    domains round-robin).  See NACA66mod_Ovst_case/README.md.

15.5 Two bugs found and fixed by the parallel / self-test checks (2026-09-06)
  * face_valid read the face-coverage mask out of bounds for corner-
    ghost stencil members (two ghost directions at once) -> rank-
    dependent donor weights; corner ghosts are now invalid (as for
    cells).  ovs_selftest (poi_verbose = 1: analytic linear fields
    interpolated through the exchange) gives identical checksums for
    np = 1, 2, 3.
  * an overset-only "early exit" of the FGMRES Arnoldi cycle on the
    2-norm estimate made every cycle stop after ~3 steps (GMRES(3):
    20100 iterations to 2.7e-6); removed.  np = 2 now reproduces np = 1
    (5550 iterations, 3.6e-8, |div|max 9e-7 1/s).
  * hole cutting uses Surfovl_tol as the inside margin (points within
    the tolerance of the surface are on it): needed for bodies made of
    overlapping patches whose surfaces deviate by up to the tolerance.
15.6 hydrofoil closed-overlap test (tools/hydrofoil_assembly_test):
  6 overlapping wall patches -> 1 body; OPEN at Surfovl_tol 2 mm (14
  edge nodes 2.0-2.9 mm away, reported with coordinates), CLOSED at
  3 mm (712 edge nodes connected by overlap), 2337 solid cells blanked;
  orphans on this stride-3 grid explained by the coarsening (README).

15.7 SPSCFoil case (2026-09-06): surface-piercing strut + hydrofoil,
  two-phase, 7-block overset grid SPSCFoil.xyz (binary, 11.4 M cells).
  Added for it: binary PLOT3D reader (read_grid_binary in read_input,
  read_plot3d_binary in bcs_b2b; record-length based precision /
  iblank detection, 4-byte markers with wrap for > 2 GB records) and the
  grid-check fast path (Ovst_only = 1: ovst_setup right after
  fill_ghost_coords, write_grid_plot3d + write_iblank_plot3d, stop -
  no field arrays; SPSCFoil: 56 s + ~7 min hole cutting, ~1.5 GB).
  Hole cutting accelerated (node tests skipped when the centre is
  farther than half the cell diagonal), per-block orphan counts.
  bcs: j = 1 walls of blocks 1-6 (like the hydrofoil example), block 7
  inlet / outflow / free-slip sides+bottom / open top p = 0, alpha flag
  2 at the inlet; nml from the Cylinder two-phase case (Char_L 0.15,
  dt 0.5 ms, upwind, cfl_switch 1).  Assembly: 30 uncovered windows;
  6 wall patches -> 1 body; OPEN at Surfovl_tol 0.5 mm (114 strut top
  edge nodes 1.1-1.6 mm from the top cap: gap in the supplied surfaces,
  in air), CLOSED at 2 mm (2112 edge nodes by overlap), 9408 solid
  holes; 330k fringe cells + 977k fringe faces; orphans (7252 cells /
  2932 faces) at the strut / cap junction, documented in the case
  README.  Validation of the assembly against SUGGAR on the reference
  NACA66mod grid (read directly as .xyz): node iblank identical on
  99.6 / 97.6 / 97.0 % of the nodes of blocks 1 / 2 / 3.

15.8 (2026-09-07) The binary PLOT3D (.xyz) grid input added in 15.7 was
  REMOVED again on request: the grid is read as PLOT3D-ASCII only
  (read_grid in read_input.F90, read_plot3d in bcs_b2b.F90).  The
  grid-check fast path (Ovst_only = 1 before the field allocation)
  and write_iblank_plot3d stay.  tools/cfdship2probcfd.py converts
  binary CFDShip-Iowa grids to the ASCII format when needed.

16.5 (2026-09-08) Outlet pressure offset / wiggle - SOLVER FIX.  The
  user's 50-step results of the four NACA setups showed a uniform
  pressure offset in every case with a p = 0 pressure outlet (overset:
  -6.5 kPa over the whole domain, dropping to 0 in the last outlet
  cell; O-grid: -1 kPa; C-grid with p = 0 outlets: -80 kPa) plus an
  odd-even pressure oscillation in the last radial cells of the O-grid,
  while the pinned / mass-flux-corrected setups were at the physical
  level.  A mass-balance diagnostic (poi_verbose, correct_velocity ->
  outlet_balance) showed the STAR outlet flux 3-5 % above the inflow
  every step, corrected exactly by the projection.  Cause: apply_bc
  (boundary.F90) overwrote the boundary face of the star normal
  velocity with the first interior face at zero-gradient / outflow
  boundaries (f(nx) = f(nx-1)), although the momentum equation is
  solved at that face; the projection then removed the difference
  between the extrapolated flux and the momentum-consistent one every
  step: dp = rho du dx / (2 dt) in the outlet cell -> the level scales
  with the outlet cell size (22 kPa per metre of outlet cell here) and
  1/dt.  Fix: never overwrite the boundary face (momentum solution +
  projection on pressure outlets, + global scaling on outflow faces),
  fill the two ghost faces by copying it (also for the corrected field,
  which previously filled one ghost only).  Result: star/inflow
  mismatch 3e-4 instead of 3-5 %; O-grid p at the upstream boundary
  136 Pa (was -700), outlet cells ~0 (was -1 kPa with oscillation).
  BFS regression: CFL identical, |div|max changes in the 3rd digit,
  iteration counts within a few percent (new reference in README).

===========================================================================
17. probCFD_부산 UPDATE 5 (2026-09-13): OVERSET ORPHAN REDUCTION
===========================================================================
Task: with the 12 overset case folders prepared so far (KCS, DTMB, GPPH,
DCWJ, ONRT, JUSKFoil, JBC, KP505, NUMRT, LCS, NSWC15E, SPSCFoil), reduce
the number of orphan points by changing the overset code; the refinement
and background resolutions may be changed, the no-slip geometry not.

17.1 Diagnosis.  In every case the orphans have the same origin: a cell
  (or fringe face) of a body-fitted block A that lies next to a hole cut
  by ANOTHER body patch B, whose containing cell in every other block is
  itself a SOLID hole (its centre lies inside B, e.g. the coarse
  background cell that straddles B's surface, or a neighbouring body
  block's cell inside B).  find_cell_donors rejected such donors
  outright, so the cell had no donor at all; the safety pass then left
  it active next to zeroed hole cells (an "orphan"), and the same
  happened to fringe faces whose staggered donor lattice contained hole
  faces.  Changing the background / refinement resolution cannot cure
  this (refining moves the domination boundary and creates more holes,
  coarsening enlarges the solid cells) - the earlier per-case tuning
  moved the counts by +-30 % only.

17.2 Fix (overset.F90): nearest-valid-donor fallback.
  * find_cell_donors: when the containing cell in a candidate block is
    T_SOLID, nearest_nonsolid() searches its +-2 index layers for the
    non-solid cell nearest to the point; if no regular (containing,
    non-solid) donor exists in any block, the finest such cell becomes
    the donor with dnr = 0.5 (the point placed at its centre -> the
    cell-lattice weights are 1 on that cell: a nearest-neighbour, 0th
    order donor).  Fallback donors are flagged (dfb) and NEVER used for
    the domination test, so the hole pattern of the regular assembly is
    unchanged; they only turn orphan cells into fringe cells.
  * build_entries: a fringe face whose face centre lies in no non-hole
    cell of another block, or whose staggered lattices are invalid,
    gets nearest_nonhole() applied to the containing cell of every
    other block; the finest non-hole cell within +-2 layers is used and
    the face lattices are evaluated at its centre (the two faces of
    that cell, 0.5 / 0.5).  face_weights() factors the lattice/validity
    logic used by both paths.
  * Reported as "fallback (nearest-neighbour) donors: n fringe cells,
    m fringe faces"; the remaining orphans (cells that no other block
    contains at all) are still listed.
  Results (orphan cells / faces, old -> new) are tabulated in
  OVERSET_CASES.md.  No change for cases without such junctions
  (JUSKFoil, NUMRT, the tutorial cases: bit-identical assembly).

17.3 Case-folder notes.  The JBC.grd inside the returned JBC.zip was
  truncated (1.44 M of 6.26 M values); the intact grid of the JBC
  delivery is used.  All 12 folders now carry Input_polyMesh_probCFD.csv
  (1 CPU per block; adapt to the machine).

17.4 Start-up check of the 12 cases (solver mode) and a metric fix
  (curvi.F90).  Every ship case aborted in build_metrics with "zero/
  negative Jacobian" at a corner position: (i) the extrapolated corner-
  ghost coordinates below a wall at a block edge collapse, and (ii) the
  block corner NODE / EDGE positions of wedge-shaped blocks whose two
  faces lie on the symmetry plane (NSWC15E, NUMRT, JBC, ... block 1) are
  exactly degenerate - legitimately, the interior cells are fine.  The
  check now aborts only for a degenerate interior CELL centre (a real
  grid defect); degenerate ghost / node / edge / face positions are
  flagged and, after the loop, filled with the metric of the nearest
  interior cell (clamped half index).  Those positions are never solved;
  the count is printed ("n degenerate ghost / wedge-edge metrics
  replaced").  With this all 12 cases pass the metric stage; the wall-
  normal CFL of the ship cases at dt = 1 ms is 30-50 in the first cell
  (uniform 1 m/s initial field) - handled implicitly, but a smaller dt
  or a longer wall ramp is advisable for production runs.

===========================================================================
18. probCFD_부산 UPDATE 6 (2026-09-15): SOLID-BODY DETECTION FOR OVERLAPPING
    SURFACES, APPENDAGE / DUCT CUT-OUT, B2B REPAIR, BLOCK ORDER
===========================================================================
18.1 Hole cutting rebuilt (overset.F90 cut_holes, ovs_geom.F90):
  1. Surfovl_tol enters only the closure test and the overlap trimming;
     the inside test has no distance margin any more (the previous
     "shrink the solid by Surfovl_tol" left 3-6 mm shells of unblanked
     cells inside SPSCFoil / DTMB).
  2. Generalized winding number w(p) = sum of the solid angles of the
     body triangles / 4 pi, signed by the fluid-side normal (surf_wind_
     prep / surf_wind): exact van Oosterom-Strackee sum for the
     triangles of the fine buckets within +-2 of the point, dipole
     (area vector at the area-weighted centroid) for the other fine
     buckets of the 27 surrounding coarse cells and for all other coarse
     cells; per body; a triangle is exact iff its centroid bucket is
     near (no double counting).  Verified on a unit cube (1 inside, 0
     outside, +-0.49 on the two sides of an open sheet).  Inside: w >
     0.5 at the cell centre or w > 0.75 at a node (a node ON the surface
     gives 0.5, more at concave vertices; 0.6 cut 4 wall cells of the
     JUSKFoil pressure-side cusp).
  3. Surface overset (trimming): a triangle whose centroid lies within
     Surfovl_tol of a STRICTLY finer triangle of another patch of the
     same body, at a point farther than Surfovl_tol from that patch's
     boundary nodes, is excluded from the hole-cutting surface (keep =
     .false.).  Mirrored patches never trim each other.
  4. The winding number is evaluated only for the cut layer (cell
     centre within 1.5 cell diagonals of the surface, surf_closest mode 2
     per body) - evaluating every cell of the body bounding box took
     6 min per body on DCWJ.  Then every block is FLOOD-FILLED (6-
     neighbour BFS) from its non-wall boundary faces; cells not reached
     are holes: body interiors, appendage roots enclosed by the
     appendage sheet and the hull wall, duct penetrations.  Boundary
     cells whose centre has w > 0.5 for any closed body are holes, not
     seeds (the symmetry-plane face of a half body runs through the
     hull; without this the fill leaked into the hull interior).
  5. Half bodies closed by a planar domain boundary: hsurf = kept wall
     triangles + their mirror image across every DISTINCT planar
     domain-boundary patch the body touches (NSWC15E: 1 plane; before
     the deduplication the 8 block faces on y = 0 produced 8 images and
     371 k holes).  Normals mirrored; used for all hole-cutting queries.
  6. Reserved region of another wall block (figure "duct pipe": the
     whole mouth region must be removed regardless of resolution, hull
     may be several blocks): a cell of wall block b inside the grid of
     wall block ob is dominated (T_DOM -> fringe next to active cells,
     hole elsewhere) when ob's wall is a different surface there, i.e.
     the closest point of b's own wall(s) to the cell is farther than
     Surfovl_tol from every wall patch of ob (reserved_by_other).
     Coincident-wall overlaps (hull patches of one surface) keep the
     normal rules.  Needs the walls to form one closed body (surface
     overlap at the penetration), otherwise nothing is cut.
  7. Cells within Ovst_fringe layers of a block's own no-slip wall are
     never dominated (near_own_wall).
  8. face_coverage prints "WARNING: block b face ... is only PARTLY
     block-to-block" when an interface window does not cover the face.
18.2 Block-to-block interfaces (NSWC15E finding): the detector handles
  arbitrary axis permutations / reversals (block 2 i-max <-> block 3
  k-min, map 3 2 -1); the windows were j = 1..2 only because the wall-
  normal redistribution of the earlier rounds was done per block with
  different point counts, so the interface nodes matched only at the
  wall.  All j-min wall blocks of every case are now resampled to one
  common point count per case with the same geometric law along the
  existing lines (tools/fixcase.py); identical lines give identical
  nodes and the interfaces are full-face again.
18.3 Case folders: block order restored to the original one (background
  last); the cap / lid blocks of KCS, GPPH, DCWJ, ONRT, KP505 now
  precede the refinement and background blocks; every folder holds the
  grid check outputs <case>.x and <case>_iblank.q of this code.

18.4 Symmetry work (2026-09-16).  DCWJ (exact mirror grid, k reversed):
  * locate_point tie-break made mirror-invariant (larger x, then larger
    |y|, then larger z of the cell centre) - no effect on DCWJ, so the
    containing-cell search was not the source.
  * cells within Ovst_fringe layers of their own wall whose nearest
    surface point lies on that wall (or on a coincident patch within
    Surfovl_tol) are never cut by the winding number (own_surface_at):
    removed ~7 400 micron-level decisions per hull.
  * exact duplicates (equal-area coincident triangles) are trimmed too
    (keep the patch with more triangles, then the lower index).
  * closure test (d): an edge node that coincides with a node of another
    edge of the SAME patch (degenerate pole line, O-cut seam) counts as
    connected -> the duplicate cap / closing-strip blocks of KP505, DCWJ,
    GPPH, ONRT are no longer needed and were removed (tools/dropblocks.py);
    the "partly block-to-block" warnings of DCWJ disappeared with them.
  Status: background and refinement blocks are mirror-symmetric to the
  node; the mirrored hull / duct / nozzle blocks still differ by 1-3 %
  of their near-body nodes (13 488 mismatched nodes of ~1.2 M, was
  19 125).  The remaining mismatches sit in the duct / nozzle overlap
  bands where two wall sheets of one surface coincide within Surfovl_tol
  and the surviving sheet decides micron-level cut / no-cut for the
  cells of the other block beyond its first two layers.  Candidate fix:
  extend own_surface_at to all cells of a block (not only the first
  Ovst_fringe layers), or trim by "any finer triangle within Surfovl_tol"
  instead of the single closest one.  Open.

18.5 Three further requests (2026-09-16):
  1. Inside a solid body every grid POINT is a hole: (a) the node test of
     the cut layer uses w > 0.5 + 1e-3 again (a node inside the body,
     however close to the surface; own-wall cells are exempted by
     own_surface_at instead of the former 0.75 threshold); (b) the node
     iblank written to <case>_iblank.q is -1 for every node whose
     winding number is > 0.5 for a closed body, even where a neighbouring
     cell (centre outside) is active or fringe - only nodes at the hole
     boundary are tested; hsurf is kept alive for this.
  2. Boundary-layer cells stay ACTIVE: cells within Ovst_wall_layers
     (new name-list variable, default 8) of a no-slip wall of their own
     block are never dominated, never reserved (their own wall is the
     nearest one) and are not turned into fringe cells by the safety pass
     when the holes next to them are SOLID holes (cut by another body's
     surface: those holes act as that body's wall).  They are turned
     into fringes only next to dominated holes, and they are of course
     cut when they lie inside another body.
  3. KCS deck region: refinement block 72 refined back to 79 x 33 x 81
     and 2x more in z for -0.02 < z < 0.13 (79 x 33 x 141, dz 7 mm at the
     deck); background 76 refined 2x in z for -0.02 < z < 0.16 and 2x in
     x, y over the ship footprint (91 x 47 x 94, dz 12 mm at the deck).
     Result: see OVERSET_CASES.md (KCS row).

18.6 NSWC15E block 2 i-max <-> block 3 (i-min + k-min) interface
  (figure "B2B_condition_Fringe_Point_issue"): region A had fringe cells
  because the interface was not detected at all: the two sides had been
  resampled from lines with different point counts (block 3 from 19,
  block 2 from 20 points), so the interface nodes differed by up to
  2e-5 m (> b2b_tol = 1e-8) and both faces fell back to overset
  boundaries.  Blocks 1-4 were rebuilt in ONE step from the original
  grid with a common 20-point, 10x first-spacing distribution; the
  interface nodes now match to 3e-9 and the detector finds block 2
  i-max <-> block 3 i-min (k 63..111) and <-> block 3 k-min (k 1..63,
  map 3 2 -1) over the full j range.  Region A is active like region B
  (block 3 corner cells iblank 1; the L-shaped two-window connection
  needs no code change).  New diagnostic in bcs_b2b.F90: faces whose
  nodes coincide only within 1000 x b2b_tol print a "nearly coincides
  ... NOT connected" warning at start-up.  Rule for grid preparation:
  never resample the two sides of a block-to-block face from different
  source resolutions.

18.7 Holes inside another block's boundary layer (always on; replaces the
  Ovst_strict_cut = 1 option, which was removed).  Cells of any block
  whose centre lies inside the near-wall band (the Ovst_wall_layers cells
  next to a no-slip wall) of ANOTHER wall block are holes and never
  fringe cells (new cell type T_BLH, status -1; excluded from donors and
  donor stencils; cells next to their own wall are exempt).  That band is
  computed by the grid carrying the wall; blanking instead of fringing
  there keeps the fringe boundary of coarse Cartesian grids one layer
  away from the body surface, so no background / refinement fringe cell
  touches the near-wall band (the "fringe cells inside the propeller"
  picture of KP505: those cells were outside the solid - winding number
  <= 0.21 at every fringe node of the refinement block near the x = 0
  slice - but inside the blades' near-wall band; a blanked cell is up to
  one cell larger than the solid, which makes them look inside).
  KP505: 1070 such cells; the flood fill / cut counts are unchanged.
  Note: KP505 block 1 is a 7 mm long hub CAP (r <= 3.2 mm); the hub /
  shaft cylinder is not part of this grid, the region inside r < 10 mm
  around the axis is fluid.

18.8 Orphan reduction in the current algorithm (2026-09-17).
  * Active cells without donor whose surrounding holes (within Ovst_fringe)
    are SOLID holes only are not orphans: no other grid covers them (the
    background above a deck, cells at a symmetry plane next to the hull);
    they stay active with the solid holes acting as the wall of that body
    ("solid-bounded" in the report).  True orphans are now only cells next
    to dominated holes without donor.
  * Fringe faces without donor that are adjacent to a solid hole become
    wall faces (status -1, zero velocity) after the donor search; doing
    this before the donor search (as a general rule) removes valid lattice
    members near the body and creates new orphan faces (SPSCFoil 18 ->
    1669), so the conversion is applied only to the donor-less faces.
  * Fringe faces lying on a covered block face (physical boundary patch or
    block-to-block window) get their value from the boundary condition /
    exchange and need no donor ("boundary faces" in the report).
  * poi_verbose = 1 lists the first 40 orphan faces (q, block, i, j, k,
    xyz, cell status) - used to find that NSWC15E's remaining 81 orphan
    faces were on block 4's i-max face, which lies on the symmetry plane
    (|y| < 2e-6) but had no symmetry patch; the patch was added to
    NSWC15E.bcs.  Result NSWC15E: 36 / 103 -> 0 / 0, SPSCFoil 18 / 18 ->
    0 / 2.

18.9 Blank value in the PLOT3D iblank file (2026-09-20).  The node blank
  value written to <case>_iblank.q is 1 (active), -1 (hole), 0 (fringe)
  - except that a fringe node DIRECTLY ADJACENT to an active node (one of
  its six index neighbours within the block is active) is written as 0.5,
  marking the interpolation boundary that the active region touches.
  Implemented in output.F90 write_iblank_plot3d (real(4) array built from
  the integer node iblank of ovs_node_iblank; the assembly itself is
  unchanged).  Verified bit-for-bit against the equivalent post-processor
  tools/iblank_half_fringe.py on JUSKFoil and NSWC15E; the iblank files
  of the other 18 case folders were converted with that script (same
  result as re-running the grid check).  Adjacency is evaluated inside
  each block; across block-to-block interfaces the neighbour of the other
  block is not consulted.

18.10 Second fringe layer in the iblank file (2026-09-20): a fringe node
  with blank value 0 that is directly adjacent (six index neighbours) to
  a 0.5 node is written as 0.25.  Values: 1 active, 0.5 fringe next to
  active, 0.25 fringe next to 0.5, 0 other fringe, -1 hole.  Implemented
  as a second pass in write_iblank_plot3d (on a copy of the 0.5 field, so
  the 0.25 layer does not propagate); tools/iblank_half_fringe.py applies
  the same two rules and gives bit-identical files (checked on KP505,
  whose case folder holds the grid check output of this binary).

===========================================================================
19. probCFD_부산 UPDATE 7 (2026-09-20): GEOMETRIC VoF (vof_scheme = 3)
===========================================================================
19.1 New transport option in vof.F90 (advect_alpha_plic / plic_sweep /
  liquid_through / plic_plane / cube_frac / plane_const / slab_volume),
  selected by  vof_scheme = 3  in the $VoF group (1 explicit CICSAM,
  2 implicit CICSAM unchanged):
  * PLIC reconstruction in every mixed cell (0 < alpha < 1) in the cell's
    computational unit cube: Youngs normal from the 3x3x3 index stencil
    of alpha (weights 1/4, 1/2, 1/4), axes mirrored so that all normal
    components are >= 0, normalized to sum 1; plane constant from the
    analytical volume relation of Scardovelli & Zaleski (JCP 164, 2000)
    inverted by bisection (48 steps).  The reconstruction lives in index
    space, i.e. the same approximation as other curvilinear PLIC codes:
    exact on Cartesian / parallelepiped cells, first-order on skewed
    cells.
  * Direction-split advection: three sweeps per sub-step (order rotated
    x-y-z / y-z-x / z-x-y), per face the liquid volume of the flux slab
    of the donor cell (slab width = face volume flux x dts / cell
    volume, evaluated with the same volume relation on the rescaled
    cube), and the Weymouth-Yue divergence correction (JCP 229, 2010):
    alpha <- alpha - (L+ - L-)/V + c (F+ - F-) dts/V with c = 1 where
    alpha > 1/2 at the start of the sweep.  With the discretely
    divergence-free face fluxes of the projection the three sweeps sum
    to a conservative update; boundedness needs face Courant <= 0.5,
    enforced by the sub-cycling (nsub = ceil(Co_max / min(vof_cfl_max,
    0.5))).  Ghost / interface / overset fringe values are refreshed
    (refresh_all) after every sweep, so multi-block and overset faces
    see the intermediate field of the neighbour block.
  * Only active cells (cstat = 1) are updated; fringe cells take the
    interpolated donor value as before.  alpha is clipped to [0,1]
    (volume reported in the monitor as before) and values within 1e-8
    of 0 / 1 are snapped.
  * The monitor prints itVoF = 0 and res(VoF) = 0 for this scheme.
19.2 Test (Cylinder_Overset_2phase, 8 steps, dt 5 ms, 4 blocks incl.
  overset): liquid volume 191.132 m3 conserved to the printed digits for
  both scheme 2 and scheme 3, alpha in [0,1], clipped volume 1.3e-4 m3
  cumulative (7e-7 of the liquid volume; scheme 2: 8e-8 m3), same
  interface band width at this early time; 4-15 sub-cycles per step at
  the initial transient.  Cost: about 30 s per step on one core for this
  case (3 sweeps x nsub, each with a ghost refresh).
19.3 Not included: the level-set coupling (curvature from a distance
  function) and a geometric overset fringe reconstruction; fringe cells
  receive alpha by trilinear interpolation as for the algebraic schemes.

19.4 Geometric overset fringe reconstruction of alpha (2026-09-21).  For
  vof_scheme = 3 the overset fringe cells no longer receive alpha by
  trilinear interpolation: at assembly (build_entries) every regular
  cell entry stores the donor cell that contains the receiver centre
  (e_dcel) and the receiver cell's node box in that donor cell's natural
  coordinates (e_box, 8 locate_point calls per fringe cell); at run time
  the donor side (exchange, kd = 0, iv = 8) evaluates alpha_geometric:
  Youngs normal from the donor's 3x3x3 stencil, plane constant from the
  donor's alpha, and the liquid fraction of the receiver box on the
  liquid side of that plane (plic_box_frac, the plane extended beyond the
  donor cell).  The geometric value is used when the donor cell is mixed
  (0 < alpha < 1) and the box spans at most two donor cells per axis;
  otherwise (pure donor cell, donor much finer than the receiver,
  fallback donor, receiver node outside the donor block) the trilinear
  value is kept.  p and the pressure rows keep the linear weights.  The
  PLIC kernels (plic_cube_frac, plic_plane_const, plic_box_frac,
  plic_normal27) moved to ovs_geom.F90 and are shared with vof.F90.
19.5 Cylinder overset grid: the uniform z distribution (35 nodes,
  dz = 0.103 m) of Cylinder_Overset_1phase / _2phase was replaced by a
  two-sided sinh stretching with a node exactly at z = 0, dz = 0.040 m at
  the free surface on both sides, 0.262 m at the bottom and 0.239 m at
  the top (19 intervals below, 15 above; growth ratio <= 1.18); node
  count and block dimensions unchanged, all three blocks share the
  distribution so the interfaces and the overset assembly are unaffected
  (1130 cut cells, 0 orphans; .x / _iblank.q regenerated).  The monitor
  reference files of the two folders are for the old uniform grid.
19.6 Check-up (Cylinder_Overset_2phase, vof_scheme = 3, 20 steps of
  5 ms, new grid): alpha in [0, 1] throughout, clipped volume 3.9e-4 m3
  cumulative, liquid volume 191.132 -> 191.129 m3 (1.6e-5 relative) -
  the residual drift is the non-conservative fringe coupling where the
  free surface crosses the overset boundary; the interface band stays
  at two node layers (2409 / 2345 / 4881 mixed nodes of 38 500 / 38 500 /
  128 030 in the cylinder blocks / background), 1-15 sub-cycles per
  step.  The implicit CICSAM run (scheme 2) on the same grid keeps
  191.132 m3 over its first steps but is several times slower here.

19.7 Third fringe layer in the iblank file (2026-09-21): a fringe node
  with blank value 0 directly adjacent (six index neighbours) to a 0.25
  node is written as 0.125.  Values: 1 active, 0.5 / 0.25 / 0.125 the
  first / second / third fringe layer counted from the active region, 0
  other fringe, -1 hole.  Third pass on a copy of the 0.25 field in
  write_iblank_plot3d; tools/iblank_half_fringe.py applies the same
  three rules (bit-identical, checked on KP505; the other 19 folders
  converted with the script).

19.8 Log monitor label (2026-09-21): the per-step VoF monitor block in
  Main.F90 printed the hard-coded title "< CICSAM-VoF >" for every scheme
  although the transport already branched on vof_scheme (advect_alpha ->
  advect_alpha_plic for scheme 3).  It now prints "< CICSAM-VoF (explicit)
  >", "< CICSAM-VoF (implicit) >" or "< Geometric PLIC-VoF (direction-
  split, Weymouth-Yue) >" with the sub-cycle count and, for scheme 3, the
  effective face CFL limit min(vof_cfl_max, 0.5); the start-up banner
  prints the same limit.  Label only - no change of the transport.

===========================================================================
20. probCFD_부산 UPDATE 8 (2026-09-21): DIVERGENCE GUARDS, ACTIVE OVERLAP,
    GEOMETRIC FRINGE RECONSTRUCTION WITHOUT PLANE EXTENSION
===========================================================================
20.1 Divergence at t ~ 1 s (cylinder, both VoF schemes).  The logs show
  a sudden NaN in the pressure solve (20 000 iterations, residual NaN)
  one step after a perfectly normal step (residuals 1e-8, div 1e-7):
  a breakdown of the Krylov solver, not a physical instability.  Guards
  added in poisson.F90: inner BiCGSTAB stops on relative breakdown
  (|r0.v| < 1e-30 rho, |rho_new| < 1e-30 |r0|^2, |omega| < 1e-30,
  |beta| > 1e30) or NaN and returns its last finite iterate; FGMRES
  drops a Krylov direction whose norm or Givens rotation is degenerate
  or NaN, guards the back substitution, and skips a cycle whose Krylov
  solution contains NaN (reported once); a zero pressure diagonal
  (isolated cell) gets no update instead of 1/0.  A long run of the
  cylinder case with the geometric VoF is in progress to confirm.
20.2 Blanked gaps in the plots (empty ring at the O-grid's outer face):
  (a) a FORCED fringe cell (the layers behind an overset face) no longer
  dominates the cells of the other grid under it - domination is decided
  by the finest containing cell that can be active - so the background
  stays active under the O-grid's forced layers; (b) Ovst_active_overlap
  (new name-list variable, default 1): after the fixed point of rules A
  and B a dominated fringe cell next to an active cell of its own block
  becomes active as well (rule C, applied once, then A/B iterated
  again), so that the active regions of two overlapping grids overlap by
  one cell; the interpolated band starts one layer further out.
  Cylinder overset: 2 652 overlap cells; on z = 0 the O-grid's active
  nodes reach r = 1.3 and the background's come in to r = 0.8-1.0 -
  blanking with < 1 shows no gap.  Cells of the overlap carry weight 0.5
  in the liquid-volume monitor (ovs_overlap_weight); the reported volume
  is therefore higher than before (the ring under the forced layers was
  never counted).
20.3 Unphysical bubbles at start-up with the geometric VoF: they came
  from the geometric fringe reconstruction extending a donor cell's PLIC
  plane beyond that cell into its neighbours (a receiver box of up to
  two donor cells).  The reconstruction now splits the receiver box at
  the donor cell boundaries and takes each donor cell's own liquid
  fraction of the overlap (its plane restricted to its cell, 0 / 1 for
  pure cells); no plane is extended, so an interface cell can no longer
  put liquid into the neighbouring air cells.

===========================================================================
21. probCFD_부산 UPDATE 9 (2026-09-22/23): TURBULENCE MODELS, IDDES, WALL
    FUNCTION (specification turbulence_IDDES_wallfunction_spec.md)
===========================================================================
21.1 New module turbulence.F90 (after vof.F90 in the build).  Name list
  ($Control, before $Overset): TURB_model 0 none / 1 Spalart-Allmaras /
  2 k-epsilon / 3 k-omega (Wilcox 1988) / 4 k-omega SST (Menter 2003);
  IDDES_model 0/1; Max_Itr_Turb, Min_Itr_Turb, Rsd_Lmt_Turb,
  SOR_fct_Turb, Cnvt_scheme_Turb (1 first-order upwind, 2 second-order
  upwind as deferred correction), Cnvt_mode_Turb (1 explicit, sub-cycled
  to cfl_switch_Turb; 2 implicit), Dffs_mode_Turb (1 explicit, sub-
  cycled to vn_switch_Turb; 2 implicit).  All 20 case name lists carry
  the block with the defaults (model 0 = unchanged behaviour).
21.2 Coupling: the eddy viscosity mu_t/mu1 is a cell field (mut) and is
  added inside mu_h(), so the momentum viscous terms (and the Von Neumann
  check) see mu + mu_t without any change of the convection scheme.  mut
  is refreshed like alpha (field id 16); the model variables tk (k or
  nu_tilde) and tw (omega or epsilon) are field ids 14 / 15 with zero-
  gradient ghost values on physical faces, wall / inlet values set by the
  module (k = 0, nu_tilde = 0, Menter's omega at the first cell, free-
  stream values at velocity inlets: I = 1 %, nu_t/nu = 10, k-epsilon
  dissipation from a 0.07 Lch mixing length).
21.3 Models: standard forms (constants in the module header).  The
  production terms are treated implicitly (P = (P/phi_old) phi) with the
  growth per time step limited to 10x - the explicit form exploded in the
  start-up transient (SA mu_t/mu 1e19).  Production limiters P_k <= 10
  D_k for the two-equation models; SST with F1 / F2 / a1 limiter and
  cross-diffusion; SA with the fv1/fv2/fw functions and r <= 10.
21.4 IDDES (Shur et al. 2008): the RANS length scale of the destruction
  term (d for SA, k^1.5/eps or sqrt(k)/(beta* omega) for the two-equation
  models) is replaced by l_IDDES = fd~ (1+fe) l_RANS + (1-fd~) C_DES
  Delta_IDDES with Delta_IDDES = min(max(Cw dw, Cw hmax, hwn), hmax) and
  the standard fd~, fe, fb, fl, ft functions; C_DES 0.65 (SA), 0.78 /
  0.61 (k-omega branches).  TURB_model = 0 with IDDES_model = 1 gives the
  WMLES branch alone: nu_t = (Cs Delta_IDDES)^2 |S| (1 - exp(-y+/26))^2.
21.5 Wall function, bcs flag 11 on u, v, w: treated as the no-slip wall
  (Dirichlet 0) everywhere, recognised as a solid wall by the overset
  body detection (is_wall), and the wall shear stress of the first cell
  is taken from the law of the wall: u_tau from Spalding's profile
  (SA, none), from k (k-epsilon, viscous limit below y+ 11), from
  Menter's automatic wall treatment (k-omega, SST; omega1 = sqrt(omega_
  vis^2 + omega_log^2), k1 blended); imposed through the wall-face
  viscosity factor mu_eff = rho u_tau^2 y1 / u1 in mu_h(); first-cell
  model values nu_tilde1 = kappa u_tau y1, k1 = u_tau^2/sqrt(Cmu), eps1
  = u_tau^3/(kappa y1).  With TURB_model = 0 and IDDES_model = 0 the SST
  wall function (Spalding u_tau, no model variables) is used.  Flag 11
  on only one or two velocity components stops the run with an error.
  The monitor prints the y+ range and the mean Cf of the flag-11 faces.
21.6 Wall distance: nearest wall vertex by brute force, then the exact
  distance to the triangles that have a vertex within d_vertex + h_max
  (seconds for 1e5 cells; the bucket search with a global bound was
  minutes).  Tests (NACA66mod O-grid, Re 1e6, 3 steps each): SA, k-eps,
  k-omega, SST wall-resolved and with flag 11 (y+ 6.7-31, mean Cf
  1.27e-2 in the impulsive start-up), IDDES+SST, IDDES alone, wall
  function alone - all run stably without NaN; the flag-11 consistency
  error is triggered by a 11/1/1 surface.  Not yet done: a converged
  skin-friction comparison against the flat-plate correlation.
21.7 Also fixed: an out-of-bounds read at the outermost ghost layer in
  the velocity mirror of apply_bc (boundary.F90, cf() now clamps).

21.8 Flat-plate verification of the wall function (FlatPlate_WF, 2026-09-23):
  161x65x5 single block, Re_L 1e6, SST with flag 11 on the plate, first
  cell y+ 29-50.  Converged (t = 4 and 6 s identical): Cf within +2.5 /
  -3.1 % of Schultz-Grunow over Re_x 2.5e5-2e6 (mean -1.9 %), Prandtl
  1/7 within the same band.  Fix found with this case: the wall-function
  cells are now Dirichlet cells of the turbulence transport (their k1 /
  omega1 / nu_tilde1 / eps1 from the wall function are not overwritten by
  the transport sweeps); before, the leading-edge production entered the
  first cell and the boundary layer grew too fast (Cf -20 % downstream).
  New output <case>_wallfn_NNNNNNNN.dat with y+, Cf and u_tau of every
  wall-function cell, written with each solution file.

21.9 Two new case folders (2026-09-23): NACA66mod_3D (10 point-matched
  blocks, two-phase, wall function on the walls of types 20/24 of the
  source boundary file, symmetry 61, slip 28, inlet 9, outflow 12) and
  DelftTwist (3 overset blocks, single phase, wall function on the foil,
  slip tunnel walls); see their README files.  Interface detector
  (bcs_b2b.F90): the index map of a window is now derived from an
  interior node and a node whose matched index differs from the map is
  accepted when the node predicted by the map lies at the same position
  (folded donor faces: the C-grid wake cut of NACA66mod_3D block 4 made
  its window to block 5 fail with "inconsistent index mapping").

===========================================================================
22. probCFD_부산 UPDATE 10 (2026-09-23): CALCULATION-STOP CRITERIA
===========================================================================
  New module stop_criteria.F90 (checks called from Main.F90):
  1. NaN: after the momentum solve, the pressure solve and the projection
     / VoF / turbulence step, u, v, w, p, alpha, k / nu_tilde, omega / eps
     of every own cell are scanned; the first NaN is reported with the
     field name, grid block, global (i, j, k) index and position, then the
     run stops (MPI_Abort after the message is flushed).
  2. Residual / update norm: the momentum residual, the pressure residual
     |b-Ap|/|b| and the pressure update norm max|dp| are compared with
     Stop_Rsd_Max (name list, default 1e9); the solvers now record the
     block / index of their rank-local maximum residual (u_momentum.F90
     mom_res_loc, poisson.F90 poi_res_loc / poi_upd_loc) and the owning
     rank prints it with the position before the stop.
  3. CFL: the maximum local CFL number (cfl_numbers in Main.F90 records
     its cell, cfl_loc) is compared with Stop_CFL_Max (default 1e4); the
     block, index, position and velocity of that cell are printed before
     the stop.  This is the check that would have stopped the
     NACA66mod_3D run at step 2 (CFL 24 656) instead of sub-cycling.
  Both thresholds are in the $Control group before $Overset and echoed in
  the banner; every case name list carries the defaults.  Verified on the
  cylinder overset case by lowering the thresholds (criterion 3: CFL 7.05
  at block 1 cell (19,1,2); criterion 2: residual 3.6e-9 at block 3 cell
  (61,1,18), update at block 1 cell (1,1,1)).

===========================================================================
23. probCFD_부산 UPDATE 11 (2026-09-23): SURFACE TENSION (CSF)
===========================================================================
  New module surface_tension.F90; name list ($VoF): Surf_tension = 0 / 1,
  Sigma_ST [N/m] (default 0.0728).  Active only for Num_of_Phase = 2 (a
  single-phase run with Surf_tension = 1 prints "requested but inactive")
  and independent of the VoF scheme (the force depends on alpha only).
  * Curvature: alpha smoothed twice with a 27-point 1-2-1 kernel (field
    alps, id 17, ghosts refreshed), unit normal from the cell-centred
    gradient (metric m-matrix), kappa = -div(n) from the contravariant
    fluxes J m(:,q).n through the six faces (field kap, id 18); zero in
    pure cells.
  * Force: sigma kappa grad(alpha) / rho as an acceleration on every
    velocity face (compute_rhs_momentum, next to gravity), with the face
    gradient of alpha in the pressure-gradient stencil (direct difference
    across the face, averaged central differences tangentially) - the
    balanced-force form.  sigma_nd = Sigma_ST / (rho1 Uch^2 Lch) = 1/We.
  * Monitor: interface-cell count, max |kappa| [1/m], capillary time-step
    limit sqrt((rho1+rho2)/2 h^3 / (2 pi sigma)) on the smallest interface
    cell, with a warning when dt exceeds it.
  * Test: Cylinder_Overset_2phase, 3 steps with each VoF scheme (1, 2, 3):
    runs, max |kappa| 23-230 1/m, capillary limit 0.05 s (dt 0.005 s).
  Three latent defects found with the debug build on the way:
  - output.F90: the node-averaged solution was copied into the file array
    with a mismatched last dimension when the overset iblank is appended
    (array bound mismatch 6/5; now 1:nvar);
  - mod_global.F90 comp_axis: the new cell-centred field ids >= 14
    (turbulence 14-16, alps / kappa 17-18) were classified as face
    components in the block-to-block exchange (out-of-range xp%q index);
    they are cell fields now - this also affected the turbulence fields
    on multi-block grids;
  - grad_scalar (turbulence.F90) / grad_cell: assumed-shape dummy arrays
    lose the lower bounds, so the stencil was shifted by ng = 3 cells;
    the dummies are now declared with the domain lower bounds.

23.1 Initial drop shapes and the surface-tension test (2026-09-23): new
  name-list variables Init_shape (0 flat free surface, 1 circular
  cylinder along y, 2 sphere, 3 square prism along y), Drop_Center,
  Drop_Radius; drop_fraction in vof.F90 (5x5x5 sub-samples of the
  trilinear cell map).  Case Drop2D_SurfTension (40x2x40 cells, no
  gravity): square drop without surface tension stays square and at
  rest; with sigma = 0.0728 N/m it oscillates about the circle (period
  ~0.095 s vs 0.082 s linear theory for the n = 4 mode, coarse grid);
  circular drop Laplace pressure 9.30 Pa vs sigma/R = 9.10 Pa (+2.2 %).

===========================================================================
24. probCFD_부산 UPDATE 12 (2026-09-28): 6-DoF RIGID-BODY MOTION, DYNAMIC OVERSET
===========================================================================
  Request: "6DoF Motion Implementation Requirements" (SNAME / Fossen
  convention, > 10 bodies, object ids in row 6 of the csv, per-body
  $SixDOF group, imposed-motion file per body, forces without double
  counting on overlapping surfaces, moments about the centre of rotation,
  four force / moment files + one motion file per body every Data_Intvl,
  test cases NACA66mod_Overset (imposed pitch + heave), KP505 (imposed
  rotation), JBC (free heave + pitch), Box_6DoF (two free boxes)).

  Study of the materials (CFDShip-Iowa v5.5 cfdship_6dof.F90, Carrica et
  al. 2007 Comp. Fluids 36, Fossen 1994 ch. 2): CFDShip integrates
  Fossen's equations in the body frame with the CG offset (parallel-axes
  form 2.111 without products of inertia) in a 100-pass Jacobi loop, in a
  frame x aft / y starboard / z up (the grid frame), and mixes frames in
  its motion file.  probCFD_부산 keeps the grid frame for the flow and
  uses SNAME for everything the user sees (files, name list).

24.1 Design
  * Frames: grid G (x downstream, y starboard, z up); SNAME earth E =
    F G with F = diag(-1, 1, -1) (a proper rotation, so forces and
    moments map alike); body B aligned with E at t = 0, origin at O =
    Rot_Center.  eta = [x y z phi theta psi] (displacement of O in E,
    Z-Y-X Euler angles), nu = [u v w p q r] in B, eta_dot = J(eta) nu.
    Grid motion x = Rg (x0 - O0) + O0 + F eta(1:3), Rg = F R F.
  * Dynamics: M_RB nu_dot + C_RB(nu) nu = tau_B (Fossen 2.90/2.91 with
    the full I_O by the parallel-axes theorem 2.108, C from 2.72/2.87
    vector form), written as M_eta eta_ddot = J^-T(tau - C nu) + M_eta
    Jdot nu; fixed DOFs (eta component held) and imposed DOFs (eta from
    the file) are known rows, the free block is solved by a 6x6 Gauss
    elimination (no fixed-point iteration).  tau_B = fluid force /
    moment about O rotated to B + weight m g at the CG (+ optional
    damping).  Explicit coupling with the added-mass stabilization
    (M_eta + M_A) a^{n+1} = b + M_A a^n (Added_Mass per DOF, body frame).
    Integration: velocity Verlet on eta.  Body state broadcast from rank
    0 (replicated grids).
  * Imposed motion: <case>_body<k>_imposed.dat, cubic Hermite (Catmull-
    Rom) interpolation with analytic first / second derivatives, end
    values held outside the table.
  * Moving grids (rigid): blk%xn/yn/zn, blk%xg, dom%Xh from the stored
    originals; the metric vectors m(:,q) rotate (Rg m0), J is invariant;
    grid-velocity face components Ug1..3 = fvw m . u_g (u_g = v_O +
    omega x r) in the units of U1..3; every advecting flux is U - Ug
    (convection.F90 fflux, vof.F90 face fluxes, turbulence.F90 face_flux,
    CFL monitor); continuity unchanged (rigid motion sweeps no net
    volume); the stored face components of the rotated faces are updated
    incrementally, U += fvw (dRg m0) . u_vec (a translation leaves them
    untouched, so a body at rest reproduces the static solver exactly);
    no-slip wall velocity = grid velocity (boundary.F90 bv / mirror);
    wall-function velocity relative to the wall (turbulence.F90).
  * Dynamic overset: overset.F90 split into ovst_setup (coverage, auto
    patches, first assembly) and the re-entrant ovst_assemble /
    ovst_clear; ovst_reassemble after every displacement.  Caches of the
    motion-invariant parts (rigid objects): self hole-cut marks (the
    cells of the object's own blocks inside its own surface), trimming
    flags of the wall triangles, and per cell the finest donor / domination
    cell / fallback among the blocks of the SAME motion group (object id),
    revalidated against the current solid / BLH status.  Cylinder test:
    10.6 s -> 3.7 s per re-assembly (hole cut 7.4 -> 0.4 s).  Cells and
    faces that were holes and are not any more take the nearest
    previously valid value before the fringe interpolation.
  * Forces (sixdof_forces): every wall face cell of the object's wall
    patches: -p n A (first-cell pressure, n into the fluid) + mu_eff
    (u_t - u_wall,t)/d A (first-cell tangential velocity relative to the
    moving wall, mu_eff incl. wall-function factor); weight = fraction of
    the face's two triangles kept by the surface-overset trimming
    (ovs_wall_weights, from wsurf%keep), hole cells skipped; moments
    about the current O; MPI_Allreduce; grid -> E (F) -> B (R^T).
  * Files: <case>_force_body<k> / _moment_body<k> / _force_earth<k> /
    _moment_earth<k>.dat (friction, pressure, total; weight NOT included)
    and <case>_motion<k>.dat (k = body number; 32 columns: eta [m, deg], nu [m/s, deg/s], eta_dot,
    nu_dot, O and CG in the grid frame) every Data_Intvl steps; moving
    grids <case>_NNNNNNNN.x with every .q.
  * Inputs: csv row 6 (Object_ID, any label accepted, optional);
    $SixDOF group (Body, Mass, CG, Rot_Center, Motion_DOF, Inertia,
    Prod_Inertia, Added_Mass, Motion_Damp, Motion_File, Motion_Start_T);
    Data_Intvl ($Export_Interval); Grid_scale ($Numerical_Scheme).
  * Makefile order: ... curvi sixdof read_input vof turbulence ...

24.2 Verification
  * Static regression: Cylinder_Overset_2phase 10 steps identical to the
    update-11b binary (no object ids).
  * Zero-motion imposed body (re-assembly every step): monitor identical
    to the static run to 1e-8; drag -24177.9 N both.
  * Cylinder O-grid yawed at 0.001 rad/s: drag -24239.70 N vs static
    -24239.72 N (rotated metrics, incremental flux update, ALE fluxes).
  * Imposed sway 0.3 m amplitude: runs, added-mass reaction force in
    the opposite direction of the acceleration.
  * Box_6DoF: closed bodies (tube + 2 caps per box, 298 overlap-
    connected edge nodes), 0 orphans, 3000 one-cell-overlap cells.
    Full grid, 100 steps of 2 ms (after the fixes of 24.3): both boxes
    rise and right themselves (pitch 4.8 deg at t = 0.2 s, heave -2 mm,
    bodies 1 and 2 within 1 %), pressure force 1033 N up at t = 0;
    air-side velocities at the box edges grow (CFL 6 at step 100) ->
    the case is delivered with dt = 1 ms (see Box_6DoF/README.md).
  * NACA66mod_Overset, 60 steps: imposed heave / pitch reproduced
    exactly, lift 35.4 -> 39.8 N for pitch 0.13 -> 0.37 deg (thin-airfoil
    2 pi dalpha gives +4.1 N), CFL ~2.
  * KP505 (imposed rotation) and JBC (free heave / pitch): assembly with
    the object ids checked (bodies closed, 0 orphans; 1104 resp. 368
    orphan cells closed from same-block neighbours, 4223 resp. 8695 donor
    promotions); the runs need a proper machine.
24.3 Solver robustness fixes found on Box_6DoF (two-phase overset, fluid at rest)
  The Box_6DoF grids (O-grid tubes wrapping the box edges, caps abutting
  the tube ends partly block-to-block, refinement boxes, coarse tank)
  exposed five weaknesses that the earlier overset cases did not:
  1. Inner solver of the coupled overset pressure operator: the BiCGSTAB
     diverged (|s|/|r0| -> 1e100) on every two-phase overset system with
     the free surface present, already on a tank with one submerged
     refinement block (also with the update-11b binary).  Replaced by a
     restarted GMRES(20) (poisson.F90 gmres_c, basis dom%kW): monotone,
     no breakdown.  Cylinder_Overset_2phase: 84 instead of 237 FGMRES
     iterations per step, 7.8 instead of 14.2 s.
  2. Pressure interpolation across the free surface: a fringe cell in
     air next to a water donor cell (donor straddling the interface, or
     a nearest-neighbour fallback) received a water-like pressure; the
     active air cells beside it saw rho_w g dz / rho_a and moved at tens
     of m/s in the first step.  The fringe rows now interpolate the
     deviation from the still-water hydrostatic pressure, p_f = I(p -
     p_hyd) + p_hyd(z_f), p_hyd = rho_ref(z) |g| (z_fs - z), through the
     per-entry offset ovs_ph (overset.F90 ovs_build_phyd, rebuilt at every
     assembly; poisson.F90 cell_div_G, iv = 4 only).
  3. Orphan fringe cells / faces (no donor anywhere): they kept p = 0
     resp. a stale velocity.  build_entries now closes them with a
     zero-gradient copy of a same-block neighbour (active first, then a
     neighbour with a regular donor, then a neighbour filled in an
     earlier pass), reported as "orphans closed by a same-block
     neighbour" in <case>.ovst.
  4. Cartesian velocities of Dirichlet faces: contra_to_cart rebuilt the
     Cartesian components of the wall / inlet boundary faces and of the
     mirrored ghost faces from the transformed components with the
     forward tangents; at the wall edges of the tubes (families strongly
     non-orthogonal, degenerate ghost corners) the loop transformed ->
     Cartesian -> mirror -> transformed had a gain of ~20 per momentum
     sweep (10^12 after 14 sweeps).  curvi.F90 dirichlet_cart now sets
     those Cartesian values directly to the prescribed velocity
     (boundary face) and its mirror (ghost layers), with the moving-wall
     velocity through the procedure pointer wall_velocity_fn (set by
     sixdof_init, no module cycle).  tq() is now the exact inverse of the
     (possibly replaced) m matrix.
  5. Singular fringe cycles: chains of fringe cells interpolating each
     other across three or more grids (tube - cap - refinement -
     background) closed into cycles on which the interpolation
     identities are singular; the Krylov solve returned |dp| ~ 1e19 with
     a converged residual.  (a) classify rule A' (fixed point and on the
     final status): a DOMINATED, non-forced donor cell of a fringe cell
     with no hole within Ovst_fringe becomes active ("donor promotion");
     (b) the fringe rows of the pressure system carry a neighbour-
     average term, (1 + eps) p_f - I p - eps <p>_nb = c_f, eps =
     ovs_prelax = 0.2 (poisson.F90 fringe_row): every fringe cell is
     coupled to its own block, which reaches solved cells; for a smooth
     field <p>_nb and I p agree to O(h^2); (c) an active cell without a
     single active face (island between fringe bands) gets the same
     Laplace-fill row (isolated_active).
  With 1-5 the full Box_6DoF grid runs (2 cores: 45-55 s per step, of
  which 4.2 s re-assembly, 550-630 pressure iterations; initial pressure
  force 1033 N upward per box = rho g V/2 = 1030 N; the water stays at
  rest while the air above the surface shows the usual start-up noise
  of a light fluid next to interpolated fringes).  The coarsened grid
  (every second node) still fails after a few steps: its caps and tube
  ends are almost entirely forced fringe layers with no active cell left
  near the box edges - a grid, not a solver, limitation.

  Open: the impulsive start of an imposed motion (velocity jump at t = 0)
  makes the momentum Jacobi sweeps diverge at CFL ~ 9 in the wall cells
  (use a smooth start or the towing-tank ramp); strong (iterated)
  fluid-body coupling is not implemented (Added_Mass covers light
  bodies); the KP505 / JBC / Box_6DoF cases were only checked for the
  assembly and the first steps in the sandbox (see the case READMEs).

  Update 12b (2026-09-28): output file names carry the body number at
  the end - <case>_motion<k>.dat, <case>_force_body<k>.dat,
  <case>_moment_body<k>.dat, <case>_force_earth<k>.dat,
  <case>_moment_earth<k>.dat (was <case>_body<k>_motion.dat etc.).  The
  imposed-motion INPUT file keeps its name <case>_body<k>_imposed.dat.

24.4 Update 12c (2026-09-29): no empty spaces in the overset assembly
  Request: the z = 0 slice of Box_6DoF and the y = 0 slice of KP505
  (the user's own runs, plotted as "active cells only") showed EMPTY
  SPACES - fluid regions solved by no grid at all: the wedges between
  the corners of the O-grid tubes and the refinement boxes, the layers
  between the propeller blade grids and the refinement block, the cut
  layer where two body grids meet at the hub.  Condition given by the
  user: no grid point may ever "permeate" into a solid body - a coarse
  cell must not be kept active where it lies inside a body.

  Why the gaps existed.  Three rules of classify() together left cells
  to nobody: (1) the safety rule turned every active cell within
  Ovst_fringe of a hole into a fringe when it had a donor block at all,
  even when the whole donor stencil consisted of FORCED fringes of the
  finer grid (the cells behind an overset face, never solvable), so the
  coarse cell interpolated from cells that themselves interpolate from
  it; (2) a DOMINATED hole was decided from the finest containing cell
  before that cell's own final status was known (it could end as a hole
  or fringe of a third grid); (3) a cell whose centre lies in a solved
  cell of another grid was taken as covered even when parts of it (its
  nodes) lie outside every solved cell (background cells three times
  larger than the two forced layers of the refinement at the NSWC15E
  bow).

  Rules added (overset.F90 classify, all on the final status cs):
  a. Safety rule with stencil_active(): an active cell next to a hole
     becomes a fringe only if its donor stencil contains at least one
     ACTIVE cell of the donor block.  Otherwise it stays active, its
     dominated holes (-1) within Ovst_fringe that have a donor become
     fringes (interpolated from the finer grid), and if only SOLID
     holes (-2) remain around it the cell is "solid-bounded": the solid
     holes act as the wall (staircase).  Solid holes and boundary-layer
     holes (T_BLH) are never reopened - no cell inside a body is ever
     solved (the user's condition).  Count "solid-bounded" in the log
     and in <case>.ovst.
  b. Uncovered holes: a dominated hole whose centre lies in no active
     cell of any other grid becomes active (then rule a again); never
     applied to FORCED cells.
  c. Coverage pass: every non-active, non-solid cell must lie inside
     solved cells - its centre in an active cell of its donor block AND
     each of its 8 nodes, 12 edge midpoints and 6 face centres in an
     active cell of some other block (nodes_covered; the 8 nodes alone
     left slivers of the tube's fringe layer inside coarse cells at the
     Box_6DoF tube corners).  A cell that fails is made active itself (an
     overlap cell, weight 1/2 in the global sums), with its dominated
     holes within Ovst_fringe turned into fringes.  A FORCED cell that
     fails cannot be made active (its stencil would read the block's
     outer ghost layer, which carries no boundary condition): instead,
     for each of its 27 test points that lies in no active cell of
     another block, the finest containing cell of another block that
     can be solved (not solid, not a boundary-layer hole, not forced) is
     made active (cover_forced) - the coarse grid is solved under the
     outer layers of the fine one.  Nodes that lie in solid cells of
     the other grids need no cover (the body is there); nodes covered
     by nothing but forced cells are counted and reported.  Iterated
     until nothing changes (<= 50 sweeps).
     The first version of this pass made the forced cells themselves
     active: Box_6DoF then gave a pressure force of 988 N per box after
     the first projection instead of the hydrostatic 1033 N (the tube
     cells behind the overset faces read undefined ghost pressures).
     With cover_forced the force is 1034 N.
  d. Winding number (ovs_geom.F90 surf_wind_prep / surf_wind): a fast
     winding number in the form of Barill et al. (2018).  For every
     body a cluster tree is built over its kept triangles: a base grid
     over the body's own bounding box (about nt^(1/3) cells per axis,
     cells at least as large as the mean triangle), each triangle in the
     cell of its centroid, and levels of 2x coarsening above it; every
     cluster carries the area-weighted dipole, the area centroid, the
     bounding radius of its vertices and the second moment
     sum a_t n_t (c_t - c)^T.  A query descends from the coarsest level:
     a cluster farther than BETA = 2 radii from the point contributes
     its multipole expansion (dipole + second moment, the first two
     orders of Barill), a nearer one is opened, and at the base level
     the triangles of a near cell are summed exactly (van Oosterom -
     Strackee).  Before, one bucket grid over the whole surface served
     all bodies and the "far" set was fixed by bucket-index distance
     (dipoles beyond the 27 coarse cells / +-1 fine buckets around the
     point): on the DelftTwist foil (thin body, triangles larger than
     the buckets) that gave w = 0.4 - 0.7 at points 4 mm from the
     surface in bands one coarse cell apart, 1865 O-grid cells 1 - 9
     layers above the foil were cut as "inside", the wall cells under
     them were flood-filled, and 15720 cells around them became
     solid-bounded.  Now: 0 cells of the O-grid inside, 38 solid-bounded
     (the corners of the refinement box); the cut counts of the other
     cases change by a few cells only (GPPH body 2: 15866 instead of
     15583 - the nodes lying ON a coincident wall surface give w = 0.5
     exactly, and the node test now uses a margin of 1e-2 instead of
     1e-3 so that the ~1e-3 accuracy of the multipole sums cannot flip
     them).  The per-body trees matter for the cost: with one grid
     sized by the hull, the propeller of GPPH (26752 triangles of 1 mm
     in 46 cells of 11 mm) needed 7750 exact triangles per query, 410 s
     for that body alone; with its own tree 36 s.  GPPH assembly 243 s
     instead of 704 s (update 9c), DelftTwist 27 s instead of 37 s.
     The log now shows per body "winding number: N queries, M exact
     triangles, K clusters visited, T s".
  e. Bodies thinner than a cell (KP505 blade tips inside background /
     refinement cells; found on the y = 0 slice after a-d: the blade
     tips appeared covered by solved background cells).  The hole cut
     tests the centre and the nodes, which are all outside such a cell
     although the blade passes through it, and the blade grid's outer
     (forced) layers do not dominate it, so it stayed active - a solved
     cell reaching into the body.  Two rules: (1) cut_holes marks a
     cut-layer cell as a solid hole when the closest surface point of
     the body lies strictly inside the cell (inv_trilinear of q in
     (0.001, 0.999)^3): the surface passes through the cell ("crossed by
     the surface, thinner than the cell" in the log; KP505: 9009 of
     18162 cut cells); the body's own wall cells are excluded as before.
     (2) The coverage pass never activates a cell one of whose 27 test
     points lies inside a body (winding number > 0.51, tested for the
     uncovered points only; point_in_body), and the uncovered-hole rule
     never re-activates a hole whose centre lies inside a body: such
     cells stay holes / fringes, the body's own grid solves around the
     body ("uncovered cells crossed by a thinner body left to the finer
     grid" in the log).
  f. Permeation check (check_permeation, after classify, first
     assembly): every ACTIVE cell within 0.6 cell diagonals of a closed
     body's surface - except the body's own wall cells, outside by
     construction - is tested with the winding number at its centre and
     its 26 other test points; the log and the .ovst report "permeation
     check : N active cells with a point inside a body, M with the
     centre inside" and lists the first offending cells; grid-check
     mode only (it costs as much as the hole cut).  Results: table below.
     (The winding number is unreliable within a few cell heights of a
     hub - blade root junction, where the hub patch and the blade
     collar are two sheets a few microns apart: w = 0.6 - 1.2 at the
     hub's own wall cells there.  Those cells are protected by the
     own-wall rule; the check excludes them for the same reason.)

  g. Node blank export (ovs_node_iblank, 2026-09-30): a node adjacent
     to a SOLID hole cell is exported as -1 (before: the max over the
     adjacent cells, so a single background cell crossed by a thin blade
     tip - a hole - had all eight nodes at 1 and looked like a solved
     cell in the [Blank_var] plots; the user's KP505 slice raised this).
     The status is unchanged, only the picture: solid holes now show as
     a closed band around every body in every grid; the price is that
     the body-side half of an active cell next to a solid hole is drawn
     blanked in an "active cells only" node plot (Tecplot: use value
     blanking with the "all corners" option, then only whole hole cells
     are blanked).
  h. Exact cell blank (output.F90 write_cell_blank): beside every
     solution <case>_NNNNNNNN.q and the grid-check <case>_iblank.q the
     solver now writes <case>_NNNNNNNN_cblank.q / <case>_cblank.q - a
     cell-centred Plot3D q file (dimensions nx, ny, nz per block, one
     variable: 1 active, 0 fringe, -1 hole) with the exact cell status.
     tools/blank_slice.py uses it when present: every slice polygon is
     coloured and kept by its parent cell's status, so single hole cells
     show and solved cells next to holes stay whole (no half-blanked
     "specks").  KP505 y = 0 with it: blades and hub are clean holes in
     every grid, nothing else blanked.

  i. Self-overlap donors (2026-09-30): a block may now be its own donor
     - locate_from() searches a block in itself excluding the cells
     within 2 Ovst_fringe + 1 index layers of the receiver (locate_point
     got the optional exclusion xi / xd).  Used in the cell donor search
     (and its cache), the face donors, the uncovered-hole rule and the
     coverage tests.  Motivation: the KP505 tip-cap grids (blocks 3, 5,
     7, 9, 11) are C-type in k with their two ends overlapping by only
     two cells without a point match, so both ends are forced layers
     that could only interpolate from block 2's forced layers there.
     They now interpolate from the far side of their own block.  Cylinder
     and Box regressions unchanged.  What this cannot do: the seam band
     of those caps (about 0.3 mm at z = -0.985 on the y = 0 slice, next
     to the blade) is still covered by no ACTIVE cell - both ends of the
     cap are entirely within their forced layers there (the overlap is
     2 cells, the forced width is Ovst_fringe = 2, and Ovst_fringe must
     be >= 2 for the stencils).  The interface fix is in the grid: make
     the cap's k-ends point-matched (the seam then becomes a block-to-
     block interface like block 2's wake cut, detected automatically)
     or extend one end by at least 2 Ovst_fringe + 1 = 5 layers so that
     active cells of one end cover the forced layers of the other.
     The j = 1 wall layer of block 2 under those cap layers is not
     blanked as a boundary-layer hole any more when block 2's own cell
     there is forced (blank_in_boundary_layers).

  Log lines (also in <case>.ovst): "coverage : N fluid cells covered by
  no solved cell of any grid made active", "coverage : N cells of coarser
  grids made active under uncovered forced (overset-face) layers of
  finer grids", "coverage : N forced cells with a node covered by nothing
  but forced cells of other grids", "solid-bounded : N active cells next
  to solid holes with no solved cell of another grid there (... N of them
  covered only by forced fringes of a finer grid)".

  Checked on every overset case folder (Ovst_only = 1 with the final
  binary; <case>.x, <case>_iblank.q, <case>.ovst and ovst_check.log are
  in the folders): the y = 0 (or z = 0) slices of the node blank
  status, drawn as "active cells only" with tools/blank_slice.py
  (blank_check/*.png in every folder): no empty space left anywhere
  except the body interiors - Box_6DoF, NACA66mod_Overset(_1phase),
  Cylinder_Overset_1phase/2phase, DelftTwist, KP505, JUSKFoil, NUMRT,
  NSWC15E, DCWJ, DTMB, GPPH, SPSCFoil, KCS, LCS, ONRT, JBC.  Regression:
  Cylinder_Overset_2phase 10 steps - monitor values and iteration counts
  identical to update 12b (the classification of that case is unchanged);
  Box_6DoF 1 step - pressure force 1034 N per box (hydrostatic 1030 N).

  Overset check of every case folder with the update-12c binary (Ovst_only = 1,
  2 MPI ranks in the sandbox; the times include the permeation check, which
  costs as much as the hole cut and runs in grid-check mode only).  Coverage
  = cells made active by the coverage pass (self / coarse cells under forced
  layers / forced cells with a test point covered by nothing / uncovered
  cells crossed by a thinner body, left as holes); crossed = cut-layer cells
  cut because the surface passes through them (bodies thinner than the
  cell); solid-bounded = active cells whose holes are solid only, the holes
  act as the wall (in brackets: covered only by forced fringes of a finer
  grid); orphans = fringe cells / faces without donor (closed from a
  same-block neighbour); permeation = active cells with any of 27 test
  points inside a body / with the centre inside.

  | case | coverage self / under forced / uncoverable / thin-body | crossed | solid-bounded | orphans | permeation pt / centre | assembly |
  |---|---|---|---|---|---|---|
  | Box_6DoF | 124 / 494 / 0 / 0 | 770 | 512 (512) | 0 / 0 | 0 / 0 | 17 s |
  | NACA66mod_Overset | 0 / 0 / 0 / 0 | 1084 | 0 | 0 / 0 | 0 / 0 | 21 s |
  | NACA66mod_Overset_1phase | 0 / 0 / 0 / 0 | 1084 | 0 | 0 / 0 | 0 / 0 | 21 s |
  | Cylinder_Overset_1phase | 0 / 0 / 0 / 0 | 884 | 0 | 0 / 0 | 0 / 0 | 9 s |
  | Cylinder_Overset_2phase | 0 / 0 / 0 / 0 | 884 | 0 | 0 / 0 | 0 / 0 | 10 s |
  | DelftTwist | 29 / 107 / 0 / 0 | 1383 | 42 (42) | 0 / 0 | 0 / 0 | 56 s |
  | KP505 | 11334 / 8263 / 0 / 1544 | 9009 | 16608 (16608) | 0 / 0 | 0 / 0 | 170 s |
  | JUSKFoil | 0 / 0 / 0 / 0 | 284 | 0 | 0 / 0 | 0 / 0 | 34 s |
  | NUMRT | 44 / 1177 / 0 / 0 | 4352 | 6822 (6822) | 0 / 0 | 13 / 0 | 304 s |
  | NSWC15E | 3198 / 3991 / 64 / 588 | 4220 | 6792 (6756) | 0 / 0 | 10 / 0 | 226 s |
  | DCWJ | 28523 / 3771 / 0 / 4927 | 28757 | 83486 (83486) | 0 / 2 | 577 / 0 | 863 s |
  | DTMB | 22893 / 4837 / 0 / 2882 | 14081 | 84110 (83618) | 0 / 0 | 136 / 0 | 449 s |
  | GPPH | 48779 / 5225 / 0 / 2650 | 19713 | 71228 (71228) | 0 / 0 | 154 / 0 | 574 s |
  | SPSCFoil | 5143 / 1559 / 0 / 1858 | 2095 | 20906 (20852) | 0 / 0 | 126 / 0 | 408 s |
  | KCS | 67884 / 21254 / 0 / 5716 | 47931 | 102904 (102904) | 1 / 0 | 106 / 0 | 2130 s |
  | LCS | 62500 / 7922 / 0 / 3641 | 13348 | 45572 (45572) | 0 / 0 | 137 / 0 | 936 s |
  | ONRT | 88103 / 12931 / 0 / 6671 | 54654 | 135194 (135188) | 30 / 0 | 411 / 0 | 2886 s |
  | JBC | 775 / 189 / 0 / 444 | 9117 | 4658 (4658) | 0 / 0 | 10 / 0 | 872 s |

  Permeation: no case has an active cell with its CENTRE inside a body.
  The cells with an edge / face / node test point inside (0.01 - 0.03 %
  of the cells of the hull cases) sit where two wall patches of the same
  body meet imperfectly: the centre-plane joint of the two hull halves of
  JBC and ONRT (ONRT's log warns that this joint is only partly
  block-to-block: 190 + 57 face cells mismatch), the waterjet-nozzle
  junction of DCWJ, the strut - foil junction of SPSCFoil, the appendage
  roots of DTMB / GPPH / KCS / LCS: the winding number of a point between
  two nearly coincident sheets is ~1, so these are cells of the second
  and third layer above one patch whose test points fall behind the
  other, slightly protruding patch - grid joints to review, not fluid
  cells inside a body.
  The large solid-bounded counts of the hull cases (DCWJ, DTMB, GPPH, KCS,
  ONRT) are the cells of the coarser hull / refinement blocks next to the
  hull surface where the finer hull block that carries the wall has only
  its forced (overset-face) layers: before update 12c those cells were
  fringes interpolated from the forced layers (which interpolate back from
  them), now they are solved with the solid holes as a staircase wall.  A
  point to review on the hull grids: extending the wall blocks by the
  fringe width past their neighbours removes them.


25. probCFD_부산 UPDATE 13 (2026-09-30): SSP TIME SCHEME, DYNAMIC PRESSURE, RUNTIME NAME LIST,
    STOP REPORT, LOG LAYOUT, OUTPUT FOLDER
===========================================================================
  Request (six items): a 2nd-order SSP time integration behind a
  Time_scheme switch; the dynamic pressure in the two-phase volume files;
  name-list changes applied to the running calculation; a markdown stop
  report; the 6-DoF log layout of the user's figure and a 3x larger gap
  between step groups; step number on the initial grid, extensions .grd /
  .p3d and a Volume_solutions_grids sub-folder.  Delivery: sources only.

25.1 Time_scheme = 2 (2nd-order SSP, Heun)  [Main.F90, mod_global.F90]
  The fractional step S (VoF transport from alpha^n, momentum predictor
  from u^n, pressure projection, velocity correction, turbulence
  transport; Iner_Itn inner iterations) is the first-order forward
  operator of the existing scheme.  With Time_scheme = 2 every step
  computes U(1) = S(U^n) and U(2) = S(U(1)) (the stage-2 base copies un /
  vn / wn / alpn are the stage-1 results; the turbulence transport takes
  its old level at solve time, so it follows) and sets
      U^(n+1) = (U^n + U(2)) / 2
  for the Cartesian face velocities u, v, w, the contravariant fluxes
  U1..U3, alpha and k / omega / mu_t (ssp_save / ssp_average, copies
  s_u ... in the domain type, allocated on first use).  Both U^n and U(2)
  are discretely divergence-free on the same faces, so is their average
  (linear); halo, boundary, ghost and overset-fringe values are linear
  too, so no extra exchange is needed; the pressure is that of the
  second stage.  Moving grids: the body motion and the re-assembly are
  done once at the start of the step (stage 2 runs on the n+1 grid), the
  averaging then mixes the two metrics - second order in time for the
  fluid on static grids, first order in the grid motion.  Cost: two
  stages per step.  Test: Cylinder_Overset_2phase, 10 steps, Time_scheme
  1 reproduces the update-12c monitor to all digits; Time_scheme 2 runs
  with the same iteration counts, CFL 6.72 instead of 6.38 at step 1
  (different velocity history), divergence 1e-8.

25.2 Dynamic pressure  [output.F90]
  Two-phase solution files: variable 6 = p - p_hyd(z), p_hyd = rho_ref(z)
  |g| (z_fs - z) with rho_ref = rho_w below the still free surface and
  rho_a above - the still-water hydrostatic part that the overset
  interpolation subtracts (24.3 item 2).  Variables: u, v, w, p, alpha,
  p_dyn, (iblank).  The log lists them at the first export.

25.3 Runtime name-list updates  [read_input.F90 reread_namelist, Main.F90]
  Every step rank 0 reads probCFD.nml and hashes its content; when the
  hash differs from the last read, the lines are broadcast (MPI_Bcast,
  all ranks parse the same text - a file edited while several ranks
  read it independently could give different values) and each rank
  applies the runtime-changeable keys (list in README.md), with range
  checks (out-of-range values are ignored and reported).  Changed values
  are written to the log ("[ name list ] probCFD.nml changed - re-read
  at step N: End_T 5.0E-2 -> 3.0E-2, time steps 10 -> 6").  The physical
  time is now accumulated (t_now = t_now + Del_T, mod_global) so that
  Del_T may change; the step count becomes (istep - 1) + (End_T -
  t_now)/Del_T, and the marching loop is a while loop over nstep.  Static
  inputs are ignored on purpose (they would need re-allocation or a new
  assembly).  Test: End_T 0.05 -> 0.03 written at step 3 of the cylinder
  run: the run ended at step 6 (t = 0.03 s) and reported the change.

25.4 Reason_for_stop.md  [stop_criteria.F90]
  Module stop_criteria (now compiled right after mod_mpi so that
  read_input can use it) keeps a ring buffer of the last 40 completed
  steps (record_step from Main: CFL, momentum / pressure iterations and
  residuals, update norm, liquid volume, wall time) and the detail lines
  of every stop criterion (the reporting rank adds them; stop_run gathers
  them on rank 0 with MPI_Send / Recv before MPI_Abort).  The report has:
  status (diverged / input error / finished normally), the stop message,
  step and time, the detail lines (field, block, index, position,
  velocity), the history table, a diagnosis section per kind of stop
  (NaN, residual, CFL) naming the likely causes and the name-list keys
  to change, and the settings of the run.  Input errors of probCFD.nml
  and the csv (nml_err, csv_err) write the same file with kind "input
  error"; the normal end writes it with "finished normally".

25.5 Log layout  [sixdof.F90 sixdof_monitor, Main.F90]
  Per body: force (E, SNAME) [N], moment about O [N m] (+ moving /
  static), surge/sway/heave [m], roll/pitch/yaw [deg], u v w [m/s],
  p q r [deg/s], one per line with the labels right-aligned to the
  colon, an empty line after each body; then "overset re-assemblies",
  "cells active / fringe / hole", "orphans / fringe entries" flush left.
  21 empty lines between step groups (was 7).

25.6 Volume files  [output.F90, Main.F90]
  Volume_solutions_grids/ (created by rank 0 with mkdir -p on first use)
  holds <case>_NNNNNNNN.grd (grid; the initial grid is _00000000, moving
  cases write one per export), <case>_NNNNNNNN.p3d (solution) and
  <case>_NNNNNNNN_cblank.p3d (exact cell status, overset).  The grid
  check writes <case>_00000000.grd, _00000000_iblank.p3d,
  _00000000_cblank.p3d there.  The force / moment / motion / monitor /
  wall-function data files stay in the working folder.
  tools/blank_slice.py reads the new layout first and falls back to the
  old <case>.x / _iblank.q / _cblank.q files.

  Every probCFD.nml of the case folders got "Time_scheme = 1" after
  Iner_Itn (the key is optional: default 1).

25.7 Tecplot variable names (update 13b, 2026-10-01)  [output.F90]
  A PLOT3D function file carries no names; Tecplot shows F1V1 ... unless
  a PLOT3D "name file" is loaded with it.  write_solution now writes the
  variables in the order u, v, w, TTL_Prs, Dym_Prs, VoF, Blank (Dym_Prs /
  VoF only for two-phase runs, Blank only for overset) and the name file
  Volume_solutions_grids/<case>.nam (one name per line) at the first
  export.  tools/blank_slice.py keeps using the last variable as the
  blank.

25.8 Restart (update 13b)  [restart.F90, Main.F90, sixdof.F90, read_input.F90]
  Name list: Mode = 0 / 1 (first line of $Case_title; 1 = restart from
  the latest <case>_NNNNNNNN.rsto), Restart_Intvl (Export_Interval,
  default 0 = off, runtime-changeable).  write_restart (every
  Restart_Intvl steps and at the last step) gathers on rank 0 and writes
  one binary file: tag, case name, step, time, nblk, ndom, phases,
  turbulence flag, number of bodies, Del_T, the monitor sums (initial
  liquid volume, clipped volume), the 42-value state of every body (eta,
  etad, etadd, nu, nud, fp, ff, mp, mf) and, per sub-domain, U1 U2 U3,
  u v w, p, [alpha], [tk tw mut] with their ghost bounds.  The file is
  listed in <case>_restart.log.  read_restart (Mode = 1, after the
  initial assembly and sixdof_init) checks the header against the run,
  scatters the arrays to their owners (bounds checked), sets the body
  states and calls sixdof_restore_pose: update_transform, move_grids and
  after_motion with the module flag "restoring" so that neither
  rotate_fluxes nor fill_uncovered touches the restored fields (they
  were saved on the moved grid); the overset is re-assembled, halos /
  fringes refreshed, the star field set.  The monitor and the 6-DoF data
  files are opened in append mode, the restored grid pose and field are
  exported, the marching loop starts at istep0 + 1 with the restored
  t_now.  Tests: Cylinder_Overset_2phase 10 steps with Restart_Intvl 5,
  then Mode 1 from step 5: steps 6-10 identical to the uninterrupted run
  in every monitor digit (|div|min at 1e-16 differs); Box_6DoF 3 steps
  with Restart_Intvl 1, restart from step 1: the pressure forces of steps
  2-3 agree to 1e-8 relative (the re-assembled interpolation and the
  Krylov path differ at round-off).  Not restored: the SSP stage copies
  (rebuilt every step), un / alpn (set every step), wall distances
  (recomputed at start-up).  A restart with a different Del_T is allowed
  (noted in the log); a different grid, csv splitting, phase count,
  turbulence flag or body count is refused.

## 26. Update 14 — child and sub-child bodies of the 6-DoF motion (2026-10-01)

Request (6DoF Motion Subroutine: Child Body and Sub-Child Body
Implementation): a three-level hierarchy body > child body > sub-child
body (e.g. hull > stern drive > propeller), children with the same
properties as a body plus the mother ids, inertia in their own axes,
imposed-motion files with the mother ids in the name, name-list order
[body1] [child1] [sub1] ... [child2] ... [body2], csv rows
Child_Object_ID / Sub_Child_Object_ID, outputs relative to the mother in
a separate folder Output_data, motion coupling (imposed or fluid-driven)
with the mother's 6-DoF motion, verification and a GPPH case.

26.1 Data model  [read_input.F90 read_split_csv, mod_global.F90, sixdof.F90]
  Every rigid part is a NODE bd(k) of the sixdof module; blk_obj(b) is
  the node id of block b (0 background) and is what the overset uses as
  the rigid motion group (mgroup) - a child is its own group, so the
  donor / self-cut caches stay valid for every part.  read_split_csv
  reads rows 6-8 (Object_ID, Child_Object_ID, Sub_Child_Object_ID; 7-8
  optional, any label accepted), checks them (a child id needs a body id,
  a sub-child id a child id, ids 1..n without gaps per parent) and builds
  the node table (global n6node, node_bid/cid/sid, node_lev, node_par):
  the bodies first (node = body id), then the children sorted by (body,
  child), then the sub-children sorted by (body, child, sub) - so the
  parent of a node always has a smaller index and loops over k = 1..n6node
  update the mothers first.  Name list: sixdof_nml keeps the entries in
  ent(:) keyed by (level, body, child, sub); "Body = k" starts a body,
  "Child_body = c" a child of the current body, "Sub_child_body = s" a
  sub-child of the current child; "Mother_body" / "Mother_child_body"
  are required for children / sub-children and must equal the enclosing
  ids (error otherwise); the property keys apply to the current entry.
  sixdof_init maps the csv nodes to the entries (errors: a node without
  name-list data, an entry without blocks).  body_t gained lev, bid, cid,
  sid, par, own_moving (any own DOF), moving (own or a moving mother),
  r0 (initial offset O_k - O_m in the mother frame), the absolute pose /
  velocity / acceleration (R, pE, vE, OmE, aE, alB), the grid-frame
  velocities vG, omG used by move_grids and sixdof_wall_velocity, the
  composite properties (Mc, rcgc, IOc, MRBc) and the subtree wrench
  (Fs, Ms).  MAXB = 64 parts in total.

26.2 Kinematics  [sixdof.F90 update_transform]
  A body keeps the absolute eta as before.  A child / sub-child keeps the
  RELATIVE state: eta(1:3) = displacement of O_k from its initial
  position measured in the mother frame B_m, eta(4:6) = Euler angles of
  B_k relative to B_m (R_k = R_m R_rel); nu(1:3) = R_rel^T etad(1:3),
  nu(4:6) = omega_rel = J2^-1 etad(4:6).  Composition (mothers first):
    R_k = R_m R_rel,  p_k = p_m + R_m r,  r = r0 + eta(1:3)
    v_k = v_m + R_m (omega_m x r + etad(1:3))
    Omega_k = Omega_m + R_k omega_rel
    a_k = a_m + R_m (alpha_m x r + omega_m x (omega_m x r)
                     + 2 omega_m x etad(1:3) + etadd(1:3))
    alpha_k (B_k) = R_rel^T alpha_m - omega_rel x (R_rel^T omega_m) + nud(4:6)
  then Rg = F R_k F, dg = F p_k - O0, vG = F v_k, omG = F Omega_k; the
  grid moves as x = Rg (x0 - O0) + O0 + dg with the grid velocity
  vG + omG x (x - O0 - dg) (move_grids, sixdof_wall_velocity).  A fixed
  child (DOF 0) has R_rel = I, eta = 0 and moves rigidly with its mother.

26.3 Dynamics  [sixdof.F90 composite_properties, transfer_wrenches, body_step]
  Before the step (at the poses of the previous step): the composite
  mass, first moment and inertia about O_k in B_k of every part with all
  its descendants (parallel axes, R_j^k Icg_j R_j^kT), the Fossen mass
  matrix MRBc of the composite, and the subtree wrench bottom-up: Fs, Ms
  (E, about O_k) = own wall forces + the children's Fs, Ms transferred to
  O_k (lever arm p_j - p_k) WITHOUT the components along the child's free
  (released) DOFs: translations along the mother-frame axes e_i = R_m(:,i),
  rotations about the Euler-angle axes R_j J2^-1(:,i).  Fixed and imposed
  DOFs transmit everything (an imposed propeller = motor-driven shaft:
  the torque is reacted by the mother; a free one spins and passes only
  the thrust).  The step of part k (rank 0, then broadcast) is the
  Newton-Euler system about O_k in B_k written for the RELATIVE
  acceleration x = [R_rel^T a_rel ; omegadot_rel]:
    MRBc x = b,
    b(1:3) = F - M (a_fr + alpha_fr x r_G + Omega x (Omega x r_G))
    b(4:6) = tau - Omega x (IOc Omega) - IOc alpha_fr - M r_G x a_fr
  with F, tau = R_k^T (Fs, Ms) + the composite weight at r_G - the
  artificial damping, Omega = R_rel^T omega_m + omega_rel (absolute, in
  B_k), a_fr = R_rel^T (R_m^T a_m + alpha_m x r + omega_m x (omega_m x r)
  + 2 omega_m x v_rel) and alpha_fr = R_rel^T alpha_m - omega_rel x
  (R_rel^T omega_m) the frame terms of the mother's motion of the CURRENT
  step (zero for a body).  In the eta representation x = Ji etadd - g,
  g = [0 ; J2^-1 J2dot omega_rel], Ji = diag(R_rel^T, J2^-1):
  Meta = Ji^T MRBc Ji, bvec = Ji^T (b + MRBc g) (+ the added-mass
  relaxation as before), the free rows solved with the fixed / imposed
  accelerations known, velocity Verlet on eta.  For a body (no frame
  terms) this is algebraically the former body_step (checked: the Box
  single-body run of the new code equals the 13b binary to 1e-7, the
  level of the Krylov convergence noise).  Approximation stated in the
  header: the inertial reaction of the RELATIVE motion of a child on its
  mother (I_child alpha_rel, gyroscopic terms) is neglected - valid for
  appendages light compared with the mother; the composite inertia of a
  free rotational child DOF includes the child's inertia about that axis.
  Order per step: composite_properties, transfer_wrenches, then the nodes
  k = 1..n (mothers first): own step if own_moving, Rg_prev, update_
  transform (composition with the mother's NEW pose); move_grids,
  after_motion as before.

26.4 Files, log, restart  [sixdof.F90 sixdof_init, sixdof_export, sixdof_monitor]
  All 6-DoF data files go to Output_data/ (created by rank 0): bodies
  <case>_{force,moment}_{body,earth}<k>.dat and <case>_motion<k>.dat
  (unchanged content, own walls); children <case>_body<k>_child<c>_
  {force,moment,motion}.dat and sub-children <case>_body<k>_child<c>_
  sub<s>_{force,moment,motion}.dat: force / moment of the own walls in
  the MOTHER's frame (R_m^T), moment about the own O; motion file:
  relative displacement in the mother frame, relative Euler angles,
  relative velocities / rates / accelerations, O and CG in the grid
  frame, and (columns 33-38) the absolute Euler angles and angular
  velocity.  Imposed files: <case>_body<k>_child<c>_imposed.dat,
  <case>_body<k>_child<c>_sub<s>_imposed.dat (or Motion_File).  The
  screen monitor indents the children under their body and prints their
  relative motion; the body line adds "with the children, force" when a
  child exists.  Restart: the 42-value state per node (relative for the
  children) as before; sixdof_restore_pose composes mothers first.

26.5 Verification and validation
  (a) Box_6DoF, 3 free steps (Motion_Start_T = 0), body 1 split into the
      body (blocks 1-3, 60 kg, 0.6 I) and a FIXED child (block 4, 40 kg,
      0.4 I, CG at the body's CG, own rotation centre offset by (0, 0.5,
      0.2) m): the composite reproduces mass 100 kg and I; the motion of
      step 1 is identical to the single-body run (8 digits), the child's
      pose / velocity / angular velocity composed from the mother agree
      with the rigid motion (grid velocity sums of block 4 equal to
      1e-14).  The forces of step 1 differ by 0.04 N of 1000 N (4e-5):
      traced to round-off - with the cell donors, fringe entries and
      weights of the re-assembly identical to 1e-13 and the solvers
      converged to 1e-12, the difference stays; with the child's rotation
      centre equal to the mother's (bit-identical grid) the runs agree to
      2e-9; a 1e-10 m shift of the single body's rotation centre changes
      the step-1 force by the same 0.025 N.  The Box case is sensitive
      at round-off level (non-smooth switches: VoF limiter / clipping);
      the donor caches of the dynamic overset show the same sensitivity
      (cache vs. fresh assembly: 0.24 N) - noted for a later look, not a
      6-DoF matter.
  (b) Composed imposed motion: body 1 static, child (block 3) imposed
      surge 5 m/s + yaw 2000 deg/s, sub-child (block 4) imposed sway +
      roll 3000 deg/s, all rotation centres at (3, 7, 0), versus a single
      body with the analytically composed table (phi, 0, psi; p = p_c +
      Rz(psi) r_rel): the block-4 grid of step 1 agrees to 5e-7 m (the
      single-precision grid file; displacement 0.11 m) and to the
      Rodrigues formula; the absolute angular velocity (3000, 104.67,
      1997.26 deg/s) and the absolute Euler angles are identical.
  (c) Restart of the body + child case from step 2: step 3 identical
      (motion to all digits, forces to 2e-9).
  (d) GPPH (26.6): the start-up on one sandbox rank - csv rows 6-8 and
      the name-list hierarchy mapped (nodes 1 hull / 2 drive / 3
      propeller, 3 + 19 + 9 blocks), composite mass 6.71 kg, imposed
      table of 8401 rows read, initial assembly, step-0 forces exported
      for the three parts in Output_data/ (hull pressure force 59.9 N
      upward = the buoyancy of the 6.21 kg displacement).  The first
      moving step (re-assembly with the propeller turned 3.6 deg) exceeded
      the sandbox memory limit (6 GB RSS on one rank; 2 ranks x 3.3 GB
      before), so the time marching of GPPH has to be done on the user's
      machine; the moving hierarchy itself is validated by (a)-(c).

26.6 GPPH case folder  [cases/GPPH]
  Rows 6-8 of the csv: body 1 = hull blocks 1-2 + refinement 31; child 1
  = stern-drive blocks 3-20 + refinement 30; sub-child 1 = propeller +
  hub blocks 21-29; background 32.  Hull data from the grid: displaced
  volume below z = 0 of the hull surface 6.21e-3 m3 (mass 6.21 kg), LCB
  x = 2.145 m (transom at x = 2.439), CG 2 cm above the water line, radii
  of gyration 0.35 B / 0.25 L; free heave + pitch from Motion_Start_T =
  0.5 s with an added-mass relaxation of the order of the mass.  Stern
  drive rigid with the hull, rotation centre on the vertical steering
  axis at the transom (the annotation says how to impose a steering
  yaw).  Propeller: hub centre (2.5186, 0, -0.0989) from the block-21
  hub surface, shaft axis (-0.99921, 0, 0.03964) in the grid frame
  (2.27 deg nose-up forward, fitted to the hub rings), 10 rev/s right-
  handed (J = 1.1); the table GPPH_body1_child1_sub1_imposed.dat (0.25 ms
  rows, 2.1 s) was written by tools/make_axis_rotation.py, which converts
  a constant-speed rotation about an arbitrary axis (Rodrigues) into the
  unwrapped Z-Y-X Euler angles the solver interpolates (checked: the
  Euler matrix equals the axis rotation to 1e-10).  Del_T = 0.001 s is
  3.6 deg per step (0.0005 s recommended for the blade passage).  Every
  step re-assembles the overset (400-570 s on 1-2 sandbox cores for this
  0.94 M-cell grid, and more memory than the sandbox has; the 39-CPU
  splitting of the csv is meant for the cluster).  An input error of the
  6-DoF data (missing imposed file, missing name-list part, ...) now
  stops the run BEFORE the assembly (sixdof_precheck, called right after
  the csv), not 7 minutes later.

26.7 Update 14b (2026-10-02): Mode = 2 grid-movement check, missing-file stops
  [Main.F90, sixdof.F90 sixdof_check_movement / move_nodes_only,
   output.F90 write_grid_to / make_dir, stop_criteria.F90 stop_clean,
   read_input.F90]
  Request: a Mode = 2 that only checks the mesh movement of the impose
  motion files before the CFD run (PLOT3D grids at the time-step / export
  interval of the name list, in a separate folder Check_grid_movement),
  explicit stops "There is no name list file" / "There is no grid file" /
  "There is no impose motion", and the GPPH case with BOTH the stern drive
  (steering about the given axis to +45 deg, then fixed) and the propeller
  (continuous) rotating, delivered in Mode = 2.
  Implementation: Mode accepts 0..2.  Main calls sixdof_check_movement
  right after sixdof_precheck (i.e. after the name list, the b2b file, the
  grid and the csv; before the domains, BCs, metrics and the assembly).
  It keeps the original node coordinates of the moving blocks, loops
  istep = 0..nint(End_T/Del_T): imposed_at for the imposed parts,
  update_transform of every part (mothers first: a child of a steering
  drive turns with it), and at every Output_Intvl step (+ the last) moves
  the node coordinates (move_nodes_only: x = Rg (x0 - O0) + O0 + dg, no
  half-index / metric / grid-velocity work) and writes
  Check_grid_movement/<case>_NNNNNNNN.grd through write_grid_to(istep,
  dir) (the former write_grid_plot3d body; write_grid_plot3d(istep) now
  calls it with Volume_solutions_grids/), and a line of
  <case>_movement.dat (12 columns per imposed part: relative x y z, roll
  pitch yaw, absolute roll pitch yaw, O in the grid frame).  Free DOFs
  stay zero, fixed children ride with their mothers.  stop_clean(kind,
  msg) (all ranks): message to the log, Reason_for_stop.md (kind 5
  "missing input", kind 6 "Mode = 2 finished"), MPI_Finalize, stop.  The
  name-list check (inquire probCFD.nml) is at the top of read_namelist,
  the grid check (inquire <case>.grd) at its end; sixdof_precheck in Mode 2
  stops when no part has an imposed DOF, and read_imposed reports a
  missing table as "There is no impose motion file: <name>".
  GPPH: the steering axis given in the design frame of the hull (Point1 =
  (1.0295788, 0, 0.0003549465), Point2 = (1.0226202, 0, 0.0498683500)) and
  the propeller axis (horizontal at z = -0.086811873) were mapped into the
  grid: the grid hull is the design hull trimmed 2.272 deg bow-up (the keel
  slope and the hub-axis tilt both give 0.0396) and shifted; with the bow
  tip at design x = 0 and the propeller axis through the measured hub
  centre the shift is (1.4442, 0, 0.0306), so the steering axis passes
  through (2.4730, 0, -0.0099) with the direction (-0.0998, 0, 0.9950) and
  the propeller axis is the hub line (2.5186, 0, -0.0989) / (-0.99921, 0,
  0.03964) as before (the shift only moves the steering point along the
  nearly vertical axis, so the uncertainty of the shift hardly matters).
  Tables: drive 0.125 rev/s (45 deg/s) to +45 deg then held
  (--stop_deg), propeller 10 rev/s; Del_T = 0.0005 s, End_T = 1.2 s,
  Output_Intvl = 20 -> 121 grids.  Sandbox run (2 ranks, 4 min): the
  drive turns about the raked axis with the propeller end swinging to
  starboard, reaching 44.74 deg SNAME yaw (+ 4.05 roll, 1.67 pitch = the
  Euler decomposition of the 45-deg rotation about the raked axis) at
  t = 1.0 s and holding; the propeller keeps turning on the turned shaft
  (relative roll 4320 deg at 1.2 s); movement_check.png in
  cases/GPPH/sandbox_mode2_check.  Regression: Box_6DoF body + child
  1-step run bit-identical to before the change; the three stop messages
  verified (empty folder, name list only, Box with Mode = 2 and no
  imposed DOF).
  Addendum (Tecplot solution time): PLOT3D grid files carry no time, so
  Tecplot showed all 121 grids at once with the solution-time bar
  inactive.  Mode = 2 now writes beside every grid a PLOT3D solution file
  <case>_NNNNNNNN.q (write_check_q: 5 single-precision variables Part_ID,
  dX, dY, dZ, Disp = 6-DoF part of the block and the node displacement
  from the initial grid; block header FSMACH 0, ALPHA 0, RE = step, TIME =
  t; names in <case>.nam) - Tecplot's PLOT3D loader takes the solution
  time of a grid / solution pair from the TIME header, so the pairs loaded
  together animate.  Per export 11.6 + 24.5 MB for GPPH (4.5 GB for the
  121 exports of the delivered setting; raise Output_Intvl to thin it).

## 27. Update 15 — cavitation models (2026-10-02)

Specification agreed in the discussion: transport-equation models
(Schnerr-Sauer, Zwart-Gerber-Belamri, Kunz) selected by `Cav_model`, the
vapour as a third phase next to the water / air VoF, the model constants
in the name list with the literature values as defaults, the absolute
pressure level taken from the boundary conditions / `RefPrs` (no
`Prs_Atm`, no `Cav_Number`), fixed internal numerics (ramp, implicit
sink, air cut, bounding), `Vapor_scheme` (CICSAM explicit / implicit /
ENO with `Cnvt_mode`), `Reboud` among the turbulence variables, `Vapor`
in the solution files, cases NACA66mod and DelftTwist.

27.1 Formulation  [cavitation.F90, mod_global.F90, vof.F90, poisson.F90]
  Cell arrays vap (vapour fraction of the liquid-vapour part), vapn,
  cavsc / cavsp (explicit source, implicit sink coefficient, 1/tch),
  cavs (continuity source), s_vap (SSP); field id 19 for the b2b /
  overset exchange.  Mixture properties at every half-index location
  (rho_h, mu_h with the new vap_h, and rho_lam / mu_lam of the
  turbulence module):
     rho* = alpha [(1 - av) + av rho_v/rho_l] + (1 - alpha) rho_a/rho_l
  so the momentum, Poisson and turbulence routines see the vapour
  through the same functions as the air.  Vapour equation
  d(av)/dt + div(av u) = mdot/rho_v solved by advect_alpha of vof.F90,
  generalised to an advected field id (module pointer avf, helper sbv):
  the explicit CICSAM update becomes (av - dt flux + dt cavsc)/(1 + dt
  cavsp), the implicit CICSAM system gets dt cavsp on the diagonal and
  dt cavsc on the right-hand side; Vapor_scheme = 3 uses minmod-limited
  2nd-order upwind (ENO) face values (as weights beta in the implicit
  machinery), integrated explicitly (sub-cycled) for Cnvt_mode = 1,
  implicitly for 2, and for the hybrid 3 explicitly when the maximum
  Courant number is below cfl_switch, implicitly otherwise.  Final
  clipping to [0,1] with the clipped volume reported.  Continuity:
  div u = S with S = (mdot/rho_v)(1 - rho_v/rho_l) evaluated from the
  final vapour fraction of the inner iteration (cavs = (cavsc - cavsp av)
  (1 - rv)); compute_rhs_poisson uses (div u* - S)/dt, divergence_stats
  reports |div u - S|, outflow_correction targets qout = qin + sum S V.
  The sources are evaluated once per inner iteration from the latest
  pressure (cav_sources), before the vapour transport; the step order is
  alpha transport, vapour sources + transport, momentum, Poisson,
  correction, turbulence.  SSP (Time_scheme = 2) averages vap like
  alpha; the restart file (tag v2, cavitation flag in the header,
  v1 files still read) carries vap; fill_uncovered copies vap; the
  overset fringe cells receive vap (ovs_update_cells) - and, found on
  the way, also k / omega / mu_t (ids 14-16), which before kept their
  initial values in the fringe cells.
27.2 Models (rates mdot/rho_v in 1/s, written non-dimensionally with
  dp* = (p - p_v)/(rho_l U^2), rv = rho_v/rho_l; the (1 - av) factor of
  the evaporation and the av factor of the condensation are implicit):
  Schnerr-Sauer: R_B = [(av + a_nuc)/(1 - av) 3/(4 pi n0)]^(1/3) (>= R0),
  a_nuc = V/(1 + V), V = n0 4/3 pi R0^3, C = 3/(rho_m* R_B/L) sqrt(2|dp*|/3),
  evaporation C (av + a_nuc)(1 - av), condensation C av (1 - av) - the
  OpenFOAM form, with the nuclei fraction seeding the evaporation from
  av = 0.  Zwart: evaporation F_vap 3 r_nuc (1 - av)/(R_b/L)
  sqrt(2|dp*|/3), condensation F_cond 3 av/(R_b/L) sqrt(2 dp*/3).  Kunz
  (t_inf = tch): evaporation C_dest (1 - av) (-2 dp*), condensation
  C_prod av (1 - av)^2.  Common: ramp 0 -> 1 over 50 steps from
  Cav_Start_T (cav_ramp, printed), no mass transfer where 1 - alpha >
  0.5, only in active cells.  Reboud (turbulence.F90 rho_reboud): the
  density in mu_t = Re rho nu_t replaced by rho_v + (rho - rho_v)^10/
  (rho_l - rho_v)^9 when Reboud = 1.
27.3 Pressure level [cavitation.F90 init_cavitation, output.F90]
  p_level = RefPrs (SW_Ref_Prs = 1), else the first Dirichlet pressure of
  the .bcs file, else p_cel_initial (with a warning that the level is
  not fixed); sigma = (p_level - p_v)/(1/2 rho U^2) printed, warnings
  when sigma <= 0.  Dym_Prs of two-phase files subtracts p_level.  The
  cases set the Dirichlet outlet pressure AND p_cel_initial to
  p_v + sigma q (NACA66mod: 2794.2 Pa for sigma 0.91 at U = 1 m/s,
  rho 998.2; DelftTwist: 2875 Pa for sigma 1.07).
27.4 Name list / output
  $Cavitation: SW_Cav, Cav_model, Vapor_scheme, Cav_Start_T (runtime-
  changeable), Prs_Vapor, Lou_v, mu_v, SS_Nuclei, SS_Radius, ZGB_Fvap,
  ZGB_Fcond, ZGB_Rnuc, ZGB_Rb, Kunz_Cdest, Kunz_Cprod; Reboud with the
  turbulence variables.  Start-up summary lines "Cavitation : ON, <model>
  (constants)", "vapour transport", "mass transfer from t = ...", and
  the cavitation block (pressure level, p_v, q, sigma).  Monitor line:
  Vvap[m3] and avmax appended; log block "[ Cell ] Center Monitor :
  < Cavitation: model >" with the vapour volume, max fraction,
  evaporating / condensing cell counts, max |div u| source, ramp,
  iterations / residual of the implicit transport, sub-cycles, clipped
  volume.  .p3d variables u, v, w, TTL_Prs, [Dym_Prs, VoF], Vapor, [Blank].
27.5 Stability of the pressure / mass-transfer coupling (what the tests taught)
  First version: the mass transfer evaluated from the pressure of the
  previous iterate, transported, and its volume source put explicitly
  into the Poisson right-hand side.  On the cylinder (impulsive start,
  sigma 0.7) and on NACA66mod (sigma 0.3) the cavity FLASHED: in the
  onset step 10-15 % of the cells (every cell with p < p_v) evaporated at
  the full rate (Schnerr-Sauer: dav/dt of order 100 1/s at 500 Pa of
  tension - the 1-micron nuclei grow to cell size in a millisecond, far
  inside one 4-5 ms time step), the incompressible response to a volume
  source switched on over a large region in one step is a pressure rise
  of order rho S L^2/dt (thousands of Pa), the next step condensed
  everything, and so on (vapour volume 1.4 -> 2e-4 -> 3.1 m3 ... on the
  cylinder; the NACA run diverged at the third cavitating step).  Three
  measures, all internal, cured it:
  (1) the volume source is PRESSURE-IMPLICIT in the Poisson equation:
      S(p) = S0 + S'(p - p0), S' <= 0, with S' the secant to p_v
      (S0/(p0 - p_v): the source vanishes when p reaches p_v) or the
      tangent of the rate law when steeper; -J S'/dt is added to the
      operator diagonal (cell_div_G, apply_M7 / apply_Mfull, the Jacobi
      and line-SGS diagonals) and S' p0 to the right-hand side.  This
      is the "compressibility" every cavitating projection solver
      carries; it pins the cavity pressure at p_v instead of letting it
      overshoot;
  (2) consistency: after the pressure solve the vapour fraction receives
      the increment the implicit part produced, dt S'(p - p0)/(1 - rv),
      and the continuity source is updated (cav_pressure_update) - the
      vapour created and the dilatation of the velocity field are the
      same quantity;
  (3) the explicit production is capped at a fraction change of 0.2 per
      step (CAPDAV): without the cap the Zwart model (F_vap = 50,
      R_b = 1 micron: 6e4 1/s) filled cells to av = 1 in one step before
      the pressure could react and the predictor blew up in the
      near-zero-density cells (CFL 276).  With the cap the implicit part
      supplies the rest and the pressure settles at p_v.
  With these the cylinder runs (12 steps, Iner_Itn 2, verification/
  update15_cavitation) grow the cavity smoothly for all three models
  and all three vapour schemes, the SSP scheme and a restart from step
  12; abs(div u - S) stays at 1e-3 1/s (SSP: the averaged field's target
  is the averaged source, s_cavs).  Recommended with cavitation:
  Iner_Itn = 2 (the sources are re-evaluated from the new pressure,
  a Picard iteration of the coupling) and Del_T small enough for the
  cavity dynamics (the cavity pressure is held at p_v whatever dt, but
  the shedding needs resolving).  NACA66mod, sigma 0.3 (sandbox: 100
  steps of 4 ms, 2 cores, 45 s per step; the delivered case is sigma
  0.91 for the user's 5-s run since the suction peak of the slow
  towing-tank start is still Cp -0.4 at 0.4 chord passages): the vapour
  appears at step 26 with a smooth pressure response (Poisson update
  norm 5e-3 Pa, before: 10 Pa and divergence), a thin sheet (a few
  cells, av up to 0.44) along the suction side from x/c 0.15 to 0.6-0.75
  where Cp < -sigma, the pressure minimum pinned at p_v (2302-2338 Pa
  vs 2340), cases/NACA66mod_cav/sandbox_sigma0.3/naca_cav.png.
  DelftTwist, sigma 1.07, SST + Reboud (sandbox: 3 steps on one core,
  20 min each): start-up, pressure level 2875 Pa / sigma 1.07 printed,
  vapour forming around the foil in the impulsive-start suction
  (151 564 evaporating cells at the onset, 5e-4 m3 after 3 steps), log
  in cases/DelftTwist_cav/sandbox_startup.  Found on the way: the
  turbulence wall-distance set-up segfaulted because init_cavitation
  ran after init_turbulence while rho_lam already read the vapour
  array - init_cavitation now follows allocate_fields.  Regressions:
  Cylinder_Overset_2phase 10 steps identical to the 12b reference
  (|div|min at 1e-17 differs), NACA66mod_1phase 10 steps unchanged.

## 28. Update 16 — heat / energy equation (2026-10-03)

Specification agreed in the discussion: `$Energy` with `SW_Energy`,
`T_init`, the specific heats and conductivities of liquid / air /
vapour (read only when the phase exists), `Beta_l` / `Beta_a`
(0 = no buoyancy), optional `Pr_t`; no `T_ref` (= `T_init`), no
`Energy_Start_T`, no `Buoyancy` switch; `Latent_Heat` and `Pv_of_T`
moved to `$Cavitation`; separate solver controls `Max_Itr_Engy`,
`Min_Itr_Engy`, `Rsd_Lmt_Engy`, `SOR_fct_Engy`, `Cnvt_scheme_Engy`,
`Cnvt_mode_Engy`, `Dffs_mode_Engy`, `cfl_switch_Engy`, `vn_switch_Engy`;
`.bcs` T line with flags 0 / 1 / 3 / 4.

28.1 Formulation  [energy.F90, mod_global.F90]
  T [K] dimensional at the cell centres (temp, tempn, s_temp; field id
  20 exchanged like alpha; `temp_bc_hook` in refresh_all(20) completes
  the physical-boundary ghosts after the generic zero-gradient
  apply_bc).  Non-dimensional equation (x/Lch, t/tch, u/Uch)
     dT/dt + u.grad T = div(Gam grad T) + S,
     Gam = k_eff/((rho c_p) Uch Lch) = (k_eff/k_l)/(Pe_l rcp),
     Pe_l = rho_l c_p,l Uch Lch/k_l,  rcp = (rho c_p)/(rho_l c_p,l)
  with the mixture k and rho c_p volume-weighted (rcp_cell, kr_cell),
  k_t/k_l = (mu_t/mu_l)(Pr_l/Pr_t)(c_p,mix/c_p,l) (kef_cell), S = -mdot
  L/(rho c_p) from the continuity source cavs of the cavitation model
  (lat_nd = rho_v L/(rho_l c_p,l)).  Viscous dissipation and pressure
  work neglected.
  Discretisation: contravariant relative face fluxes (ALE), convection
  in the NON-CONSERVATIVE form sum_f F_f (T_f - T_c) (first-order upwind
  implicit or explicit, the higher-order face value of Cnvt_scheme_Engy
  as a deferred correction) - the update is a convex combination of the
  neighbours whatever the discrete divergence (cavitation dilatation,
  overset fringes), so T stays bounded; diffusion: orthogonal part with
  the harmonic-mean face conductivity divided by the cell's rho c_p,
  implicit / explicit per direction, cross-derivative part explicit
  (face average of the central differences); masks per cell and
  direction like the momentum hybrid (cell CFL dt |F|/V, VN 2 dt Gam
  g^qq); explicit modes sub-cycled to the switches (sub-cycle count
  global).  SOR sweeps over all owned domains between two exchanges
  with one global residual per sweep (relative to the first sweep).
  Boundary faces are kept in per-domain records (ebc: cell, family,
  side, patch; map(6,i,j,k)): Dirichlet -> implicit 2 Gam_f (T_w - T_c)
  on the half-cell distance (cfb = Gam A/dn), an inflow face carries
  T_w; flux -> explicit source q_w A/((rho c_p)_c Uch); convective ->
  implicit U (T_inf - T_c), U = 1/(1/h + dn Lch/k_f); zero gradient ->
  nothing (also for the convective flux, which cancels in the
  non-conservative form).  Thermal wall function on flag-11 walls:
  q_w = rho c_p u_tau (T_w - T_c)/T+, Kader's T+(y+, Pr) with the log
  slope Pr_t/kappa (cfb = u_tau A/T+); u_tau from the wall factor of
  the momentum wall function (turbulence: wall_utau).  The ghost values
  reproduce the boundary flux with the interior stencil (energy_bc).
  Buoyancy: buoyancy_force(c, pt, H) = -beta_f (T_f - T_init) g_c added
  to the momentum right-hand side at every face (beta_f volume-weighted
  by alpha at the face), active when Beta_l / Beta_a /= 0 and g /= 0.
  Output: 'Temperature' after 'Vapor' in the .p3d files (8-cell node
  average, overset: non-hole cells), monitor block < Energy >, monitor
  columns itT res(T) Tmin Tmax, Output_data/<case>_heat.dat (step, t,
  Tmin, Tmax, Tavg, then Q [W] into the fluid and the mean T_w of every
  Dirichlet / flux / convective surface, header with the surface list).
  Restart tag probCFD_rsto_v3 (energy flag, T as the last array; v1 /
  v2 readable, missing T -> T_init).  check_nan includes T.

28.2 Cavitation coupling  [cavitation.F90]
  pv_cell(i,j,k) = Prs_Vapor or, with Pv_of_T = 1 and the energy
  equation, pv_sat(T) (Buck 1996: 611.21 exp((18.678 - t/234.5)
  t/(257.14 + t)) Pa, t in C; 2339 Pa at 20 C) in the driving pressure
  of all three models; Latent_Heat in the energy source.  Cylinder
  sigma = 0.7 (verification case of update 15) with the energy
  equation, latent heat and p_v(T): the cavity cools by 0.13 K and the
  vapour volume after 12 steps is 2.78e-2 m3 instead of 2.82e-2 m3
  (thermodynamic delay).

28.3 Incremental pressure form  [poisson.F90 add_pgrad, Main.F90]
  Prs_form = 2: compute_rhs_momentum is followed by add_pgrad(1), which
  subtracts (1/rho) grad p^n at every solved face (the corrector's own
  operators: G^q fluxes / balanced chain-rule differences Q_l), and
  solve_momentum by add_pgrad(2), which adds dt (1/rho) grad p^n back to
  us/vs/ws and Us1..3 - the unchanged Poisson equation and corrector
  (full pressure p^(n+1), all BCs, pinning, fringe interpolation, the
  cavitation terms) then realise u = u* - dt grad(p^(n+1) - p^n)/rho.
  Fixed-flux faces untouched.  Single-phase runs with gravity start from
  the hydrostatic pressure.  Finding (28.6): with form 1 the predictor
  diffuses the body-force velocity dt g at no-slip walls wherever the
  wall-normal diffusion is implicit, and the projection cannot remove
  the resulting rotational part; a cavity of air at rest on the
  wall-clustered grid reached CFL 22 in 5 steps (0.19 m/s) - exactly
  0 with form 2 (and with form 1 on the uniform grid, where the wall
  diffusion is explicit).  Form 1 remains the default (every earlier
  result unchanged); form 2 is recommended for gravity / buoyancy
  driven flows.

28.4 Name list / .bcs  [read_input.F90, overset.F90]
  $Energy keys as above (+ runtime: the _Engy controls and Pr_t),
  Prs_form in $Numerical_Scheme, Latent_Heat / Pv_of_T in $Cavitation;
  bc_flag / bc_val now 6 rows, bc_tinf for the second value of flag 4;
  the reader takes up to two extra "flag value [value]" lines after p
  (two-phase: alpha and T both required with SW_Energy = 1; single
  phase: T, or alpha (ignored) + T); check_flags validates the T flags
  and values; the overset auto-patches carry the 6 rows and bc_tinf.
  Summary prints Pe, Pr of every phase, buoyancy, solver settings.

28.5 Verification  [verification/update16_energy, tools/energy_cases.py]
  (a) 1-D transient conduction, rho = c_p = k = 1 (alpha = 1 m2/s),
      T(0) = 310, T(1) = 290, T_init 300, 40 cells, dt 0.001 (VN 1.6,
      implicit): max |T - T_series| = 0.036 K at t = 0.05, 0.011 K at
      t = 0.10 (first-order time error; SSP Time_scheme = 2 with dt
      0.004: 0.09 K).  Restart at step 50 reproduces the continuous run
      exactly; a restart from a file without T starts from T_init.
  (b) 1-D steady slab, flux 100 W/m2 at x = 0, convective h = 10,
      T_inf = 300 at x = 1 (exact: T = 310 + 100 (1 - x)): T(0) =
      409.97, T(1) = 309.996, Q1 = 1.0000 W, Q2 = -0.9996 W, T_w2 =
      310.00 (heat file) after 4 s.
  (c) step advection u = 1 m/s, inlet 310 into 300 K, 100 cells, CFL
      0.2: scheme 1 bounded (front smeared over 10 cells), 2 undershoot
      -5.6 K, 3 explicit overshoot (+12 K, unstable like the momentum
      QUICK; implicit mode: +0.4 K), 4 (ENO) bounded and sharp (3
      cells).
  (d) two-phase overset cylinder (energy, cylinder wall 303 K, inlet 293
      K, 10 steps): bounded, heat flow 3e4 W; NACA66mod with cavitation
      + energy + SST on 3 domains / 2 ranks: runs (the deadlock fixed);
      flat plate SST + wall function with the plate at 310 K: runs,
      Nu / Nu_corr follows Cf / Cf_corr during the start-up.
  (e) natural-convection cavity: cases/NaturalConvection_Cavity README.

28.6 Regressions
  Cylinder_Overset_2phase 10 steps, NACA66mod single phase 3 steps,
  flat plate SST + wall function 20 steps: bitwise identical to update
  15 with SW_Energy = 0 (Prs_form = 1).  NACA66mod_cav (SW_Cav = 1):
  differences at the 1e-4 level from the vapour-ghost fix (the old
  binary read uninitialised memory at the outlet ghost cells).

28.7 Blank output (16b)
  ovs_node_iblank now copies the exact cell status cs(b)%a(min(i,nx),
  min(j,ny), min(k,nz)) to node (i,j,k) (1 active, 0 fringe, -1 hole;
  the max-face nodes repeat the last cell) - no max / solid rule, no
  averaging; the <case>_NNNNNNNN_cblank.p3d file (ovs_cell_blank,
  write_cell_blank) is removed.  The grid-check file
  <case>_00000000_iblank.p3d (Ovst_only = 1) keeps its 0.5 / 0.25 /
  0.125 fringe-layer shading, now built on the exact node status.

## 29. Update 17 — compressible mixture, thermal cavitation (2026-10-04)

29.1 Compressible liquid-vapour mixture  [cavitation.F90 cav_compress, poisson.F90, Main.F90]
  Cav_compressible = 1.  Wood / Wallis law of the liquid-vapour part of a
  cell: 1/(rho c^2) = av/(rho_v c_v^2) + (1 - av)/(rho_l c_l^2), the liquid
  with c_l = Snd_liquid, the vapour as an isothermal ideal gas
  rho_v c_v^2 = max(p, p_v) (c_v = sqrt(p_v/rho_v) = 368 m/s for water at
  20 C); the air of a two-phase cell stays incompressible (factor alpha).
  Non-dimensional coefficient beta_c = rho_l U^2 / (rho c^2); cavbc =
  beta_c/dt and presn = p at the start of the step (copied with vapn,
  also at the second SSP stage).  Pressure equation
     (A + (cavdp - cavbc)/dt) p = (div u* - cavs + cavdp p0 - cavbc p^n)/dt
  i.e. div u^(n+1) = S_cav + dS/dp (p - p0) - beta_c (p - p^n)/dt: the
  compressibility enters exactly like the pressure-implicit mass-transfer
  term of update 15 (cell_div_G, dcav of the preconditioners, the Jacobi
  diagonal, compute_rhs_poisson); the divergence monitor uses
  cav_div_target = cavs - cavbc (p - p^n).  All-Neumann problems: no
  pinning and no compatibility correction (the level follows the
  compression).  After the pressure solve the vapour fraction follows the
  isothermal gas law between max(p^n, p_v) and max(p, p_v):
  av' = av r/(1 - av + av r), r = p_old/p_new (bounded to [0.5, 2]); the
  conservative vapour transport in the dilating field carries the
  mixture-volume part, so the compression of the vapour is over-counted by
  av beta_c dp (1 % of av in the 30 % / 3 % test below, negligible in a
  cavity at p ~ p_v).  Monitor: max acoustic CFL (|u| + c) dt/h and the
  minimum mixture sound speed.  The pressure term is implicit, so the
  acoustic CFL does not limit the stability; dt < h/c_mix is needed only
  where the waves themselves matter.
  Verification (verification/update17_compressible_thermal, tools/
  energy_cases.py comp_box / comp_wave; water with a uniform vapour
  fraction 0.3 at 1 bar, mass transfer off, vap_cel_initial):
  (a) closed tube 1 m, piston 0.01 m/s: Wood law rho_m c^2 = 3.333e5 Pa,
      dp/dt = 3333 Pa/s at 1 bar rising to 3447 Pa/s at 1.034 bar (the
      gas stiffens); solver mean dp/dt 3364 (0-0.4 s), 3412 (0-1 s) Pa/s;
      the piston wave rho c u = 152.7 Pa bounces in the tube (p max
      matches); av 0.300 -> 0.296 (gas law 0.293, see above).
  (b) tube 2 m, 200 cells, dt 0.5 ms (acoustic CFL 1.09): the pressure
      step is 152.8 Pa (rho c u = 152.7) and its front travels at the
      Wood speed 21.82 m/s (0.430 / 0.870 / 1.310 / 1.740 m at 20 / 40 /
      60 / 80 ms vs 0.436 / 0.873 / 1.309 / 1.746).
  (c) cylinder sigma = 0.7 (update 15 verification case) with
      Cav_compressible = 1: see 29.4.

29.2 Thermal cavitation  [cavitation.F90, energy.F90]
  Cav_thermal = 1 (needs SW_Energy = 1; sets Pv_of_T = 1):
  - p_v(T) by Clausius-Clapeyron through the reference point (Prs_Vapor,
    Lou_v at T_init) with the constant Latent_Heat: M = Lou_v R T_init /
    Prs_Vapor, p_v = Prs_Vapor exp(L M/R (1/T_init - 1/T)), rho_v(T) =
    p_v M/(R T) (pv_T, rhov_T) - any fluid; Pv_of_T = 1 alone keeps the
    Buck law of water (update 16).  rv_cell(i,j,k) = rho_v(T)/rho_l is
    used in the rates (rhom of Schnerr-Sauer), in the continuity source
    factor (1 - rv) and in the vapour increment of cav_pressure_update;
    the mixture density of the momentum keeps the reference rho_v.
  - thermodynamic delay (B-factor, Plesset-Zwick): where p < p_v(T) the
    evaporation is driven by p_v(T - dT*),
       dT* = 1/2 (rho_v L/(rho_l c_p,l)) R_B / sqrt(alpha_l t)  (<= 0.1 T),
    the bubble of radius R_B (Schnerr-Sauer radius of the cell's vapour
    fraction) drawing its heat of evaporation through the thermal layer
    sqrt(alpha_l t); t = Cav_thermal_t, the bubble-growth / cavity
    residence time (default 1 ms; the modelling constant, runtime-
    changeable); monitor: max dT*.  Brennen's Sigma(T_init) and the
    scale rho_v L/(rho_l c_p) are printed at start-up (a first attempt
    with dT* = Sigma sqrt(t) gave hundreds of K for LN2: Sigma sqrt(t) is
    a velocity scale, not a temperature).
  - energy equation: sink -mdot L/(rho c_p) with rv_cell, and its
    derivative dS/dT = -dS/dp dp_v/dT (cavdp dpvdT_cell) implicit on the
    diagonal (stiff for cryogens).
  - vap_cel_initial: initial vapour fraction of the liquid cells
    (cav_init_field after initialize_fields).

29.3 Name list  [$Cavitation]
  Cav_compressible, Snd_liquid, Cav_thermal, Cav_thermal_t (runtime),
  vap_cel_initial; validation in read_namelist, summary lines.

29.4 Cylinder runs (Schnerr-Sauer, Iner_Itn 2)  [verification/update17_compressible_thermal]
  Compressible mixture from a restart at step 12 of the incompressible
  run: vapour volume 4.30e-2 vs 4.33e-2 m3 at step 24, within 1-5 % all
  along, min mixture sound speed 3.06 m/s.  From an impulsive start the
  far field dips below p_v in the transient (also incompressible: 19278
  evaporating cells at step 11) and the compressible fog expands -
  vapour volume 0.2 m3 at step 12: switch the compressibility on after
  the transient.  Thermal: LN2 at 20 m/s, sigma 0.7, 12 steps: cavity
  -9 % vapour, dT* up to 0.85 K, cavity cooled to 77.01 K; water at
  1 m/s: dT* 0.012 K, +0.7 % vapour (negligible, as expected).

29.5 Regression
  Cylinder sigma = 0.7 (SW_Energy = 0, new switches off): steps 1-11
  bitwise identical to update 16b.

## 30. Update 18 — incident waves (2026-10-04)

30.1 Design  [waves.F90, wave_fields.F90]
  Regular and irregular incident waves for the two-phase VoF solver,
  after the CFDShip-Iowa v5.5 treatment (exact potential-flow wave on
  the far-field boundary: type 17 there) with three changes agreed in
  the discussion:
  - the wave is defined ONCE in the earth-fixed frame, eta(X,Y,t),
    u_w(X,Y,Z,t), p_w, and the grid position is converted with one shift
    X = x + X_shift(t), X_shift = -(u_init t + (U_mean - u_init) R(t)) -
    wave_shift_extra, R = integral of bc_ramp (t^2/2t_ramp, then t -
    t_ramp/2), U_mean = the values of the u/v/w lines of the wave patches.
    A body held in a stream therefore sees the encounter frequency
    omega + k cos(theta) U automatically (CFDShip switches the
    frequency with earth_system instead), a time-varying stream is
    exact, and the later moving-background mode only adds its
    displacement to wave_shift_extra [m] (public hook, default 0).
  - the boundary condition is a far-field incident-wave condition,
    .bcs flag 12 on the u, v, w, p and alpha lines together (checked;
    two-phase only), usable on any face: velocity = ramped mean + wave
    orbital velocity at the face position (bv / mirror of apply_bc and
    dirichlet_cart treat 12 as a Dirichlet value), alpha ghost = fill
    fraction of the ghost cell below eta at the cell centre
    (wave_fill), pressure ghost = interior + (p_w(ghost) - p_w(interior))
    with the local density ratio (cellval; the Krylov vectors stay
    zero-gradient: the face is a fixed-flux face, type 1, so the
    projection does not touch the prescribed boundary flux - CFDShip
    needs the exact gradient in its Neumann row, here it only completes
    the ghost pressure for the cross terms and the output).  The
    turbulence inlet detection treats flag 12 as an inlet; the wall
    detections (flags 1 with zero values) are not triggered.
  - relaxation zones (Jacobsen et al. 2012, waves2Foam): after the time
    step (after the SSP average) phi <- (1 - w) phi + w phi_wave for
    alpha at the active cells and for the velocity at the active
    interior faces (transformed and Cartesian component with the same
    weight, then refresh_uvw / contra_to_cart ghosts), w = (exp(chi^3.5)
    - 1)/(e - 1), chi = 1 at the boundary and 0 at the inner edge; zones
    Relax_In_L / Relax_Out_L / Relax_Side_L measured from the x-min /
    x-max / y faces of the bounding box of the blocks that carry a
    flag-12 patch (so body blocks are excluded and the box moves with
    the background later).  The velocity weight is multiplied by the
    liquid fraction of the face: the orbital velocity frozen at the
    surface value in the air is not divergence-free, and blending it
    there put a divergence of 20/s into the zone cells (first run);
    with the water-only blend the divergence left after the relaxation
    is 0.1-0.2 (non-dim) and is removed by the next projection.  The
    first weight shape (mirror image, steep at the inner edge) made a
    kink at the zone edge that radiated waves; the Jacobsen shape is
    gentle over the inner half.
  Wave models (SI inside the module, non-dimensional values returned):
  regular wave (Wave_Amp / Wave_Height, Wave_Period, Wave_Phase,
  Wave_Heading from +x, Water_Depth with the finite-depth dispersion
  relation by Newton, 0 = deep; Wave_Order 2 = Stokes second order,
  CFDShip's finite-depth forms); irregular sea (Spectrum 1
  Bretschneider / ITTC two-parameter, 2 JONSWAP with Gamma_JONSWAP, 3
  Pierson-Moskowitz with the peak from Hs; N_Wave_Comp components with
  equal d(omega) over 0.5..3.5 omega_p, a_l = sqrt(2 S d(omega));
  N_Wave_Dir > 1: cos^2 spreading over +-90 deg; phases from an own
  linear congruential generator seeded by Wave_Seed - reproducible on
  every compiler; or a component file Wave_File "amp omega heading
  phase").  Above the surface z' = min(z, eta) (velocity frozen at the
  surface value, as CFDShip), below the bottom z' >= -d.  The deep-water
  pressure keeps CFDShip's second-order mean term a e^{kz}(cos - a k
  e^{kz}/2).  Wave_Ramp_T: half-cosine amplitude ramp.  Wave_Init = 1:
  alpha and the velocity (all face arrays, ghost layers included, then
  cart_to_contra) from the exact wave at t = 0, ramp off.
  Wave probes: Wave_Probe_N points; eta = z_bottom + column integral of
  alpha over the active owned cells whose x-y footprint contains the
  probe (MPI sum / min), with the exact incident elevation next to it,
  Output_data/<case>_wave.dat every step (append on restart); the
  component list to Output_data/<case>_components.wav; monitor line
  "[ Wave ] ramp, X_shift, probe elevations".  Runtime-changeable:
  Wave_Ramp_T and the three zone lengths.  Restart: nothing to add (the
  wave is a function of the restored t_now).  Compile order: waves.F90
  after b2b_exchange.F90 (boundary.F90 / curvi.F90 use it), the field
  routines in wave_fields.F90 after poisson.F90 (mom_ops, curvi).

30.2 Name list / bcs
  $Wave: SW_Wave, Wave_Amp (Wave_Height), Wave_Period, Wave_Phase,
  Wave_Heading, Water_Depth, Wave_Order, Spectrum, Hs, Tp,
  Gamma_JONSWAP, N_Wave_Comp, N_Wave_Dir, Wave_Seed, Wave_File,
  Wave_Ramp_T, Wave_Init, Relax_In_L, Relax_Out_L, Relax_Side_L,
  Wave_Probe_N, Wave_Probe_X(:), Wave_Probe_Y(:) (wave_nml hook of
  read_namelist, like sixdof_nml).  .bcs: flag 12 with the mean velocity
  as the value of the u/v/w lines (identical on all wave patches,
  checked), 0 on p and alpha; the T line is independent.

30.3 Validation (tools/wave_cases.py, verification/update18_waves)
  Wave tank 5 x 1.25 m (2-D), depth 1 m, lambda = 1 m, a = 0.01 m (ka =
  0.063), T = 0.8003 s, 32 cells per wavelength, 2 mm cells in |z| <=
  0.03 m, dt = 0.004 s (200 steps per period), geometric VoF, Prs_form
  2, zones 1 lambda (inlet) and 1.5 lambda (outlet), flag 12 on both
  ends, probes at x = 1.5, 2, 3 m; compare_wave.py fits a sinusoid at
  the (encounter) frequency to the probe histories and compares the
  snapshot eta(x) with the exact wave.  Results (details in
  verification/update18_waves/README.md): flat start - the front travels
  at the group velocity, amplitude / a = 1.025 / 1.024 / 1.022 and phase
  lag -16.5 / -21.2 / -29.2 deg at x = 1.5 / 2 / 3 m after 12 periods,
  crests 1.00 / 1.97 / 2.94 m (phase speed 2-3 % low, no reflection
  pattern); exact initial field kept over six periods (amplitude within
  2 %); convergence of the lag after two periods 20.4 / 13.4 / 9.6 deg on
  16 / 32 / 64 cells per wavelength, 13.1 deg with dt/2 on 32 cells
  (spatial, first order; steady lag at x = 1.5 m 16.5 -> 9.9 deg on the
  doubled grid); Time_scheme = 2 and Iner_Itn = 2 both damp the wave by
  12 % in two periods with the geometric VoF (SSP leaves flotsam) -
  Time_scheme 1 / Iner_Itn 1 for waves; mean current 0.5 m/s - encounter
  period 0.5714 s reproduced, amplitude grows +5..13 % along two
  wavelengths (+4..8 % with ENO; open item); JONSWAP Hs 2 cm, Tp 0.8 s,
  60 components, 8 s - stable, Hs ratio 0.88-0.94 over the short windows.

30.4 Regression
  Cylinder_Overset_2phase, 10 steps, SW_Wave = 0: monitor and .p3d
  bitwise identical to update 17.

## 31. Update 19 — force models: Hough-Ordway disk, BEM disk, point force (2026-10-06)

31.1 Design  [force_model.F90, sixdof.F90, convection.F90, read_input.F90, Main.F90]
  Force models are "bodies without a grid" in the 6-DoF hierarchy: a
  Force_model = n ... Force_model_End = -1 group of $SixDOF hangs on a
  Body / Child_body / Sub_child_body (Mother_body, Mother_child_body,
  Mother_sub_child_body; the node is found with sixdof_node_of) or on
  nothing (Mother_body = 0: fixed in space).  Body_End / Child_body_End /
  Sub_child_body_End = -1 close the parts (sixdof_nml resets its context;
  optional for old files).  read_namelist routes every line between
  Force_model and Force_model_End to force_nml (force_in_group), the rest
  to sixdof_nml as before.  Up to 64 force models (MAXF), no children.
  Geometry (Prop_P1, Prop_P2, Force_Point; grid coordinates [m] at t = 0)
  moves with the mother part through sixdof_node_pose (x = Rg (x0 - O0)
  + O0 + dg; point velocity vG + omG x (x - O)).
  - Force_type 1, Hough-Ordway: J = U_a/(n D) from the liquid-weighted
    disk-averaged axial inflow relative to the disk (pass 1 over the
    cells of all owned domains, MPI sums), T = K_T(J) rho n^2 D^4, Q =
    K_Q(J) rho n^2 D^5 with the quadratic curves of Prop_KT / Prop_KQ;
    f_x = A_x b(r*), f_t = A_t a(r*), b = r* sqrt(1 - r*), a = b/(r* (1 -
    r_h) + r_h), A_x = T/sum(b V), A_t = Q/sum(a r V) over the discrete
    disk cells (exact totals on the actual cells, CFDShip-Iowa's rescaling
    in one pass), times the liquid fraction of the cell.
  - Force_type 2, pure BEM: per disk cell the relative inflow (U_a axial,
    U_t along the rotation direction), the blade speed omega r, beta =
    atan2(U_a, omega r - U_t), the station data interpolated in r/R
    (Prop_rR, Prop_cD, Prop_PD, Prop_fc, Prop_CD0, Prop_CLa, Prop_CDk),
    phi = atan((P/D)/(pi r/R)), alpha = phi - beta + 2 f/c, C_L = C_La
    alpha (|C_L| <= 1.2), C_D = C_D0 + C_Dk C_L^2, q = 1/2 rho W^2 c, the
    Prandtl tip and hub loss F, dT/dr = Z F q (C_L cos beta - C_D sin
    beta), dF_t/dr = Z F q (C_L sin beta + C_D cos beta), f_x = (dT/dr)/
    (2 pi r dx), f_t = (dF_t/dr)/(2 pi r dx), times the liquid fraction;
    K_T, K_Q of the result reported.  No calibration (the user's choice).
  - Force_type 3, point force: Force_Fx..Fz and the optional couples
    Force_Mx..Mz are text formulas compiled at start-up into reverse
    Polish (shunting yard; numbers, the variables x y z roll pitch yaw =
    eta of the mother and u v w p q r = nu, t, constants pi g rho, the
    functions sin cos tan asin acos atan atan2 sqrt exp log abs min max
    sign step, + - * / ^, unary minus, parentheses; a trial evaluation
    checks the operand count) and evaluated every step.  The force is in
    the mother's body frame B_m, applied at Force_Point: wrench on the
    mother F_G = Fm R_m F_B, M_O = Fm R_m M_B + (x_pt - O) x F_G.  With
    no mother the formulas see x y z = Force_Point (the rest 0) and the
    force goes into the fluid at the cell containing the point (F/V).
  - Fluid coupling: the disk forces (and the fixed point force) are
    stored per owned domain as src(d)%f(3,cells) = f/(rho1 Uch^2/Lch);
    force_source(c, pt, H) gives the face value (average of the two
    cells)/rho* that compute_rhs_momentum adds next to gnd(c) and the
    buoyancy term.  The wrench on the mother (sum of -f V and -(x - O_m)
    x f V over the disk cells, i.e. thrust and torque together) goes to
    bd(k)%fext / mext (sixdof_set_ext) and transfer_wrenches adds it to
    the fluid wrench, so it climbs the hierarchy with the same free-joint
    rules.  force_compute runs after sixdof_forces every step (and once
    after the initial condition): explicit coupling like the fluid forces.
  - Output: Output_data/<case>_force_model<n>.dat (Data_Intvl) - F, M
    about O_m and M about the own centre in B_m (grid frame without a
    mother), for the disks n, J, K_T, K_Q, T, Q, liquid fraction, cells;
    "[ Force model n ]" monitor line.  Restart: no state (files appended).

31.2 Verification  [verification/update19_force_models, tools/force_cases.py]
  Actuator disk in a uniform stream (box 2 x 1 x 1 m, U = 1 m/s, D =
  0.25 m, n = 8 1/s, fixed in space, 96 disk cells, 600 steps): the
  Hough-Ordway disk converges to J = 0.667, K_T = 0.172, T = 43.0 N against
  43.7 N of momentum theory at that inflow (ratio 0.985; C_T = 1.75); the
  pure BEM (4 blades, P/D = 1, B4-70-like chords, camber 2 %, C_D0 0.01)
  gives K_T = 0.330 at J = 0.794 with C_La = 2 pi and 0.276 at J = 0.750
  with C_La = 3.8 (momentum-theory ratios 0.90 / 0.94 at C_T = 3.4 / 2.8;
  about twice the B4-70 series K_T: the face pitch + thin-airfoil polar
  over-predict the angle of attack, the user must set the section data
  against open-water K_T, K_Q, which the output reports).  Point-force
  spring on the NACA66mod_Overset foil (m = 1e4 kg, k = 3.95e7 N/m, period
  0.1 s, 50 steps per period, 100 steps at 20-45 s each): amplitude 9.97 mm
  (exact 10.0), period 0.1007 s (+0.7 %, the added mass of the foil), max
  deviation 0.41 mm over two periods.  The first attempt grew the
  amplitude by 5 % per half period: body_step advanced the velocity with
  a^n only (first order in energy); the velocity-Verlet half-kick 1/2 dt
  (a^n - a^(n-1)) is now applied when the next acceleration is known
  (vhalf flag per DOF) - moving-body results of earlier updates change at
  O(dt).  g in the formulas = |Grvty| (0 in the foil case).

31.3 Regression
  Cylinder_Overset_2phase (10 steps, no force model): bitwise identical to
  update 18 (the hooks are inactive without a Force_model group).

## 32. Update 20 — moving background (2026-10-06, validation 2026-10-11)

32.1 Design  [sixdof.F90 bgr_setup / bgr_follow, mod_global.F90, read_input.F90, waves.F90]
  $Case_title: Moving_BGR (0/1), Moving_BGR_DoF (six 0/1 flags, SNAME
  order), BGR_Following_Body (a Body, level 1), Moving_BGR_Filter_T [s].
  With Moving_BGR = 1 the blocks with Object_ID = 0 become the node
  kbg = nbd + 1 of the 6-DoF module (bgr_setup, after the node table):
  lev 0, no mass, no DOFs of its own, moving, rot0 = the followed body's
  Rot_Center; blk_obj(b) = kbg for those blocks, so they are their own
  rigid motion group of the overset (mgroup) with their own donor cache,
  move_grids / the ALE grid fluxes / Xh0 copies / the moved-grid output
  treat them like a body's blocks.  Every step after the body loop
  (sixdof_advance) bgr_follow sets eta_bg(c) = eta_body(c) for the
  permitted c, 0 otherwise (yaw only -> R = Rz(psi) about the body's
  rotation centre, which itself follows the permitted translations),
  low-pass filtered when Moving_BGR_Filter_T > 0 (f <- a f + (1 - a) raw,
  a = exp(-dt/T), etad = (f_new - f_old)/dt), then update_transform(kbg)
  (par = 0) gives Rg / dg / vG / omG; the grids move only when the pose
  changed (bgr_changed).  The pose is published to mod_global (bgr_on,
  bgr_Rg, bgr_O0, bgr_dg) for the modules compiled before sixdof.
  Conventions: Dirichlet velocities on the background blocks are
  EARTH-frame values - sixdof_wall_velocity returns 0 for kbg (the
  far field sees still water / the wave; a free-running ship has
  U_inf = 0); the background's tank walls carry no body
  (sixdof_forces skips kbg); the moved grid coordinates are the earth
  positions, so the wave module needs no shift (wave_shift_extra of
  update 18 removed) - only the relaxation zones, attached to the grid
  faces, take the position back to the original grid coordinates with
  the published pose (wave_weight), and the probes are earth-fixed.
  Output: Output_data/<case>_bgr_motion.dat (eta, etad, O of the
  background every Data_Intvl steps).  Restart: the background node is
  one more state record (sixdof_nbody = nbd + 1; eta = the filtered
  pose, so the filter continues), hence a restart needs the same
  Moving_BGR as the file.  Relaxation zones per face: Relax_L(1:6)
  (x-min x-max y-min y-max z-min z-max), the update-18 keys as aliases.
  Re-assembly tolerance ($Overset Ovst_Reassemble_Tol [m], default 0 =
  every moving step): after_motion skips ovst_reassemble (and the fill
  of uncovered cells) while every moving part's rotation centre and its
  three unit offsets have moved less than the tolerance RELATIVE to the
  background (to the fixed grid without a moving background) since the
  last assembly (reassembly_needed / remember_asm_pose); the summary
  line counts the re-used steps.
  Two fixes found by the validation: (a) sixdof_nbody returns nbd + 1
  with a moving background, so the restart file carries the background
  record (the filtered pose continues; the header check of restart.F90
  rejects a file written with another Moving_BGR); (b) the wall velocity
  of the force integration (sixdof_forces) and of the wall-function
  relative velocity (turbulence.F90, both places) now contains the ramped
  Dirichlet value u_init + bc_ramp (u_bc - u_init) exactly as
  boundary.F90 applies it - before, only bc_ramp u_bc was used, so during
  a towing-tank start (u_cel_initial = U, t_ramp > 0: the wall moves with
  the fluid at U (1 - ramp)) the reported friction and the wall-function
  u_tau treated the wall as fixed.  Only runs with t_ramp > 0 and
  u_cel_initial /= 0 are affected, and only during the ramp.  (c) A
  restart applied the Dirichlet boundary values of t = 0 (bc_ramp = 0:
  inlet = u_init, moving wall = u_init + body velocity) to the restored
  field, so the first step after a restart of a case with t_ramp > 0
  started from wrong level-n boundary fluxes (towed cylinder: the wall
  flux +0.5 instead of -0.5 m/s, a 1 % pressure blip and a one-step lag
  of the boundary layer); read_restart now sets bc_ramp = min(t/t_ramp, 1)
  and the wave time (wave_update(t)) before the restored pose is applied.
  Restarts of cases without a ramp were exact before (25.8) and remain so.
  (d) The grid velocity of the ALE fluxes (Ug1..3, move_grids) and the
  wall velocity of the no-slip patches (sixdof_wall_velocity) are now the
  DISCRETE velocity of the grid point, (x^(n+1) - x^n)/dt from the poses
  Rg_prev / dg_prev and Rg / dg (new function grid_velocity), instead of
  v_O + omega x r evaluated at t^(n+1): the relative flux of a face
  moving with the grid is then exactly zero and the convection is
  consistent with the displacement the time derivative sees.  Identical
  for a uniform translation or rotation; for an accelerating body the old
  value convected the near field with the spurious velocity a dt / 2 (a
  towed cylinder accelerating at 5 m/s2 with dt = 5 ms: 12.5 mm/s through
  the wall).  Moving-body results of earlier updates change at O(a dt).
  (e) sixdof_nml: Body_End / Child_body_End keep the current Body /
  Child_body as the mother of the groups written after them, so the
  order of probCFD_template.nml (Body ... Body_End, then Child_body with
  Mother_body, then Sub_child_body) is accepted - before, Child_body after
  Body_End stopped with "must follow the Body it belongs to"; property
  lines after an end mark are still refused.  The whole annotated
  template now parses (checked: the run stops only at the missing grid).
  (f) Wave probes (waves.F90 wave_export): the cell whose x-y footprint
  contains the probe is found with a point-in-cell test (inverse
  bilinear map of the four bottom nodes, half-open) instead of the
  bounding box of the cell, which on a yawed grid counted the probe in
  two neighbouring columns (elevation + 1 m); a probe exactly on a node
  belongs to the cell on its + side (margin 1e-9 in xi, eta).

32.2 Validation  [verification/update20_moving_background, tools/bgr_cases.py]
  (a) Towed cylinder (Cylinder_Overset_1phase grid, cases/Cylinder_MovingBGR):
  fluid at rest, inlet ramped over 0.1 s to 1.0 m/s, body 1 towed
  upstream at 0.5 min(t/0.1, 1) m/s with the background following it,
  against the fixed cylinder with the inlet ramped to 1.5 m/s - the same
  relative motion; the frames differ by a uniform acceleration during the
  ramp.  Pressure force during the ramp: fixed -82 191 N = rho V (1 + C_a)
  a with C_a = 1.00 (theory -82 466); towed -68 776 N = rho V a_fluid +
  rho V C_a a_rel (theory -68 722); difference = Froude-Krylov rho V
  (15 - 10) to 2.4 %.  The body moving through the static background
  (Moving_BGR = 0) gives the same forces (0.1 % during the ramp, 0.1 % of
  1/2 rho U^2 D L up to t = 0.2 s); Ovst_Reassemble_Tol = 0.02 m re-uses
  the assembly in 59 of 60 steps and reproduces the every-step run to
  3e-6.  Open item (moving-grid path of update 12, not the background):
  after the ramp the towed pressure force settles 235 N (6 % of 1/2 rho
  U^2 D L) below the fixed one, also without overset (Cylinder_1phase as
  one body), proportional to the towing speed, in the separated region;
  the impulsive start and the towing-tank ramp are not frame-equivalent
  set-ups (verification README, section 1) - the ramped inflow is.
  (b) Restart (towing-tank variant, constant surge): the background
  record restored, steps 31-33 identical to the continuous run, <= 1.2e-5
  afterwards; Cylinder_Overset_2phase restarted at step 5 identical in
  every monitor digit.
  (c) Waves (cases/WaveTank_MovingBGR): the update-18 tank with a dummy
  overset block as body 1, the tank following its imposed surge (0.15
  m/s, 1.2 m in 10 periods) or yaw (2 deg about x = 2.5 m, tank 0.1 m
  wide); earth-fixed probes at x = 2.0, 2.5, 3.0 m, Relax_L = 1.0 1.5 0 0
  0 0 moving with the tank.  Last two periods, crest-to-trough height /
  2a at x = 2.0 / 2.5 / 3.0 m: surging tank 1.018 / 0.966 / 0.942 (phase
  lag -6 / -13 / -17 deg: the inlet zone has moved towards the probes);
  the same tank and block with nothing moving 1.034 / 1.022 / 1.047
  (-21 / -26 / -29 deg, the update-18 numbers); the yawed tank 1.037 /
  (0.91-1.00 of a at the trough instants) / 1.041 with lags -21 / - /
  -23 deg (the first run's probe column counted a yawed cell twice, fix
  (f); the repeated 9-period run missed the probe on the yaw-centre node
  by round-off, now margined).  The translating grid adds 2-10 % damping
  to the wave with the first-order time scheme; the yawed tank shows
  none; the background pose follows the body exactly in both runs.

32.3 Regression
  Cylinder_Overset_2phase (10 steps, Moving_BGR = 0): monitor and .p3d
  bitwise identical to update 19.  Changed on purpose: the reported
  friction / wall-function u_tau during a towing-tank ramp (32.1 b), the
  first step after a restart of a ramped case (32.1 c), moving-body
  results at O(a dt) during accelerations (32.1 d).

## 33. Update 21 — PID controllers (2026-10-12)

33.1 Design  [pid_control.F90 (new, before sixdof.F90), sixdof.F90, force_model.F90, restart.F90]
  A controller is a block of PID_* lines inside a Body / Child_body /
  Sub_child_body group or inside a Force_model group of $SixDOF, opened
  by PID_use = 1 (several blocks per force model, up to 4; one per grid
  part).  sixdof_nml / force_nml forward every pid_* key to pid_nml of
  the new module, which keeps the controllers in pidc(1:npid) (MAXPID =
  64) with their owner (kind 1 = grid part b/c/s, kind 2 = force model
  id).  Three types (PID_type): 1 throttle, output = rotation rate
  [rev/s]; 2 steering, output = angle [deg]; 3 direct force, output =
  one component [N or N m] of a Force_type-3 point force (PID_Component
  = Force_Fx .. Force_Mz), ADDED to the formula value (a point force now
  needs a formula or a controller on at least one component).
  Actuation:
  - grid part + type 1/2 (pid_setup_parts, called in sixdof_precheck
    before the moving flags): the part is rotated about PID_Axis (grid
    frame at t = 0, converted to the E-aligned mother frame by Fm) through
    its Rot_Center by the Rodrigues rotation R(axis, angle); the angle
    is integrated from the rate (throttle: angd = 2 pi n) or set by the
    output (steering: ang = delta, angd = d(delta)/dt).  Its Motion_DOF is
    ignored (set to 0 with a NOTE), own_moving = .true., pid_rot = .true.;
    update_transform uses Rrel = R, om_rel = axis * angd in place of the
    Euler-angle pose, and pid_rotate (in sixdof_advance in place of
    body_step) writes the equivalent Euler angles / rates into eta, nu,
    etad so that the motion files and the coupling (a child's wrench to
    its mother, a mother's motion to its children) stay as in update 14.
    A controlled Body therefore has no free DOF; a controlled child is a
    steered rudder / drive or a motor-driven shaft whose speed is set by
    the controller.
  - Hough-Ordway / BEM disk + type 1: the output replaces Prop_RPS
    (Prop_RPS = the initial value, PID_Out_init is ignored with a NOTE).
  - disk + type 2: PID_Axis is required; the disk points P1, P2 are
    rotated about PID_Axis through PID_Axis_Center (default: the disk
    centre (P1 + P2)/2) by the output angle before the mother's pose is
    applied (move_points); the rps stays Prop_RPS (or the throttle's
    output - a disk may carry both).
  - point force + type 3: pid_add(comp) is added to the force / couple of
    point_force in the mother's body frame.
  Measurement (pid_measure / sixdof_node_pose): the controlled variable
  is a 6-DoF quantity of ANY part, PID_Target_body [, _child_body [,
  _sub_child_body]] (a chain: the sub-child id needs the child id),
  PID_Target_type 1 = location eta (x y z [m], roll pitch yaw [deg]; a
  Body in the earth frame, a child / sub-child relative to its mother)
  or 2 = velocity nu (u v w [m/s], roll_rate pitch_rate yaw_rate = p q r
  [deg/s], body frame), PID_Target_var by name.  PID_Target_value lists
  up to 32 targets reached in sequence: when |e| <= PID_Target_tol for
  PID_Target_hold_T seconds the next value becomes the target (the last
  one is held).  Angle errors are wrapped to (-180, 180].
  Law (pid_step, once per step on every rank with the replicated 6-DoF
  state, before the body / force update):
    out = out0 + P e + I int(e dt) - D d(meas)/dt,   e = target - meas
  (derivative on the measurement: no kick at a target change), clipped
  to PID_Out_limit and to +- PID_Rate_limit dt per step, with
  conditional integration as the anti-windup: the integral used in a
  step is the one accumulated up to the previous step, and the step's
  error is added for the next step unless a limiter - output or rate
  limit - cut the output in the error's direction (found with the
  steered foil: the integral wound up during the rate-limited ramp and
  overshot the target; a back-calculation variant was tried and
  rejected, since pulling the integral back to a rate-limited output
  cancels the proportional action); before PID_Start_T and on the first call the output is the
  initial value and the derivative memory is primed (a grid part acts
  from the second step, a force model from step 1).  The gains carry the
  sign of the actuator (positive gains when a positive output drives the
  measured variable up).
  Output: Output_data/<case>_<stem>_pid<i>.dat (stem = body1,
  body1_child1, body1_child1_sub1, force_model<n>; i = the controller's
  index within its owner): step, time, target index, target, measured,
  error, P / I / D terms, output, angle [deg] (the grid-part angle or
  the disk steering angle); the summary prints every controller and the
  parts rotated by one.  Restart: the states (target index, hold timer,
  integral, previous measurement, output, angle, started) are appended
  to the restart file after the field arrays; an older file without
  them (or with another number of controllers) starts the controllers
  fresh with a NOTE; the data file is appended.
  Checks (pid_check, pid_setup_parts, force_init): type in 1..3, target
  part given and existing (a force model needs the 6-DoF module on),
  type 1/2 of a force model needs a disk, type 3 needs Force_type 3 and
  PID_Component, PID_Axis for every grid part and for a disk steering,
  limits ordered, rate limit / tolerance / hold >= 0, one controller per
  grid part, 4 per force model.
  Fixes found by the restart test of the controllers (33.2 b):
  (a) several force models on one mother: sixdof_set_ext overwrote the
  mother's external wrench, so only the last model's force acted on the
  body (update 19; the fluid sources were complete) - force_compute now
  clears the wrenches (sixdof_clear_ext) and the models ADD theirs;
  (b) restart with force models: force_compute(0, 0) ran before the
  restart file was read, so the first restarted step used the momentum
  sources and wrenches of the INITIAL field and poses - Main.F90 now
  re-evaluates the force models after read_restart with the controllers
  held (pid_hold: pid_step returns the restored output without touching
  the state; no duplicate export rows); (c) the velocity-Verlet half-kick
  flag (vhalf, update 19) was not restored, so the first restarted step
  of a free DOF skipped the half-kick 1/2 dt (a^n - a^(n-1)) - a 0.3 %
  position offset in the 100-kg box tests; sixdof_restore_pose(t) sets
  the flag for every free DOF released at the restored time; (d) the PID
  states were read after the restart file had been closed (the "0
  controller states" NOTE) - read before the close; the PID-rotated
  part's angle and rate are set from the restored controller state.
  With (b)-(d) the restarted box reproduces the continuous run to all
  printed digits (12 significant); before, the position differed by
  0.35 % from the first restarted step on (also without a controller:
  the formula spring of Foil_Spring restarted the same way).

33.2 Verification  [verification/update21_pid, tools/pid_cases.py, cases/PID_*]
  (a) Steering of a grid part (PID_Steer): the NACA66mod_Overset foil
  as body 1 rotated about its span axis through the quarter chord by a
  steering controller on its own pitch, targets +4 / -4 deg (tol 0.1,
  hold 0.01 s), P 0.2, I 300, rate limit 300 deg/s, 50 steps of 2 ms:
  ramp at the rate limit (0.6 deg per step), +4 deg within 0.1 deg at
  t = 0.020 s (overshoot 4.18), switch at 0.036 s, monotone ramp down
  (max change per step exactly 0.6 deg), -4.000 deg from t = 0.092 s;
  pitch and q of the motion file equal the controller's angle and rate
  to all digits.  The loop of a kinematic rotation is algebraic (angle
  = output), stable for P < 0.7; the first anti-windup (integral frozen
  only at the output limits) wound up 4 deg during the ramp - now
  conditional integration at both limits.  (b) Throttle of a grid part
  (PID_Throttle): the Cylinder_MovingBGR cylinder spun about the
  vertical axis to its own yaw rate 60 deg/s (P 0.001, I 0.3; the
  algebraic loop's eigenvalues 0.65 / -0.55): the output starts at P e = 0.06 rev/s, the yaw rate follows
  kinematically (360 x output), overshoots to 61.3 deg/s at step 6 and
  is exactly 60.000 deg/s (output 1/6 rev/s) from step 14 of 40; with
  the first gains (P 0.002, I 0.05) the slow eigenvalue 0.95 left a
  4.5 deg/s error after 40 steps.  A side force of -210..-225 N (port)
  appears with the spin at 0.2 D/U - the early boundary-layer
  displacement effect, opposite to the steady Magnus lift (several D/U
  to build up; not the subject of the test).  (c) Force
  models on a damped 1-DOF mass (the ActuatorDisk box declared body 1,
  no wall: m = 100 kg, Motion_Damp 2 or 20 /s): a direct-force
  controller on Force_Fz (P 4e4 N/m, D 2500 N s/m) brings the heave to
  +5 mm at 0.14 s and, after the hold, to -5 mm at 0.34 s (-5.000 mm with
  a zero output at 1 s); a disk throttle (rps in 0..12, T = rho n^2 D^4
  K_T(J)) pushes the box to +5 mm at 0.89 s and +10 mm at 1.47 s with n
  of the force-model file equal to the output; a disk steering about
  -y through the disk centre (n = 10 fixed, +-40 deg, 200 deg/s) gives
  Fz = T sin(delta) and reaches +5 mm at 0.81 s, +10 mm at 1.44 s.
  (d) Restart (direct-force box at step 100 of 200): controller and
  motion files identical to the continuous run in every printed digit
  after the fixes 33.1 (a)-(d).

33.3 Regression
  Cylinder_Overset_2phase (10 steps) and ActuatorDisk (5 steps, a disk
  without a controller): monitor, .p3d and force-model file bitwise
  identical to update 20.  Changed on purpose: several force models on
  one mother all act (33.1 a); restarts of cases with force models or
  free DOFs continue exactly (33.1 b-d).

## 34. Update 22 — probe points, shear-stress field, wall surface files (2026-10-13)

34.1 Probe points  [probes.F90 (new, after surface_output.F90), sixdof.F90 (Probe_point lines of the parts),
                    read_input.F90 ($Probe_background), Main.F90]
  Probe_point<n> = x y z [m] (grid frame at t = 0): in the new group
  $Probe_background the point follows the background (fixed, or moving
  with the background blocks through the published pose bgr_Rg / bgr_O0
  / bgr_dg when Moving_BGR = 1); inside a Body / Child_body /
  Sub_child_body group of $SixDOF (sixdof_nml stores them in prdef(:),
  probe_init resolves the node) the point follows the part's 6-DoF
  motion, x = Rg (x0 - O0) + O0 + dg of the part - the formula of the
  force-model points, so a point on a propeller rotates with it.  Every
  Data_Intvl steps (and step 0) each point is located in the current,
  moved grid: the +-2-cell neighbourhood of its last cell first
  (inv_trilinear), then the bucket search of every block (ovs_geom
  build_buckets / locate_point; the bucket of a moving block is rebuilt
  when needed, that of a static block once).  The grids are replicated,
  so every rank finds the same candidate cell per block; the owner rank
  of each candidate evaluates it when it is not a hole and the best key
  wins (active cell before fringe cell, then the lower block, then the
  lower rank) through an Allreduce; the values are the trilinear
  interpolation (the inverse map's r s t) of the eight node values
  formed as in the volume files (face averages of u v w, cell averages
  of p, alpha, vapour, T and of the cell shear stress, holes excluded).
  No candidate (a hole - inside a solid, behind a moving body -, outside
  every block, or a cell no rank holds) writes "NoData" for the step.
  Files Output_data/<case>_<stem>probe<n>.dat (stem = '' / body1_ /
  body1_child1_ / body1_child1_sub1_; a point number used twice in one
  owner is refused): step, time, x y z (current position, grid frame),
  u v w [m/s], TTL_Prs, Dym_Prs = p - p_hyd(z) [Pa] (the Dym_Prs rule
  of the volume files, also in single-phase runs), [VoF] [Vapor]
  [Temperature], Shear_Stress [Pa], block i j k, status.  Appended on a
  restart.

34.2 Shear_Stress field of the volume files  [surface_output.F90 cell_shear, output.F90]
  New variable directly after TTL_Prs in every .p3d (the .nam file lists
  it; the variable count of the files grows by one): at the cells
  tau = mu_eff sqrt(2 S_ij S_ij) with S from the cell-centred velocity
  gradient of turbulence.F90 (velocity_gradient, now public) and mu_eff =
  laminar mixture viscosity (mu_lam: phase / vapour weighted) + the eddy
  viscosity mut of the turbulence or IDDES model (zero without a model),
  dimensional with mu1 Uch / Lch - for simple shear this is mu_eff du/dy;
  node values = averages of the surrounding cells whose gradient stencil
  lies in the domain's ghost range (holes excluded).  At every solid-wall
  node (no-slip flag-1 patches with zero values and wall-function flag-11
  patches, with or without a body) the value is replaced by the friction
  stress magnitude |tau_w| of the adjacent wall faces (the node average
  of 34.3), i.e. the wall-function tau_w where one is used.  The field
  and the wall value differ in a thin, under-resolved boundary layer
  (the wall value is u_1 / dn of the first cell, the field the central
  difference across two cells: e.g. 2.0 Pa at the wall against 0.5 Pa at
  the first interior node of the impulsively started cylinder at
  t = 0.025 s) - the optional smoothing of the request was NOT applied:
  the field shows the discrete gradients as they are and the wall value
  stays untouched.

34.3 Solid-wall surface files  [surface_output.F90, sixdof.F90 wall_face_stress / wall_patch_range]
  $Export_Interval: SW_Prs_Frc_surf = 1 writes, with every volume
  solution (Output_Intvl, step 0, the last step), Volume_solutions_grids/
  <case>_surf_NNNNNNNN.grd / .p3d: one PLOT3D block of ni x nj x 1 nodes
  per wall patch (every block, bodies or not), six variables [Pa] Prs_x
  Prs_y Prs_z = -p n (n into the fluid: the pressure force per unit area
  on the wall), Frc_x Frc_y Frc_z = the friction stress vector on the
  wall; <case>_surf.nam and <case>_surf_patches.txt (zone -> block,
  patch, face, ranges, body) once.  The per-face values are those of the
  force integration: sixdof_forces' face loop was factored into
  wall_face_stress (pressure of the first cell, mu_eff u_t / dn with
  the ramped Dirichlet + body wall velocity and the wall-function
  mu_eff of mu_h) and wall_patch_range (the domain's part of a patch),
  used by both - the force files are bitwise unchanged.  Node values:
  v_n = sum_{f at n} w_f v_f A_f/4 / sum_{f at n} A_f/4 (w_f = overset
  wall weight, A_f = face area; a hole face counts with zero stress and
  its area, as the force integration skips it), so that the trapezoidal
  rule over the patch faces, sum_f A_f (mean of the 4 node values),
  equals sum_f w_f v_f A_f exactly: integrating the surface file gives
  the force-file values (check_surface.py: towed cylinder, relative
  difference 9e-8 for the pressure, 2e-9 for the friction - the
  single-precision storage of the PLOT3D file).  Frame_Prs_Frc_surf = 0:
  the grid axes of the volume files (x y z of the grid, the moved
  coordinates when the background moves) - note the force files use the
  SNAME earth axes (X = -x, Y = y, Z = -z); = 1: the background frame
  (coordinates and vectors turned back by the background pose, so the
  surface sits at its original grid position; identical to 0 unless
  Moving_BGR = 1).  Both switches can be changed at run time (name-list
  re-read).  The surface data of a step are computed once (wall_node_data)
  and shared by the volume file's wall-node Shear_Stress and the surface
  file.

34.4 Verification  [verification/update22_probes_surface]
  (a) Probe points: on the towed cylinder (cases/Probes_TowedCylinder,
  20 steps) the body-following point 0.1 m behind the cylinder moves as
  x = 0.6 - 0.5 t^2/0.2 (0.57500 at t = 0.1 s) and stays in the same
  O-grid cell; the background points follow the moving background
  (x = 1.975 at 0.1 s); the point at the cylinder's centre writes NoData.
  On the steered foil (cases/Probes_SteeredFoil, 5 steps, 1.88 deg) the
  trailing-edge point follows z = -0.97 - 0.8 sin(theta) to 1e-5 and the
  point above the quarter chord x = 0.25 + 0.07 sin(theta).  (b)
  Shear_Stress: equals |Frc| of the surface file at the wall nodes
  (single-precision round-off); the first interior node reads 0.5 Pa
  against 2.0 Pa at the wall in the impulsively started cylinder's
  0.16-mm boundary layer (34.2).  (c) Surface files: the trapezoidal
  integral of the surface file equals the pressure and friction force of
  the force files to 1.9e-8 / 3.5e-9 (cylinder, step 20), 1.4e-8 / 2.1e-8
  (overset foil with the wall-trimming weights, step 5) - the
  single-precision storage; Frame_Prs_Frc_surf = 1 keeps the towed
  surface at x = -0.5..0.5 m while the frame-0 file moves with the body.

34.5 Regression
  Cylinder_MovingBGR (5 steps, towed cylinder with the background
  following, surface files on): monitor, force and motion files bitwise
  identical to update 21; the volume file has the new Shear_Stress
  variable and all other variables bitwise identical.

## 35. Update 23 — on-the-fly turbulence statistics (2026-10-13, 23b face-difference gradient 2026-10-18)

35.1 Design  [statistics.F90 (new, after output.F90), probes.F90 (shared sampler), read_input.F90 ($Statistics),
              restart.F90, Main.F90, mod_global.F90]
  $Statistics with SW_statistics = 1: raw moments accumulated every
  Statistics_Sample_Intvl steps inside [Statistics_Start_T,
  Statistics_End_T] (0 = End_T), the fluctuation statistics formed at
  output: mean = S/n, rms = sqrt(max(Q/n - mean^2, 0)), <x'y'> = C/n -
  mean_x mean_y, with per-point sample counters n1 (the point valid: not
  a hole) and n2 (the gradient stencil valid) so that blanked steps of a
  moving-grid point are left out; a point never sampled is written as 0
  with Coverage (= n1/N_samples) = 0 - the PLOT3D files hold numbers
  only.  Sums (not running means) so that partial sums of different
  ranks and of a restart combine.  Two output groups: Statistics_first_
  group (Coverage, mean and rms of u v w TTL_Prs Dym_Prs Shear_Stress
  [VoF] [Vapor] [Temperature] [nu_t] [k_model, omega/eps_model or
  nu_tilde], Reynolds stresses uu vv ww uv vw wu, k_res) and
  Statistics_second_group (velocity only: Coverage_grad, 9 mean
  gradients, 9 rms fluctuating gradients, 3 cross terms, eps_res = nu
  <(du'_i/dx_j + du'_j/dx_i)^2>/2 = nu [2 sum_i <(g'_ii)^2> + sum_{i/=j}
  <(g'_ij)^2> + 2 (<g'_12 g'_21> + <g'_23 g'_32> + <g'_31 g'_13>)] with
  nu = the mean laminar kinematic viscosity at the point; with a
  turbulence / IDDES model also eps_model, eps_total, k_model, k_total;
  eta = (nu^3/eps_total)^(1/4), lambda_g = sqrt(10 nu k_total/eps_total),
  Re_lambda = sqrt(2k/3) lambda/nu, Delta/eta with Delta = the mean
  local cell size, and the mean vorticity).  Modelled terms per sample:
  k-epsilon k_model = k, eps_model = epsilon; k-omega / SST (also with
  IDDES) k_model = k, eps_model = beta* k omega; SA and the WMLES branch
  eps_model = nu_t |S|^2, k_model = nu_t |S| / sqrt(beta*) (equilibrium
  estimates); none: 0.
  Where: Statistics_type1 = the original grid blocks, NODE based (the
  node values of the volume files: face / cell averages with the holes
  excluded, the gradient averaged from the cells whose stencil lies in
  the ghost range), written on the t = 0 grid (<case>_stat_grid.grd
  stored at init); Statistics_Type1_Frame = 1 turns the vectors and
  the gradient tensor of a moving block into its part's frame (R^T u,
  R^T G R with R the part's Rg; the background blocks of Moving_BGR
  through the background node) so that a propeller block gives
  blade-fixed statistics; coincident nodes of different blocks (block
  cuts, O-grid seams) get the average of their statistics at output,
  the rms columns averaged as variances so that k_res and eps_res of
  the file stay consistent with the exported components (23b fix: the
  first version averaged the rms linearly there, which made eps_res
  recomputed from the components differ by up to 7 % at the O-grid
  seam nodes of the towed cylinder).  Statistics_type2/3/4 = boxes
  / cylinders / planes: structured lattices (Statistics_box/cyl/plane =
  n ... _End = -1 groups: Origin / Facecenter, Axis1 / Axis, Axis2 /
  Ref, Length, nPoints, Radius, Theta, Polar_Vel, Mother_body /
  Mother_child_body / Mother_sub_child_body) whose points are defined
  in the grid frame at t = 0 and move with their frame (x = R (x0 - O0)
  + O0 + dg of the part or the background); every sample step each
  point is located in the moved grid with the probe machinery of update
  22 (pt_locate: last cell neighbourhood, then the bucket search per
  block; pt_owner_key: the active cell before a fringe cell; one
  Allreduce of the keys per lattice) and sampled by trilinear
  interpolation (pt_sample); the vectors are turned into the frame
  (and for a cylinder with Polar_Vel = 1 into axial / radial /
  tangential components, the gradient tensor P^T G P).  The sampler of
  probes.F90 was generalised to a fixed layout of NSAMP = 24 values (u
  v w, p, p - p_hyd, tau, alpha, vapour, T, nu_t, k, omega/eps/nu_tilde,
  the 9 gradients, nu, Delta) used by the probe files (unchanged
  output) and the statistics.
  Output (rank 0) into Statistics_volume_solutions/: <case>_stat_grid_g1
  / _g2_NNNNNNNN.p3d on <case>_stat_grid.grd; <case>_stat_box<n>_g1/
  _g2_NNNNNNNN.p3d on <case>_stat_box<n>.grd (the lattice at t = 0),
  the same for cyl<n> / plane<n>; a .nam name file per family with the
  SI units ("u_mean [m/s]"); every Statistics_Output_Intvl steps (0 =
  end only; changeable at run time) and at the end.  The volume-file and
  surface-file .nam files now carry the units too ("u [m/s]", "TTL_Prs
  [Pa]", "Prs_x [Pa]").  Restart: the type-1 node sums per domain and
  the reduced lattice sums are appended after the PID states (a file
  without them, or with another accumulator layout, starts fresh with a
  NOTE); the restored lattice sums stay on rank 0 and combine with the
  new partial sums at output.  Memory: per type-1 node 2 nf1 + 3 + 25
  doubles (<= 52) plus two counters.

  Velocity gradient of the statistics (23b, replacing the central
  difference of the first version): velocity_gradient_compact
  (turbulence.F90) builds the cell-centred tensor from the FACE
  DIFFERENCES of the staggered velocities - d u/d xi, d v/d eta,
  d w/d zeta are the compact differences of the cell's own two faces
  (the differences the solver's viscous term and divergence use; they
  resolve the 2-cell waves the central difference loses: a 4-cell wave
  keeps 0.90 of its gradient against 0.45), the other six derivatives
  the central difference of the cell-centred values across the two
  neighbouring cells (= the average of the compact differences at the
  four edges of the cell in that plane; the cell-centred value of a
  transverse derivative cannot be sharper without forming the squares
  at the edges, which would defeat the part-frame rotation of moving
  blocks).  With it the type-1 statistics of group 2 are accumulated
  AT THE CELLS (the cell's own tensor, no 8-cell node averaging of the
  signal, cell_grad_sample of probes.F90; the node values of the files
  are the averages of the per-cell statistics over the node's sampled
  cells, the rms columns averaged as variances so that eps_res of a
  node recomputed from its exported components equals the average of
  its cells' eps_res; coverage of the gradient group = the cell
  counter); the lattices sample the interpolated compact tensor.  The
  Shear_Stress field of update 22 and the turbulence models keep the
  central gradient.  The cell accumulators (6 velocity sums for k_res
  + the 25 of group 2) are in the restart file; the first version's
  files are refused with a NOTE (layout code nacc + 1000 na2c).

35.2 Verification  [verification/update23_statistics, tools/stat_cases.py]
  (a) cases/Statistics_ActuatorDisk (static Cartesian box, 10 steps, the
  solution written every step, statistics of steps 4-10): the type-1
  mean / rms of u v w TTL_Prs Shear_Stress, Rss_uv and k_res agree with
  the statistics of the per-step volume files to the single precision of
  the files (max difference 8.5e-8 of 1.29 m/s in u_mean, 1.7e-5 of 282
  Pa in TTL_Prs_mean, 1.2e-9 of 5e-3 in Rss_uv); eps_res recomputed from
  the exported 12 gradient components agrees to 3e-13 (type 1 with the
  nodes' averaged variances, and the lattices); the face-difference
  gradient raises the peak eps_res of the starting disk from 3.2e-6
  (central, first version) to 5.8e-6 m2/s3 - the 2-cell content of
  the velocity jump across the disk that the central difference
  cannot see.  Restarted at step 7:
  the type-1 files of step 10 bitwise identical to the continuous run,
  the lattice files to the last bit of the single-precision files
  (<= 1e-7 relative; summation order of the restored and the new partial
  sums).  k-omega SST on the same box: nu_t, k_model,
  omega_model means and eps_model = beta* k omega (1.97e-4 m2/s3 =
  0.09 x 1.5e-4 x 14.6), eps_total, k_total, eta, lambda, Re_lambda,
  Delta/eta all finite and consistent.
  (b) cases/Statistics_TowedCylinder (towed cylinder with the moving
  background, 20 steps, statistics of steps 4-20): the lattice point of
  the body-following box that coincides with the body-following probe
  point and the background-box point that coincides with the
  background probe point reproduce the mean / rms of the probe time
  histories to 1e-7 relative (the same sampler); the background plane
  across the cylinder shows Coverage 0 inside the solid (the hole of
  the background block, 9 x 9 points of the 17 x 17 lattice at 0.1 m)
  and 1 outside; the ring following the cylinder with polar velocities
  gives the antisymmetric u_tan (-1.15 / +1.15 m/s at theta = 60 / 300
  deg at r = 0.55 m) and the symmetric u_rad of the flow around the
  cylinder, u_axial = 0 (2-D); the type-1 coverage is 1 on the O-grid
  blocks and 0..1 on the background block (hole region).  The probe
  files and the force file are bitwise identical to update 22 (the
  sampler refactoring changed no value).

35.3 Regression
  Probes_TowedCylinder (20 steps, statistics added): probe files, force
  and motion files bitwise identical to update 22.  The volume .nam file
  changed on purpose (units appended to the names).

## 36. Update 24 — fluid-structure interaction of flexible parts (2026-10-19)

36.1 Design  [fsi_shell.F90 (new, structural solver, no MPI), fsi.F90 (new,
     coupling), sixdof.F90 (FSI keys of the part groups, morph / grid-velocity /
     flux hooks, "moving" flag of flexible parts), mod_global.F90 (hook
     pointers, fsi_blk, fsi_reasm_request), curvi.F90 (half-index wall
     velocity), convection.F90 (geometric-conservation term of morphing
     blocks), boundary.F90 / turbulence.F90 (wall velocity of morphing
     blocks), Main.F90 (predictor, coupling loop in the inner iterations,
     finish), restart.F90, surface_output.F90 (Disp_x/y/z, vonMises),
     fsi_test.F90 (stand-alone structural tests, make fsi_test)]
  SW_FSI = 1 in a Body / Child_body / Sub_child_body group makes the part
  FLEXIBLE.  Inputs (sixdof.F90 fsi_nml): Youngs_modulus, Poisson_ratio,
  Density (Mass = Density x enclosed volume) or Mass (Density = Mass /
  volume), Num_NoFSI_Group + NoFSI_block_num<n> / NoFSI_i_range<n> /
  NoFSI_j_range<n> / NoFSI_k_range<n> (clamped wall-node ranges, one value
  = that index, two = a range), Num_FSI_Submat_regions + Submat_block_num<n>
  / Submat_i/j/k_range<n> / Submat_Density<n> / Submat_Youngs_modulus<n> /
  Submat_Poisson_ratio<n>, and the optional FSI_Start_T, Damping_ratio,
  FSI_Min_Thickness, FSI_Max_Iter, FSI_Tol, FSI_Relax, FSI_Mode,
  FSI_Monitor_node<n> (block i j k).  Time_scheme = 1 is required.
  * STRUCTURE FROM THE SURFACE GRID (fsi.F90 build_structure).  Every
    no-slip wall patch of the part's blocks is a structured sheet.  Its
    AROUND direction is the periodic one (O-grid) or the one whose two end
    lines have opposed normals (C-grid), the other is the SPAN.  At every
    span station the chord is the principal axis of the section's nodes,
    n_up = chord x e_span (sign continued along the span), and a node is
    upper (n . n_up > 0.5), lower (< -0.5) or cap.  The longest run of
    upper nodes at the middle station fixes the upper index range used at
    all stations.  From every upper node a ray along the inward normal is
    intersected with the lower / cap faces (Moeller-Trumbore): mid-surface
    node = midpoint, thickness = ray length (a miss or a non-opposed hit
    falls back to the lower node nearest to the ray; floor
    FSI_Min_Thickness x the largest thickness outside the NoFSI groups).
    Shell elements = the quads of the nup x nspan lattice; coincident mid
    nodes of different sheets of one part (block cuts) are merged.  A
    sheet that is NOT sheet-like (no periodic / folded direction, fewer
    than 2 upper nodes) or that OVERLAPS a larger sheet of the same part
    (more than 40 % of its wall nodes within one node spacing of the
    larger sheet's wall nodes: the tip-cap blocks of the KP505 blades)
    is PASSIVE: it adds no structural nodes; its wall nodes are attached
    to the mid-surface of the active sheets.  Every wall node (upper,
    lower, cap, passive) is attached to the mid-surface by the bilinear
    weights of its projection onto the nearest element (nearest node
    beyond the edges) and the rigid offset r0 = X_wall - X_node of each
    attached node; x_wall = sum_a w_a (X_a + u_a + R_a r0_a).  NoFSI
    groups clamp the mid nodes whose upper or lower source node lies in
    the range; a Submat region owns the elements with >= 2 nodes in its
    range.  Enclosed volume by the divergence theorem over the wall faces
    (the median of the three int x_c n_c dA estimates, which disagree
    only where the surface is open: a span end, a blade root).
  * SHELL (fsi_shell.F90): 4-node flat shell in a co-rotational frame -
    bilinear plane-stress membrane + MITC4 Mindlin-Reissner plate
    (assumed transverse shear, no locking) + Hughes-Brezzi drilling term;
    constant local stiffness, the nonlinearity in the kinematics (element
    frame from the current corner positions, local displacements =
    current - initial local coordinates, local rotations = log(E^T R_a
    E0) with multiplicatively updated nodal triads); f = T f_l, K = T K_l
    T^T (simplified co-rotational tangent, symmetric; Newton converges
    linearly - 2-3 iterations per step in the benchmarks).  Per-element
    E, nu, rho, t.  Lumped mass (rotational m (t^2/12 + A/48)), Rayleigh
    stiffness damping a_K = 2 zeta / omega_1 (omega_1 by inverse
    iteration at set-up).  Generalized-alpha (rho_inf = 0.8) with Newton
    on the nonlinear residual; shell_step repeats the step from t_n with
    new loads (sub-iterations), shell_accept commits; shell_static for
    FSI_Mode = 0.  Banded Cholesky after a reverse Cuthill-McKee
    numbering (bandwidth ~ 6 x the shorter lattice side).
  * LOADS AND MOTION (energy-conserving transfer): the wall node forces
    (-p n + tau) A/4 of the faces around a node, overset wall weights
    included (the same values the force files integrate), go to the
    attached mid nodes a as f_a = w_a f, m_a = w_a (R_a r0_a) x f.  The
    structure lives in the part's t = 0 frame: forces are turned by
    Rg^T, displacements by Rg; gravity and the inertial loads of a moving
    frame (-m [a_O + alpha x r + omega x (omega x r) + 2 omega x v]) are
    added (a rotating propeller blade sees its centrifugal load).
  * MESH MORPHING (fsi_morph, called at the end of sixdof.move_grids):
    along every grid line normal to the wall (the column of a wall node)
    x = X0 + d_w + (R(theta_w) - I)(X0 - X_w0): the whole column moves
    rigidly with its wall node (translation and rotation), the block
    bends with the wall like a thick shell and the overset outer boundary
    follows (24b; the first version decayed the rotation with a cosine
    from 0.25 of the arc length to the outer boundary, which tangled the
    fan of columns around the plate's rounded end once the tip rotation
    reached 0.3 rad - the columns of the cap crossed those of the flat
    sides).  The rotation of a column is that of its wall node only near
    the wall: level by level away from the wall the wall rotation is
    diffused along the wall (around / span) by nint(100 sbar^2) Jacobi
    sweeps, sbar = the mean arc-length fraction of the level (24c) - rigid
    columns with the wall rotation itself sheared the outer cells between
    the clamped columns and the first rotating ones of the 20-mm fan at
    the plate root of the one-body grid: at a 22 mm tip deflection the
    cell area there fell to 6 % and the imposed-motion test diverged at
    full amplitude (with the diffusion 80 %, 64 % at 30 mm).  The
    displacement of the outer levels is then diffused along the wall as
    well, by nint(200 sbar^2) sweeps (24c): on the concave side of the
    swiveling plate the columns converge - the free benchmark run
    diverged at t = 0.65 s with 6 % of the cell area left there; the same
    state re-morphed: 37 %.  Both diffusions act per tangential direction
    (a field uniform along the span stays uniform: the k planes of a
    quasi-2-D grid stay identical - averaged over all neighbours, the end
    nodes of a 3-node span got different weights and the k planes drifted
    0.03 mm apart); then the rigid pose of the part.  The metrics of the
    morphed sub-domains are rebuilt, the discrete grid velocity (x^(n+1)
    - x^n)/dt of every half-index position feeds the ALE grid fluxes and
    the wall velocity (curvi / boundary / turbulence take the half index
    for morphing blocks), the stored face fluxes are re-projected with the
    new metrics (U_new = U_old + fvw_new (m_new . u) - fvw_old (m_old . u),
    the Cartesian face velocity carried; fvw of the old metrics follows
    curvi_form), and the momentum equation carries the geometric-
    conservation term phi div(u_grid) on morphing blocks.
    THE OLD METRICS OF THE RE-PROJECTION (24c fix): move_grids rebuilt
    Xh / mh of every moving sub-domain by the rigid rotation of the
    UNDEFORMED Xh0 / mh0 before calling fsi_morph, which then saved those
    as the "level-n" metrics m_old.  The re-projection therefore added
    fvw (m_deformed - m_undeformed) . u to the fluxes in every step - zero
    while the plate is straight, growing with the deflection, and present
    even on a grid held still (FSI_Test_Motion mode 1: a static 15 mm
    deflection blew up 18 steps after the hold began, while the same
    deformed grid run as a static case was stable; the imposed 22 mm
    bending mode blew up at the trailing end at full amplitude).
    move_grids now skips the morphing sub-domains (fsi_morph builds their
    coordinates, metrics and grid fluxes), and the hold test then follows
    the static run step for step.
  * COUPLING (Fabbri et al. 2023 Fig. 7): at the step start the interface
    displacement is predicted, d0 = d_n + dt (1.5 v_n - 0.5 v_n-1) with the
    DISCRETE velocities of the last two steps (displacement / rotation
    increment / dt, 24c; the generalized-alpha velocity can be far from
    the increment - at the release of the pluck it was 0.9 m/s while the
    tip moved 0.08 mm, and the predicted grid motion blew up the next
    fluid step; the discrete velocity is in the restart file), the
    grid morphed and (when the wall moved more than Ovst_Reassemble_Tol
    since the last assembly) the overset re-assembled; every inner
    iteration of the fractional step (momentum, Poisson, correction from
    the level-n field) ends with wall forces -> structural step from t_n
    with the new loads -> interface residual r = d~ - d_k-1 -> the first
    sub-iteration relaxes with FSI_Relax, the later ones take the
    interface quasi-Newton step (IQN-ILS, Degroote et al. 2009, on the
    residual / output history of the step, QR with column filtering) ->
    convergence ||r|| / ||d~ - d_n|| <= FSI_Tol or ||r|| <= FSI_Tol x
    shell thickness (the pressure-solver tolerance leaves load noise of
    that order on a very compliant structure) -> morph + re-assembly ->
    next inner iteration, at most FSI_Max_Iter.  Every sub-iteration
    after the first restarts the fluid step from the LEVEL-n STATE (24b):
    the Cartesian field and the face fluxes of every owned sub-domain
    (saved at the step start) are restored, the fluxes of the morphed
    blocks re-projected from the level-n grid (its metrics are kept
    through the step) and those of rigid moving blocks rotated again -
    the fluid step is then the same function of the interface
    displacement in every sub-iteration, which the quasi-Newton needs
    (fsi_restore_mode = 1).  The first version reset only the Cartesian
    field and re-projected the latest iterate's fluxes incrementally: the
    fluid solution drifted between sub-iterations and the residual stalled
    at 1e-2 of the step increment (6-12 sub-iterations); with the level-n
    state 4-6 and a residual that falls.  24c tried a mode 2 (level-n
    fluxes only, the Cartesian field and pressure of the latest iterate)
    against the breakdown of the benchmark: the sub-iterations converged
    less well and the run broke down earlier; the breakdown had another
    cause (the old metrics, above), mode 1 stays the default.  With the
    full-pressure form (Prs_form = 1) the pressure is not restored: it is
    only the initial guess of the Poisson solve (the predictor has no
    pressure gradient), and the latest sub-iterate is a far better guess
    than p^n - the benchmark's repeated steps take ~45 % fewer Krylov
    iterations with coupling residuals unchanged to 3 digits (24c; with
    Prs_form = 2 p^n is restored, the predictor uses grad p^n).  The
    quasi-Newton works in a scaled space (the rotation entries times half
    the nodal thickness = the wall offset, commensurate with the
    translations; unscaled, the least squares was dominated by the
    rotations, 24b).  FSI_Start_T holds the structure rigid
    (the flow develops first); FSI_Mode = 0 replaces the dynamic step by
    the quasi-static equilibrium.
  * STRUCTURAL SAFEGUARDS (24b): when the Newton iterations of a whole
    structural step fail (a coupling sub-iteration far from the converged
    interface position can change the load abruptly), shell_step repeats
    the step in 2, 4 and 8 sub-steps with the load interpolated linearly
    from t_n (a NOTE in the log); the result is again a trial state
    relative to t_n.
  * QUASI-2-D CASES (24c).  SW_Quasi_2D = 3 ($Numerical_Scheme): after
    every velocity correction w and U3 are set to zero and the fields are
    averaged along the grid's k (the span; extruded grid, k_split = 1,
    checked): u / U1 and v / U2 in face columns whose faces are all
    solved, p and the scalars in all-solved cell columns.  On an extruded
    grid with slip span walls the averaged fluxes keep the discrete
    continuity (the U3 differences of a column sum to zero).  Without it
    a laminar thin-span run of the benchmark grew w from the cylinder's
    separation line from step ~140 (with 10 mm span and 2 cells a
    spanwise sawtooth of the wake; with SST a spanwise mode at the
    separation plus twice the cost).  FSI_Plane_Strain = 1/2/3 (inside the
    part's group): the span axis - the structure bends cylindrically: the
    displacement along that axis and the rotations about the two other
    axes are held at every node (restored states are projected on
    restart); without it the plate twisted along the span and bent the
    span walls of the grid.  The trial states of the nodes that differ
    only along the span axis (span groups) are averaged after every
    structural solve and on restart (solved independently, the span nodes
    of the free benchmark drifted apart by 0.1 mm in 0.3 s and the fluid
    grid was no longer extruded); a restored grid that is not exactly
    extruded only gives a NOTE (SW_Quasi_2D checks the first step of a
    fresh run strictly).
  * HOLE CUTTING OF SLAB BODIES (overset.F90 build_hsurf / slab_caps,
    24c): a body touching two PARALLEL planar domain-boundary patches with
    the same touching bodies (the two span walls of an extruded grid) is
    closed by end caps on both planes instead of being mirrored: the
    in-plane edges of the body's kept wall triangles are chained by
    coincidence into the cross-section loops, every loop triangulated by
    ear clipping in plane coordinates (collinear vertices dropped), the
    cap normal pointing out of the body ("bodies closed by N end-cap
    triangles on M slab planes" in the log; a loop that does not close is
    reported).  Before, the mirror images across one span wall overlapped
    the domain across the other on a thin span, the body stayed open and
    148 active cells remained inside the plate of the 0.4-mm grid; now 0.
    Bodies touching a single plane (half bodies on a symmetry plane) are
    mirrored as before; the blanking of the other cases is unchanged
    (36.2).
  * THE FSI-PfS-2a GRID (tools/fsi_cases.py pfs2a, 24c): ONE elliptic
    O-grid around cylinder + plate as a single wall patch (the cylinder
    part clamped) instead of a cylinder O-grid with a plate O-grid whose
    root lay inside the cylinder (the staircase cut of the junction
    carried spurious velocities).  Three grid defects were found with the
    imposed-motion test and fixed: (1) the outer boundary followed the
    outline at a constant distance - concave at the junction, a fringe
    cell there had no valid donor (190 m/s): now the convex hull of the
    circles R + 12 mm and h/2 + 12 mm; (2) the seam (i = 1 = ni) next to
    the holes of the plate: now at the cylinder front; (3) a duplicate
    outline node at the end of the rounded tip made a zero-width sliver
    column - the root of the first-row and fringe velocity spikes (the
    "thin layer" of the earlier runs); the generator asserts a minimum
    wall spacing.  The interior: the Winslow equations on an 81-level
    computational grid, every i-line then re-sampled with a geometric
    distribution of first cell 0.08 mm.  Span 0.4 mm (2 cells), laminar,
    SW_Quasi_2D = 3, FSI_Plane_Strain = 3; 201 x 21 + 107 x 61 + 98 x 70
    nodes x 3 (34 k cells).
  * FSI_Test_Motion = A, f [, mode [, shape]] (verification and start-up
    aid, 24c): the structure is not solved, its nodes follow the imposed
    bending mode d = A r xi^2 sin(w tau) e_n (xi = distance from the clamp
    along the chord / free length, r = one-period ramp, the rotation the
    slope 2 A xi / L) - the fluid sees a known wall motion of the
    benchmark's amplitude and speed without structural feedback (a test
    of the moving-grid treatment of the fluid alone).  Needs
    FSI_Plane_Strain (the span axis).  Mode 1: ramp to the static
    deflection A xi^2 over one period, then hold.  Mode 2 = PLUCK, a
    start-up aid of a real FSI run: the shape phi is the linear static
    deflection of the structure under a tip load (shape 0) or its m-th
    vibration mode in vacuum (shape m, inverse iteration with deflation,
    fsi_shell.shell_mode), normalized to a unit tip deflection along e_n
    (a mode: to a unit LARGEST deflection, the tip moving along +e_n -
    normalized by the tip, the vacuum mode 2 of the benchmark (33.9 Hz)
    bent the rubber into an S of 1.4 x the tip amplitude and the
    morphing columns on its concave side crossed at t = 0.21 s;
    the shape is written to Output_data/<case>_fsi_pluck_shape.dat),
    applied INEXTENSIBLY (its rotation c theta_phi(s) integrated along the
    chord; the node at its offset along the rotated normal), oscillated as
    in mode 0 up to the maximum at tau = 2.25 / f (zero velocity), with
    the imposed acceleration in the structural state; then the structure
    is released and solved.  The reasons for each choice: the mode-0
    shape A xi^2 also bends the steel rear mass, and any linear shape
    stretches the chord by w'^2 / 2 (2 % at 11 mm: 4e9 Pa in the steel -
    the released structure snapped, its tip velocity jumped to 1.7 m/s
    and the structural Newton failed two steps later); the inextensible
    shape has 3e4 Pa at the same deflection.  Modes 0 / 1 keep the linear
    xi^2 shape (the structure is never solved there, its stresses are not
    meaningful).  Debugging aid: PROBCFD_DUMP_GRID=<step> (environment)
    writes the grid of that step in double precision
    (<case>_grid_NNNNNNNN.grd).
  * FILES: Output_data/<case>_fsi_<part>_node<n>.dat (step, t, d_x d_y
    d_z |d| in the part frame, velocity, von Mises, sub-iterations,
    residual; every Data_Intvl steps), <case>_fsi_<part>_shell.grd/.p3d/
    .nam (the mid-surface with Thickness, Submat_id, Clamped; active
    sheets only); the wall surface files gain Disp_x Disp_y Disp_z [m]
    and vonMises [Pa]; the log prints the volume, mass, mid-surface size,
    thickness range, clamped count, Submat regions, shell mass, first
    natural frequency in vacuum, and per step the sub-iterations,
    residual, max |d|, max von Mises, Newton iterations, kinetic / strain
    energy.  Restart: u, v, a, f_ext, previous velocity and the nodal
    triads of every part after the statistics accumulators (an older
    file: the structures start undeformed with a NOTE); the restoring
    grid motion morphs the blocks with the restored displacement before
    the first step.
  * Also fixed on the way: overset.F90 coincident-wall trimming read one
    vertex before the patch (index v0 + m1 - 1 instead of v0 + m1 in the
    boundary-distance test) - a 1-based / 0-based slip present since the
    overset update; affects only the trimming test of overlapping wall
    patches (KP505 hub / blade roots).

36.2 Verification  [verification/update24_fsi, tools/fsi_cases.py, fsi_test.F90]
  (a) Structure alone (make fsi_test): clamped plate with a tip load =
  beam theory to 1.2 % (shear-deformable shell, 20 x 2 elements), root
  von Mises 5.4e7 vs 6.0e7 Pa (the nodal average of the element
  stresses); free-vibration frequency 8.232 vs 8.180 Hz Euler-Bernoulli;
  Turek-Hron CSM3 (rod under gravity 2 m/s2, plane strain, 40 x 2 shell,
  dt 5 ms): d_y mean -63.47 / amplitude 65.11 mm, d_x -14.28 +- 14.28
  mm, f 1.096 Hz against -63.61 +- 65.16 mm, -14.31 +- 14.31 mm, 1.0995 Hz
  (Fabbri et al. 2023 with 70 quadratic elements: -64.16 +- 64.87, 1.097
  Hz) - a 130-mm deflection of a 350-mm rod, 2.9 Newton iterations per
  step.  The large-deflection static cantilever (Mattiasson) stalls at
  P L^2/EI ~ 3 with the simplified tangent (w/L 0.79 vs 0.81 at 10): a
  static limitation only, the dynamic steps take small increments.
  (b) Coupling machinery on the FSI-PfS-2a grid: a flexible part held
  rigid (FSI_Start_T > End_T) goes through morphing / discrete grid
  velocity / flux re-projection with zero deformation and reproduces the
  rigid-child run to 5e-5 relative in the forces of the first steps (the
  rebuilt metrics equal the originals to 4e-13; the same 5e-5 level
  separates the moving-grid path with a zero imposed table from the
  static path: the pressure-solver tolerance).  Restart (6 steps with
  the structure active vs 3 steps + restart, re-assembly every step):
  the first restarted step agrees with the continuous run to 2e-4 in
  the monitor-node displacement and 7e-5 in the plate force - the
  coupling tolerance (the sub-iteration history restarts, so bitwise
  agreement is not expected); the relaxed interface displacement the
  grid was morphed with is stored in the restart file (24b), so that the
  restored grid is the one the fluid fields were saved on.
  Regression: Statistics_ActuatorDisk bitwise, Probes_TowedCylinder
  (moving grids) probe / force / motion / surface files bitwise identical
  to the update-22 delivery.  The passive-patch classification of the
  KP505 tip caps gives five identical blade structures (1050 nodes, 984
  elements, 336 clamped).
  (c) FSI-PfS-2a (cases/FSI_PfS2a, quasi-2-D overset grid of 32 k cells,
  laminar, dt 0.5 ms, released at 0.15 s): rebuilt in 24c on the
  one-body grid (34 k cells, span 0.4 mm, SW_Quasi_2D, FSI_Plane_Strain):
  rigid to 0.15 s at 0.5 ms, then dt = 1 ms with the pluck (the static
  tip-load shape, 11 mm at the tip, 11.25 Hz for 2.25 periods) and the free
  strongly coupled response from 0.35 s to 0.85 s (5 - 6 sub-iterations,
  ~85 s per step on 2 cores).  Over t = 0.60 - 0.85 s: f = 11.7 Hz (zero
  crossings; spectrum peak 11.45 Hz) against 11.25 Hz measured and 11.73 Hz
  (Fabbri et al.); d_y/D amplitude 0.54 (the six extremes 0.46 - 0.65) at the paper's probe (x = 59
  mm) against 0.648 / 0.629, 0.87 (0.78 - 1.01) at the trailing end of the rear
  mass, which moves 5 ms behind the rubber end (the swiveling mode).
  Released from rest (no pluck) the plate of the laminar quasi-2-D model
  oscillated at ~12 Hz with 0.5 mm and did not build up within 0.1 s; with
  the pluck the oscillation is self-sustained and settles at the measured
  frequency.  The remaining amplitude deficit (-17 %) is that of a laminar
  2-D model at Re_D = 30 470 on a 34 k-cell grid.
  (d) KP505_FSI: initialisation verified (mid-surfaces, passive tip caps,
  clamped roots, 2454 Hz first mode for E = 3 GPa); the time stepping of
  the 900 k cells needs a proper machine.  The first momentum step of the
  KP505 cases diverges with the towing-tank ramp (t_ramp = 0.02 of
  cases/KP505: fluid passes through the ramped blade-root walls and the
  implicit upwind of that inflow with the mirrored ghost value is not
  diagonally dominant on 0.1-mm cells - also without the FSI, also for
  the single-body case, which had never been time-stepped in the
  sandbox); cases/KP505_FSI uses t_ramp = 0, with which the momentum step
  passes (open item: keep the walls impermeable during the ramp).
