# Verification of update 24 — fluid-structure interaction

Cases written by `tools/fsi_cases.py` (`pfs2a` → `cases/FSI_PfS2a`, `kp505` →
`cases/KP505_FSI`); their sandbox results are in the case folders
(`sandbox_results/`).  `plot_fsi.py` evaluates a case (`tip`: monitor-node
histories, spectrum, mean / amplitude / frequency of the last periods; `grid`:
the morphed overset grid of a step; `shell`: the mid-surface with the measured
thickness and the clamped / steel nodes).  `make fsi_test` builds the
stand-alone structural tests of `fsi_shell.F90` (`fsi_test cant | free | large
| csm3 [dt]`); their output is in `csm3/`.

## 1. Structure alone (`fsi_test`, folder `csm3/`)

| test | result | reference |
|---|---|---|
| clamped plate, tip load (20 x 2 shell) | tip deflection 19.76 mm, root von Mises 5.40e7 Pa | beam theory P L³/3EI = 20.0 mm (0.988), 6 P L/(W t²) = 6.0e7 Pa |
| free vibration of the plate | 8.232 Hz | Euler-Bernoulli 8.180 Hz (0.6 %) |
| Turek-Hron CSM3 (40 x 2 shell, plane strain, gravity 2 m/s², dt 5 ms, 10 s) | d_y mean −63.47 mm, amplitude 65.11 mm; d_x −14.28 ± 14.28 mm; f = 1.096 Hz | Turek & Hron: d_y −63.61 ± 65.16 mm, d_x −14.31 ± 14.31 mm, f = 1.0995 Hz |
| CSM3 with dt 2.5 ms | d_y −64.06 ± 64.93 mm, f = 1.096 Hz | (Fabbri et al. 2023, 70 elements: −64.16 ± 64.87 mm, 1.097 Hz) |

The CSM3 rod deflects 130 mm on 350 mm of length (rotations of 40° at the
tip): the co-rotational shell reproduces the large-deformation reference to
0.2 % in the mean, 0.1 % in the amplitude and 0.3 % in the frequency with 2.9
Newton iterations per step (`csm3_tip.png`).  The large-deflection static test
(Mattiasson's cantilever, P L²/EI up to 10, `fsi_test large`) converges to
w/L = 0.79 / u/L = 0.49 against 0.81 / 0.56: the simplified co-rotational
tangent stalls beyond P L²/EI ≈ 3 in a static solve (the dynamic steps of an
FSI run take hundreds of small increments per period and are not affected).

## 2. FSI-PfS-2a (`cases/FSI_PfS2a`)

The benchmark of Kalmbach & Breuer (2013) in the set-up of Fabbri et al.
(2023): cylinder D = 22 mm, para-rubber plate 50 x 2 mm (E 4.1 MPa, nu 0.48,
rho 1090) with a steel rear mass 10 x 2 mm, u = 1.385 m/s, Re_D 30 470;
quasi-2-D overset grid of 34 k cells (one elliptic O-grid around cylinder +
plate, a refinement box and the background; span 0.4 mm, 2 cells; the paper:
40 M cells, LES), laminar with `SW_Quasi_2D = 3`, `FSI_Plane_Strain = 3`.
Rigid to t = 0.15 s at 0.5 ms; from there **Del_T = 1 ms**: the pluck
(`FSI_Test_Motion = 0.011, 11.25, 2, 0`: the static tip-load shape, 11 mm at
the tip, imposed at 11.25 Hz for 2.25 periods) and from t = 0.35 s the free,
strongly coupled response (5 - 6 sub-iterations per step) to t = 0.85 s
(`fsi_pfs2a_result.png`, `cases/FSI_PfS2a/sandbox_results`).

| t = 0.60 - 0.85 s (3 periods) | probCFD | experiment (Kalmbach & Breuer) | Fabbri et al. 2023 (LES) |
|---|---|---|---|
| frequency (zero crossings / spectrum peak) | **11.7 Hz** / 11.45 Hz | 11.25 Hz | 11.73 Hz |
| d_y/D amplitude, node 1 (upper wall, x = 59 mm) | **0.54** (the six extremes 0.46 - 0.65) | 0.648 | 0.629 |
| d_y/D amplitude, trailing end of the rear mass | 0.87 (0.78 - 1.01) | - | - |
| mean d_y/D | -0.07 / -0.09 | ~0 | ~0 |

* Released at a finite amplitude, the oscillation is self-sustained: after a
  transient of ~0.25 s with a slower swing (6 - 9 Hz) it locks at 11.7 Hz -
  4 % above the experiment and equal to the LES of Fabbri et al.; the
  trailing end moves 1.6 x as far as the rubber end and 5 ms (20 deg) later:
  the rear mass swivels (the experiment's second swiveling mode).
* The amplitude at the rubber end is 17 % below the experiment; its extremes still
  varied by +-17 % from half-period to half-period at the end of the run.  The model is
  laminar and two-dimensional at Re_D = 30 470 (the experiment and the LES are
  turbulent and 3-D) and the grid has 34 k cells; the start needed the pluck
  (with the rigid splitter plate the laminar wake hardly sheds - lift ~1e-4 N -
  and released from rest at t = 0.15 s the plate oscillated at ~12 Hz with
  0.5 mm, neither growing nor decaying within 0.1 s).
* Breakdowns found and fixed on the way (all in PROJECT_RECORD 36.1): the
  old metrics of the flux re-projection; the shear of the outer cells at the
  plate root and their crushing on the concave side of the swiveling plate
  (column rotation and displacement diffused along the wall in the outer
  levels); the predictor with the generalized-alpha velocity at the release;
  the second-mode pluck normalized by the tip (S-bend, crossing columns); the
  span drift of the structure (span groups).

## 3. Consistency checks of the coupling machinery

* A flexible part held rigid (`FSI_Start_T` beyond `End_T`) runs through the
  morphing / grid-velocity / flux re-projection path with zero deformation: the
  rebuilt metrics equal the original ones to 4e-13 (max |dmh| / max |mh|), and
  the forces of the first 5 steps agree with the rigid-child run to 5e-5
  relative (the pressure solver's tolerance; the same level at which the
  moving-grid path without motion, an imposed zero table, differs from the
  static path).
* Imposed bending mode on the one-body grid (`FSI_Test_Motion = 0.022, 11.25`,
  from the developed rigid flow at t = 0.05 s, 0.5 ms; `fsi_test_motion.png`):
  the fluid alone on a grid bent by 1 D at the tip at the benchmark frequency.
  With rigid morphing columns the run diverged at full amplitude (step 317:
  the outer cells between the clamped columns of the cylinder and the first
  rotating columns of the plate root sheared to 6 % of their area); with the
  column rotation diffused along the wall in the outer levels (24c) it runs
  through the full swing (max CFL 6.1, max |V| of the active cells 3.5 m/s,
  2.5 U).  Before the 24c fix of the old metrics of the flux re-projection
  the same test (and a static 15 mm deflection held still) blew up at the
  trailing end; that fix (move_grids no longer overwrites the level-n metrics
  of the morphing blocks) was verified with the hold test: morphing path and
  static path from the same restart agreed step for step.
* Hole cutting with slab end caps: 0 active cells inside the body of the
  0.4-mm-span grid (148 before); the node iblank of the update-24 code and of
  the new code are identical for Statistics_TowedCylinder,
  Probes_TowedCylinder, Probes_SteeredFoil, Cylinder_MovingBGR and
  Foil_Spring (bodies now closed by 72 - 682 end-cap triangles) and for
  WaveTank_MovingBGR (no slab body).
* Restart with the structure active (6 steps continuous vs 3 steps + restart,
  re-assembly every step): the first restarted step agrees to 2·10⁻⁴ in the
  monitor-node displacement and 7·10⁻⁵ in the plate force (the coupling
  tolerance `FSI_Tol`; the sub-iteration history starts afresh).  The
  interface displacement of the grid is in the restart file (24b), so the
  restored grid is the one the fluid fields were saved on.
* Regression: `cases/Statistics_ActuatorDisk` (static grids, force model,
  statistics) is bitwise identical to the update-23 delivery; the towed
  cylinder of `cases/Probes_TowedCylinder` (moving grids, probes, surface
  files): probe, force, moment and motion files and the wall surface files of step 20 bitwise identical to the update-22 delivery.
* Non-sheet-like / overlapping wall patches: the tip-cap blocks of the KP505
  blades (35 x 46 nodes, 72 % of them within one node spacing of the root
  block's wall) are classified passive; the five blades give identical
  structures (1050 nodes, 984 elements, 336 clamped, f₁ = 2454 Hz in vacuum
  for E = 3 GPa).
