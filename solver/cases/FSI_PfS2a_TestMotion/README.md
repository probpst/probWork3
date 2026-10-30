# FSI_PfS2a_TestMotion — the moving-grid fluid alone (update 24c)

The FSI-PfS-2a grid and flow of `cases/FSI_PfS2a` with `FSI_Test_Motion = 0.022, 11.25`:
the structure is NOT solved, its nodes follow an imposed bending mode
d = A r xi^2 sin(2 pi f tau) (A = 22 mm = 1 D at the tip, f = 11.25 Hz, ramped over
one period from `FSI_Start_T = 0.05 s`).  The fluid sees a wall motion of the
benchmark's amplitude and frequency without structural feedback - a test of the
morphing, the discrete grid velocity, the flux re-projection, the overset
re-assembly and the ALE momentum equation.  0.5 ms, to t = 0.175 s (the full
amplitude is reached at t = 0.139 s).

Run: `Mode = 0` from scratch (the sandbox run restarted from the 50-ms restart file
of the benchmark run, which is the same flow).

`sandbox_results/`: monitor and node histories, the log; `fsi_test_motion.png`.
With rigid morphing columns the run diverged at full amplitude (step 317;
`*_rigid_columns_diverged.*`, `Reason_for_stop_rigid_columns.md`: the outer cells
at the plate root of the one-body grid sheared to 6 % of their area); with the
column rotation diffused along the wall in the outer levels (the delivered code)
it passes the full swing with max CFL 6.1 and max |V| 3.5 m/s.
