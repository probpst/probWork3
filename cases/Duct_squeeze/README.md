# Duct_squeeze — verification of the prescribed grid deformation (morph_frames.F90, H1)

Closed-end square duct (2 blocks, 16,000 cells) whose side walls squeeze periodically:
y' = y f(x,t), z' = z f(x,t), f = 1 - A sin²(πx/L) sin²(πt/T) with A = 0.3, T = 0.5 s, L = 0.1 m, h0 = 0.02 m.
The outlet face (x = L) is a rigid pressure boundary p = 0; the fluid pushed out must equal -dV/dt.

Regenerate everything (grid, bcs, csv, frame file, namelist) with
    python3 tools/make_duct_case.py cases/Duct_squeeze
(the 11 MB frame file `Duct_frames.bin` is not committed).  Run:
    Mode = 2 in probCFD.nml -> grid-movement check (Check_grid_movement/), then Mode = 0:
    mpirun -np 2 ./probCFD > log.txt 2>&1 < /dev/null
    python3 tools/check_duct.py cases/Duct_squeeze
check_duct.py compares the grid volume of every written grid with V(t), the outflow through the open
end with -dV/dt, and prints the monitor extremes (CFL, |div|max, Poisson iterations).
