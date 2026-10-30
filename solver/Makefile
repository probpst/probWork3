# probCFD - simple Makefile
# (modules first so their .mod files exist for the rest; the order below
#  is the dependency order - keep it when adding files)
# NOTE: this Makefile does not track modification times; after editing a
#       source run  make clean && make
FC   = mpif90
FFLAGS = -O2 -ffree-line-length-none
SRCS = mod_global.F90 mod_mpi.F90 stop_criteria.F90 bcs_b2b.F90 ovs_geom.F90 overset.F90 b2b_exchange.F90 waves.F90 curvi.F90 valves.F90 pid_control.F90 sixdof.F90 force_model.F90 read_input.F90 \
       vof.F90 cavitation.F90 turbulence.F90 energy.F90 surface_tension.F90 boundary.F90 convection.F90 u_momentum.F90 poisson.F90 wave_fields.F90 \
       fsi_shell.F90 morph_frames.F90 fsi.F90 patch_flux.F90 surface_output.F90 probes.F90 output.F90 statistics.F90 restart.F90 Main.F90

probCFD:
	$(FC) $(FFLAGS) -c $(SRCS)
	$(FC) *.o -o probCFD

# debug build with run-time checks (slower)
debug:
	$(FC) -O0 -g -ffree-line-length-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -c $(SRCS)
	$(FC) *.o -o probCFD

# stand-alone checks of the structural solver (update 24)
fsi_test:
	gfortran -O2 -ffree-line-length-none -o fsi_test fsi_shell.F90 fsi_test.F90

clean:
	rm -f *.o *.mod probCFD fsi_test
