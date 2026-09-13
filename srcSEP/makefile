SHELL =/bin/sh

DEFAULT_TARGET : lib  

# The parser/registry unit target is intentionally dependency-light and can run
# from a source-only archive.  Every production build target still requires the
# enclosing AMPS configuration exactly as before.
ifeq ($(MAKECMDGOALS),test-cli-unit)
-include ../../Makefile.conf
else
include ../../Makefile.conf
endif

MAINLIBOBJ = exosphere_module.o  main_lib.o inner_boundary_injection.o outer_boundary_injection.o mesh.o mover.o parker.o scattering.o drift.o sep.o field_line.o sampling.o diffusion.o output.o shock_injection.o BoundaryInjection.o BoundaryInjection_GCR.o GCR_Badavi2011ASR.o parser.o diffusion_dxx.o tests.o component_tests.o util/sep_test_registry.o shock_analytical_model.o diffusion_dmumu.o DLT.o QLT.o cooling.o QLT1.o sampling_output.o output_fl_background.o shock_analytical_model2D.o sample3d_init.o sample3d_sampling.o sample3d_output.o solar_wind.o   alfven_turbulence.o parker_streaming_calculator.o turbulence_advection_kolmogorov.o wave_energy_initialization.o test_wave_energy_initialization.o wave_particle_coupling_kolmogorov.o fte_mover.o growth_rate_validation_test.o kolmogorov.o kolmogorov_scatter.o fte_mover_dmumu.o  shock1d_vertex_distance.o sw1d.o turbulence_cascade_kolmogorov.o turbulence_reflection_kolmogorov.o turbulence_wave_number_resolved/turbulence_wave_number_resolved.o
# sep_cli.o is kept under srcSEP/util to separate command-line parsing from
# the physics driver.  It is archived into main.a because only the standalone
# executable uses these run-time switches; the AMPS library objects remain
# independent of the CLI.
MAINOBJ = main.o util/sep_cli.o debug/sep_debug_fieldline_datum.o  

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)

# The linked executable is produced by the enclosing AMPS application build,
# not by this archive-only subdirectory makefile.  Override SEP_EXECUTABLE when
# that build uses a different path.  Every integration target below invokes the
# real production CLI; selector semantics are never reimplemented in Make.
SEP_EXECUTABLE ?= ../amps

.PHONY: test-cli-unit require-sep-executable test-cli test-list test-case \
	test-group test-turbulence test

test-cli-unit:
	./test/run_step1_tests.sh

require-sep-executable:
	@test -x "$(SEP_EXECUTABLE)" || { \
		echo "ERROR: linked srcSEP executable not found at $(SEP_EXECUTABLE)." >&2; \
		echo "       Build the enclosing AMPS application or set SEP_EXECUTABLE=/path/to/executable." >&2; \
		exit 2; \
	}

test-cli: test-cli-unit require-sep-executable
	@"$(SEP_EXECUTABLE)" --help >/dev/null
	@"$(SEP_EXECUTABLE)" --list-tests >/dev/null

test-list: require-sep-executable
	@"$(SEP_EXECUTABLE)" --list-tests

test-case: require-sep-executable
	@test -n "$(CASE)" || { echo "ERROR: use make test-case CASE=<stable-ID>." >&2; exit 2; }
	@"$(SEP_EXECUTABLE)" --test "$(CASE)"

test-group: require-sep-executable
	@test -n "$(GROUP)" || { echo "ERROR: use make test-group GROUP=<name>." >&2; exit 2; }
	@"$(SEP_EXECUTABLE)" --test-group "$(GROUP)"

test-turbulence: require-sep-executable
	@"$(SEP_EXECUTABLE)" --test-group turbulence

# Keep the top-level sequence explicit and serial even when invoked as
# `make -j test`: component tests share initialized model/output state and must
# not be launched concurrently.  The embedded SWCME suite retains its own
# internal parallel build through its native target.
test: test-cli-unit require-sep-executable
	@"$(SEP_EXECUTABLE)" --list-tests >/dev/null
	@"$(SEP_EXECUTABLE)" --all-tests
	@$(MAKE) -C swcme/test -j test
