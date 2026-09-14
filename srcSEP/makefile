SHELL =/bin/sh

DEFAULT_TARGET : lib

# The parser/registry unit target is intentionally dependency-light and can run
# from a source-only archive.  Every production build target still requires the
# enclosing AMPS configuration exactly as before.
ifneq ($(filter test-cli-unit test-state-unit test-geometry-source-unit test-mover-api-unit test-field-line-scope-unit test-transport-common-unit test-parker-unit test-fte-dmumu-unit test-fte-mfp-unit test-coefficients-unit test-turbulence-core-unit test-reproducibility-unit test-acceptance-unit test-documentation-unit test-controlled-analytical test-scientific-validation test-native-amps-validation test-swmf-validation test-observational-validation test-sanitizer test-stochastic-repeat,$(MAKECMDGOALS)),)
-include ../../Makefile.conf
else
include ../../Makefile.conf
endif

MAINLIBOBJ = exosphere_module.o  main_lib.o inner_boundary_injection.o outer_boundary_injection.o mesh.o mover_state.o production_mover_runtime.o transport_common.o coefficient_providers.o parker_mover.o focused_transport_dmumu.o focused_transport_mfp.o parker.o scattering.o sep.o field_line.o flux_tube_geometry.o sampling.o diffusion.o shock_injection.o BoundaryInjection.o BoundaryInjection_GCR.o GCR_Badavi2011ASR.o parser.o diffusion_dxx.o tests.o component_tests.o util/sep_test_registry.o util/sep_acceptance_cases.o util/sep_scientific_validation.o util/sep_mover_validation.o util/sep_turbulence_validation.o util/sep_background_snapshot.o util/sep_background_runtime.o util/sep_flux_tube_geometry_core.o util/sep_production_mover.o util/sep_transport_common.o util/sep_parker_core.o util/sep_focused_transport_core.o util/sep_focused_transport_mfp_core.o util/sep_coefficient_registry.o util/sep_turbulence_core.o util/sep_reproducible_reduction.o shock_analytical_model.o diffusion_dmumu.o DLT.o QLT.o cooling.o QLT1.o sampling_output.o output_fl_background.o shock_analytical_model2D.o solar_wind.o alfven_turbulence.o parker_streaming_calculator.o turbulence_advection_kolmogorov.o wave_energy_initialization.o test_wave_energy_initialization.o wave_particle_coupling_kolmogorov.o growth_rate_validation_test.o kolmogorov.o kolmogorov_scatter.o shock1d_vertex_distance.o sw1d.o turbulence_cascade_kolmogorov.o turbulence_reflection_kolmogorov.o turbulence_wave_number_resolved/turbulence_wave_number_resolved.o
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

.PHONY: test-cli-unit test-state-unit test-geometry-source-unit test-mover-api-unit test-field-line-scope-unit test-transport-common-unit test-parker-unit test-fte-dmumu-unit test-fte-mfp-unit test-coefficients-unit test-turbulence-core-unit test-reproducibility-unit test-acceptance-unit test-documentation-unit test-controlled-analytical test-scientific-validation test-native-amps-validation test-swmf-validation test-observational-validation test-sanitizer test-stochastic-repeat require-sep-executable test-cli test-list test-case \
	test-group test-turbulence test-stress test-mpi test

test-cli-unit:
	./test/run_step1_tests.sh

test-state-unit:
	./test/run_step2_tests.sh

test-geometry-source-unit:
	./test/run_step3_tests.sh

test-mover-api-unit:
	./test/run_step4_tests.sh

test-field-line-scope-unit:
	./test/run_step5_tests.sh

test-transport-common-unit:
	./test/run_step6_tests.sh

test-parker-unit:
	./test/run_step7_tests.sh

test-fte-dmumu-unit:
	./test/run_step8_tests.sh

test-fte-mfp-unit:
	./test/run_step9_tests.sh

test-coefficients-unit:
	./test/run_step10_tests.sh

test-turbulence-core-unit:
	./test/run_step11_tests.sh

test-reproducibility-unit:
	./test/run_step12_tests.sh

test-acceptance-unit:
	./test/run_step13_tests.sh

test-documentation-unit:
	./test/run_step14_tests.sh

# This bounded source-only gate contains the complete controlled analytical
# catalog and its registry/refinement infrastructure.  It deliberately does
# not invoke a linked AMPS executable, ingest SWMF products, or read spacecraft
# observations; those evidence classes have separate targets below.
test-controlled-analytical:
	@set -e; for step in 1 6 7 8 9 11; do \
		./test/run_step$${step}_tests.sh; \
	done

# Step 15 is intentionally opt-in because VAL02 and VAL03 use statistically
# meaningful ensembles.  The target also verifies that absent real SWMF and
# spacecraft evidence produces INCOMPLETE and fails the release-mode gate.
test-scientific-validation:
	./test/run_step15_tests.sh

# Native, SWMF, and observational checks remain distinct so a PASS in one
# class cannot satisfy another.  The external-data targets authenticate all
# referenced files through validation/check_external_gate.py before reporting
# success; unset inputs fail with usage guidance rather than being skipped.
test-native-amps-validation: require-sep-executable
	@"$(SEP_EXECUTABLE)" --all-tests $(if $(JSON),--test-json "$(JSON)",) $(if $(JUNIT),--test-junit "$(JUNIT)",)

test-swmf-validation:
	@test -n "$(SWMF_MANIFEST)" || { echo "ERROR: set SWMF_MANIFEST=/path/to/real-swmf-manifest.json." >&2; exit 2; }
	@PYTHONDONTWRITEBYTECODE=1 python3 validation/check_external_gate.py \
		--evidence-class COUPLED_INTEGRATION "$(SWMF_MANIFEST)"

test-observational-validation:
	@test -n "$(OBSERVATIONAL_MANIFESTS)" || { echo "ERROR: set OBSERVATIONAL_MANIFESTS to one or more held-out event manifests." >&2; exit 2; }
	@PYTHONDONTWRITEBYTECODE=1 python3 validation/check_external_gate.py \
		--evidence-class OBSERVATIONAL_VALIDATION $(foreach manifest,$(OBSERVATIONAL_MANIFESTS),"$(manifest)")

# Runs every dependency-light sanitizer target sequentially.  Each child uses a
# disposable build directory and propagates its nonzero status immediately.
test-sanitizer:
	@set -e; for step in 3 6 7 8 9 10 11 12 13; do \
		./test/run_step$${step}_tests.sh; \
	done

test-stochastic-repeat:
	./test/run_stochastic_repeat_tests.sh

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

# Extended diagnostics never enter --all-tests implicitly.  This target makes
# the opt-in stress campaign visible and allows the same report destinations as
# the routine acceptance target.
test-stress: require-sep-executable
	@"$(SEP_EXECUTABLE)" --test PARKER02 --test SCAT01 $(if $(JSON),--test-json "$(JSON)",) $(if $(JUNIT),--test-junit "$(JUNIT)",)

MPIEXEC ?= mpiexec
MPI_NP ?= 2
test-mpi: require-sep-executable
	@command -v "$(MPIEXEC)" >/dev/null || { echo "ERROR: MPI launcher '$(MPIEXEC)' not found." >&2; exit 2; }
	@"$(MPIEXEC)" -n "$(MPI_NP)" "$(SEP_EXECUTABLE)" --all-tests $(if $(JSON),--test-json "$(JSON)",) $(if $(JUNIT),--test-junit "$(JUNIT)",)

# Keep the top-level sequence explicit and serial even when invoked as
# `make -j test`: component tests share initialized model/output state and must
# not be launched concurrently.  The embedded SWCME suite retains its own
# internal parallel build through its native target.
test: test-cli-unit test-state-unit test-geometry-source-unit test-mover-api-unit test-field-line-scope-unit test-transport-common-unit test-parker-unit test-fte-dmumu-unit test-fte-mfp-unit test-coefficients-unit test-turbulence-core-unit test-reproducibility-unit test-acceptance-unit test-documentation-unit require-sep-executable
	@"$(SEP_EXECUTABLE)" --list-tests >/dev/null
	@"$(SEP_EXECUTABLE)" --all-tests $(if $(JSON),--test-json "$(JSON)",) $(if $(JUNIT),--test-junit "$(JUNIT)",)
	@$(MAKE) -C swcme/test -j test
