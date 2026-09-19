SHELL := /bin/sh

DEFAULT_TARGET: lib

# Match srcSEP3D's relocation-safe path discovery. AMPS compiles an
# application after copying it from AMPS/srcSEP to AMPS/build/main, and it may
# invoke this makefile with -f from a third working directory. Therefore no
# shared-model or configuration path may be relative to the process CWD.
SEP_MAKEFILE_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

ifneq ($(wildcard $(SEP_MAKEFILE_DIR)/../Makefile.conf),)
AMPS_ROOT ?= $(abspath $(SEP_MAKEFILE_DIR)/..)
else ifneq ($(wildcard $(SEP_MAKEFILE_DIR)/../../Makefile.conf),)
AMPS_ROOT ?= $(abspath $(SEP_MAKEFILE_DIR)/../..)
else
AMPS_ROOT ?= $(abspath $(SEP_MAKEFILE_DIR)/..)
endif

AMPS_CONFIG ?= $(AMPS_ROOT)/Makefile.conf

# Both shared model directories are resolved from the enclosing AMPS root.
# This is the same source/build-main contract used by srcSEP3D: a normal AMPS
# build needs no SEP_COMMON_DIR or SWCME_DIR environment variable, while a
# deliberately detached component checkout may still override either path.
SEP_COMMON_DIR ?= $(AMPS_ROOT)/src/models/sep_common
SEP_COMMON_ARCHIVE ?= $(SEP_COMMON_DIR)/sep_common.a
SEP_COMMON_OBJECTS := \
    $(SEP_COMMON_DIR)/build/sep_transport_common.o \
    $(SEP_COMMON_DIR)/build/sep_coefficient_physics.o \
    $(SEP_COMMON_DIR)/build/sep_coefficient_registry.o \
    $(SEP_COMMON_DIR)/build/sep_background_snapshot.o \
    $(SEP_COMMON_DIR)/build/sep_test_registry.o \
    $(SEP_COMMON_DIR)/build/sep_injection_spectrum.o \
    $(SEP_COMMON_DIR)/build/sep_species_source.o

# SWCME has one source and binary owner under src/models. srcSEP consumes the
# canonical header-only 1-D interface and inserts the canonical 3-D object into
# mainlib.a, exactly as srcSEP3D does.
SWCME_DIR ?= $(AMPS_ROOT)/src/models/swcme
SWCME_ARCHIVE ?= $(SWCME_DIR)/swcme.a
SWCME_OBJECTS := $(SWCME_DIR)/build/swcme3d.o

# The parser/registry unit target is intentionally dependency-light and can run
# from a source-only archive.  Every production build target still requires the
# enclosing AMPS configuration exactly as before.
STANDALONE_GOALS := test-cli-unit test-state-unit test-geometry-source-unit \
                    test-initialization-unit \
                    test-mover-api-unit test-field-line-scope-unit \
                    test-transport-common-unit test-parker-unit \
                    test-fte-dmumu-unit test-fte-mfp-unit \
                    test-coefficients-unit test-turbulence-core-unit \
                    test-reproducibility-unit test-wp01-wp10-unit \
                    test-wp11-wp20-unit test-wp21-wp30-unit \
                    test-wp31-wp41-unit test-wp42-wp64-experimental \
                    test-wp59-wp64-native test-python-runner-unit \
                    test-cv01-unit test-cv02-cv05-unit test-cv06-cv12-unit \
                    test-iv01-iv06-unit test-xm01-xm03-unit \
                    print-configuration-matrix test-acceptance-unit \
                    test-documentation-unit test-controlled-analytical \
                    test-scientific-validation test-native-amps-validation \
                    test-swmf-validation test-observational-validation \
                    test-sanitizer test-stochastic-repeat \
                    test-swcme-relocation-unit test-sep-common-ownership-unit \
                    test-swcme-fail-closed-unit \
                    test-swcme-configuration-unit \
                    test-d03-native-integration \
                    test-d03-runner-unit \
                    test-package-hygiene-unit \
                    print-layout-paths
ifneq ($(filter $(STANDALONE_GOALS),$(MAKECMDGOALS)),)
-include $(AMPS_CONFIG)
else
include $(AMPS_CONFIG)
endif

# Makefile.conf revisions in supported AMPS trees use one or more of these
# variables in their generic .cpp -> .o rule. Populate all three, as srcSEP3D
# does, so canonical SEP-common and SWCME headers resolve identically in srcSEP
# and in the copied build/main tree.  Application-local duplicate headers are
# deliberately absent, making include-order ambiguity impossible.
MODEL_INCLUDE_FLAGS := -I$(SEP_COMMON_DIR) -I$(SWCME_DIR)
CPPFLAGS += $(MODEL_INCLUDE_FLAGS)
CXXFLAGS += $(MODEL_INCLUDE_FLAGS)
INCLUDE += $(MODEL_INCLUDE_FLAGS)

# Compile only the two SWCME-facing translation units with the canonical model
# directory on their private include path. Historical AMPS Makefile.conf
# generic rules do not consistently consume application-appended CPPFLAGS,
# CXXFLAGS, or INCLUDE. In particular, a copied build/main must not rely on
# those variables to find swcme1d_input.hpp. Keeping this list private also
# leaves sep.h and every unrelated AMPS translation unit independent of SWCME.
SWCME_CONSUMER_CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -O3 \
                           -Iadapters -Iutil -I$(SEP_COMMON_DIR) -I$(SWCME_DIR)

MAINLIBOBJ = exosphere_module.o  main_lib.o inner_boundary_injection.o outer_boundary_injection.o mesh.o mover_state.o production_mover_runtime.o transport_common.o coefficient_providers.o parker_mover.o focused_transport_dmumu.o focused_transport_mfp.o turbulence_production_adapter.o parker.o scattering.o sep.o field_line.o flux_tube_geometry.o sampling.o diffusion.o shock_injection.o BoundaryInjection.o BoundaryInjection_GCR.o GCR_Badavi2011ASR.o parser.o diffusion_dxx.o tests.o component_tests.o validation/cases/CV01/cv01_model.o validation/cases/controlled_transport_models.o validation/cases/advanced_validation_models.o util/sep_acceptance_cases.o util/sep_scientific_validation.o util/sep_mover_validation.o util/sep_swcme_validation.o util/sep_turbulence_validation.o util/sep_background_runtime.o util/sep_flux_tube_geometry_core.o util/sep_sampling_products.o util/sep_shock_source_core.o util/sep_transactional_output.o util/sep_run_configuration.o util/sep_population_control.o util/sep_configuration_matrix.o util/sep_evidence.o util/sep_system_ledger.o util/sep_validation_tools.o util/sep_observation_forward_model.o util/sep_production_mover.o util/sep_parker_core.o util/sep_focused_transport_core.o util/sep_focused_transport_mfp_core.o util/sep_turbulence_core.o util/sep_reproducible_reduction.o util/sep_initialization.o util/sep_initialization_validation.o shock_analytical_model.o diffusion_dmumu.o DLT.o QLT.o cooling.o QLT1.o sampling_output.o output_fl_background.o shock_analytical_model2D.o solar_wind.o alfven_turbulence.o parker_streaming_calculator.o turbulence_advection_kolmogorov.o wave_energy_initialization.o test_wave_energy_initialization.o wave_particle_coupling_kolmogorov.o growth_rate_validation_test.o kolmogorov.o kolmogorov_scatter.o shock1d_vertex_distance.o adapters/swcme1d_adapter.o turbulence_cascade_kolmogorov.o turbulence_reflection_kolmogorov.o turbulence_wave_number_resolved/turbulence_wave_number_resolved.o
# Integrated IV01-IV06 model stages are linked into the same application-level
# registry as CV01-CV12; they are not standalone substitutes for AMPS.
MAINLIBOBJ += validation/cases/integrated_validation_models.o
# XM01-XM03 callbacks must be linked into the production executable so their
# native artifacts cannot be generated by an untracked Python replacement.
MAINLIBOBJ += validation/cases/cross_model_validation_models.o

# sep_cli.o is kept under srcSEP/util to separate command-line parsing from
# the physics driver.  It is archived into main.a because only the standalone
# executable uses these run-time switches; the AMPS library objects remain
# independent of the CLI.
MAINOBJ = main.o util/sep_cli.o debug/sep_debug_fieldline_datum.o

lib: $(SEP_COMMON_ARCHIVE) $(SWCME_ARCHIVE) $(MAINLIBOBJ)
	# Recreate the application archive. Plain `ar -r` would retain a stale
	# shared member after a manifest change and could silently introduce a
	# second shared-kernel or SWCME definition.
	rm -f mainlib.a
	ar -rc mainlib.a $(MAINLIBOBJ) $(SEP_COMMON_OBJECTS) $(SWCME_OBJECTS)
	@members=$$(ar t mainlib.a); \
	for member in $(notdir $(SEP_COMMON_OBJECTS) $(SWCME_OBJECTS)); do \
	  count=$$(printf '%s\n' "$$members" | awk -v expected="$$member" \
	    '$$0 == expected { found++ } END { print found + 0 }'); \
	  test "$$count" -eq 1 || { \
	    echo "ERROR: canonical shared member $$member occurs $$count times" >&2; \
	    exit 1; \
	  }; \
	done
	@echo "Target 'main' is done."

# Phony prerequisites are intentional: the application archive consumes the
# canonical objects rather than nesting .a files, so both build directories
# must exist after an independent clean of either shared model.
.PHONY: sep-common swcme
sep-common:
	$(MAKE) -C $(SEP_COMMON_DIR)

$(SEP_COMMON_ARCHIVE): sep-common

swcme:
	$(MAKE) -C $(SWCME_DIR)

$(SWCME_ARCHIVE): swcme

amps : $(MAINOBJ)
	ar -rc main.a $(MAINOBJ)

adapters/%.o: adapters/%.cpp
	$(CXX) $(SWCME_CONSUMER_CXXFLAGS) -c -o $@ $<

# D01-D03 registry callbacks consume the canonical SWCME input schema just as
# the runtime adapter consumes the canonical model. Give this object an exact
# rule so it cannot fall through to an AMPS generic compile recipe that omits
# the application include flags. The rule works unchanged in srcSEP and in the
# copied AMPS/build/main tree because SWCME_DIR is AMPS-root-relative.
util/sep_swcme_validation.o: util/sep_swcme_validation.cpp
	$(CXX) $(SWCME_CONSUMER_CXXFLAGS) -c -o $@ $<

# As in srcSEP3D, a production verification must enter the enclosing AMPS
# build. Compiling AMPS-dependent sources directly in srcSEP would omit the
# generated build/pic headers and configuration assembled by the parent build.
.PHONY: strict-production audit-production-shared-models audit-production-swcme
strict-production:
	@$(MAKE) -C "$(AMPS_ROOT)" amps
	@$(MAKE) --no-print-directory -f "$(SEP_MAKEFILE_DIR)/makefile" \
	  audit-production-shared-models AMPS_ROOT="$(AMPS_ROOT)"

PRODUCTION_ARCHIVE_DIR ?= $(AMPS_ROOT)/build/main
PRODUCTION_MAINLIB := $(PRODUCTION_ARCHIVE_DIR)/mainlib.a

audit-production-shared-models:
	@command -v ar >/dev/null 2>&1 || { \
	  echo "ERROR: ar is required for the production archive audit." >&2; \
	  exit 2; \
	}
	@test -f "$(PRODUCTION_MAINLIB)" || { \
	  echo "ERROR: enclosing AMPS build did not produce $(PRODUCTION_MAINLIB)" >&2; \
	  exit 1; \
	}
	@members=$$(ar t "$(PRODUCTION_MAINLIB)"); \
	for member in $(notdir $(SEP_COMMON_OBJECTS) $(SWCME_OBJECTS)); do \
	  count=$$(printf '%s\n' "$$members" | awk -v expected="$$member" \
	    '$$0 == expected { found++ } END { print found + 0 }'); \
	  test "$$count" -eq 1 || { \
	    echo "ERROR: canonical shared member $$member occurs $$count times in $(PRODUCTION_MAINLIB)" >&2; \
	    exit 1; \
	  }; \
	done
	@duplicates=$$(nm -g --defined-only "$(PRODUCTION_MAINLIB)" | \
	  awk '$$2 ~ /^[TDBRS]$$/ { count[$$3]++ } END { \
	    for (symbol in count) if (count[symbol] > 1) print symbol }'); \
	test -z "$$duplicates" || { \
	  echo "ERROR: duplicate strong definitions in $(PRODUCTION_MAINLIB):" >&2; \
	  printf '%s\n' "$$duplicates" >&2; exit 1; \
	}
	@echo "[srcSEP] enclosing AMPS build contains one canonical copy of every shared member"

# Compatibility target retained for existing build scripts.  Its stronger B02
# replacement audits both sep_common and SWCME, so the old name cannot weaken
# release evidence.
audit-production-swcme: audit-production-shared-models

# The linked executable is produced by the enclosing AMPS application build,
# not by this archive-only subdirectory makefile.  Override SEP_EXECUTABLE when
# that build uses a different path.  Every integration target below invokes the
# real production CLI; selector semantics are never reimplemented in Make.
SEP_EXECUTABLE ?= ../amps

.PHONY: test-cli-unit test-initialization-unit test-state-unit test-geometry-source-unit test-mover-api-unit test-field-line-scope-unit test-transport-common-unit test-parker-unit test-fte-dmumu-unit test-fte-mfp-unit test-coefficients-unit test-turbulence-core-unit test-reproducibility-unit test-wp01-wp10-unit test-wp11-wp20-unit test-wp21-wp30-unit test-wp31-wp41-unit test-wp42-wp64-experimental test-wp59-wp64-native test-python-runner-unit test-cv01-unit test-cv02-cv05-unit test-cv06-cv12-unit test-iv01-iv06-unit test-xm01-xm03-unit print-configuration-matrix test-acceptance-unit test-documentation-unit test-controlled-analytical test-scientific-validation test-native-amps-validation test-swmf-validation test-observational-validation test-sanitizer test-stochastic-repeat test-swcme-relocation-unit test-sep-common-ownership-unit test-swcme-fail-closed-unit test-swcme-configuration-unit test-d03-native-integration test-d03-runner-unit test-package-hygiene-unit require-sep-executable test-cli test-list test-case \
	test-group test-turbulence test-stress test-mpi test

# B01 is a release-tree check, not a srcSEP dependency.  Resolve the tool from
# the enclosing AMPS root so the same target works when this makefile is copied
# to build/main; the tool itself never imports either application.
test-package-hygiene-unit:
	@PYTHONDONTWRITEBYTECODE=1 python3 \
	  "$(AMPS_ROOT)/tools/sep_package_hygiene.py" \
	  --root "$(AMPS_ROOT)" --self-test

# Fast structural test for source-tree/build-main path equivalence and the
# absence of application-local SWCME headers or sources.
test-swcme-relocation-unit:
	./test/run_swcme_relocation_tests.sh

test-swcme-fail-closed-unit:
	./test/run_d01_tests.sh

test-swcme-configuration-unit:
	./test/run_d02_tests.sh

# D03PRE first compiles and executes the exact extended descriptor registered
# in the linked application. The Python gate then owns external build/run
# policy. Environment-backed arguments let test/run_tests.py retain the same
# outer JSON/JUnit wrapper as every source suite while direct Make users can
# set the variables explicitly.
test-d03-runner-unit:
	./test/run_d03_preflight_tests.sh
	@PYTHONDONTWRITEBYTECODE=1 python3 test/test_native_integration_runner.py

test-d03-native-integration: test-d03-runner-unit
	@PYTHONDONTWRITEBYTECODE=1 python3 test/run_native_integration.py

test-sep-common-ownership-unit:
	./test/run_sep_common_ownership_tests.sh

test-cli-unit:
	./test/run_step1_tests.sh

test-initialization-unit:
	./test/run_initialization_tests.sh

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

test-wp01-wp10-unit:
	./test/run_wp01_wp10_tests.sh

# WP11-WP20 combine numerical boundary/refinement checks with pure coefficient,
# species-source, and compatibility tests.  The runner also verifies that the
# repaired pure kernels are reached by production adapter source paths.
test-wp11-wp20-unit:
	./test/run_wp11_wp20_tests.sh

test-wp21-wp30-unit:
	./test/run_wp21_wp30_tests.sh

# WP31-WP41 provide bounded contracts for operator limits, population
# invariants, compatibility, ledgers, statistics, fuzzing, forward modeling,
# performance, and evidence governance. External evidence stays separate.
test-wp31-wp41-unit:
	./test/run_wp31_wp41_tests.sh

# B05 keeps these prototypes out of MAINLIBOBJ until their disposition manifest
# names and satisfies a native promotion gate.  The suite still compiles every
# implementation with sanitizers and verifies all 23 component contracts.
test-wp42-wp64-experimental:
	./test/run_wp42_wp64_tests.sh

# Site-owned linked/MPI/replay/scaling evidence is intentionally separate from
# the dependency-light prototype suite.  With no command it emits the exact
# runner SKIP marker; a reviewed workflow supplies SRCSEP_NATIVE_GATE.
test-wp59-wp64-native:
	./test/run_wp59_wp64_native_gates.sh

# The runner test uses a temporary synthetic registry report to exercise both
# metric-level and pointwise analytical plots.  It neither links AMPS nor
# writes generated figures into the source tree.
test-python-runner-unit:
	@PYTHONDONTWRITEBYTECODE=1 python3 test/test_python_test_runner.py

test-xm01-xm03-unit:
	./test/run_xm01_xm03_tests.sh

# CV01 establishes the reusable linked-application case pattern. The target
# always compiles the native model source strictly and runs the full case only
# when SEP_EXECUTABLE names the linked srcSEP/AMPS binary; it never substitutes
# a standalone driver for missing application evidence.
test-cv01-unit:
	./test/run_cv01_tests.sh

# CV02-CV05 extend the same linked-application contract to diffusion,
# nonuniform stochastic drift, adiabatic cooling, and magnetic focusing.
test-cv02-cv05-unit:
	./test/run_cv02_cv05_tests.sh

# CV06-CV12 add modal scattering, telegraph/first-passage transport, DSA,
# spectral turbulence evolution, and a closed particle-wave energy ledger.
test-cv06-cv12-unit:
	./test/run_cv06_cv12_tests.sh

# Integrated manufactured verification is a distinct gate because these cases
# couple geometry/operators and remain required nightly in the campaign plan.
test-iv01-iv06-unit:
	./test/run_iv01_iv06_tests.sh

print-configuration-matrix:
	./test/print_configuration_matrix.sh

test-acceptance-unit:
	./test/run_step13_tests.sh

test-documentation-unit:
	./test/run_step14_tests.sh

# This bounded gate contains the dependency-light controlled analytical catalog
# and its registry/refinement infrastructure. CV01's compile contract always
# runs; its application-level comparison runs only when SEP_EXECUTABLE is
# available. SWMF and observational evidence remain separate targets below.
test-controlled-analytical:
	@set -e; for step in 1 6 7 8 9 11; do \
		./test/run_step$${step}_tests.sh; \
	done
	./test/run_cv01_tests.sh
	./test/run_cv02_cv05_tests.sh
	./test/run_cv06_cv12_tests.sh
	./test/run_iv01_iv06_tests.sh

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
	done; \
	./test/run_wp11_wp20_tests.sh
	./test/run_wp21_wp30_tests.sh
	./test/run_wp31_wp41_tests.sh

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
test: test-swcme-relocation-unit test-cli-unit test-initialization-unit test-state-unit test-geometry-source-unit test-mover-api-unit test-field-line-scope-unit test-transport-common-unit test-parker-unit test-fte-dmumu-unit test-fte-mfp-unit test-coefficients-unit test-turbulence-core-unit test-reproducibility-unit test-wp01-wp10-unit test-wp11-wp20-unit test-wp21-wp30-unit test-wp31-wp41-unit test-python-runner-unit test-cv01-unit test-acceptance-unit test-documentation-unit require-sep-executable
	@"$(SEP_EXECUTABLE)" --list-tests >/dev/null
	@"$(SEP_EXECUTABLE)" --all-tests $(if $(JSON),--test-json "$(JSON)",) $(if $(JUNIT),--test-junit "$(JUNIT)",)
	@PYTHONDONTWRITEBYTECODE=1 python3 "$(SWCME_DIR)/test/run_tests.py" \
	  --routine --output-dir "$(abspath test_output/swcme-r1)"

# Diagnostic target copied from srcSEP3D's relocation contract. It is used by
# the focused test and is also useful when diagnosing an installed AMPS tree.
.PHONY: print-layout-paths
print-layout-paths:
	@echo "SEP_MAKEFILE_DIR=$(SEP_MAKEFILE_DIR)"
	@echo "AMPS_ROOT=$(AMPS_ROOT)"
	@echo "AMPS_CONFIG=$(AMPS_CONFIG)"
	@echo "SEP_COMMON_DIR=$(SEP_COMMON_DIR)"
	@echo "SWCME_DIR=$(SWCME_DIR)"
	@echo "MODEL_INCLUDE_FLAGS=$(MODEL_INCLUDE_FLAGS)"
