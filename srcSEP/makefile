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

# The shared parallel-diffusion coefficient library (revision 1.4) is built by
# its own makefile.  Its two objects are flattened into mainlib.a exactly like
# the sep_common kernels, so srcSEP and srcSEP3D link the same canonical
# objects.  srcSEP reaches it only through adapters/parallel_diffusion_adapter.
PARALLEL_DIFFUSION_DIR ?= $(SEP_COMMON_DIR)/parallel_diffusion
PARALLEL_DIFFUSION_ARCHIVE ?= $(PARALLEL_DIFFUSION_DIR)/libparallel_diffusion.a
PARALLEL_DIFFUSION_OBJECTS := \
    $(PARALLEL_DIFFUSION_DIR)/build/parallel_diffusion.o \
    $(PARALLEL_DIFFUSION_DIR)/build/parallel_diffusion_advanced.o

# SWCME has one source and binary owner under src/models. srcSEP consumes the
# canonical header-only 1-D interface and inserts the canonical 3-D object into
# mainlib.a, exactly as srcSEP3D does.
SWCME_DIR ?= $(AMPS_ROOT)/src/models/swcme
SWCME_ARCHIVE ?= $(SWCME_DIR)/swcme.a
SWCME_OBJECTS := $(SWCME_DIR)/build/swcme3d.o

# The reduced srcSEP adapter consumes the same canonical coronal ambient and
# corona/SWCME composite as srcSEP3D.  AMPS links application archives but does
# not recursively extract nested archives, so enumerate the maintained model
# objects and flatten exactly one copy into mainlib.a.
SEP_CORONAL_CME_DIR ?= $(AMPS_ROOT)/src/models/sep_coronal_cme
SEP_CORONAL_CME_ARCHIVE ?= $(SEP_CORONAL_CME_DIR)/build/libsep_coronal_cme.a
SEP_CORONAL_CME_OBJECTS := $(patsubst \
    $(SEP_CORONAL_CME_DIR)/src/%.cpp,$(SEP_CORONAL_CME_DIR)/build/%.o,\
    $(wildcard $(SEP_CORONAL_CME_DIR)/src/*.cpp)) \
    $(SEP_CORONAL_CME_DIR)/build/common_sep_field_line_exchange.o \
    $(SEP_CORONAL_CME_DIR)/build/common_sep_field_line_bundle_io.o \
    $(SEP_CORONAL_CME_DIR)/build/common_sep_coherent_transport.o
SEP_CORONA_SWCME_DIR ?= $(AMPS_ROOT)/src/models/sep_corona_swcme
SEP_CORONA_SWCME_ARCHIVE ?= $(SEP_CORONA_SWCME_DIR)/build/libsep_corona_swcme.a
SEP_CORONA_SWCME_OBJECTS := $(patsubst \
    $(SEP_CORONA_SWCME_DIR)/src/%.cpp,$(SEP_CORONA_SWCME_DIR)/build/%.o,\
    $(wildcard $(SEP_CORONA_SWCME_DIR)/src/*.cpp)) \
    $(patsubst \
    $(SEP_CORONA_SWCME_DIR)/shock_front/%.cpp,$(SEP_CORONA_SWCME_DIR)/build/shock_front_%.o,\
    $(wildcard $(SEP_CORONA_SWCME_DIR)/shock_front/*.cpp))

# The parser/registry unit target is intentionally dependency-light and can run
# from a source-only archive.  Every production build target still requires the
# enclosing AMPS configuration exactly as before.
STANDALONE_GOALS := test-cli-unit test-state-unit test-geometry-source-unit \
                    test-backend-admission-unit test-backend-admission-native test-material-background-unit \
                    test-standalone-application-runner-unit test-standalone-application-native \
                    test-standalone-integration-native \
                    test-material-checkpoint-native test-material-checkpoint-unit \
                    test-stage7-checkpoint-native test-stage7-checkpoint-unit \
                    test-standalone-integration-application \
                    test-rng-lifecycle-unit \
                    test-background-derivatives-unit \
                    test-analytic-radial-background-unit \
                    test-analytic-shock-epoch-unit \
                    test-wave-normalization-unit \
                    test-real-wave-normalization-unit \
                    test-polarized-wave-state-unit \
                    test-birth-window-unit \
                    test-warm-parallel-study-unit \
                    test-finite-amplitude-study-unit \
                    test-frozen-collision-transaction-unit \
                    test-flux-tube-state-unit \
                    test-polarized-wkb-unit \
                    test-maxwell-juttner-flux-unit \
                    test-maxwell-juttner-volume-unit \
                    test-material-tube-kinematics-unit \
                    test-standard-batsrus-wave-unit \
                    test-conservation-resampling-unit \
                    test-resonant-collision-unit \
                    test-initialization-unit \
                    test-parallel-diffusion-binding-unit \
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
                    test-ov01-ov05-unit test-ev01-ev02-unit \
                    print-configuration-matrix test-acceptance-unit \
                    test-documentation-unit test-controlled-analytical \
                    test-scientific-validation test-native-amps-validation \
                    test-swmf-validation test-observational-validation \
                    test-sanitizer test-stochastic-repeat \
                    test-swcme-relocation-unit test-sep-common-ownership-unit \
                    test-swcme-fail-closed-unit \
                    test-swcme-configuration-unit \
                    test-reduced-shock-background-unit \
                    test-standalone-integration-unit \
                    test-d03-native-integration \
                    test-d03-runner-unit \
                    test-stage1-stage2-contracts \
                    test-stage3-contracts \
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
# $(AMPS_ROOT)/src/general/input holds the header-only generic runtime-input
# reader (amps_input.h, src/INPUT_FILE.md) used by util/sep_initialization.cpp
# and by the shared models' event readers.
MODEL_INCLUDE_FLAGS := -I$(SEP_COMMON_DIR) -I$(SWCME_DIR) \
                       -I$(SEP_CORONAL_CME_DIR)/include \
                       -I$(SEP_CORONA_SWCME_DIR)/include \
                       -I$(SEP_CORONA_SWCME_DIR)/shock_front \
                       -I$(AMPS_ROOT)/src/general/input
CPPFLAGS += $(MODEL_INCLUDE_FLAGS)
CXXFLAGS += $(MODEL_INCLUDE_FLAGS)
INCLUDE += $(MODEL_INCLUDE_FLAGS)

# The enclosing AMPS suffix rule consumes FLAGC_EXTRA, rather than CXXFLAGS.
# Track every transitive header: changing a stored physics struct must rebuild
# all consumers, otherwise a successfully linked archive can contain mixed ABIs.
FLAGC_EXTRA += -MMD -MP

# Compile only the two SWCME-facing translation units with the canonical model
# directory on their private include path. Historical AMPS Makefile.conf
# generic rules do not consistently consume application-appended CPPFLAGS,
# CXXFLAGS, or INCLUDE. In particular, a copied build/main must not rely on
# those variables to find swcme1d_input.hpp. Keeping this list private also
# leaves sep.h and every unrelated AMPS translation unit independent of SWCME.
SWCME_CONSUMER_CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -O3 \
                           -MMD -MP -Iadapters -Iutil $(MODEL_INCLUDE_FLAGS)

MAINLIBOBJ = exosphere_module.o  main_lib.o inner_boundary_injection.o outer_boundary_injection.o mesh.o mover_state.o production_mover_runtime.o transport_common.o coefficient_providers.o parker_mover.o focused_transport_dmumu.o focused_transport_mfp.o turbulence_production_adapter.o parker.o scattering.o sep.o field_line.o flux_tube_geometry.o sampling.o diffusion.o shock_injection.o BoundaryInjection.o BoundaryInjection_GCR.o GCR_Badavi2011ASR.o parser.o diffusion_dxx.o tests.o component_tests.o validation/cases/CV01/cv01_model.o validation/cases/controlled_transport_models.o validation/cases/advanced_validation_models.o util/sep_acceptance_cases.o util/sep_scientific_validation.o util/sep_mover_validation.o util/sep_swcme_validation.o util/sep_turbulence_validation.o util/sep_background_runtime.o util/sep_flux_tube_geometry_core.o util/sep_sampling_products.o util/sep_shock_source_core.o util/sep_transactional_output.o util/sep_run_configuration.o util/sep_population_control.o util/sep_configuration_matrix.o util/sep_evidence.o util/sep_system_ledger.o util/sep_validation_tools.o util/sep_observation_forward_model.o util/sep_production_mover.o util/sep_parker_core.o util/sep_focused_transport_core.o util/sep_focused_transport_mfp_core.o util/sep_turbulence_core.o util/sep_reproducible_reduction.o util/sep_initialization.o util/sep_initialization_validation.o shock_analytical_model.o diffusion_dmumu.o DLT.o QLT.o cooling.o QLT1.o sampling_output.o output_fl_background.o shock_analytical_model2D.o solar_wind.o alfven_turbulence.o parker_streaming_calculator.o turbulence_advection_kolmogorov.o wave_energy_initialization.o test_wave_energy_initialization.o wave_particle_coupling_kolmogorov.o growth_rate_validation_test.o kolmogorov.o kolmogorov_scatter.o shock1d_vertex_distance.o adapters/swcme1d_adapter.o adapters/reduced_shock_background_adapter.o adapters/parallel_diffusion_adapter.o turbulence_cascade_kolmogorov.o turbulence_reflection_kolmogorov.o turbulence_wave_number_resolved/turbulence_wave_number_resolved.o
# Integrated IV01-IV06 model stages are linked into the same application-level
# registry as CV01-CV12; they are not standalone substitutes for AMPS.
MAINLIBOBJ += validation/cases/integrated_validation_models.o
# The native registry and the focused runner share the stochastic lifecycle
# descriptors, including complete RNG/cache checkpoint and repeated-call tests.
MAINLIBOBJ += util/sep_rng_validation.o
MAINLIBOBJ += util/sep_background_derivative_validation.o
MAINLIBOBJ += util/sep_analytic_radial_native_validation.o
MAINLIBOBJ += util/sep_analytic_radial_background.o util/sep_analytic_radial_validation.o
MAINLIBOBJ += reduced_shock_geometry.o util/sep_reduced_geometry.o
MAINLIBOBJ += reduced_shock_particles.o util/sep_reduced_source.o util/sep_absorbing_bridge.o
MAINLIBOBJ += reduced_shock_boundary.o
MAINLIBOBJ += util/sep_host_step_native_validation.o
MAINLIBOBJ += util/sep_analytic_shock_epoch.o util/sep_analytic_shock_epoch_validation.o
MAINLIBOBJ += util/sep_coefficient_radius_native_validation.o
MAINLIBOBJ += util/sep_wave_normalization_validation.o
MAINLIBOBJ += util/sep_real_wave_normalization_validation.o
MAINLIBOBJ += util/sep_polarized_wave_state.o util/sep_polarized_wave_state_validation.o
MAINLIBOBJ += polarized_wave_native_adapter.o util/sep_polarized_native_validation.o
MAINLIBOBJ += util/sep_birth_window.o util/sep_birth_window_validation.o
MAINLIBOBJ += birth_window_native.o
MAINLIBOBJ += util/sep_birth_native_validation.o
MAINLIBOBJ += util/sep_frozen_collision_transaction.o util/sep_frozen_collision_transaction_validation.o
MAINLIBOBJ += util/sep_flux_tube_state.o util/sep_flux_tube_state_validation.o
MAINLIBOBJ += util/sep_exact_nonnegative_sum.o util/sep_polarized_wkb_transport.o util/sep_polarized_wkb_validation.o
MAINLIBOBJ += util/sep_maxwell_juttner_flux.o util/sep_maxwell_juttner_flux_validation.o
# Volume initialization has a different phase measure from boundary flux.
# Link the production sampler, fresh SWMF source adapter and the same registered
# verification callbacks into both standalone and coupled application archives.
# Omitting these objects can compile callers successfully yet fail SWMF linking.
MAINLIBOBJ += util/sep_maxwell_juttner_volume.o util/sep_maxwell_juttner_volume_validation.o thermal_population.o
MAINLIBOBJ += util/sep_material_tube_kinematics.o util/sep_material_tube_kinematics_validation.o
MAINLIBOBJ += util/sep_standard_batsrus_waves.o standard_batsrus_wave_adapter.o
MAINLIBOBJ += util/sep_standard_batsrus_wave_validation.o util/sep_standard_batsrus_wave_native_validation.o
MAINLIBOBJ += util/sep_population_native_validation.o
MAINLIBOBJ += util/sep_source_lifecycle.o
MAINLIBOBJ += util/sep_conservative_resampling.o
MAINLIBOBJ += util/sep_conservation_resampling_validation.o
MAINLIBOBJ += util/sep_population_production_adapter.o
MAINLIBOBJ += util/sep_resonant_collision.o util/sep_resonant_collision_validation.o
MAINLIBOBJ += physics_checkpoint.o util/sep_checkpoint_native_validation.o
MAINLIBOBJ += source_capabilities.o util/sep_source_native_validation.o
MAINLIBOBJ += util/sep_maxwellian_native_validation.o
MAINLIBOBJ += util/sep_cli_campaign_validation.o
MAINLIBOBJ += util/sep_production_cli_validation.o
MAINLIBOBJ += adapters/swcme3d_spherical_material_adapter.o test/stage13_material_adapter_validation.o
# XM01-XM03 callbacks must be linked into the production executable so their
# native artifacts cannot be generated by an untracked Python replacement.
MAINLIBOBJ += validation/cases/cross_model_validation_models.o
# The command-line parser is library code: component_tests.o and util/sep_production_cli_validation.o
# call it, and the SWMF-coupled build (PT/AMPS) links mainlib.a without main.o.
MAINLIBOBJ += util/sep_cli.o

# Define executable objects BEFORE the immediate dependency-list expansion.
# Otherwise := sees an empty MAINOBJ, omits main/debug .d files, and a changed
# PIC ABI or worker-count header can leave old template instantiations in
# main.a even though every library/mover object has rebuilt successfully.
MAINOBJ = main.o debug/sep_debug_fieldline_datum.o util/sep_stage13_validation.o

# Bootstrap dependency files once for objects left by an older configured tree.
# This preserves generated configuration while preventing stale header consumers
# from surviving the first build after dependency tracking is introduced.
SEP_APP_OBJECTS := $(MAINOBJ) $(MAINLIBOBJ)
SEP_DEPENDENCY_OBJECTS := $(patsubst %.d,%.o,$(wildcard $(SEP_APP_OBJECTS:.o=.d)))
.PHONY: sep-header-dependency-bootstrap
$(filter-out $(SEP_DEPENDENCY_OBJECTS),$(SEP_APP_OBJECTS)): sep-header-dependency-bootstrap
-include $(SEP_APP_OBJECTS:.o=.d)

# sep_cli.o (srcSEP/util, command-line parsing apart from the physics driver) is in mainlib.a: library objects
# use it (see MAINLIBOBJ above), and main.a, linked first, takes it from there.

lib: $(SEP_COMMON_ARCHIVE) $(PARALLEL_DIFFUSION_ARCHIVE) $(SWCME_ARCHIVE) \
     $(SEP_CORONAL_CME_ARCHIVE) $(SEP_CORONA_SWCME_ARCHIVE) $(MAINLIBOBJ)
	# Recreate the application archive. Plain `ar -r` would retain a stale
	# shared member after a manifest change and could silently introduce a
	# second shared-kernel or SWCME definition.
	rm -f mainlib.a
	ar -rc mainlib.a $(MAINLIBOBJ) $(SEP_COMMON_OBJECTS) \
	  $(PARALLEL_DIFFUSION_OBJECTS) $(SWCME_OBJECTS) \
	  $(SEP_CORONAL_CME_OBJECTS) $(SEP_CORONA_SWCME_OBJECTS)
	@members=$$(ar t mainlib.a); \
	for member in $(notdir $(SEP_COMMON_OBJECTS) $(PARALLEL_DIFFUSION_OBJECTS) \
	    $(SWCME_OBJECTS) \
	    $(SEP_CORONAL_CME_OBJECTS) $(SEP_CORONA_SWCME_OBJECTS)); do \
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
.PHONY: sep-common parallel-diffusion swcme sep-coronal-cme sep-corona-swcme
sep-common:
	$(MAKE) -C $(SEP_COMMON_DIR)

$(SEP_COMMON_ARCHIVE): sep-common

parallel-diffusion:
	$(MAKE) -C $(PARALLEL_DIFFUSION_DIR) -f makefile all

$(PARALLEL_DIFFUSION_ARCHIVE): parallel-diffusion

swcme:
	$(MAKE) -C $(SWCME_DIR)

$(SWCME_ARCHIVE): swcme

sep-coronal-cme:
	$(MAKE) -C $(SEP_CORONAL_CME_DIR)

$(SEP_CORONAL_CME_ARCHIVE): sep-coronal-cme

sep-corona-swcme:
	$(MAKE) -C $(SEP_CORONA_SWCME_DIR)

$(SEP_CORONA_SWCME_ARCHIVE): sep-corona-swcme

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

# Use the canonical public leaf headers from their owning include directories
# in both source and copied build/main layouts.
adapters/swcme3d_spherical_material_adapter.o: adapters/swcme3d_spherical_material_adapter.cpp
	$(CXX) $(SWCME_CONSUMER_CXXFLAGS) -I$(SEP_CORONA_SWCME_DIR)/include/sep_corona_swcme -c -o $@ $<

util/sep_stage13_validation.o: util/sep_stage13_validation.cpp
	$(CXX) $(SWCME_CONSUMER_CXXFLAGS) -c -o $@ $<

# The linked bridge test owns no MPI/PIC state. Compile with the same private
# model headers as its adapter and the canonical manufactured test fixture.
test/stage13_material_adapter_validation.o: test/stage13_material_adapter_validation.cpp
	$(CXX) $(SWCME_CONSUMER_CXXFLAGS) -I$(SEP_CORONA_SWCME_DIR)/test -c -o $@ $<

# This native consumer needs both PIC/MPI and the installed coronal source
# headers. Preserve the configured native compiler/FLAGCC contract (including
# its MPI definitions and application language mode), adding only model header
# paths missing from the generic suffix rule. The dependency-light physics
# kernels have their separate strict ISO/Werror verification build.
reduced_shock_particles.o: reduced_shock_particles.cpp
	$(COMPILE.mpicxx) $(FLAGCC) $(MODEL_INCLUDE_FLAGS) $< -o $@

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
	for member in $(notdir $(SEP_COMMON_OBJECTS) $(PARALLEL_DIFFUSION_OBJECTS) \
	    $(SWCME_OBJECTS) \
	    $(SEP_CORONAL_CME_OBJECTS) $(SEP_CORONA_SWCME_OBJECTS)); do \
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

.PHONY: test-cli-unit test-initialization-unit test-state-unit test-geometry-source-unit test-mover-api-unit test-field-line-scope-unit test-transport-common-unit test-parker-unit test-fte-dmumu-unit test-fte-mfp-unit test-coefficients-unit test-turbulence-core-unit test-reproducibility-unit test-wp01-wp10-unit test-wp11-wp20-unit test-wp21-wp30-unit test-wp31-wp41-unit test-wp42-wp64-experimental test-wp59-wp64-native test-python-runner-unit test-cv01-unit test-cv02-cv05-unit test-cv06-cv12-unit test-iv01-iv06-unit test-xm01-xm03-unit test-ov01-ov05-unit test-ev01-ev02-unit print-configuration-matrix test-acceptance-unit test-documentation-unit test-controlled-analytical test-scientific-validation test-native-amps-validation test-swmf-validation test-observational-validation test-sanitizer test-stochastic-repeat test-swcme-relocation-unit test-sep-common-ownership-unit test-swcme-fail-closed-unit test-swcme-configuration-unit test-reduced-shock-background-unit test-d03-native-integration test-d03-runner-unit test-stage1-stage2-contracts test-stage3-contracts test-package-hygiene-unit require-sep-executable test-cli test-list test-case \
	test-group test-turbulence test-stress test-mpi test

# Source-only audit for the Stage 1–2 integration boundaries.  It stays
# separate from linked native tests so a detached review archive can still
# prove single-owner calls, provider-aware initialization, and fixture guards.
test-stage1-stage2-contracts:
	@PYTHONDONTWRITEBYTECODE=1 python3 test/check_stage1_stage2_contracts.py

# Stage 3 is a source/release-structure gate.  It deliberately runs before
# expensive scientific cases and proves that retired code, case-registration
# policy, and the proton-only srcSEP3D contract cannot drift independently.
test-stage3-contracts:
	@PYTHONDONTWRITEBYTECODE=1 python3 test/check_stage3_contracts.py

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

test-reduced-shock-background-unit:
	./test/run_reduced_shock_background_tests.sh

.PHONY: test-standalone-integration-unit
test-standalone-integration-unit:
	python3 test/run_standalone_integration_tests.py

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

# PDB01-PDB04: schema-4 [parallel_diffusion] location, CLI/section/mover
# cross-check, host-capability gate, and kappa through the shared library's
# ActiveParallelDiffusion pointer (no AMPS/MPI; the PIC sampling wrapper in
# coefficient_providers.cpp is covered only by a native build).
test-parallel-diffusion-binding-unit:
	./test/run_parallel_diffusion_binding_tests.sh

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

.PHONY: test-rng-lifecycle-unit
test-rng-lifecycle-unit:
	./test/run_rng_lifecycle_tests.sh

# The same bounded state-continuity gate belongs to both routine development
# aggregates; the native registry separately executes its shared descriptors.
test test-controlled-analytical: test-rng-lifecycle-unit

.PHONY: test-background-derivatives-unit
test-background-derivatives-unit:
	./test/run_background_derivative_tests.sh

test test-controlled-analytical: test-background-derivatives-unit

.PHONY: test-wave-normalization-unit
test-wave-normalization-unit:
	./test/run_wave_normalization_tests.sh

test test-controlled-analytical: test-wave-normalization-unit

# An independent real circular-wave Lorentz integration fixes the one-sided
# polarization convention; this is deliberately distinct from unit conversion.
.PHONY: test-real-wave-normalization-unit
test-real-wave-normalization-unit:
	./test/run_real_wave_normalization_tests.sh

test test-controlled-analytical: test-real-wave-normalization-unit

.PHONY: test-polarized-wave-state-unit
test-polarized-wave-state-unit:
	./test/run_polarized_wave_state_tests.sh

test test-controlled-analytical: test-polarized-wave-state-unit

.PHONY: test-birth-window-unit
test-birth-window-unit:
	./test/run_birth_window_tests.sh

test test-controlled-analytical: test-birth-window-unit

# Validation-only parallel warm-plasma reference; no production dependency or
# inferred physical cutoff is introduced by this numerical verification suite.
.PHONY: test-warm-parallel-study-unit
test-warm-parallel-study-unit:
	./test/run_warm_parallel_study_tests.sh

test test-controlled-analytical: test-warm-parallel-study-unit

# Bounded mathematical/protocol verification of the independent prescribed-wave
# Lorentz study. Production limits and the opt-in orbit campaigns are separate.
.PHONY: test-finite-amplitude-study-unit
test-finite-amplitude-study-unit:
	./test/run_finite_amplitude_study_tests.sh

test test-controlled-analytical: test-finite-amplitude-study-unit

# Prepare-only coupled bookkeeping: accepted frozen-wave proposals must pass
# provenance, representable transfer and all-bin budgets before a candidate
# exists. This target does not claim native publication or exact coupled rates.
.PHONY: test-frozen-collision-transaction-unit
test-frozen-collision-transaction-unit:
	./test/run_frozen_collision_transaction_tests.sh

test test-controlled-analytical: test-frozen-collision-transaction-unit

.PHONY: test-flux-tube-state-unit
test-flux-tube-state-unit:
	./test/run_flux_tube_state_tests.sh

test test-controlled-analytical: test-flux-tube-state-unit

# Controlled ideal-MHD WKB equations and exact external-account arithmetic;
# publication into evolving native material geometry is a separate gate.
.PHONY: test-polarized-wkb-unit
test-polarized-wkb-unit:
	./test/run_polarized_wkb_tests.sh

test test-controlled-analytical: test-polarized-wkb-unit

# These retained POSIX test scripts have source mode 0644. Invoke their
# declared shell explicitly so a source checkout/archive runs the exact tests
# without depending on executable-bit repair or bypassing any sanitizer gate.
# The exact relativistic material-plane flux and its independently evaluated
# cold limit are tested without enabling a native guiding-center injection.
.PHONY: test-maxwell-juttner-flux-unit
test-maxwell-juttner-flux-unit:
	$(SHELL) ./test/run_maxwell_juttner_flux_tests.sh

test test-controlled-analytical: test-maxwell-juttner-flux-unit

# A spatial-volume initial condition requires its own Lorentz measure and
# normalization checks; boundary-flux tests cannot stand in for this gate.
.PHONY: test-maxwell-juttner-volume-unit
test-maxwell-juttner-volume-unit:
	bash ./test/run_maxwell_juttner_volume_tests.sh

test test-controlled-analytical: test-maxwell-juttner-volume-unit

# Cauchy induction, true material mass and cap/side geometric conservation are
# qualified here only for the explicit neutral deformation/map capabilities.
.PHONY: test-material-tube-kinematics-unit
test-material-tube-kinematics-unit:
	./test/run_material_tube_kinematics_tests.sh

test test-controlled-analytical: test-material-tube-kinematics-unit

# Standard BATSRUS total wave energy has an explicit immutable magnetic view;
# the dependency-light conversion is separate from actual native array reads.
.PHONY: test-standard-batsrus-wave-unit
test-standard-batsrus-wave-unit:
	$(SHELL) ./test/run_standard_batsrus_wave_tests.sh

test test-controlled-analytical: test-standard-batsrus-wave-unit

# A fixed analytic provider must use its declared origin/speed and exact
# Cartesian derivative; geometry-only creation has its own validation boundary.
.PHONY: test-analytic-radial-background-unit
test-analytic-radial-background-unit:
	./test/run_analytic_radial_background_tests.sh

test test-controlled-analytical: test-analytic-radial-background-unit

# Source radius and speed must belong to one immutable departure epoch; this
# tests the existing trajectory law without reading or advancing a host clock.
.PHONY: test-analytic-shock-epoch-unit
test-analytic-shock-epoch-unit:
	./test/run_analytic_shock_epoch_tests.sh

test test-controlled-analytical: test-analytic-shock-epoch-unit

.PHONY: test-conservation-resampling-unit
test-conservation-resampling-unit:
	./test/run_conservation_resampling_tests.sh

test test-controlled-analytical: test-conservation-resampling-unit

.PHONY: test-resonant-collision-unit
test-resonant-collision-unit:
	./test/run_resonant_collision_tests.sh

test test-controlled-analytical: test-resonant-collision-unit

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
	@PYTHONDONTWRITEBYTECODE=1 python3 test/test_routine_fixtures.py
	@PYTHONDONTWRITEBYTECODE=1 python3 test/test_routine_wiring.py

test-xm01-xm03-unit:
	./test/run_xm01_xm03_tests.sh

# The shell gates authenticate the fixed inputs and provenance without
# claiming a scientific PASS when no linked AMPS executable is available.
# If SEP_EXECUTABLE is set, OV01-OV05 additionally execute through the linked
# registry; EV01/EV02 remain explicitly labelled as a nine-event pilot.
test-ov01-ov05-unit:
	./test/run_ov01_ov05_tests.sh

test-ev01-ev02-unit:
	./test/run_ev01_ev02_tests.sh

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
	@PYTHONDONTWRITEBYTECODE=1 python3 test/check_stage1_stage2_contracts.py
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
test-scientific-validation: test-ov01-ov05-unit test-ev01-ev02-unit
	./test/run_step15_tests.sh

# Native, SWMF, and observational checks remain distinct so a PASS in one
# class cannot satisfy another.  The external-data targets authenticate all
# referenced files through validation/check_external_gate.py before reporting
# success; unset inputs fail with usage guidance rather than being skipped.
# One maximum native startup cannot satisfy all Routine physics fixtures.
# Keep every compiled Routine ID, but let the maintained descriptor-aware
# runner prepare each actual protocol independently. This includes the full
# independent CV/IV comparisons and all Birth control roles; no SKIP/missing
# fixture can qualify. JSON/JUNIT retain their existing Make destinations.
SEP_ROUTINE_OUTPUT ?= $(if $(SEP_TEST_OUTPUT_DIR),$(abspath $(SEP_TEST_OUTPUT_DIR))/routine,$(if $(SRCSEP_REPORT_DIR),$(abspath $(SRCSEP_REPORT_DIR))/routine,$(abspath test_output/routine)))
SEP_ROUTINE_COMMAND = PYTHONDONTWRITEBYTECODE=1 python3 "$(SEP_MAKEFILE_DIR)/test/run_tests.py" --amps "$(SEP_EXECUTABLE)" --routine --output-dir "$(SEP_ROUTINE_OUTPUT)" --plot none $(if $(JSON),--report-json "$(JSON)",) $(if $(JUNIT),--report-junit "$(JUNIT)",)

test-native-amps-validation: require-sep-executable
	@$(SEP_ROUTINE_COMMAND)

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
	./test/run_rng_lifecycle_tests.sh
	./test/run_background_derivative_tests.sh
	./test/run_analytic_radial_background_tests.sh
	./test/run_analytic_shock_epoch_tests.sh
	./test/run_wave_normalization_tests.sh
	./test/run_real_wave_normalization_tests.sh
	./test/run_polarized_wave_state_tests.sh
	./test/run_birth_window_tests.sh
	./test/run_frozen_collision_transaction_tests.sh
	./test/run_flux_tube_state_tests.sh
	./test/run_polarized_wkb_tests.sh
	$(SHELL) ./test/run_maxwell_juttner_flux_tests.sh
	./test/run_material_tube_kinematics_tests.sh
	$(SHELL) ./test/run_standard_batsrus_wave_tests.sh
	./test/run_conservation_resampling_tests.sh
	./test/run_resonant_collision_tests.sh
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

.PHONY: test-cli-native test-cli-campaign-native
# These execute the linked program, including applied model globals and MPI
# failure propagation. Their reports retain unqualified production/SWMF gates;
# the software compatibility target does not turn those SKIPs into validation.
test-cli-native: require-sep-executable
	@PYTHONDONTWRITEBYTECODE=1 python3 test/run_native_production_cli_tests.py \
	  --amps "$(SEP_EXECUTABLE)"

test-cli-campaign-native: require-sep-executable
	@PYTHONDONTWRITEBYTECODE=1 python3 test/run_native_cli_campaign_tests.py \
	  --amps "$(SEP_EXECUTABLE)"

test-cli: test-cli-unit test-cli-native test-cli-campaign-native require-sep-executable
	@"$(SEP_EXECUTABLE)" --help >/dev/null
	@"$(SEP_EXECUTABLE)" --list-tests >/dev/null

test-list: require-sep-executable
	@"$(SEP_EXECUTABLE)" --list-tests

test-case: require-sep-executable
	@test -n "$(CASE)" || { echo "ERROR: use make test-case CASE=<stable-ID>." >&2; exit 2; }
	@PYTHONDONTWRITEBYTECODE=1 python3 "$(SEP_MAKEFILE_DIR)/test/run_tests.py" --amps "$(SEP_EXECUTABLE)" --test "$(CASE)" --output-dir "$(SEP_ROUTINE_OUTPUT)/selected" --plot none $(if $(JSON),--report-json "$(JSON)",) $(if $(JUNIT),--report-junit "$(JUNIT)",)

test-group: require-sep-executable
	@test -n "$(GROUP)" || { echo "ERROR: use make test-group GROUP=<name>." >&2; exit 2; }
	@PYTHONDONTWRITEBYTECODE=1 python3 "$(SEP_MAKEFILE_DIR)/test/run_tests.py" --amps "$(SEP_EXECUTABLE)" --group "$(GROUP)" --output-dir "$(SEP_ROUTINE_OUTPUT)/selected" --plot none $(if $(JSON),--report-json "$(JSON)",) $(if $(JUNIT),--report-junit "$(JUNIT)",)

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
	@$(SEP_ROUTINE_COMMAND) --mpi-np "$(MPI_NP)" --mpiexec "$(MPIEXEC)"

.PHONY: test-stage13-cases
test-stage13-cases: require-sep-executable
	@$(MAKE) -C "$(SEP_CORONAL_CME_DIR)" test-tube-isopoly
	@$(MAKE) -C "$(SEP_CORONA_SWCME_DIR)" test-spherical-shock-only
	@python3 -B "$(SEP_MAKEFILE_DIR)/test/test_stage13_validation_cli.py" --amps "$(SEP_EXECUTABLE)"

# Keep the top-level sequence explicit and serial even when invoked as
# `make -j test`: component tests share initialized model/output state and must
# not be launched concurrently.  The embedded SWCME suite retains its own
# internal parallel build through its native target.
test: test-swcme-relocation-unit test-cli-unit test-reduced-shock-background-unit test-initialization-unit test-state-unit test-geometry-source-unit test-mover-api-unit test-field-line-scope-unit test-transport-common-unit test-parker-unit test-fte-dmumu-unit test-fte-mfp-unit test-coefficients-unit test-turbulence-core-unit test-reproducibility-unit test-wp01-wp10-unit test-wp11-wp20-unit test-wp21-wp30-unit test-wp31-wp41-unit test-python-runner-unit test-cv01-unit test-acceptance-unit test-documentation-unit require-sep-executable
	@"$(SEP_EXECUTABLE)" --list-tests >/dev/null
	@$(MAKE) test-stage13-cases SEP_EXECUTABLE="$(SEP_EXECUTABLE)"
	@$(SEP_ROUTINE_COMMAND)
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

# Real native qualification requires an explicitly matching configured archive.
# The existing Python source-suite interface propagates every child ERROR/FAIL;
# it never substitutes a serial dependency implementation for actual MPI.
.PHONY: test-standalone-integration-native
test-standalone-integration-native:
	python3 test/run_native_standalone_integration.py

# Moving restart is a real-MPI gate, separate from source-only codec controls.
# Preparing the same fixture records the complete PIC dependency closure; it
# does not execute or certify transport. The second command executes all paired
# roles (including 1<->8 repartition) and propagates missing/failed jobs. Never
# substitute the serial MPI dependency adapter when this native gate is chosen.
MATERIAL_CHECKPOINT_OUTPUT ?= $(abspath test_output/material-checkpoint)
NATIVE_BUILD_ROOT ?= $(if $(AMPS_FIELD_LINE_NATIVE_BUILD),$(AMPS_FIELD_LINE_NATIVE_BUILD),$(AMPS_ROOT))
.PHONY: test-material-checkpoint-native test-material-checkpoint-unit
test-material-checkpoint-native:
	python3 test/run_native_standalone_integration.py --build-root "$(NATIVE_BUILD_ROOT)" --output-dir "$(MATERIAL_CHECKPOINT_OUTPUT)/native" --ranks 1 2 4 8 --prepare-only
	python3 test/run_material_checkpoint_tests.py --native-prepared "$(MATERIAL_CHECKPOINT_OUTPUT)/native/prepared.json" --output "$(MATERIAL_CHECKPOINT_OUTPUT)/campaign"

test-material-checkpoint-unit:
	PYTHONDONTWRITEBYTECODE=1 python3 test/test_material_checkpoint_runner.py
	PYTHONDONTWRITEBYTECODE=1 python3 test/test_native_dependency_closure.py

# The combined gate requires separately configured, dependency-audited worker
# builds. Explicit paths prevent silently running a single-thread executable
# while claiming two-/four-thread restart coverage. All outputs must be fresh.
STAGE7_CHECKPOINT_OUTPUT ?= $(if $(SRCSEP_REPORT_DIR),$(abspath $(SRCSEP_REPORT_DIR))/stage7-native,$(abspath test_output/stage7-native))
.PHONY: test-stage7-checkpoint-native test-stage7-checkpoint-unit
test-stage7-checkpoint-native:
	python3 test/run_stage7_checkpoint_tests.py --build-root "$(NATIVE_BUILD_ROOT)" --output "$(STAGE7_CHECKPOINT_OUTPUT)"

test-stage7-checkpoint-unit:
	PYTHONDONTWRITEBYTECODE=1 python3 test/test_stage7_checkpoint_runner.py

# These bounded core/adapter cases are routine code-validation prerequisites.
# Native MPI qualification remains an explicit configured-build step above.
test test-controlled-analytical: test-standalone-integration-unit test-material-checkpoint-unit test-stage7-checkpoint-unit

# Backend admission is routine verification. Unsupported core modes retain
# their diagnostic APIs but cannot silently enter moving material transport.
.PHONY: test-backend-admission-unit test-backend-admission-native test-standalone-application-runner-unit test-standalone-application-native
test-backend-admission-unit:
	python3 -B test/run_backend_admission_tests.py
test-backend-admission-native:
	python3 -B test/run_backend_admission_tests.py --native --amps "$(NATIVE_BUILD_ROOT)/amps"
test-standalone-application-runner-unit:
	python3 -B test/test_standalone_application_runner.py
test-standalone-application-native:
	python3 -B test/run_standalone_application.py --build-root "$(NATIVE_BUILD_ROOT)" --output-dir "$(APPLICATION_OUTPUT)"
APPLICATION_OUTPUT ?= $(abspath test_output/standalone-application)
test test-controlled-analytical: test-backend-admission-unit test-standalone-application-runner-unit

# Audit completed real-MPI application receipts and a complete physical output
# window. Explicit evidence is required: missing/crashed/partial campaigns are
# errors, never a reason to substitute the shorter native/component fixtures.
.PHONY: test-standalone-integration-application
test-standalone-integration-application:
	@test -n "$(EVIDENCE_DIR)" || { echo "Set EVIDENCE_DIR to the completed stage-6 application campaign"; exit 1; }
	python3 test/verify_standalone_application.py --campaign-dir "$(EVIDENCE_DIR)"
	python3 test/test_standalone_application_verifier.py --campaign-dir "$(EVIDENCE_DIR)"

# Stage 11 shares actual canonical adapters with production. Missing physics
# stays an explicit capability error; these tests never admit another backend.
.PHONY: test-material-background-unit
test-material-background-unit:
	python3 -B test/run_material_background_tests.py
test test-controlled-analytical: test-material-background-unit
