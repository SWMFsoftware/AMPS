SHELL := /bin/sh

CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O2

CPPFLAGS += -isystem .. -Icore
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

TARGET := output/test_swcme
PST08_TARGET := output/test_swcme_pst08_performance
OUT08_DEBUG_DIR := output/out08-debug
OUT08_OPTIMIZED_DIR := output/out08-optimized
SEP_REFERENCE_TARGET := output/sep_reference
DEMO_TARGETS := output/demo1d output/demo3d_1 output/demo3d_2
# Every validation remains in one deterministic executable. The common tests
# use separate translation units so each physical contract has an auditable
# implementation without introducing another executable or build system.
SOURCES := core/test_swcme.cpp core/test_state_immutability.cpp core/test_prepared_state.cpp core/test_configuration_state.cpp core/test_state_integrity.cpp core/test_concurrent_state.cpp core/test_output.cpp core/test_output_parsing.cpp core/test_demo_execution.cpp core/test_adapter_equivalence.cpp core/test_state_lifetime.cpp core/test_state_performance.cpp core/test_strict_writer.cpp core/test_smoothing_policy.cpp core/test_radius_domain.cpp core/test_config.cpp core/test_units.cpp core/test_density.cpp core/test_kinematics.cpp core/test_shock.cpp core/test_regions.cpp core/test_acceleration.cpp core/test_status.cpp core/test_divergence.cpp core/test_defaults.cpp core/test_sep_interface.cpp 3d/test_parker.cpp 3d/test_geometry.cpp 3d/test_mesh.cpp 3d/test_connectivity.cpp 3d/test_shock_surface.cpp 1d/test_ambient.cpp test_1d3d.cpp
OBJECTS := output/swcme3d.o output/test_swcme.o output/test_state_immutability.o output/test_prepared_state.o output/test_configuration_state.o output/test_state_integrity.o output/test_concurrent_state.o output/test_output.o output/test_output_parsing.o output/test_demo_execution.o output/test_adapter_equivalence.o output/test_state_lifetime.o output/test_state_performance.o output/test_strict_writer.o output/test_smoothing_policy.o output/test_radius_domain.o output/test_config.o output/test_units.o output/test_density.o output/test_kinematics.o output/test_shock.o output/test_regions.o output/test_acceleration.o output/test_status.o output/test_divergence.o output/test_defaults.o output/test_sep_interface.o output/test_parker.o output/test_geometry.o output/test_mesh.o output/test_connectivity.o output/test_shock_surface.o output/test_ambient.o output/test_1d3d.o

# POSIX thread support is isolated behind one overridable variable.  PST04's
# translation unit needs the compile flag on implementations that expose
# std::thread through pthreads, and the final executable needs the link flag.
THREAD_FLAGS ?= -pthread

# LeakSanitizer cannot inspect threads when the test process is itself observed
# through ptrace (as in several CI/container runners).  Keep the mandatory
# address/use-after-free gate portable; an untraced host can opt into leak
# scanning with `make PST05_ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 ...`.
PST05_ASAN_OPTIONS ?= detect_leaks=0:halt_on_error=1

# OUT08 deliberately does not inherit CXXFLAGS: both variants must use the
# documented CI warning policy even when a caller's normal development build
# is permissive.  -Werror turns every selected diagnostic into a build failure;
# the explicit -Werror=format-security keeps the security requirement visible
# in logs and review even though it is also covered by the global promotion.
OUT08_WARNING_FLAGS := -std=c++17 -Wall -Wextra -Wpedantic -Wformat=2 \
	-Wformat-security -Wconversion -Wsign-conversion -Wshadow \
	-Werror=format-security -Werror

.PHONY: all test tools demos demo-run campaign-test pst05-sanitize pst08-performance out08-strict clean

all: $(TARGET) demos

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(THREAD_FLAGS) $(OBJECTS) -o $@

# Both dimensional model headers include the shared explicit-field integrity
# serializer and checked output layer.  Make every validation object observe
# both headers so an incremental build cannot link stale inline seal or output
# logic from an older test.
$(OBJECTS): ../swcme_prepared_integrity.hpp ../swcme_output.hpp

# Compile the production 3-D implementation exactly once and link it normally.
# Earlier validation revisions included swcme3d.cpp from test_config.cpp merely
# to obtain its symbols, which coupled the entire production implementation to
# one unrelated test translation unit.  A normal object-library style link is
# clearer, scales to additional tests, and prevents accidental duplicate code.
output/swcme3d.o: ../swcme3d.cpp ../swcme3d.hpp ../swcme_constants.hpp ../swcme_units.hpp ../swcme_defaults.hpp ../swcme_status.hpp ../swcme_output.hpp ../swcme_divergence.hpp ../swcme_config.hpp ../swcme_kinematics.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_shock.hpp ../swcme_regions.hpp ../swcme_acceleration.hpp ../swcme_sep_source.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Explicit rules keep source ownership obvious and let `make -j` build the
# independent validation units concurrently without introducing another build
# system. Header dependencies force the relevant objects to refresh.
output/test_swcme.o: core/test_swcme.cpp core/test_framework.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# PST01 freezes model configuration after the first successful prepared state,
# exercises every legacy mutation entry point, and verifies bitwise-stable
# reevaluation plus the replacement-model builder in both dimensions.
output/test_state_immutability.o: core/test_state_immutability.cpp core/test_framework.hpp ../swcme_status.hpp ../swcme1d.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# PST02 verifies the model-instance ownership token at direct 1-D/3-D APIs,
# AMPS-facing adapters, legacy geometry wrappers, and checked file writers.
# Keeping it in a dedicated translation unit makes this cross-cutting safety
# contract easy to audit without coupling it to unrelated physics fixtures.
output/test_prepared_state.o: core/test_prepared_state.cpp core/test_framework.hpp ../swcme_status.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_sep_interface.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# PST03 verifies deterministic complete configuration digests, one-field and
# multi-field cross-configuration rejection, and stale same-model 1-D states.
# The separate unit keeps configuration provenance independent of PST02's
# exhaustive API-boundary matrix.
output/test_configuration_state.o: core/test_configuration_state.cpp core/test_framework.hpp ../swcme_status.hpp ../swcme1d.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# PST06 protects canonical caches and public compatibility mirrors with a
# private explicit-field seal, verifies valid copy/move behavior, and rejects
# controlled record corruption before any checked output is modified.
output/test_state_integrity.o: core/test_state_integrity.cpp core/test_framework.hpp ../swcme_status.hpp ../swcme_prepared_integrity.hpp ../swcme1d.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# PST04 shares the exact same prepared state among 1, 2, 4, and 8 threads and
# compares scalar background, batch field, shock, source, connectivity, and
# failure-diagnostic records with an exact serial oracle under three schedules.
output/test_concurrent_state.o: core/test_concurrent_state.cpp core/test_framework.hpp ../swcme_status.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_sep_source.hpp ../swcme_sep_interface.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(THREAD_FLAGS) -c $< -o $@

# OUT02 injects deterministic stream faults, OUT03 checks staging/commit, OUT05
# proves domain preflight, OUT04 validates BoxSpec, and OUT06 seals mesh/metric
# records before output.  Keeping these failure contracts in one unit lets them
# cover header-only 1-D and compiled 3-D paths without adding test-only state to
# either physics model.
output/test_output.o: core/test_output.cpp core/test_framework.hpp ../swcme_status.hpp ../swcme_output.hpp ../swcme1d.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# OUT01 is intentionally a separate translation unit with a validation-owned
# Tecplot grammar and unit schema.  It calls public writers but imports no
# production parsing/format constants, so a writer defect cannot update both
# sides of the comparison through one shared helper.
output/test_output_parsing.o: core/test_output_parsing.cpp core/tecplot_parser.hpp core/test_framework.hpp ../swcme_status.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# OUT07 launches the real user-facing executables in isolated empty
# directories, validates exact artifact manifests, parses every Tecplot/CSV
# product, and cleans successful runs.  Keeping process orchestration outside
# OUT01 leaves the parser reusable while preserving its writer independence.
output/test_demo_execution.o: core/test_demo_execution.cpp core/tecplot_parser.hpp core/test_framework.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# PST07 treats the AMPS adapter as a transparent integration boundary.  The
# test compares direct and adapter status/provenance, background fields,
# pressure/focusing, shock/source records, spectra, and connectivity length in
# both dimensions without importing any adapter-private implementation detail.
output/test_adapter_equivalence.o: core/test_adapter_equivalence.cpp core/test_framework.hpp ../swcme_status.hpp ../swcme_solarwind.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_sep_source.hpp ../swcme_sep_interface.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# PST05 exercises destruction, model/adapter moves, rejected assignment into a
# frozen destination, vector relocation, and asynchronous evaluation.  It uses
# pthread at compile time because std::async is an explicit part of the test.
output/test_state_lifetime.o: core/test_state_lifetime.cpp core/test_framework.hpp ../swcme_status.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_sep_interface.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(THREAD_FLAGS) -c $< -o $@

# PST08 measures the full state guard separately from direct and AMPS scalar/
# batch queries, emulates the removed duplicate-validation topology, and counts
# every C++17 allocation form in warmed hot paths.
output/test_state_performance.o: core/test_state_performance.cpp core/test_framework.hpp ../swcme_status.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_sep_interface.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(THREAD_FLAGS) -c $< -o $@

# OUT08's registered runtime half exercises literal and formatted record
# semantics.  Its compile-time half is the isolated out08-strict target below.
output/test_strict_writer.o: core/test_strict_writer.cpp core/test_framework.hpp ../swcme_status.hpp ../swcme_output.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# CFG03 verifies the shared reject-before-prepare smoothing policy at zero,
# ordinary, exact-limit, one-ULP-invalid, non-finite, and multi-field inputs in
# both dimensions, then audits the effective local widths and blend geometry.
output/test_smoothing_policy.o: core/test_smoothing_policy.cpp core/test_framework.hpp ../swcme_config.hpp ../swcme_regions.hpp ../swcme1d.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# CFG04 checks every configured/public radius at the exact 1.05-Rsun domain
# boundary and adjacent ULPs, including indexed PCHIP diagnostics and the
# source-surface/observer ordering used by 3-D connectivity.
output/test_radius_domain.o: core/test_radius_domain.cpp core/test_framework.hpp ../swcme_config.hpp ../swcme_solarwind.hpp ../swcme1d.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output/test_config.o: core/test_config.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme_units.hpp ../swcme_defaults.hpp ../swcme_status.hpp ../swcme_divergence.hpp ../swcme_config.hpp ../swcme_kinematics.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_regions.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output/test_units.o: core/test_units.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme_units.hpp ../swcme_defaults.hpp ../swcme_status.hpp ../swcme_divergence.hpp ../swcme_config.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_regions.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# DEN01 calls the production 3-D implementation linked through swcme3d.o and
# keeps its analytical Leblanc reference independent of production code.
output/test_density.o: core/test_density.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_regions.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# KIN01-KIN09 validate the common apex-kinematics engine directly and through
# both dimensional wrappers.  The test contains independent closed-form DBM
# references, dense PCHIP monotonicity checks, and explicit extrapolation-domain
# stress at radial crossings, turning points, poles, and floating-point limits;
# production kinematics remain shared so the dimensional models cannot diverge.
output/test_kinematics.o: core/test_kinematics.cpp core/test_framework.hpp ../swcme_kinematics.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_regions.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# SHK01-SHK12 exercise the shared production ideal-MHD jump solver and its
# integration into the 1-D/3-D models.  The test file contains independent
# limiting-case references and conservation residual checks; it does not
# duplicate the production nonlinear solver.
output/test_shock.o: core/test_shock.cpp core/test_framework.hpp ../swcme_shock.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_regions.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


# REG01-REG05 validate the shared SHOCK_ONLY/FULL_ICME region contract,
# including exact RH downstream boundaries, ejecta factors below unity,
# self-similar local layer geometry, and C1 artificial interface blends.
output/test_regions.o: core/test_regions.cpp core/test_framework.hpp ../swcme_regions.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_shock.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# ACC01-ACC05 enforce the single shock-acceleration representation, verify the
# SOURCE no-double-counting gate, and compare the shared resolved-compression
# shock smoothing used by the 1-D and 3-D transport-facing fields.
output/test_acceleration.o: core/test_acceleration.cpp core/test_framework.hpp ../swcme_acceleration.hpp ../swcme_regions.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_shock.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


# ERR01-ERR05 verify explicit numerical-status propagation, outside-domain
# handling, shock-solver outcome codes, bad-vector rejection, and writer
# refusal to sanitize non-finite physics values.
output/test_status.o: core/test_status.cpp core/test_framework.hpp ../swcme_status.hpp ../swcme_divergence.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_shock.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# DIV01-DIV03 validate exact radial divergence, analytical radial-profile
# derivatives, and the shared full Cartesian second-order divergence operator.
output/test_divergence.o: core/test_divergence.cpp core/test_framework.hpp ../swcme_divergence.hpp ../swcme_regions.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_status.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# DEF01-DEF04 guard the canonical shared baseline, explicit model-scope
# convention, observer-local pre-shock validity gate, and deterministic
# resolved-configuration metadata consumed by the campaign/SEP integration layer.
output/test_defaults.o: core/test_defaults.cpp core/test_framework.hpp ../swcme_defaults.hpp ../swcme1d.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


# SEP01-SEP06 validate the stable SWCME-to-AMPS contract: direct background
# identity, unit/spectrum conversion, 1-D/3-D source identity, area-weighted
# surface patches, cobpoint source reuse, and SOURCE disablement in resolved mode.
output/test_sep_interface.o: core/test_sep_interface.cpp core/test_framework.hpp ../swcme_sep_source.hpp ../swcme_sep_interface.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_status.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# PAR01-PAR03 exercise the public production 3-D field evaluator linked through
# swcme3d.o; the test translation unit contains only independent references.
output/test_parker.o: 3d/test_parker.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme3d.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


# GEO01-GEO08 validate finite shock geometry using the production swcme3d.o
# implementation and independent geometric references.
output/test_geometry.o: 3d/test_geometry.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme3d.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


# MSH01-MSH05 validate the explicit triangular-manifold topology, cell
# orientation/area convergence, and canonical area-weighted source sampler.
output/test_mesh.o: 3d/test_mesh.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme3d.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# CON01-CON10 validate analytical Parker-line mapping, cobpoint behavior,
# explicit resolution-budget qualification, and observer-domain classification
# against the production swcme3d.o implementation and operational SEP adapter.
output/test_connectivity.o: 3d/test_connectivity.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme3d.hpp ../swcme_sep_interface.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# SHK13-SHK14 are 3-D integration guards for the surface-owned shock-state
# contract.  They prove that arbitrary Cartesian query radii cannot leak into
# local shock strength and that diagnostics/mesh nodes reuse the canonical
# shock_state_direction() result rather than maintaining parallel calculations.
output/test_shock_surface.o: 3d/test_shock_surface.cpp core/test_framework.hpp ../swcme3d.hpp ../swcme_shock.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_regions.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output/test_ambient.o: 1d/test_ambient.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_regions.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# 1D3D01-02 validate that the two public dimensional wrappers reduce to the
# same common solar-wind, apex-kinematics, and shock state in an exactly
# equivalent spherical/equatorial geometry.
output/test_1d3d.o: test_1d3d.cpp core/test_framework.hpp ../swcme_core.hpp ../swcme_solarwind.hpp ../swcme_shock.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_regions.hpp ../swcme_acceleration.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


tools: $(SEP_REFERENCE_TARGET)

$(SEP_REFERENCE_TARGET): output/swcme3d.o ../tools/sep_reference.cpp ../swcme_sep_interface.hpp ../swcme_sep_source.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) ../tools/sep_reference.cpp output/swcme3d.o -o $@

# Demos are first-class warning-checked build products.  The 3-D executables
# link the same single production object as the validation binary, while the
# header-only 1-D example needs no duplicate compiled implementation.
demos: $(DEMO_TARGETS)

output/demo1d: ../demo1d.cpp ../swcme1d.hpp ../swcme_output.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

output/demo3d_1: ../demo3d_1.cpp output/swcme3d.o ../swcme3d.hpp ../swcme_output.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) ../demo3d_1.cpp output/swcme3d.o -o $@

output/demo3d_2: ../demo3d_2.cpp output/swcme3d.o ../swcme3d.hpp ../swcme_output.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) ../demo3d_2.cpp output/swcme3d.o -o $@

# This target is the documented one-command acceptance gate: build all three
# examples, execute them through OUT07, parse every declared product, and fail
# make if any example or independent reader rejects the result.
demo-run: all
	./$(TARGET) --test OUT07

# Campaign-manager unit tests exercise profile expansion, convergence fitting,
# sweep/comparison execution, exit codes, manifest generation, and machine-
# readable reporting without duplicating production physics.
campaign-test: all tools
	python3 python/test_run_tests.py

# Build a separate fully instrumented validation executable so AddressSanitizer
# never contaminates normal incremental objects.  Running only PST05 keeps the
# dedicated lifetime gate fast while every production/test translation unit is
# compiled with identical ASan instrumentation. Leak scanning is opt-in because
# LeakSanitizer is incompatible with ptrace-based CI runners.
pst05-sanitize: | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -O1 -g -fsanitize=address -fno-omit-frame-pointer $(THREAD_FLAGS) ../swcme3d.cpp $(SOURCES) -o output/test_swcme_pst05_asan
	ASAN_OPTIONS=$(PST05_ASAN_OPTIONS) ./output/test_swcme_pst05_asan --test PST05

# Rebuild only production 3-D code and PST08 with pinned optimized flags.  This
# avoids stale normal objects when CXXFLAGS changes and produces a compact,
# reproducible performance log without compiling unrelated validation tests.
pst08-performance: | output
	$(CXX) $(CPPFLAGS) -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic $(THREAD_FLAGS) ../swcme3d.cpp core/test_state_performance.cpp core/test_state_performance_main.cpp -o $(PST08_TARGET)
	./$(PST08_TARGET)

# Compile fresh debug and optimized objects so incremental products built with
# ordinary flags cannot make the strict gate pass accidentally.  The matrix
# covers the compiled 3-D writers, header-only 1-D writer/status stack, all
# three user demonstrations, and a dedicated literal/format runtime probe.
# CXX remains overridable: CI runs this target once with each supported compiler
# (for example `make out08-strict CXX=g++` and `CXX=clang++`).
out08-strict: | output
	mkdir -p $(OUT08_DEBUG_DIR) $(OUT08_OPTIMIZED_DIR)
	$(CXX) $(CPPFLAGS) $(OUT08_WARNING_FLAGS) -O0 -g3 -c ../swcme3d.cpp -o $(OUT08_DEBUG_DIR)/swcme3d.o
	$(CXX) $(CPPFLAGS) $(OUT08_WARNING_FLAGS) -O0 -g3 ../demo1d.cpp -o $(OUT08_DEBUG_DIR)/demo1d
	$(CXX) $(CPPFLAGS) $(OUT08_WARNING_FLAGS) -O0 -g3 ../demo3d_1.cpp $(OUT08_DEBUG_DIR)/swcme3d.o -o $(OUT08_DEBUG_DIR)/demo3d_1
	$(CXX) $(CPPFLAGS) $(OUT08_WARNING_FLAGS) -O0 -g3 ../demo3d_2.cpp $(OUT08_DEBUG_DIR)/swcme3d.o -o $(OUT08_DEBUG_DIR)/demo3d_2
	$(CXX) $(CPPFLAGS) $(OUT08_WARNING_FLAGS) -O0 -g3 core/test_strict_writer.cpp core/test_strict_writer_main.cpp -o $(OUT08_DEBUG_DIR)/test_out08
	$(CXX) $(CPPFLAGS) $(OUT08_WARNING_FLAGS) -O3 -DNDEBUG -c ../swcme3d.cpp -o $(OUT08_OPTIMIZED_DIR)/swcme3d.o
	$(CXX) $(CPPFLAGS) $(OUT08_WARNING_FLAGS) -O3 -DNDEBUG ../demo1d.cpp -o $(OUT08_OPTIMIZED_DIR)/demo1d
	$(CXX) $(CPPFLAGS) $(OUT08_WARNING_FLAGS) -O3 -DNDEBUG ../demo3d_1.cpp $(OUT08_OPTIMIZED_DIR)/swcme3d.o -o $(OUT08_OPTIMIZED_DIR)/demo3d_1
	$(CXX) $(CPPFLAGS) $(OUT08_WARNING_FLAGS) -O3 -DNDEBUG ../demo3d_2.cpp $(OUT08_OPTIMIZED_DIR)/swcme3d.o -o $(OUT08_OPTIMIZED_DIR)/demo3d_2
	$(CXX) $(CPPFLAGS) $(OUT08_WARNING_FLAGS) -O3 -DNDEBUG core/test_strict_writer.cpp core/test_strict_writer_main.cpp -o $(OUT08_OPTIMIZED_DIR)/test_out08
	./$(OUT08_DEBUG_DIR)/test_out08
	./$(OUT08_OPTIMIZED_DIR)/test_out08

output:
	mkdir -p $@

test: all
	./$(TARGET)

clean:
	rm -f $(TARGET) $(SEP_REFERENCE_TARGET) $(DEMO_TARGETS) $(OBJECTS) output/test_swcme_pst05_asan $(PST08_TARGET) output/bench_tmp \
		$(OUT08_DEBUG_DIR)/swcme3d.o $(OUT08_DEBUG_DIR)/demo1d $(OUT08_DEBUG_DIR)/demo3d_1 $(OUT08_DEBUG_DIR)/demo3d_2 $(OUT08_DEBUG_DIR)/test_out08 \
		$(OUT08_OPTIMIZED_DIR)/swcme3d.o $(OUT08_OPTIMIZED_DIR)/demo1d $(OUT08_OPTIMIZED_DIR)/demo3d_1 $(OUT08_OPTIMIZED_DIR)/demo3d_2 $(OUT08_OPTIMIZED_DIR)/test_out08
