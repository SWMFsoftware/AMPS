SHELL := /bin/sh

CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O2

CPPFLAGS += -isystem .. -Icore
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

TARGET := output/test_swcme
SEP_REFERENCE_TARGET := output/sep_reference
# Every validation remains in one deterministic executable. The common tests
# use separate translation units so each physical contract has an auditable
# implementation without introducing another executable or build system.
SOURCES := core/test_swcme.cpp core/test_state_immutability.cpp core/test_prepared_state.cpp core/test_configuration_state.cpp core/test_state_integrity.cpp core/test_concurrent_state.cpp core/test_config.cpp core/test_units.cpp core/test_density.cpp core/test_kinematics.cpp core/test_shock.cpp core/test_regions.cpp core/test_acceleration.cpp core/test_status.cpp core/test_divergence.cpp core/test_defaults.cpp core/test_sep_interface.cpp 3d/test_parker.cpp 3d/test_geometry.cpp 3d/test_mesh.cpp 3d/test_connectivity.cpp 3d/test_shock_surface.cpp 1d/test_ambient.cpp test_1d3d.cpp
OBJECTS := output/swcme3d.o output/test_swcme.o output/test_state_immutability.o output/test_prepared_state.o output/test_configuration_state.o output/test_state_integrity.o output/test_concurrent_state.o output/test_config.o output/test_units.o output/test_density.o output/test_kinematics.o output/test_shock.o output/test_regions.o output/test_acceleration.o output/test_status.o output/test_divergence.o output/test_defaults.o output/test_sep_interface.o output/test_parker.o output/test_geometry.o output/test_mesh.o output/test_connectivity.o output/test_shock_surface.o output/test_ambient.o output/test_1d3d.o

# POSIX thread support is isolated behind one overridable variable.  PST04's
# translation unit needs the compile flag on implementations that expose
# std::thread through pthreads, and the final executable needs the link flag.
THREAD_FLAGS ?= -pthread

.PHONY: all test tools campaign-test clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(THREAD_FLAGS) $(OBJECTS) -o $@

# Both dimensional model headers include the shared explicit-field integrity
# serializer.  Make every validation object observe changes to that header so
# an incremental build cannot link stale inline seal logic from an older test.
$(OBJECTS): ../swcme_prepared_integrity.hpp

# Compile the production 3-D implementation exactly once and link it normally.
# Earlier validation revisions included swcme3d.cpp from test_config.cpp merely
# to obtain its symbols, which coupled the entire production implementation to
# one unrelated test translation unit.  A normal object-library style link is
# clearer, scales to additional tests, and prevents accidental duplicate code.
output/swcme3d.o: ../swcme3d.cpp ../swcme3d.hpp ../swcme_constants.hpp ../swcme_units.hpp ../swcme_defaults.hpp ../swcme_status.hpp ../swcme_divergence.hpp ../swcme_config.hpp ../swcme_kinematics.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_shock.hpp ../swcme_regions.hpp ../swcme_acceleration.hpp ../swcme_sep_source.hpp | output
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

output/test_config.o: core/test_config.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme_units.hpp ../swcme_defaults.hpp ../swcme_status.hpp ../swcme_divergence.hpp ../swcme_config.hpp ../swcme_kinematics.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_regions.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output/test_units.o: core/test_units.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme_units.hpp ../swcme_defaults.hpp ../swcme_status.hpp ../swcme_divergence.hpp ../swcme_config.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_regions.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# DEN01 calls the production 3-D implementation linked through swcme3d.o and
# keeps its analytical Leblanc reference independent of production code.
output/test_density.o: core/test_density.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_regions.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# KIN01-KIN08 validate the common apex-kinematics engine directly and through
# both dimensional wrappers.  The test contains an independent closed-form DBM
# reference and dense PCHIP monotonicity checks; production kinematics remain in
# swcme_kinematics.hpp so the 1-D and 3-D implementations cannot diverge.
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

# CON01-CON08 validate analytical Parker-line mapping and cobpoint behavior
# against the production swcme3d.o implementation.
output/test_connectivity.o: 3d/test_connectivity.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme3d.hpp ../swcme_status.hpp ../swcme_divergence.hpp | output
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

# Campaign-manager unit tests exercise profile expansion, convergence fitting,
# sweep/comparison execution, exit codes, manifest generation, and machine-
# readable reporting without duplicating production physics.
campaign-test: all tools
	python3 python/test_run_tests.py

output:
	mkdir -p $@

test: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(SEP_REFERENCE_TARGET) $(OBJECTS)
