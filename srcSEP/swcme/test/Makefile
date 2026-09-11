SHELL := /bin/sh

CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O2

CPPFLAGS += -isystem .. -Icore
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

TARGET := output/test_swcme
# Every validation remains in one deterministic executable. The common tests
# use separate translation units so each physical contract has an auditable
# implementation without introducing another executable or build system.
SOURCES := core/test_swcme.cpp core/test_config.cpp core/test_units.cpp core/test_density.cpp core/test_kinematics.cpp core/test_shock.cpp 3d/test_parker.cpp 3d/test_geometry.cpp 3d/test_connectivity.cpp 1d/test_ambient.cpp test_1d3d.cpp
OBJECTS := output/swcme3d.o output/test_swcme.o output/test_config.o output/test_units.o output/test_density.o output/test_kinematics.o output/test_shock.o output/test_parker.o output/test_geometry.o output/test_connectivity.o output/test_ambient.o output/test_1d3d.o

.PHONY: all test clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

# Compile the production 3-D implementation exactly once and link it normally.
# Earlier validation revisions included swcme3d.cpp from test_config.cpp merely
# to obtain its symbols, which coupled the entire production implementation to
# one unrelated test translation unit.  A normal object-library style link is
# clearer, scales to additional tests, and prevents accidental duplicate code.
output/swcme3d.o: ../swcme3d.cpp ../swcme3d.hpp ../swcme_constants.hpp ../swcme_units.hpp ../swcme_config.hpp ../swcme_kinematics.hpp ../swcme_solarwind.hpp ../swcme_core.hpp ../swcme_shock.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Explicit rules keep source ownership obvious and let `make -j` build the
# independent validation units concurrently without introducing another build
# system. Header dependencies force the relevant objects to refresh.
output/test_swcme.o: core/test_swcme.cpp core/test_framework.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output/test_config.o: core/test_config.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme_units.hpp ../swcme_config.hpp ../swcme_kinematics.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output/test_units.o: core/test_units.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme_units.hpp ../swcme_config.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# DEN01 calls the production 3-D implementation linked through swcme3d.o and
# keeps its analytical Leblanc reference independent of production code.
output/test_density.o: core/test_density.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# KIN01-KIN08 validate the common apex-kinematics engine directly and through
# both dimensional wrappers.  The test contains an independent closed-form DBM
# reference and dense PCHIP monotonicity checks; production kinematics remain in
# swcme_kinematics.hpp so the 1-D and 3-D implementations cannot diverge.
output/test_kinematics.o: core/test_kinematics.cpp core/test_framework.hpp ../swcme_kinematics.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# SHK01-SHK12 exercise the shared production ideal-MHD jump solver and its
# integration into the 1-D/3-D models.  The test file contains independent
# limiting-case references and conservation residual checks; it does not
# duplicate the production nonlinear solver.
output/test_shock.o: core/test_shock.cpp core/test_framework.hpp ../swcme_shock.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# PAR01-PAR03 exercise the public production 3-D field evaluator linked through
# swcme3d.o; the test translation unit contains only independent references.
output/test_parker.o: 3d/test_parker.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


# GEO01-GEO08 validate finite shock geometry using the production swcme3d.o
# implementation and independent geometric references.
output/test_geometry.o: 3d/test_geometry.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


# CON01-CON08 validate analytical Parker-line mapping and cobpoint behavior
# against the production swcme3d.o implementation.
output/test_connectivity.o: 3d/test_connectivity.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output/test_ambient.o: 1d/test_ambient.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# 1D3D01-02 validate that the two public dimensional wrappers reduce to the
# same common solar-wind, apex-kinematics, and shock state in an exactly
# equivalent spherical/equatorial geometry.
output/test_1d3d.o: test_1d3d.cpp core/test_framework.hpp ../swcme_core.hpp ../swcme_solarwind.hpp ../swcme_shock.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme_solarwind.hpp ../swcme_core.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output:
	mkdir -p $@

test: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(OBJECTS)
