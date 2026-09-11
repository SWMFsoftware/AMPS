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
SOURCES := core/test_swcme.cpp core/test_config.cpp core/test_units.cpp core/test_density.cpp core/test_shock.cpp 3d/test_parker.cpp 3d/test_geometry.cpp 1d/test_ambient.cpp
OBJECTS := output/test_swcme.o output/test_config.o output/test_units.o output/test_density.o output/test_shock.o output/test_parker.o output/test_geometry.o output/test_ambient.o

.PHONY: all test clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

# Explicit rules keep source ownership obvious and let `make -j` build the
# independent validation units concurrently without introducing another build
# system. Header dependencies force the relevant objects to refresh.
output/test_swcme.o: core/test_swcme.cpp core/test_framework.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output/test_config.o: core/test_config.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme3d.cpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output/test_units.o: core/test_units.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# DEN01 links against the production 3-D implementation already compiled into
# test_config.o; this rule compiles only its independent reference and calls to
# the public 1-D/3-D APIs, preventing duplicate production symbol definitions.
output/test_density.o: core/test_density.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# SHK01-SHK12 exercise the shared production ideal-MHD jump solver and its
# integration into the 1-D/3-D models.  The test file contains independent
# limiting-case references and conservation residual checks; it does not
# duplicate the production nonlinear solver.
output/test_shock.o: core/test_shock.cpp core/test_framework.hpp ../swcme_shock.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# PAR01-PAR03 exercise the public production 3-D field evaluator.  The current
# test harness obtains the swcme3d.cpp production symbols from test_config.o;
# this separate test translation unit contains only independent references and
# calls through swcme3d.hpp, so the production implementation is not duplicated.
output/test_parker.o: 3d/test_parker.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


# GEO01-GEO08 validate the corrected finite shock geometry, exact analytical
# normals, strict angular support, and self-similar normal speed.  As with the
# Parker tests, the production swcme3d.cpp symbols are linked through the
# existing test_config.o arrangement; this file contains only independent
# geometric references and public-API calls.
output/test_geometry.o: 3d/test_geometry.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme3d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output/test_ambient.o: 1d/test_ambient.cpp core/test_framework.hpp ../swcme_constants.hpp ../swcme1d.hpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

output:
	mkdir -p $@

test: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(OBJECTS)
