# Pitch-angle diffusion coefficient library (Release A).
# Standalone C++17 build; no PIC, MPI or application headers.  srcSEP3D
# flattens the library objects into its mainlib.a (srcSEP3D/makefile).
#
#   make            build libpitch_angle_diffusion.a
#   make tests      build the unit, acceptance and driver binaries
#   make verify     python3 run_all_tests.py (every library test)
#   make clean
SHELL := /bin/sh

CXX ?= c++
AR ?= ar
CPPFLAGS ?=
CXXFLAGS ?= -O2
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

BUILD_DIR := build
LIBRARY := libpitch_angle_diffusion.a
OBJECTS := $(BUILD_DIR)/pitch_angle_diffusion.o $(BUILD_DIR)/pitch_angle_kernels.o \
           $(BUILD_DIR)/pitch_angle_certification.o
UNIT_TEST := $(BUILD_DIR)/test_pitch_angle_diffusion
ACCEPTANCE_TEST := $(BUILD_DIR)/test_pitch_angle_montecarlo
DRIVER := $(BUILD_DIR)/pitch_angle_driver
HEADERS := pitch_angle_diffusion.h pitch_angle_kernels.h pitch_angle_quadrature.h \
           pitch_angle_registry.inc pitch_angle_certified_pairs.inc

.PHONY: all tests verify clean

all: $(LIBRARY)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/%.o: %.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIBRARY): $(OBJECTS)
	rm -f $@
	$(AR) rcs $@ $^

$(UNIT_TEST): test_pitch_angle_diffusion.cpp $(LIBRARY) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -L. -lpitch_angle_diffusion -o $@

# The acceptance campaign compiles the library sources into the test with
# -O3 -flto (no -march/-ffast-math, so the floating-point semantics equal the
# production build); it is a throughput optimization of the inner loop only.
$(ACCEPTANCE_TEST): test_pitch_angle_montecarlo.cpp pitch_angle_diffusion.cpp \
                    pitch_angle_kernels.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -O3 -flto -pthread test_pitch_angle_montecarlo.cpp \
	    pitch_angle_diffusion.cpp pitch_angle_kernels.cpp -o $@

$(DRIVER): pitch_angle_driver.cpp $(LIBRARY) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -L. -lpitch_angle_diffusion -o $@

tests: $(UNIT_TEST) $(ACCEPTANCE_TEST) $(DRIVER)

verify: tests
	python3 run_all_tests.py

clean:
	rm -rf $(BUILD_DIR) $(LIBRARY)
