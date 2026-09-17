# ==============================================================================
# srcSEP3D/makefile — R0-R2 plus Phase-M/B/T/P/A/O production and tests
#
# Production and dependency-light tests are intentionally separate:
#
#   lib / amps
#       These archive targets are called by the enclosing AMPS build after the
#       application has been copied to AMPS/build/main.  They are not complete
#       standalone source-directory build commands.
#
#   strict-production
#       Delegates to the enclosing `make amps` target so compilation receives
#       the generated headers, include paths, definitions, and link settings
#       owned by AMPS.  It then audits the resulting build/main archives.
#
#   test-stage1 / check-layering
#       Compile the AMPS-independent mesh, background, turbulence, transport,
#       Runtime, neutral adapters, output/restart, and shared SEP test registry.
#       These targets work without AMPS or MPI and are driven by the runner.
#
# The old SEP3D.cpp object is absent by design. main_lib.cpp and the explicit
# AMPS particle adapter are the only AMPS-dependent application-library
# sources; main.cpp supplies the driver archive. R1 shared code is inserted
# from canonical archives, while M/B/T/P/A/O model sources stay AMPS-free.
# ==============================================================================

SHELL := /bin/sh

DEFAULT_TARGET: lib

# The source application lives at AMPS/srcSEP3D, but AMPS copies it to
# AMPS/build/main for compilation.  The build may also invoke this makefile
# through -f from a different current working directory.  Resolve the absolute
# directory of the active makefile first; no configuration or shared-library
# path below depends on AMPS/build/main being the process working directory.
SEP3D_MAKEFILE_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

ifneq ($(wildcard $(SEP3D_MAKEFILE_DIR)/../Makefile.conf),)
AMPS_ROOT ?= $(abspath $(SEP3D_MAKEFILE_DIR)/..)
else ifneq ($(wildcard $(SEP3D_MAKEFILE_DIR)/../../Makefile.conf),)
AMPS_ROOT ?= $(abspath $(SEP3D_MAKEFILE_DIR)/../..)
else
AMPS_ROOT ?= $(abspath $(SEP3D_MAKEFILE_DIR)/..)
endif

AMPS_CONFIG ?= $(AMPS_ROOT)/Makefile.conf
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
SWCME_DIR ?= $(AMPS_ROOT)/src/models/swcme
SWCME_ARCHIVE ?= $(SWCME_DIR)/swcme.a
SWCME_OBJECTS := $(SWCME_DIR)/build/swcme3d.o

# Dependency-light goals must remain usable in a source-only archive.  All
# production goals deliberately fail during makefile parsing if AMPS_CONFIG is
# missing, because a local mock configuration is not production evidence.
STANDALONE_GOALS := check-layering test/stage1 test-stage1 test-harn-unit \
                    test-lay-unit test-bld-unit test-util-unit test-life-unit \
                    test-mesh-unit test-background-unit test-turbulence-unit \
                    test-transport-unit test-adapter-unit test-output-unit \
                    clean clean-standalone help print-layout-paths
ifneq ($(filter $(STANDALONE_GOALS),$(MAKECMDGOALS)),)
-include $(AMPS_CONFIG)
else
include $(AMPS_CONFIG)
endif

CXX ?= g++

# Makefile.conf owns most production flags, but shared model headers are not
# part of every historic AMPS include list.  Use absolute paths resolved from
# the active makefile so these work identically from srcSEP3D and build/main.
MODEL_INCLUDE_FLAGS := -I$(SEP_COMMON_DIR) -I$(SWCME_DIR)
CPPFLAGS += $(MODEL_INCLUDE_FLAGS)
CXXFLAGS += $(MODEL_INCLUDE_FLAGS)
# Historic AMPS Makefile.conf revisions name their include accumulator INCLUDE
# rather than CPPFLAGS. Supplying both variables keeps the application portable
# across those revisions; duplicate -I arguments are harmless and deterministic.
INCLUDE += $(MODEL_INCLUDE_FLAGS)

# L0 numerical core.  Headers are not listed here; normal compiler dependency
# tracking in the enclosing AMPS build determines when they trigger rebuilds.
CORE_SRC := mesh/mesh_model.cpp
BG_SRC := \
    background/bg_provider.cpp \
    background/bg_parker.cpp \
    background/bg_swmf.cpp \
    background/background_snapshot.cpp
TURBULENCE_SRC := turbulence/turbulence_models.cpp
TRANSPORT_SRC := \
    transport/keyed_random.cpp \
    transport/time_step.cpp \
    transport/parker_transport.cpp \
    transport/focused_transport.cpp
RUNTIME_SRC := \
    runtime/run_configuration.cpp \
    runtime/runtime.cpp \
    runtime/runtime_adapters.cpp
ADAPTER_SRC := \
    adapters/transport_adapter.cpp \
    adapters/particle_ledger.cpp \
    adapters/swcme_source_adapter.cpp
OUTPUT_SRC := \
    output/sampling.cpp \
    output/publication.cpp \
    output/restart.cpp \
    output/output_coordinator.cpp

# L3 AMPS application sources.  Retired prototype sources must never be added
# here.  BLDL3D02 also rejects them if they reappear elsewhere in production.
AMPS_SRC := main_lib.cpp amps/amps_particle_adapter.cpp
DRIVER_SRC := main.cpp

CORE_OBJ := $(CORE_SRC:.cpp=.o)
BG_OBJ := $(BG_SRC:.cpp=.o)
TURBULENCE_OBJ := $(TURBULENCE_SRC:.cpp=.o)
TRANSPORT_OBJ := $(TRANSPORT_SRC:.cpp=.o)
RUNTIME_OBJ := $(RUNTIME_SRC:.cpp=.o)
ADAPTER_OBJ := $(ADAPTER_SRC:.cpp=.o)
OUTPUT_OBJ := $(OUTPUT_SRC:.cpp=.o)
AMPS_OBJ := $(AMPS_SRC:.cpp=.o)
DRIVER_OBJ := $(DRIVER_SRC:.cpp=.o)

MAINLIBOBJ := $(CORE_OBJ) $(BG_OBJ) $(TURBULENCE_OBJ) $(TRANSPORT_OBJ) \
              $(RUNTIME_OBJ) $(ADAPTER_OBJ) $(OUTPUT_OBJ) $(AMPS_OBJ)
MAINOBJ := $(DRIVER_OBJ)

STAGE1_TEST_SRC := \
    test/stage1.cpp \
    test/individual-test/test_harness.cpp \
    test/individual-test/test_layering.cpp \
    test/individual-test/test_build.cpp \
    test/individual-test/test_kernels.cpp \
    test/individual-test/test_runtime.cpp \
    test/individual-test/test_mesh.cpp \
    test/individual-test/test_background.cpp \
    test/individual-test/test_turbulence.cpp \
    test/individual-test/test_transport.cpp \
    test/individual-test/test_adapters.cpp \
    test/individual-test/test_output.cpp
STAGE1_TEST_OBJ := $(STAGE1_TEST_SRC:.cpp=.o)

SEP_COMMON_INCLUDES := -I$(SEP_COMMON_DIR)
STAGE1_CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -Werror -O2 \
                   -Icore -Ibackground -Iruntime -Imesh -Iturbulence -Itransport \
                   -Iadapters -Ioutput -I$(SWCME_DIR) $(SEP_COMMON_INCLUDES)

.PHONY: lib amps strict-production audit-production-symbols

lib: $(SEP_COMMON_ARCHIVE) $(SWCME_ARCHIVE) $(MAINLIBOBJ)
	# Keep the enclosing AMPS link contract unchanged: common objects are built
	# once in their canonical model directories and inserted into this
	# application's archive. Recreate the archive so ar cannot retain a retired
	# member from a previous application version.
	rm -f mainlib.a
	ar -rc mainlib.a $(MAINLIBOBJ) $(SEP_COMMON_OBJECTS) $(SWCME_OBJECTS)
	@members=$$(ar t mainlib.a); \
	for member in $(notdir $(SEP_COMMON_OBJECTS) $(SWCME_OBJECTS)); do \
	  count=$$(printf '%s\n' "$$members" | awk -v expected="$$member" \
	    '$$0 == expected { found++ } END { print found + 0 }'); \
	  test "$$count" -eq 1 || { \
	    echo "ERROR: shared member $$member occurs $$count times" >&2; exit 1; }; \
	done
	@echo "[srcSEP3D] mainlib.a built from application + shared model objects"

amps: $(MAINOBJ)
	ar -rc main.a $(MAINOBJ)
	@echo "[srcSEP3D] main.a built from: $(MAINOBJ)"

# A source-directory submake does not inherit the include list assembled by the
# top-level AMPS build.  Compiling main_lib.cpp here would therefore test an
# artificial command and fail to find generated build/pic/pic.h.  Delegate to
# the enclosing build instead; this exercises the same build/main copy and
# compiler command used to create the real AMPS executable.
strict-production:
	@$(MAKE) -C "$(AMPS_ROOT)" amps
	@$(MAKE) --no-print-directory -f "$(SEP3D_MAKEFILE_DIR)/makefile" \
	  audit-production-symbols AMPS_ROOT="$(AMPS_ROOT)"

# The enclosing build writes these archives in AMPS/build/main.  Keeping their
# location absolute prevents an audit launched from srcSEP3D, AMPS, or another
# working directory from accidentally inspecting stale local archives.
PRODUCTION_ARCHIVE_DIR ?= $(AMPS_ROOT)/build/main
PRODUCTION_ARCHIVES := $(PRODUCTION_ARCHIVE_DIR)/mainlib.a \
                       $(PRODUCTION_ARCHIVE_DIR)/main.a

audit-production-symbols:
	@command -v nm >/dev/null 2>&1 && command -v ar >/dev/null 2>&1 || { \
	  echo "ERROR: nm and ar are required for the production archive audit." >&2; \
	  exit 2; \
	}
	@for archive in $(PRODUCTION_ARCHIVES); do \
	  test -f "$$archive" || { \
	    echo "ERROR: enclosing AMPS build did not produce $$archive" >&2; \
	    exit 1; \
	  }; \
	done
	@if nm -A $(PRODUCTION_ARCHIVES) | grep -E \
	  'Mover_Axisymmetric_SecondOrder|TotalParticleAcceleration|GlobalEnergyDistribution|inject_particle_onto_field_line' \
	  >/dev/null; then \
	  echo "ERROR: a retired prototype symbol is present in a production archive." >&2; \
	  exit 1; \
	fi
	@members=$$(ar t $(PRODUCTION_ARCHIVE_DIR)/mainlib.a); \
	for member in $(notdir $(SEP_COMMON_OBJECTS) $(SWCME_OBJECTS)); do \
	  count=$$(printf '%s\n' "$$members" | \
	    awk -v expected="$$member" '$$0 == expected { found++ } \
	      END { print found + 0 }'); \
	  if [ "$$count" -ne 1 ]; then \
	    echo "ERROR: shared archive member $$member occurs $$count times in mainlib.a" >&2; \
	    exit 1; \
	  fi; \
	done
	@echo "[srcSEP3D] enclosing AMPS build and production archive audit passed."

.PHONY: sep-common
sep-common:
	$(MAKE) -C $(SEP_COMMON_DIR)

$(SEP_COMMON_ARCHIVE): sep-common

.PHONY: swcme
swcme:
	$(MAKE) -C $(SWCME_DIR)

$(SWCME_ARCHIVE): swcme

test/stage1.o: test/stage1.cpp
	$(CXX) $(STAGE1_CXXFLAGS) -c -o $@ $<

test/individual-test/%.o: test/individual-test/%.cpp
	$(CXX) $(STAGE1_CXXFLAGS) -c -o $@ $<

core/%.o: core/%.cpp
	$(CXX) $(STAGE1_CXXFLAGS) -c -o $@ $<

background/%.o: background/%.cpp
	$(CXX) $(STAGE1_CXXFLAGS) -c -o $@ $<

runtime/%.o: runtime/%.cpp
	$(CXX) $(STAGE1_CXXFLAGS) -c -o $@ $<

mesh/%.o: mesh/%.cpp
	$(CXX) $(STAGE1_CXXFLAGS) -c -o $@ $<

turbulence/%.o: turbulence/%.cpp
	$(CXX) $(STAGE1_CXXFLAGS) -c -o $@ $<

transport/%.o: transport/%.cpp
	$(CXX) $(STAGE1_CXXFLAGS) -c -o $@ $<

adapters/%.o: adapters/%.cpp
	$(CXX) $(STAGE1_CXXFLAGS) -c -o $@ $<

output/%.o: output/%.cpp
	$(CXX) $(STAGE1_CXXFLAGS) -c -o $@ $<

test/stage1: check-layering $(CORE_OBJ) $(BG_OBJ) $(TURBULENCE_OBJ) $(TRANSPORT_OBJ) $(RUNTIME_OBJ) $(ADAPTER_OBJ) $(OUTPUT_OBJ) $(STAGE1_TEST_OBJ) \
             $(SEP_COMMON_ARCHIVE)
	$(CXX) -o $@ $(CORE_OBJ) $(BG_OBJ) $(TURBULENCE_OBJ) $(TRANSPORT_OBJ) $(RUNTIME_OBJ) $(ADAPTER_OBJ) $(OUTPUT_OBJ) $(STAGE1_TEST_OBJ) \
	  $(SEP_COMMON_ARCHIVE)
	@echo "[srcSEP3D] standalone test binary built: test/stage1"

.PHONY: test-stage1 test-harn-unit test-lay-unit test-bld-unit test-util-unit \
        test-life-unit test-mesh-unit test-background-unit test-turbulence-unit \
        test-transport-unit test-adapter-unit test-output-unit
test-stage1: test/stage1
	./test/stage1 --all-tests

test-harn-unit: test/stage1
	./test/stage1 --test-group HARN

test-lay-unit: test/stage1
	./test/stage1 --test-group LAY

test-bld-unit: test/stage1
	./test/stage1 --test-group BLD

test-util-unit: test/stage1
	./test/stage1 --test-group UTIL

test-life-unit: test/stage1
	./test/stage1 --test-group LIFE3D

test-mesh-unit: test/stage1
	./test/stage1 --test-group MSH3D

test-background-unit: test/stage1
	./test/stage1 --test-group BGP3D --test-group SNAP3D

test-turbulence-unit: test/stage1
	./test/stage1 --test-group TUR3D --test COEF3D01 --test COEF3D02

test-transport-unit: test/stage1
	./test/stage1 --test-group PRK3D --test-group FTE3D --test-group RNG3D

test-adapter-unit: test/stage1
	./test/stage1 --test-group ADP3D --test NAT3D04 --test NAT3D05 --test NAT3D08 --test-group SHK3D

test-output-unit: test/stage1
	./test/stage1 --test NAT3D06 --test NAT3D07 --test-group RST3D

.PHONY: check-layering
check-layering:
	@echo "[srcSEP3D] checking that core/background/runtime/mesh/turbulence/transport/adapters/output do not depend on AMPS or MPI"
	@FOUND=0; \
	for dir in core background runtime mesh turbulence transport adapters output; do \
	  if grep -rn -E --include='*.h' --include='*.cpp' \
	    -e '#include.*pic\.h' -e '#include.*mpi\.h' -e '\bPIC::' \
	    $$dir 2>/dev/null; then FOUND=1; fi; \
	done; \
	if [ $$FOUND -ne 0 ]; then \
	  echo "ERROR: srcSEP3D AMPS-independent layering violation." >&2; exit 1; \
	fi

.PHONY: clean clean-production clean-standalone
clean: clean-production clean-standalone

clean-production:
	rm -f $(MAINLIBOBJ) $(MAINOBJ) mainlib.a main.a

clean-standalone:
	rm -f $(CORE_OBJ) $(BG_OBJ) $(TURBULENCE_OBJ) $(TRANSPORT_OBJ) $(RUNTIME_OBJ) $(ADAPTER_OBJ) $(OUTPUT_OBJ) $(STAGE1_TEST_OBJ) test/stage1

.PHONY: help
help:
	@echo "Production:  make lib | make amps | make strict-production"
	@echo "Standalone:  make test-stage1 | make test-life-unit | make test-mesh-unit | make test-background-unit | make test-turbulence-unit | make test-transport-unit | make test-adapter-unit | make test-output-unit | make check-layering"
	@echo "Preferred:   python3 test/run_tests.py --help"

# Diagnostic target used by BLDL3D05 and useful when inspecting a configured
# AMPS build tree.  In particular, it proves that a makefile copied into
# AMPS/build/main resolves the same canonical files as AMPS/srcSEP3D/makefile.
.PHONY: print-layout-paths
print-layout-paths:
	@echo "SEP3D_MAKEFILE_DIR=$(SEP3D_MAKEFILE_DIR)"
	@echo "AMPS_ROOT=$(AMPS_ROOT)"
	@echo "AMPS_CONFIG=$(AMPS_CONFIG)"
	@echo "SEP_COMMON_DIR=$(SEP_COMMON_DIR)"
	@echo "SWCME_DIR=$(SWCME_DIR)"
