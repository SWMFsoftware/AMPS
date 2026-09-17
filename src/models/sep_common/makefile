# ==============================================================================
# sep_common/makefile
#
# Builds the dependency-free SEP kernel library shared by srcSEP and srcSEP3D.
# The sources live here, beside both applications, so neither model owns or
# privately recompiles the canonical production implementation.
#
# Public targets:
#   all      build sep_common.a (default)
#   verify   check layering and exact archive membership
#   clean    remove generated objects/dependency files/archive
# ==============================================================================

SHELL := /bin/sh

CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O2

# C++17 is the common-library contract.  Avoid -ffast-math: srcSEP3D's UTIL02
# record intentionally detects last-bit changes in the numerical kernels.
CPPFLAGS += -I.
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

KERNELS := \
  sep_transport_common \
  sep_coefficient_physics \
  sep_coefficient_registry \
  sep_background_snapshot \
  sep_test_registry \
  sep_injection_spectrum \
  sep_species_source

BUILD_DIR := build
OBJECTS := $(addprefix $(BUILD_DIR)/,$(addsuffix .o,$(KERNELS)))
DEPENDENCIES := $(OBJECTS:.o=.d)
ARCHIVE := sep_common.a

.PHONY: all verify clean

# Every application build receives the membership/layering audit, including an
# incremental checkout where the archive itself appears newer than its inputs.
all: verify

$(ARCHIVE): $(OBJECTS)
	# Recreate rather than update in place: ar never removes a member whose
	# source left KERNELS, which would violate the one-definition R1 contract.
	rm -f $@
	$(AR) rcs $@ $(OBJECTS)
	@echo "[sep_common] built $@ from $(words $(OBJECTS)) canonical objects"

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/%.o: %.cpp %.h | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

# SEP_COMMON01 combines two architectural checks.  First, the library must
# remain independent of AMPS/PIC/MPI so source-only model tests can link it.
# Second, stale objects must never survive a kernel-list change in the archive.
verify: $(ARCHIVE)
	@echo "[sep_common] SEP_COMMON01: checking dependency boundary"
	@if grep -nE '#[[:space:]]*include[[:space:]]*[<\"](pic|mpi)\.h[>\"]' \
	    *.h *.cpp; then \
	  echo "ERROR: sep_common must not depend on PIC or MPI" >&2; exit 1; \
	fi
	@if command -v nm >/dev/null 2>&1 && \
	    nm -C $(ARCHIVE) | grep -E ' U (PIC::|MPI_)'; then \
	  echo "ERROR: sep_common contains an unresolved PIC/MPI reference" >&2; \
	  exit 1; \
	fi
	@if command -v nm >/dev/null 2>&1; then \
	  duplicates=$$(nm -g --defined-only $(ARCHIVE) | \
	    awk '$$2 ~ /^[TDBRS]$$/ { count[$$3]++ } END { \
	      for (symbol in count) if (count[symbol] > 1) print symbol }'); \
	  if [ -n "$$duplicates" ]; then \
	    echo "ERROR: sep_common has multiple strong definitions:" >&2; \
	    printf '%s\n' "$$duplicates" >&2; exit 1; \
	  fi; \
	fi
	@members=$$($(AR) t $(ARCHIVE) | sort); \
	expected=$$(printf '%s\n' $(notdir $(OBJECTS)) | sort); \
	if [ "$$members" != "$$expected" ]; then \
	  echo "ERROR: $(ARCHIVE) membership differs from KERNELS" >&2; \
	  echo "archive:" >&2; printf '%s\n' "$$members" >&2; \
	  echo "expected:" >&2; printf '%s\n' "$$expected" >&2; exit 1; \
	fi
	@echo "[sep_common] SEP_COMMON01 PASS: PIC/MPI-free, exact members, unique strong symbols"

clean:
	rm -rf $(BUILD_DIR) $(ARCHIVE)

-include $(DEPENDENCIES)
