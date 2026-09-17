# ==============================================================================
# Canonical SWCME model-library build
#
# SWCME is shared by srcSEP and srcSEP3D.  This makefile is deliberately kept
# beside the canonical sources so application makefiles consume one compiled
# object instead of rebuilding swcme3d.cpp with application-specific flags.
# The 1-D interface is header-only; swcme3d.cpp is the only production
# translation unit and therefore the only member of swcme.a.
# ==============================================================================

SHELL := /bin/sh

CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O2
AR ?= ar

CPPFLAGS += -I.
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

BUILD_DIR := build
OBJECT := $(BUILD_DIR)/swcme3d.o
DEPENDENCY := $(OBJECT:.o=.d)
ARCHIVE := swcme.a

.PHONY: all verify test clean

# Run the exact-membership audit on incremental application builds as well as
# explicit verification requests.
all: verify

$(BUILD_DIR):
	mkdir -p $@

$(OBJECT): swcme3d.cpp swcme3d.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(ARCHIVE): $(OBJECT)
	# Recreate the archive so a removed historical member cannot survive an
	# incremental build and reintroduce a second SWCME definition at link time.
	rm -f $@
	$(AR) rcs $@ $^
	@echo "[swcme] built $@ from the canonical 3-D implementation"

# R1SW01 guards the single-build contract.  Exact membership matters because
# ar replaces matching members but does not remove stale members after a
# manifest change; accepting one would silently retain obsolete definitions.
verify: $(ARCHIVE)
	@members=$$($(AR) t $(ARCHIVE)); \
	expected=$$(basename $(OBJECT)); \
	if [ "$$members" != "$$expected" ]; then \
	  echo "ERROR: $(ARCHIVE) contains unexpected members: $$members" >&2; \
	  exit 1; \
	fi
	@if command -v nm >/dev/null 2>&1 && \
	    nm -C $(ARCHIVE) | grep -E ' U (PIC::|MPI_)' >/dev/null; then \
	  echo "ERROR: SWCME model library depends on AMPS or MPI" >&2; exit 1; \
	fi
	@if command -v nm >/dev/null 2>&1; then \
	  duplicates=$$(nm -g --defined-only $(ARCHIVE) | \
	    awk '$$2 ~ /^[TDBRS]$$/ { count[$$3]++ } END { \
	      for (symbol in count) if (count[symbol] > 1) print symbol }'); \
	  if [ -n "$$duplicates" ]; then \
	    echo "ERROR: SWCME has multiple strong definitions:" >&2; \
	    printf '%s\n' "$$duplicates" >&2; exit 1; \
	  fi; \
	fi
	@echo "[swcme] R1SW01 PASS: exact members, unique strong symbols, no AMPS/MPI dependency"

test: verify
	@PYTHONDONTWRITEBYTECODE=1 test/run_tests.py --routine

clean:
	rm -rf $(BUILD_DIR) $(ARCHIVE)

-include $(DEPENDENCY)
