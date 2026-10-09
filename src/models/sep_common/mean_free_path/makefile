SHELL := /bin/sh

CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O2
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

BUILD_DIR := build
LIBRARY := $(BUILD_DIR)/libmean_free_path.a
OBJECT := $(BUILD_DIR)/mean_free_path.o
TEST := $(BUILD_DIR)/test_mean_free_path

.PHONY: all test verify clean

all: $(LIBRARY)

$(BUILD_DIR):
	mkdir -p $@

$(OBJECT): mean_free_path.cpp mean_free_path.h | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(LIBRARY): $(OBJECT)
	rm -f $@
	$(AR) rcs $@ $^

$(TEST): test_mean_free_path.cpp $(LIBRARY) mean_free_path.h | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) test_mean_free_path.cpp $(LIBRARY) -o $@

test: $(TEST)
	$(TEST) --bundle MEAN_FREE_PATH_MODEL_DATA
	python3 MEAN_FREE_PATH_MODEL_DATA/reference_verification.py
	python3 MEAN_FREE_PATH_MODEL_DATA/independent_physics_checks.py

verify: test
	@if grep -nE '#[[:space:]]*include[[:space:]]*[<"](pic|mpi)\.h[>"]' \
	    mean_free_path.h mean_free_path.cpp; then \
	  echo "ERROR: mean_free_path must remain independent of PIC and MPI" >&2; exit 1; \
	fi
	@if command -v nm >/dev/null 2>&1 && \
	    nm -C $(LIBRARY) | grep -E ' U (PIC::|MPI_)'; then \
	  echo "ERROR: mean_free_path has a PIC/MPI dependency" >&2; exit 1; \
	fi
	@echo "[mean_free_path] host-neutral build and tests PASS"

clean:
	rm -rf $(BUILD_DIR)

-include $(OBJECT:.o=.d)
