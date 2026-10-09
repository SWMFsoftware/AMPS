SHELL := /bin/sh

CXX ?= c++
AR ?= ar
CPPFLAGS ?=
CXXFLAGS ?= -O2
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

BUILD_DIR := build
LIBRARY := libperpendicular_diffusion.a
OBJECTS := $(BUILD_DIR)/perpendicular_diffusion.o \
           $(BUILD_DIR)/perpendicular_diffusion_models.o
TEST_BINARY := $(BUILD_DIR)/test_perpendicular_diffusion
PARALLEL_DIR := ../parallel_diffusion
PARALLEL_LIBRARY := $(PARALLEL_DIR)/libparallel_diffusion.a

.PHONY: all verify verify-high-precision clean

all: $(LIBRARY)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/%.o: %.cpp perpendicular_diffusion.h perpendicular_diffusion_internal.h | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -I$(PARALLEL_DIR) -c $< -o $@

$(LIBRARY): $(OBJECTS)
	rm -f $@
	$(AR) rcs $@ $^

$(PARALLEL_LIBRARY):
	$(MAKE) -C $(PARALLEL_DIR) all

$(TEST_BINARY): test_perpendicular_diffusion.cpp $(LIBRARY) $(PARALLEL_LIBRARY) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(PARALLEL_DIR) $< \
		-L. -lperpendicular_diffusion -L$(PARALLEL_DIR) -lparallel_diffusion -o $@

verify: $(TEST_BINARY)
	$(TEST_BINARY)
	cd perpendicular_diffusion_model_data && python3 reference_verification.py
	cd perpendicular_diffusion_model_data && python3 paired_reference_verification.py
	@if python3 -c 'import mpmath' >/dev/null 2>&1; then \
	  cd perpendicular_diffusion_model_data && python3 review_additions_verification.py; \
	else \
	  echo "SKIP review_additions_verification.py: optional Python dependency mpmath is unavailable"; \
	fi
	cd perpendicular_diffusion_model_data && python3 document_fixtures.py --document ../PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL.md
	cd perpendicular_diffusion_model_data && python3 document_audit.py --document ../PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL.md

# Release/maintenance environments that provide NumPy, SciPy and mpmath use
# this target when a skipped high-precision audit is not acceptable.
verify-high-precision:
	cd perpendicular_diffusion_model_data && python3 review_additions_verification.py

clean:
	rm -rf $(BUILD_DIR) $(LIBRARY)

-include $(OBJECTS:.o=.d)
