CXX ?= g++
.DEFAULT_GOAL := all
AR ?= ar
PYTHON ?= python3

# The shared model must remain independently compilable.  An AMPS include
# path here would defeat ARCHSCCM01 and obscure an accidental dependency.
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic -Werror
CPPFLAGS := -Iinclude -I../sep_common -Isrc
BUILD := build
LIB_SOURCES := $(wildcard src/*.cpp)
LIB_OBJECTS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(LIB_SOURCES))
COMMON_SOURCES := ../sep_common/sep_field_line_exchange.cpp \
	../sep_common/sep_field_line_bundle_io.cpp ../sep_common/sep_coherent_transport.cpp
COMMON_OBJECTS := $(BUILD)/common_sep_field_line_exchange.o \
	$(BUILD)/common_sep_field_line_bundle_io.o $(BUILD)/common_sep_coherent_transport.o
TEST_SOURCES := $(wildcard test/*.cpp)
TEST_OBJECTS := $(patsubst test/%.cpp,$(BUILD)/test_%.o,$(TEST_SOURCES))
ADAPTER_OBJECTS := $(BUILD)/adapter_field_line_export.o \
	$(BUILD)/adapter_field_line_bundle.o
# Public structs can change ABI even when a .cpp timestamp does not. Track
# transitive header dependencies so incremental builds never combine old test
# layouts with new model objects (also important for overlay source updates).
DEPENDENCY_FLAGS := -MMD -MP
-include $(LIB_OBJECTS:.o=.d) $(COMMON_OBJECTS:.o=.d) \
    $(TEST_OBJECTS:.o=.d) $(ADAPTER_OBJECTS:.o=.d)

.PHONY: all lib test test-stage0 test-stage1 test-stage2 test-stage3 \
	test-stage4 test-stage5 test-stage6 test-stage7 test-stage8 test-stage9 \
	test-stage10 test-stage11 test-stage12 test-stage13 test-stage14 \
	check-architecture check-adapters clean

all: lib

lib: $(BUILD)/libsep_coronal_cme.a

$(BUILD):
	mkdir -p $@

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPENDENCY_FLAGS) -c $< -o $@

$(BUILD)/common_%.o: ../sep_common/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPENDENCY_FLAGS) -c $< -o $@

$(BUILD)/test_%.o: test/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) -Itest $(CXXFLAGS) $(DEPENDENCY_FLAGS) -c $< -o $@

# Compile the two application boundaries without linking them into the shared
# archive.  This proves that the adapters consume only the installed public
# APIs while preserving the architecture rule that srcSEP never links SCCM.
$(BUILD)/adapter_field_line_export.o: \
		../../../srcSEP3D/adapters/field_line_export_adapter.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) -I../../../srcSEP3D/adapters $(CXXFLAGS) $(DEPENDENCY_FLAGS) -c $< -o $@

$(BUILD)/adapter_field_line_bundle.o: \
		../../../srcSEP/adapters/field_line_bundle_adapter.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) -I../../../srcSEP/adapters $(CXXFLAGS) $(DEPENDENCY_FLAGS) -c $< -o $@

$(BUILD)/libsep_coronal_cme.a: $(LIB_OBJECTS) $(COMMON_OBJECTS)
	$(AR) rcs $@ $^

$(BUILD)/sep_coronal_cme_tests: $(TEST_OBJECTS) $(BUILD)/libsep_coronal_cme.a \
		$(BUILD)/adapter_field_line_bundle.o
	$(CXX) $(CXXFLAGS) $(TEST_OBJECTS) $(BUILD)/adapter_field_line_bundle.o \
		$(BUILD)/libsep_coronal_cme.a -o $@

check-architecture: $(BUILD)/libsep_coronal_cme.a
	$(PYTHON) test/test_architecture.py

check-adapters: $(ADAPTER_OBJECTS)

test-stage0: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 0 --expect-count 13

test-stage1: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 1 --expect-count 22

test-stage2: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 2 --expect-count 53

test-stage3: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 3 --expect-count 82

test-stage4: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 4 --expect-count 96

test-stage5: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 5 --expect-count 106

test-stage6: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 6 --expect-count 117

test-stage7: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 7 --expect-count 128

test-stage8: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 8 --expect-count 153

test-stage9: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 9 --expect-count 173

test-stage10: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 10 --expect-count 196

test-stage11: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 11 --expect-count 202

test-stage12: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 12 --expect-count 206

# Passing this software gate verifies release machinery. Actual qualification
# additionally consumes the explicit production evidence profile and D1--D10.
test-stage13: $(BUILD)/sep_coronal_cme_tests check-adapters
	$(PYTHON) test/run_tests.py --stage 13 --expect-count 209

test-stage14: $(BUILD)/sep_coronal_cme_tests check-adapters
	$(PYTHON) test/run_tests.py --stage 14 --expect-count 222

test: $(BUILD)/sep_coronal_cme_tests check-adapters
	$(PYTHON) test/run_tests.py --all --expect-count 222

clean:
	rm -rf $(BUILD)
