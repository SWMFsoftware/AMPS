CXX ?= g++
AR ?= ar
PYTHON ?= python3

# The shared model must remain independently compilable.  An AMPS include
# path here would defeat ARCHSCCM01 and obscure an accidental dependency.
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic -Werror
CPPFLAGS := -Iinclude -I../sep_common -Isrc
BUILD := build
LIB_SOURCES := $(wildcard src/*.cpp)
LIB_OBJECTS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(LIB_SOURCES))
TEST_SOURCES := $(wildcard test/*.cpp)
TEST_OBJECTS := $(patsubst test/%.cpp,$(BUILD)/test_%.o,$(TEST_SOURCES))

.PHONY: all lib test test-stage0 test-stage1 test-stage2 check-architecture clean

all: lib

lib: $(BUILD)/libsep_coronal_cme.a

$(BUILD):
	mkdir -p $@

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD)/test_%.o: test/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) -Itest $(CXXFLAGS) -c $< -o $@

$(BUILD)/libsep_coronal_cme.a: $(LIB_OBJECTS)
	$(AR) rcs $@ $^

$(BUILD)/sep_coronal_cme_tests: $(TEST_OBJECTS) $(BUILD)/libsep_coronal_cme.a
	$(CXX) $(CXXFLAGS) $(TEST_OBJECTS) $(BUILD)/libsep_coronal_cme.a -o $@

check-architecture: $(BUILD)/libsep_coronal_cme.a
	$(PYTHON) tools/check_architecture.py

test-stage0: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 0

test-stage1: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 1

test-stage2: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --stage 2

test: $(BUILD)/sep_coronal_cme_tests
	$(PYTHON) test/run_tests.py --all

clean:
	rm -rf $(BUILD)
