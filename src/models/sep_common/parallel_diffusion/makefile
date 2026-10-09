SHELL := /bin/sh

CXX ?= c++
AR ?= ar
CPPFLAGS ?=
CXXFLAGS ?= -O2
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

BUILD_DIR := build
LIBRARY := libparallel_diffusion.a
OBJECTS := $(BUILD_DIR)/parallel_diffusion.o $(BUILD_DIR)/parallel_diffusion_advanced.o
TEST_OBJECT := $(BUILD_DIR)/test_parallel_diffusion.o
TEST_BINARY := $(BUILD_DIR)/test_parallel_diffusion
ADVANCED_TEST_DRIVER := $(BUILD_DIR)/test_parallel_diffusion_advanced_driver
TEST_REPORT := $(BUILD_DIR)/test-report.json

.PHONY: all verify clean

all: $(LIBRARY)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/parallel_diffusion.o: parallel_diffusion.cpp parallel_diffusion.h parallel_diffusion_advanced.h | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/parallel_diffusion_advanced.o: parallel_diffusion_advanced.cpp parallel_diffusion.h parallel_diffusion_advanced.h nlgce_f_coefficients.inc | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(TEST_OBJECT): test_parallel_diffusion.cpp parallel_diffusion.h | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(LIBRARY): $(OBJECTS)
	rm -f $@
	$(AR) rcs $@ $^

$(TEST_BINARY): $(TEST_OBJECT) $(LIBRARY)
	$(CXX) $(CXXFLAGS) $(TEST_OBJECT) -L. -lparallel_diffusion -o $@

$(ADVANCED_TEST_DRIVER): test_parallel_diffusion_advanced_driver.cpp $(LIBRARY) parallel_diffusion.h | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -pthread $< -L. -lparallel_diffusion -o $@

verify: $(TEST_BINARY) $(ADVANCED_TEST_DRIVER)
	python3 generate_nlgce_coefficients.py --check
	$(TEST_BINARY) --json $(TEST_REPORT)
	python3 test_parallel_diffusion_advanced.py --driver $(ADVANCED_TEST_DRIVER) --json $(BUILD_DIR)/advanced-test-report.json

clean:
	rm -rf $(BUILD_DIR) $(LIBRARY)

-include $(OBJECTS:.o=.d) $(TEST_OBJECT:.o=.d)
