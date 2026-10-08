SHELL := /bin/sh

CXX ?= c++
AR ?= ar
CPPFLAGS ?=
CXXFLAGS ?= -O2
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

BUILD_DIR := build
LIBRARY := libparallel_diffusion.a
OBJECT := $(BUILD_DIR)/parallel_diffusion.o
TEST_OBJECT := $(BUILD_DIR)/test_parallel_diffusion.o
TEST_BINARY := $(BUILD_DIR)/test_parallel_diffusion
TEST_REPORT := $(BUILD_DIR)/test-report.json

.PHONY: all verify clean

all: $(LIBRARY)

$(BUILD_DIR):
	mkdir -p $@

$(OBJECT): parallel_diffusion.cpp parallel_diffusion.h | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(TEST_OBJECT): test_parallel_diffusion.cpp parallel_diffusion.h | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(LIBRARY): $(OBJECT)
	rm -f $@
	$(AR) rcs $@ $^

$(TEST_BINARY): $(TEST_OBJECT) $(LIBRARY)
	$(CXX) $(CXXFLAGS) $(TEST_OBJECT) -L. -lparallel_diffusion -o $@

verify: $(TEST_BINARY)
	$(TEST_BINARY) --json $(TEST_REPORT)

clean:
	rm -rf $(BUILD_DIR) $(LIBRARY)

-include $(OBJECT:.o=.d) $(TEST_OBJECT:.o=.d)
