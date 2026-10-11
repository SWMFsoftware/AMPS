# Generic AMPS runtime input reader (header-only; see src/INPUT_FILE.md).
# This makefile only builds and runs the dependency-free grammar test; the
# header itself is compiled by whichever application or model includes it.
#
#   make -C src/general/input test
CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic -Werror
BUILD := build

.PHONY: test clean
test: $(BUILD)/test_amps_input
	./$(BUILD)/test_amps_input

$(BUILD)/test_amps_input: test_amps_input.cpp amps_input.h
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -I. test_amps_input.cpp -o $@

clean:
	rm -rf -- $(BUILD)
