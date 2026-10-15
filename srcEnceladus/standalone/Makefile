CXX ?= g++
CXXFLAGS ?= -O2 -Wall -Wextra -pedantic
CSPICE_ROOT ?= /home/vtenishe/SPICE/cspice
BUILD_DIR ?= build

.PHONY: all spice-ground-tracks spice-observer-ranges clean

all: spice-ground-tracks spice-observer-ranges

spice-ground-tracks: $(BUILD_DIR)/spice_ground_tracks

$(BUILD_DIR)/spice_ground_tracks: tools/spice_ground_tracks.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) -std=c++11 $(CXXFLAGS) -I$(CSPICE_ROOT)/include $< \
	  $(CSPICE_ROOT)/lib/cspice.a $(CSPICE_ROOT)/lib/csupport.a -lm -o $@

spice-observer-ranges: $(BUILD_DIR)/spice_observer_ranges

$(BUILD_DIR)/spice_observer_ranges: tools/spice_observer_ranges.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) -std=c++11 $(CXXFLAGS) -I$(CSPICE_ROOT)/include $< \
	  $(CSPICE_ROOT)/lib/cspice.a $(CSPICE_ROOT)/lib/csupport.a -lm -o $@

clean:
	rm -f $(BUILD_DIR)/spice_ground_tracks $(BUILD_DIR)/spice_observer_ranges
