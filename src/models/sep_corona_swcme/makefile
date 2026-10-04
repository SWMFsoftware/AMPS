CXX ?= g++
AR ?= ar
PYTHON ?= python3
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic -Werror
CPPFLAGS := -Iinclude -I../sep_coronal_cme/include -I../sep_common
BUILD := build
SOURCES := $(wildcard src/*.cpp)
OBJECTS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SOURCES))
CORONAL_ARCHIVE := ../sep_coronal_cme/build/libsep_coronal_cme.a

.PHONY: all lib test clean coronal-lib check-architecture
.DEFAULT_GOAL := all

all: lib

lib: $(BUILD)/libsep_corona_swcme.a

$(BUILD):
	mkdir -p $@

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD)/libsep_corona_swcme.a: $(OBJECTS)
	$(AR) rcs $@ $^

coronal-lib:
	$(MAKE) -C ../sep_coronal_cme -j16 lib

$(BUILD)/test_bg3d1: test/test_bg3d1.cpp $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d2: test/test_bg3d2.cpp $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d3: test/test_bg3d3.cpp $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

$(BUILD)/test_bg3d4: test/test_bg3d4.cpp $(BUILD)/libsep_corona_swcme.a coronal-lib
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD)/libsep_corona_swcme.a \
		$(CORONAL_ARCHIVE) -o $@

check-architecture: $(BUILD)/libsep_corona_swcme.a
	$(PYTHON) test/check_architecture.py

test: $(BUILD)/test_bg3d1 $(BUILD)/test_bg3d2 $(BUILD)/test_bg3d3 $(BUILD)/test_bg3d4 check-architecture
	$(BUILD)/test_bg3d1
	$(BUILD)/test_bg3d2
	$(BUILD)/test_bg3d3
	$(BUILD)/test_bg3d4

clean:
	rm -rf $(BUILD)

-include $(OBJECTS:.o=.d)
