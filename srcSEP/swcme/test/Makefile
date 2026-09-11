SHELL := /bin/sh

CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O2

CPPFLAGS += -isystem .. -Icore
CXXFLAGS += -std=c++17 -Wall -Wextra -Wpedantic

TARGET := output/test_swcme
SOURCES := core/test_swcme.cpp core/test_config.cpp 1d/test_ambient.cpp

.PHONY: all test clean

all: $(TARGET)

$(TARGET): $(SOURCES) core/test_framework.hpp ../swcme_constants.hpp ../swcme1d.hpp ../swcme3d.hpp ../swcme3d.cpp | output
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SOURCES) -o $@

output:
	mkdir -p $@

test: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
