# LiouvilleCOPS build
#   make            optimized build
#   make CSPICE=... path to the CSPICE toolkit

CSPICE   ?= /home/vtenishe/SPICE/cspice
CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -fopenmp -Wall -Wextra -Wno-unused-parameter
INCLUDES  = -I$(CSPICE)/include
LIBS      = $(CSPICE)/lib/cspice.a -lm -fopenmp

SRC = $(wildcard src/*.cpp)
OBJ = $(patsubst src/%.cpp,build/%.o,$(SRC))

liouville_cops: $(OBJ)
	$(CXX) -o $@ $(OBJ) $(LIBS)

# Small provenance helper used when a surface-slope template needs the
# centrifugal acceleration.  It shares Geometry.cpp with the forward solver
# instead of implementing a second SPICE frame-rate calculation.
spice_angular_velocity: tools/spice_angular_velocity.cpp src/Geometry.cpp src/Input.cpp $(wildcard src/*.h)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Isrc -o $@ tools/spice_angular_velocity.cpp src/Geometry.cpp src/Input.cpp $(LIBS)

build/%.o: src/%.cpp $(wildcard src/*.h)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf build liouville_cops spice_angular_velocity

.PHONY: clean
