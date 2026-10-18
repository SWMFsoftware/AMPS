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

build/%.o: src/%.cpp $(wildcard src/*.h)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf build liouville_cops

.PHONY: clean
