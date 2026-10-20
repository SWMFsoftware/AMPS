# liouville_cops_cuda: the unchanged LiouvilleCOPS forward model with its two expensive loops
# on the GPU (ROADMAP C11).  See README.md in this directory.
#
#   make                         build ./liouville_cops_cuda
#   make CSPICE=... CUDA=...     other toolkit locations
#   make clean
#
# The reference sources ../../LiouvilleCOPS/src/*.cpp are compiled UNCHANGED with the
# reference flags (../../LiouvilleCOPS/Makefile).  In two translation units the CPU
# definition of one member function is renamed so that the GPU definition in
# gpu_liouville.cu provides the symbol main() calls:
#   Liouville.cpp: cLiouville::EvaluateKernels -> EvaluateKernelsCPU
#   Source.cpp:    cSourceModel::Evaluate      -> EvaluateCPU
# (Source.cpp's TU does not include Geometry.h, so no other "Evaluate" is affected.)

REF      ?= ../../LiouvilleCOPS
CSPICE   ?= /home/vtenishe/SPICE/cspice
CUDA     ?= /usr/local/cuda-12.4
CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -fopenmp -Wall -Wextra -Wno-unused-parameter
NVCC      = $(CUDA)/bin/nvcc
ARCH     ?= sm_86
NVFLAGS   = -O3 -std=c++17 -arch=$(ARCH) --fmad=false -Xcompiler "-O3 -march=native -fopenmp"
INCLUDES  = -I$(CSPICE)/include -I$(REF)/src
LIBS      = $(CSPICE)/lib/cspice.a -L$(CUDA)/lib64 -lcudart_static -ldl -lrt -lpthread -lm -fopenmp

REF_SRC  = $(wildcard $(REF)/src/*.cpp)
REF_OBJ  = $(patsubst $(REF)/src/%.cpp,build/ref/%.o,$(REF_SRC))
OBJ      = $(REF_OBJ) build/access_shim.o build/gpu_liouville.o
HDR      = $(wildcard $(REF)/src/*.h)

liouville_cops_cuda: $(OBJ)
	$(CXX) -o $@ $(OBJ) $(LIBS)

build/ref/Liouville.o: $(REF)/src/Liouville.cpp $(HDR)
	@mkdir -p build/ref
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DEvaluateKernels=EvaluateKernelsCPU -c $< -o $@

build/ref/Source.o: $(REF)/src/Source.cpp $(HDR)
	@mkdir -p build/ref
	$(CXX) $(CXXFLAGS) $(INCLUDES) -DEvaluate=EvaluateCPU -c $< -o $@

build/ref/%.o: $(REF)/src/%.cpp $(HDR)
	@mkdir -p build/ref
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

build/access_shim.o: access_shim.cpp access_shim.h $(HDR)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

build/gpu_liouville.o: gpu_liouville.cu access_shim.h $(HDR)
	@mkdir -p build
	$(NVCC) $(NVFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf build liouville_cops_cuda

.PHONY: clean
