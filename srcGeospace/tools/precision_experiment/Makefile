# Float-vs-double backtracing experiment (see README.md and ../../doc/PRECISION.md).
#
#   make            build ./precision_experiment
#   make run        run all three movers in parallel, write $(OUT)/out_<MOVER>.txt
#   make clean
#
# Results are written outside srcGeospace: OUT defaults to
# $GEOSPACE_RUNS/results/precision_experiment/run (GEOSPACE_RUNS defaults to
# <repository>/../AMPS_runs/srcGeospace, as for scripts/build_amps.sh); override with
# "make run OUT=<directory>".
#
# -ffp-contract=off keeps the compiler from fusing a*b+c into FMA instructions, so the
# float and double variants differ only by the storage/arithmetic type (and the results
# do not depend on whether -march=native enables FMA).

CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -ffp-contract=off
MOVERS   := BORIS HC4 RK4

GEOSPACE_RUNS ?= $(abspath $(CURDIR)/../../../../AMPS_runs/srcGeospace)
OUT           ?= $(GEOSPACE_RUNS)/results/precision_experiment/run

precision_experiment: precision_experiment.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

$(OUT)/out_%.txt: precision_experiment
	@mkdir -p $(OUT)
	/usr/bin/time -f "wall %e s, max RSS %M kB" ./precision_experiment $* > $@ 2> $(OUT)/time_$*.txt

run: $(addprefix $(OUT)/out_,$(addsuffix .txt,$(MOVERS)))

clean:
	rm -f precision_experiment

.PHONY: run clean
