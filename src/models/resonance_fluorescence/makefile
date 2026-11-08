SHELL=/bin/sh

DEFAULT_TARGET: resonance_fluorescence

include ../../../Makefile.conf

RESONANCE_FLUORESCENCE_OBJ=resonance_fluorescence.o

resonance_fluorescence: $(RESONANCE_FLUORESCENCE_OBJ)
	${AR} resonance_fluorescence.a $(RESONANCE_FLUORESCENCE_OBJ)
	@echo "Target 'resonance_fluorescence' is done."

test:
	python3 tools/generate_reference_tables.py --check
	PYTHONDONTWRITEBYTECODE=1 python3 tests/test_reference_data.py
	$${CXX:-g++} -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror \
		-I. resonance_fluorescence.cpp tests/test_resonance_fluorescence.cpp \
		-o tests/test_resonance_fluorescence.tmp
	./tests/test_resonance_fluorescence.tmp
	rm -f tests/test_resonance_fluorescence.tmp

.PHONY: resonance_fluorescence test
