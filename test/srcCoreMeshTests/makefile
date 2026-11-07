SHELL = /bin/sh

DEFAULT_TARGET: lib

include ../../Makefile.conf

MAINLIBOBJ = main_lib.o
MAINOBJ = main.o mesh_test_fixture.o mesh_topology.o mesh_adaptation.o mesh_coarsen_data.o pic_interpolation.o mesh_cutcell.o

lib: $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'lib' is done."

amps: $(MAINOBJ)
	ar -rc main.a $(MAINOBJ)
	@echo "Target 'amps' is done."
