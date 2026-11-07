SHELL = /bin/sh

DEFAULT_TARGET: lib

include ../../Makefile.conf

MAINLIBOBJ = main_lib.o
MAINOBJ = main.o particle_buffer.o particle_layout_t2.o split_merge_t2.o physics_photochemistry.o physics_collisions.o physics_distribution.o physics_idf.o physics_movers.o physics_t4.o physics_t5.o species.o units.o weight_time_t2.o

lib: $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'lib' is done."

amps: $(MAINOBJ)
	ar -rc main.a $(MAINOBJ)
	@echo "Target 'amps' is done."
