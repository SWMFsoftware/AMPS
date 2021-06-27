SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

MAINLIBOBJ = exosphere_module.o  main_lib.o inner_boundary_injection.o outer_boundary_injection.o mesh.o mover.o parker.o scattering.o drift.o sep.o field_line.o sampling.o diffusion.o output.o 
MAINOBJ = main.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)

