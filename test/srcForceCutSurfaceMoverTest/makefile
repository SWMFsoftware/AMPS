SHELL = /bin/sh
DEFAULT_TARGET: lib
include ../../Makefile.conf
MAINLIBOBJ = main_lib.o
MAINOBJ = main.o
lib: $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
amps: $(MAINOBJ)
	ar -rc main.a $(MAINOBJ)
