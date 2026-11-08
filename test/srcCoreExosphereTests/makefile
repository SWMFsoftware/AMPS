SHELL = /bin/sh

DEFAULT_TARGET: lib

include ../../Makefile.conf

MAINLIBOBJ = main_lib.o
MAINOBJ = main.o

lib: $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'lib' is done."

amps: $(MAINOBJ)
	ar -rc main.a $(MAINOBJ)
	@echo "Target 'amps' is done."
