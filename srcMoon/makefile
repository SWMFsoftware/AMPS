SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

MAINLIBOBJ = Exosphere_Helium.o Exosphere_Neon.o Moon.o Moon_SampleVelocityDistribution.o Moon_SampleSubSolarLimbColumnIntegrals.o tvis.Kaguya.o main_lib.o 
MAINOBJ = main.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)

