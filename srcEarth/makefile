SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

MAINLIBOBJ = Earth_Sampling.o \
	main_lib.o main_lib_geo.o Earth.o \
	BoundaryInjection.o \
	BoundaryInjection_ELECTRONS.o \
	BoundaryInjection_SEP.o \
	BoundaryInjection_GCR.o \
	BoundaryInjection_SW.o \
	GCR_Badavi2011ASR.o \
	CutoffRigidity.o \
	ImpulseSource.o \
	DomainBoundaryParticlePropertyTable.o \
	parser.o \
	SamplingIndividualPoints.o \
	OutsideParticleFlux.o \
	EnergeticParticlePhysics.o \
	Atmosphere.o \
	pt.o 
	
	
MAINOBJ = main.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)

