SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

# Dimension=1 replaces these two object lists during input preprocessing.
MAINLIBOBJ = Gravity.o Comet.o CometData.o CometModel.o CometChemistry.o CometPhysics.o main_lib.o OrbitCalculation.o ColumnIntegral.o SurfaceDataCG.o SampleSizeDistribution.o RadialMesh.o RadialMeshProfile.o RadialSimulation.o
MAINOBJ = main.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)
