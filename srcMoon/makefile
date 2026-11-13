SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

MAINLIBOBJ = Exosphere_Helium.o Exosphere_Neon.o species/Argon.o species/SodiumCase1.o Moon.o MoonInput.o LunarSurface.o MoonSurfaceEnvironment.o MoonSurfaceTemperature.o MoonThermalModel.o MoonSurfacePhysics.o MoonDivinerObservation.o MoonIntegrationTests.o Moon_SampleVelocityDistribution.o Moon_SampleSubSolarLimbColumnIntegrals.o tvis.Kaguya.o main_lib.o
MAINLIBOBJ += MoonParticleBudget.o
MAINLIBOBJ += MoonM2Budget.o
MAINOBJ = main.o

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)
