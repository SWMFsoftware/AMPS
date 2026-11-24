SHELL =/bin/sh

DEFAULT_TARGET : lib  

include ../../Makefile.conf

MAINLIBOBJ = Exosphere_Helium.o Exosphere_Neon.o species/Argon.o species/SodiumCase1.o Moon.o MoonInput.o LunarSurface.o MoonSurfaceEnvironment.o MoonSurfaceTemperature.o MoonThermalModel.o MoonSurfacePhysics.o MoonDivinerObservation.o MoonIntegrationTests.o Moon_SampleVelocityDistribution.o Moon_SampleSubSolarLimbColumnIntegrals.o tvis.Kaguya.o main_lib.o
MAINLIBOBJ += MoonParticleBudget.o
MAINLIBOBJ += MoonM2Budget.o MoonM2Cohorts.o
MAINLIBOBJ += MoonPolarDiagnostics.o
MAINLIBOBJ += species/ElectronImpact.o
# The normal application owns these same rate/transaction/tracing services
# exercised by the thin probes. Linking their production objects here avoids
# a test-only executor or parser implementation; runtime guards still keep
# unsupported daughter kinetics and legacy receipts from enabling reactions.
MAINLIBOBJ += species/Photochemistry.o species/MoonChargeExchange.o species/ReactionExecution.o
MAINLIBOBJ += MoonReactionAdapter.o MoonReactionController.o MoonTracing.o
MAINLIBOBJ += MoonNativeReactionRegistration.o MoonSemanticRandom.o MoonSodiumCohort.o
MAINLIBOBJ += MoonTransportClock.o species/PiecewiseHazard.o species/ReactionExposure.o
# Approved P01--P09 component services. These are the same production kernels
# linked by supporting probes; linking does not select missing state/threshold
# data, activate physical reactions or install the optional ideal collector.
MAINLIBOBJ += species/ReducedKinematics.o species/SodiumPartialCapture.o
MAINLIBOBJ += species/IonCollector.o species/HeliumSilicaChargeSurvival.o
# Source-bound M3 completion services. These objects provide the same profile,
# event, exposure and keyed-random kernels used by qualification adapters.
# Linking them does not admit a missing physical source/state crosswalk.
MAINLIBOBJ += species/PhysicalStateProfile.o species/PhysicalReactionEvents.o species/HeliumElectronImpact.o
MAINLIBOBJ += species/QualificationRandom.o species/TrajectoryExposure.o MoonTrajectoryExposure.o
MAINOBJ = main.o

# libAMPS.a flattens archive member basenames. The upstream shared module also
# owns ChargeExchange.o, so this application object needs its unique name;
# otherwise the later core archive step silently replaces the Moon service.
# Preserve the normal configured compiler/MPI/optimization flags verbatim.
species/MoonChargeExchange.o: species/ChargeExchange.cpp species/ChargeExchange.h species/Photochemistry.h
	${COMPILE.mpicxx} ${FLAGCC} $< -o $@

lib : $(MAINLIBOBJ)
	ar -rc mainlib.a $(MAINLIBOBJ)
	@echo "Target 'main' is done." 

amps : $(MAINOBJ) 
	ar -rc main.a $(MAINOBJ)
