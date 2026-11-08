/*
 * Application-owned closures for the native spherical-exosphere fixture.
 *
 * The source and surface kernels are production Exosphere code.  Only the
 * target-specific quantities that the shared library intentionally leaves to
 * applications are fixed here.  SI units are used throughout: temperature is
 * kelvin, positions are metres, and the sticking callback returns a
 * dimensionless probability and represented-particle fraction.
 */
#include "pic.h"
#include "Exosphere.h"

#include <cstring>

char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_] = "TestSphere";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_] = "IAU_TEST_SPHERE";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_] = "TEST_SPHERE_SO";

double Exosphere::GetSurfaceTemperature(
    double, double*) {
  // A constant 250 K isolates the source/ledger algorithms from a target
  // thermal map.  It is positive and finite, as required by the Maxwellian
  // emission and accommodation kernels.
  return 250.0;
}

double Exosphere::SurfaceInteraction::StickingProbability(
    int, double& reemission_particle_fraction, double) {
  // Every selected impact sticks, while 75% of its represented weight enters
  // the retained inventory.  This yields exact independent closure oracles:
  // return=W, sticking=W, adsorption=0.75 W.
  reemission_particle_fraction = 0.75;
  return 1.0;
}

int Exosphere::ColumnIntegral::GetVariableList(char* variable_list) {
  if (variable_list != NULL) variable_list[0] = '\0';
  return 0;
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(
    double* result, int result_length, double*,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*) {
  if (result_length != 0)
    exit(__LINE__, __FILE__,
         "native sphere fixture declares no column-integral observables");
  for (int i = 0; i < result_length; ++i) result[i] = 0.0;
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(
    double*, int result_length) {
  // GetVariableList() declares an empty column-integral vector for this
  // source-only fixture, so there are no density/speed pairs to normalize.
  // A nonzero length contradicts that contract and therefore indicates a
  // framework/fixture integration defect rather than a value that can be
  // normalized meaningfully.
  if (result_length != 0)
    exit(__LINE__, __FILE__,
         "native sphere fixture declares no column-integral observables");
}

double Exosphere::SourceProcesses::SolarWindSputtering::
    EnergyDistributionFunction(double energy_joule, int* species) {
  // Use the shared Thompson-like spectrum exactly as supplied by the
  // production exosphere library:
  //   f(E) proportional to E/(E+E_b)^3 * [1-sqrt((E+E_b)/T_m)],
  // with E, E_b and T_m in joules.  Target applications may override this
  // callback with laboratory-calibrated material data; this generic sphere
  // has no material-specific measurements and therefore must not invent a
  // different law.  The native test qualifies dispatch, ownership and rate
  // ledgers, not the provenance of a target-specific sputter spectrum.
  return DefaultEnergyDistributionFunction(energy_joule, species);
}

double Exosphere::OrbitalMotion::GetTAA(SpiceDouble) { return 0.0; }
