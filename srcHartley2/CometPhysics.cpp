#include "CometPhysics.h"

#include <cstdio>
#include <cstdlib>

#include "pic.h"
#include "CometPhysics.generated.h"

const CometaryPhysics::cCollisionDatabase&
Comet::Physics::GetCollisionDatabase() {
  // Function-local construction avoids cross-translation-unit initialization
  // ordering.  The generated arrays have static storage for the whole run.
  static const CometaryPhysics::cCollisionDatabase database(
      Generated::Collisions,Generated::CollisionCount,PIC::nTotalSpecies);
  return database;
}

const char *Comet::Physics::GetDatabaseName() {
  return Generated::DatabaseName;
}

bool Comet::Physics::Validate(std::string *error) {
  return GetCollisionDatabase().Validate(error);
}

double Comet::Physics::GetNeutralCrossSection(double *v0,double *v1,
    int species0,int species1,PIC::Mesh::cDataBlockAMR *block,
    PIC::Mesh::cDataCenterNode *cell) {
  // The current Rubin database uses velocity-independent hard-sphere values.
  // Keep the AMPS callback arguments in this adapter so a future database can
  // add a relative-speed law without changing the core interface.
  (void)v0;
  (void)v1;
  (void)block;

  const CometaryPhysics::cCollision *entry=
      GetCollisionDatabase().Find(species0,species1);
  if (entry==NULL) {
    // Returning zero here would silently turn an undocumented pair into a
    // collisionless gas.  Stop at the first such encounter and force the case
    // input/database to supply traceable physics.
    std::fprintf(stderr,
        "Comet collision database '%s' has no cross section for species pair (%d,%d)\n",
        GetDatabaseName(),species0,species1);
    std::abort();
  }

  double temperature=1.0;
  if (entry->Law==CometaryPhysics::CollisionLaw::TemperaturePowerLaw) {
    if (cell==NULL) {
      std::fprintf(stderr,
          "temperature-dependent comet collision requested without a cell\n");
      std::abort();
    }

    // Rubin Table 3 has one temperature-dependent entry: H2O-H2O.  Its
    // Crifo law is defined in terms of the water translational temperature,
    // so no mixed-species temperature closure is needed or assumed.
    if (species0!=species1) {
      std::fprintf(stderr,
          "temperature law for mixed comet species (%d,%d) needs an explicit temperature definition\n",
          species0,species1);
      std::abort();
    }
    temperature=cell->GetTranslationalTemperature(species0);
    if (!(temperature>0.0)) {
      std::fprintf(stderr,
          "non-positive translational temperature for comet collision species %d\n",
          species0);
      std::abort();
    }
  }

  bool found=false,temperatureInRange=false;
  const double sigma=GetCollisionDatabase().GetCrossSection(
      species0,species1,temperature,&found,&temperatureInRange);
  if (!temperatureInRange) {
    std::fprintf(stderr,
        "comet collision temperature %.17g K is outside [%.17g,%.17g] K for species pair (%d,%d)\n",
        temperature,entry->MinimumTemperature,entry->MaximumTemperature,
        species0,species1);
    std::abort();
  }
  if (!found || !(sigma>0.0)) {
    std::fprintf(stderr,
        "invalid comet collision cross section for species pair (%d,%d)\n",
        species0,species1);
    std::abort();
  }
  return sigma;
}

// The legacy DIM3 Comet.cpp owns this AMPS callback on three-dimensional
// builds.  Radial builds do not link Comet.cpp, so provide the core symbol
// here only for DIM1.  Comet.cpp delegates to the same adapter for DIM3.
#if DIM == 1
double PIC::MolecularCollisions::ParticleCollisionModel::UserDefined::
GetTotalCrossSection(double *v0,double *v1,int species0,int species1,
    PIC::Mesh::cDataBlockAMR *block,PIC::Mesh::cDataCenterNode *cell) {
  return Comet::Physics::GetNeutralCrossSection(
      v0,v1,species0,species1,block,cell);
}
#endif
