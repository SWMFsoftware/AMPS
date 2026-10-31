#include "CometPhysics.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "pic.h"
#include "CometPhysics.generated.h"

double Comet::Physics::CollisionTemperatureFloor=-1.0;
double Comet::Physics::CollisionTemperatureMinimumSampleParticleNumber=0.0;

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
  const CometaryPhysics::cCollisionDatabase &database=GetCollisionDatabase();
  if (!database.Validate(error)) return false;
  for (std::size_t entry=0;entry<database.GetCollisionCount();entry++)
    if (database.GetCollision(entry).Law==
        CometaryPhysics::CollisionLaw::TemperaturePowerLaw &&
        (!(CollisionTemperatureFloor>0.0) ||
         !std::isfinite(CollisionTemperatureFloor) ||
         !(CollisionTemperatureMinimumSampleParticleNumber>=0.0) ||
         !std::isfinite(CollisionTemperatureMinimumSampleParticleNumber))) {
      if (error!=NULL) {
        if (!(CollisionTemperatureFloor>0.0) ||
            !std::isfinite(CollisionTemperatureFloor))
          *error="a temperature-dependent comet collision law requires a positive CollisionTemperatureFloor";
        else
          *error="CollisionTemperatureMinimumSampleParticleNumber must be finite and non-negative";
      }
      return false;
    }
  return true;
}

double Comet::Physics::GetNeutralCrossSection(double *v0,double *v1,
    int species0,int species1,PIC::Mesh::cDataBlockAMR *block,
    PIC::Mesh::cDataCenterNode *cell) {
  // Rubin et al. (2011), Table 3, supplies hard-sphere cross sections for the
  // configured parent neutrals.  Most pairs are constant; H2O-H2O is the one
  // exception and uses the temperature-dependent Crifo (1989) power law below.
  // Neither law currently depends on relative speed, but keep both velocity
  // arguments in this adapter so a future database can add such a law without
  // changing the AMPS core callback interface.
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
    // Both accessors read the completed AMPS sampling buffer, not the buffer
    // currently being collected.  DatumParticleNumber is the mean number of
    // model particles occupying this cell during that completed window; it is
    // also the value printed in the "Particle Number" output column.  Do not
    // interpret a noisy temperature from a case-defined undersampled cell.
    const double sampledParticleNumber=cell->GetDatumAverage(
        PIC::Mesh::DatumParticleNumber,species0);
    if (CollisionTemperatureMinimumSampleParticleNumber>0.0 &&
        sampledParticleNumber<CollisionTemperatureMinimumSampleParticleNumber)
      temperature=CollisionTemperatureFloor;
    else {
      temperature=cell->GetTranslationalTemperature(species0);
      if (!(temperature>CollisionTemperatureFloor))
        temperature=CollisionTemperatureFloor;
    }

    // The validity interval is part of the selected collision-law record.
    // Cap a high sampled temperature at its documented upper endpoint rather
    // than extrapolating the power law beyond the cited data.  The Halley
    // Crifo record therefore evaluates every temperature above 300 K at 300 K.
    if (temperature>entry->MaximumTemperature)
      temperature=entry->MaximumTemperature;
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
