// Geometry-independent neutral collision data for comet coma models.
//
// Applications generate cCollision entries after resolving species names to
// AMPS species indices.  A missing pair is reported explicitly through the
// optional `found` result; zero is not treated as a physical cross section.

#ifndef _AMPS_COMETARY_PHYSICS_COLLISION_DATABASE_H_
#define _AMPS_COMETARY_PHYSICS_COLLISION_DATABASE_H_

#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

namespace CometaryPhysics {

enum class CollisionLaw {
  // ReferenceCrossSection is independent of the gas state.  The temperature
  // fields in cCollision are unused for this law and must be zero.
  Constant,

  // sigma(T)=ReferenceCrossSection*
  //          (T/ReferenceTemperature)^TemperatureExponent.
  // The database also records the temperature interval over which the cited
  // source states that this fit is valid.
  TemperaturePowerLaw
};

struct cCollision {
  // Species0/Species1 form an unordered pair.  cCollisionDatabase::Find
  // therefore gives the same result for (a,b) and (b,a).
  int Species0;
  int Species1;
  CollisionLaw Law;

  // Cross sections are in m^2 and all temperatures are in K.  A constant-law
  // entry stores zero in MinimumTemperature and MaximumTemperature.  A
  // temperature-law entry stores the closed validity interval published with
  // the fit; evaluating outside it is rejected instead of extrapolated.
  double ReferenceCrossSection;
  double ReferenceTemperature;
  double TemperatureExponent;
  double MinimumTemperature;
  double MaximumTemperature;

  // Source is diagnostic provenance, not part of the numerical calculation.
  const char *Source;
};

class cCollisionDatabase {
 public:
  cCollisionDatabase(const cCollision *collisions,std::size_t collisionCount,
      int speciesCount) : Collisions(collisions),CollisionCount(collisionCount),
      SpeciesCount(speciesCount) {
    BuildPairIndex();
  }

  std::size_t GetCollisionCount() const {return CollisionCount;}
  int GetSpeciesCount() const {return SpeciesCount;}
  const cCollision& GetCollision(std::size_t collision) const {
    return Collisions[collision];
  }

  bool Validate(std::string *error=NULL) const {
    if (SpeciesCount<=0) return Fail("the species count must be positive",error);
    if (CollisionCount!=0 && Collisions==NULL)
      return Fail("the collision table is missing",error);

    for (std::size_t i=0;i<CollisionCount;i++) {
      const cCollision &entry=Collisions[i];
      if (entry.Species0<0 || entry.Species0>=SpeciesCount ||
          entry.Species1<0 || entry.Species1>=SpeciesCount) {
        std::ostringstream out;
        out << "collision entry " << i << " has a species index outside [0,"
            << SpeciesCount << ')';
        return Fail(out.str(),error);
      }
      if (!std::isfinite(entry.ReferenceCrossSection) ||
          entry.ReferenceCrossSection<=0.0)
        return Fail("a collision entry has a non-positive cross section",error);
      if (entry.Law==CollisionLaw::Constant) {
        if (entry.MinimumTemperature!=0.0 || entry.MaximumTemperature!=0.0)
          return Fail("a constant collision entry has a temperature range",error);
      }
      else if (entry.Law==CollisionLaw::TemperaturePowerLaw) {
        if (!std::isfinite(entry.ReferenceTemperature) ||
            entry.ReferenceTemperature<=0.0 ||
            !std::isfinite(entry.TemperatureExponent) ||
            !std::isfinite(entry.MinimumTemperature) ||
            !std::isfinite(entry.MaximumTemperature) ||
            entry.MinimumTemperature<=0.0 ||
            entry.MaximumTemperature<entry.MinimumTemperature)
          return Fail("a temperature-dependent collision entry has invalid parameters or validity range",error);
      }
      else {
        return Fail("a collision entry has an unknown collision law",error);
      }

      for (std::size_t j=0;j<i;j++)
        if (SamePair(entry,Collisions[j])) {
          std::ostringstream out;
          out << "duplicate collision pair (" << entry.Species0 << ','
              << entry.Species1 << ')';
          return Fail(out.str(),error);
        }
    }

    return true;
  }

  const cCollision *Find(int species0,int species1) const {
    if (species0<0 || species0>=SpeciesCount ||
        species1<0 || species1>=SpeciesCount || PairIndex.empty()) return NULL;
    return PairIndex[static_cast<std::size_t>(species0)*SpeciesCount+species1];
  }

  // Evaluate an entry already returned by Find().  Keeping record evaluation
  // separate from pair lookup prevents application adapters from paying for a
  // second search merely to apply the selected law.  It also leaves a single
  // documented implementation of the temperature-range contract.
  double Evaluate(const cCollision &entry,double temperature,
      bool *temperatureInRange=NULL) const {
    if (temperatureInRange!=NULL) *temperatureInRange=false;
    if (entry.Law==CollisionLaw::Constant) {
      if (temperatureInRange!=NULL) *temperatureInRange=true;
      return entry.ReferenceCrossSection;
    }

    // A zero result with temperatureInRange=false identifies a request outside
    // the cited law's domain.  Application adapters can report the offending
    // pair using the entry they already hold.
    if (!std::isfinite(temperature) ||
        temperature<entry.MinimumTemperature ||
        temperature>entry.MaximumTemperature) return 0.0;
    if (temperatureInRange!=NULL) *temperatureInRange=true;
    return entry.ReferenceCrossSection*std::pow(
        temperature/entry.ReferenceTemperature,entry.TemperatureExponent);
  }

  double GetCrossSection(int species0,int species1,double temperature,
      bool *found=NULL,bool *temperatureInRange=NULL) const {
    const cCollision *entry=Find(species0,species1);
    if (found!=NULL) *found=(entry!=NULL);
    if (temperatureInRange!=NULL) *temperatureInRange=false;
    if (entry==NULL) return 0.0;
    return Evaluate(*entry,temperature,temperatureInRange);
  }

 private:
  const cCollision *Collisions;
  std::size_t CollisionCount;
  int SpeciesCount;
  std::vector<const cCollision*> PairIndex;

  void BuildPairIndex() {
    // Like the reaction-network index, construction must be safe even for an
    // invalid generated table so Validate() can return the intended detailed
    // error.  Out-of-range entries are skipped.  For duplicates, retain the
    // first record in the index; Validate() still rejects the database before
    // it can be used by a simulation.
    if (SpeciesCount<=0) return;
    PairIndex.assign(static_cast<std::size_t>(SpeciesCount)*SpeciesCount,NULL);
    if (CollisionCount==0 || Collisions==NULL) return;

    for (std::size_t i=0;i<CollisionCount;i++) {
      const cCollision &entry=Collisions[i];
      if (entry.Species0<0 || entry.Species0>=SpeciesCount ||
          entry.Species1<0 || entry.Species1>=SpeciesCount) continue;
      const std::size_t forward=
          static_cast<std::size_t>(entry.Species0)*SpeciesCount+entry.Species1;
      const std::size_t reverse=
          static_cast<std::size_t>(entry.Species1)*SpeciesCount+entry.Species0;
      if (PairIndex[forward]==NULL) PairIndex[forward]=Collisions+i;
      if (PairIndex[reverse]==NULL) PairIndex[reverse]=Collisions+i;
    }
  }

  static bool SamePair(const cCollision &a,const cCollision &b) {
    return (a.Species0==b.Species0 && a.Species1==b.Species1) ||
        (a.Species0==b.Species1 && a.Species1==b.Species0);
  }

  static bool Fail(const std::string &message,std::string *error) {
    if (error!=NULL) *error=message;
    return false;
  }
};

} // namespace CometaryPhysics

#endif
