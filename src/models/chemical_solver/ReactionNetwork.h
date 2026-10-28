// A small, model-independent description of a chemical reaction network.
//
// AMPS species remain a compile-time list.  A model input parser can therefore
// translate species names to their integer AMPS indices and generate the flat
// term/reaction tables consumed here.  The network itself does not depend on a
// particular comet, geometry, or particle implementation.

#ifndef _AMPS_CHEMICAL_SOLVER_REACTION_NETWORK_H_
#define _AMPS_CHEMICAL_SOLVER_REACTION_NETWORK_H_

#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>

namespace ChemicalSolver {

enum class RateLaw {
  // k=A.  A carries the units required by the reaction molecularity.
  Constant,

  // k(r)=A*(ReferenceHeliocentricDistanceAU/r)^2.  This is used for
  // unshielded photoprocesses tabulated at a reference solar distance.
  HeliocentricInverseSquare,

  // k(T)=A*T^TemperatureExponent*exp(-ActivationTemperature/T).
  Arrhenius
};

struct cStoichiometricTerm {
  // Species is the compile-time AMPS species index.  Coefficient is the
  // positive stoichiometric multiplicity stored as a double so the generic
  // continuum evaluator can also represent noninteger kinetic mechanisms.
  int Species;
  double Coefficient;
};

struct cReaction {
  // Id points to static generated storage and is used in diagnostics.
  const char *Id;

  // Reactants and products occupy nonoverlapping slices of the shared Terms
  // array.  ProductCount may be zero for loss to untracked states, such as
  // photoionization in a neutral-only model.
  std::size_t ReactantOffset;
  std::size_t ReactantCount;
  std::size_t ProductOffset;
  std::size_t ProductCount;
  RateLaw Law;

  // Constant and heliocentric laws use A directly.  For an Arrhenius law,
  // k(T)=A*T^TemperatureExponent*exp(-ActivationTemperature/T).
  double A;
  double TemperatureExponent;
  double ActivationTemperature;
  double ReferenceHeliocentricDistanceAU;
  double ExcessEnergyJ;
};

class cReactionNetwork {
 public:
  // The network is a non-owning view over generated static tables.  Keeping
  // ownership outside this class makes the same data usable in application
  // callbacks without dynamic initialization or geometry dependencies.
  cReactionNetwork(const cStoichiometricTerm *terms,std::size_t termCount,
      const cReaction *reactions,std::size_t reactionCount,int speciesCount) :
      Terms(terms),TermCount(termCount),Reactions(reactions),
      ReactionCount(reactionCount),SpeciesCount(speciesCount) {}

  std::size_t GetReactionCount() const {return ReactionCount;}
  std::size_t GetTermCount() const {return TermCount;}
  int GetSpeciesCount() const {return SpeciesCount;}
  const cReaction& GetReaction(std::size_t reaction) const {
    return Reactions[reaction];
  }
  const cStoichiometricTerm& GetTerm(std::size_t term) const {
    return Terms[term];
  }
  const cStoichiometricTerm& GetReactantTerm(std::size_t reaction,
      std::size_t term) const {
    const cReaction &entry=Reactions[reaction];
    return Terms[entry.ReactantOffset+term];
  }
  const cStoichiometricTerm& GetProductTerm(std::size_t reaction,
      std::size_t term) const {
    const cReaction &entry=Reactions[reaction];
    return Terms[entry.ProductOffset+term];
  }

  // Check table bounds, species indices, coefficients, and rate-law
  // parameters before any particle or continuum solver uses the network.
  bool Validate(std::string *error=NULL) const {
    if (SpeciesCount<=0) return Fail("the species count must be positive",error);
    if (ReactionCount!=0 && Reactions==NULL)
      return Fail("the reaction table is missing",error);
    if (TermCount!=0 && Terms==NULL)
      return Fail("the stoichiometric-term table is missing",error);

    for (std::size_t i=0;i<ReactionCount;i++) {
      const cReaction &reaction=Reactions[i];
      const std::string id=(reaction.Id==NULL) ? std::string("<unnamed>") :
          std::string(reaction.Id);

      if (reaction.ReactantCount==0)
        return Fail("reaction '"+id+"' must have at least one reactant",error);
      if (reaction.ReactantOffset+reaction.ReactantCount>TermCount ||
          reaction.ProductOffset+reaction.ProductCount>TermCount)
        return Fail("reaction '"+id+"' has an invalid term range",error);
      if (!std::isfinite(reaction.A) || reaction.A<0.0)
        return Fail("reaction '"+id+"' has an invalid rate coefficient",error);
      if (!std::isfinite(reaction.ExcessEnergyJ) || reaction.ExcessEnergyJ<0.0)
        return Fail("reaction '"+id+"' has invalid excess energy",error);
      if (reaction.Law==RateLaw::HeliocentricInverseSquare &&
          (!std::isfinite(reaction.ReferenceHeliocentricDistanceAU) ||
           reaction.ReferenceHeliocentricDistanceAU<=0.0))
        return Fail("reaction '"+id+"' requires a positive reference distance",error);
      if (reaction.Law==RateLaw::Arrhenius &&
          (!std::isfinite(reaction.TemperatureExponent) ||
           !std::isfinite(reaction.ActivationTemperature)))
        return Fail("reaction '"+id+"' has invalid Arrhenius parameters",error);

      for (std::size_t j=0;j<reaction.ReactantCount;j++)
        if (!ValidateTerm(Terms[reaction.ReactantOffset+j],id,error)) return false;
      for (std::size_t j=0;j<reaction.ProductCount;j++)
        if (!ValidateTerm(Terms[reaction.ProductOffset+j],id,error)) return false;
    }

    return true;
  }

  // Evaluate k for one reaction.  Invalid nonpositive environmental inputs
  // produce zero for laws that require them; Validate covers static data.
  double GetRateCoefficient(std::size_t reaction,double temperature,
      double heliocentricDistanceAU) const {
    const cReaction &entry=Reactions[reaction];

    switch (entry.Law) {
    case RateLaw::Constant:
      return entry.A;
    case RateLaw::HeliocentricInverseSquare:
      if (!(heliocentricDistanceAU>0.0)) return 0.0;
      return entry.A*std::pow(
          entry.ReferenceHeliocentricDistanceAU/heliocentricDistanceAU,2.0);
    case RateLaw::Arrhenius:
      if (!(temperature>0.0)) return 0.0;
      return entry.A*std::pow(temperature,entry.TemperatureExponent)*
          std::exp(-entry.ActivationTemperature/temperature);
    }

    return 0.0;
  }

  // Return the mass-action event rate k*product(n_i^nu_i).  Number-density
  // units must be consistent with A and are deliberately chosen by the
  // caller because the table itself is unit-system agnostic.
  double GetReactionRate(std::size_t reaction,const double *numberDensity,
      double temperature,double heliocentricDistanceAU) const {
    const cReaction &entry=Reactions[reaction];
    double rate=GetRateCoefficient(reaction,temperature,
        heliocentricDistanceAU);

    for (std::size_t i=0;i<entry.ReactantCount;i++) {
      const cStoichiometricTerm &term=Terms[entry.ReactantOffset+i];
      const double density=numberDensity[term.Species];
      if (!(density>0.0)) return 0.0;
      rate*=std::pow(density,term.Coefficient);
    }

    return rate;
  }

  // Accumulate separate nonnegative production and loss rates for every
  // species.  A solver can combine them using its own integration method.
  void EvaluateSourceLoss(const double *numberDensity,double temperature,
      double heliocentricDistanceAU,double *source,double *loss) const {
    for (int species=0;species<SpeciesCount;species++)
      source[species]=0.0,loss[species]=0.0;

    for (std::size_t i=0;i<ReactionCount;i++) {
      const cReaction &reaction=Reactions[i];
      const double rate=GetReactionRate(i,numberDensity,temperature,
          heliocentricDistanceAU);

      for (std::size_t j=0;j<reaction.ReactantCount;j++) {
        const cStoichiometricTerm &term=Terms[reaction.ReactantOffset+j];
        loss[term.Species]+=term.Coefficient*rate;
      }
      for (std::size_t j=0;j<reaction.ProductCount;j++) {
        const cStoichiometricTerm &term=Terms[reaction.ProductOffset+j];
        source[term.Species]+=term.Coefficient*rate;
      }
    }
  }

  // Discrete particle decay callbacks can execute only one unit parent at a
  // time.  This predicate lets an application reject general mass-action
  // reactions before starting such a run.
  bool IsUnitUnimolecular(std::size_t reaction) const {
    const cReaction &entry=Reactions[reaction];
    return entry.ReactantCount==1 &&
        Terms[entry.ReactantOffset].Coefficient==1.0;
  }

 private:
  const cStoichiometricTerm *Terms;
  std::size_t TermCount;
  const cReaction *Reactions;
  std::size_t ReactionCount;
  int SpeciesCount;

  static bool Fail(const std::string &message,std::string *error) {
    if (error!=NULL) *error=message;
    return false;
  }

  bool ValidateTerm(const cStoichiometricTerm &term,const std::string &id,
      std::string *error) const {
    if (term.Species<0 || term.Species>=SpeciesCount) {
      std::ostringstream out;
      out << "reaction '" << id << "' references species index "
          << term.Species << " outside [0," << SpeciesCount << ')';
      return Fail(out.str(),error);
    }
    if (!std::isfinite(term.Coefficient) || term.Coefficient<=0.0)
      return Fail("reaction '"+id+"' has a non-positive stoichiometric coefficient",error);
    return true;
  }
};

} // namespace ChemicalSolver

#endif
