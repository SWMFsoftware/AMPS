#include "ReactionNetwork.h"

#include <cassert>
#include <cmath>
#include <string>

namespace {
bool Near(double a,double b) {
  return std::fabs(a-b)<=1.0e-12*(1.0+std::fabs(a)+std::fabs(b));
}
}

int main() {
  using namespace ChemicalSolver;

  // R0: A -> B + 2 C, k=0.5 s^-1
  // R1: A -> B, k(2 AU)=1 s^-1
  // R2: A + B -> C, k(T)=2*T with activation temperature zero
  // R3: C -> no tracked neutral product (for example, photoionization)
  const cStoichiometricTerm terms[]={
      {0,1.0},{1,1.0},{2,2.0},
      {0,1.0},{1,1.0},
      {0,1.0},{1,1.0},{2,1.0},
      {2,1.0}};
  const cReaction reactions[]={
      {"decay",0,1,1,2,RateLaw::Constant,0.5,0.0,0.0,1.0,0.0},
      {"photo",3,1,4,1,RateLaw::HeliocentricInverseSquare,
          4.0,0.0,0.0,1.0,0.0},
      {"binary",5,2,7,1,RateLaw::Arrhenius,
          2.0,1.0,0.0,1.0,0.0},
      {"ionization",8,1,9,0,RateLaw::HeliocentricInverseSquare,
          1.0e-6,0.0,0.0,1.0,0.0}};

  cReactionNetwork network(terms,sizeof(terms)/sizeof(terms[0]),reactions,
      sizeof(reactions)/sizeof(reactions[0]),3);
  std::string error;
  assert(network.Validate(&error));
  assert(network.IsUnitUnimolecular(0));
  assert(network.IsUnitUnimolecular(1));
  assert(!network.IsUnitUnimolecular(2));
  assert(network.IsUnitUnimolecular(3));
  assert(network.GetSpeciesCount()==3);
  assert(network.GetTermCount()==9);
  assert(network.GetReactantTerm(3,0).Species==2);
  assert(network.GetReaction(3).ProductCount==0);
  assert(Near(network.GetRateCoefficient(1,100.0,2.0),1.0));

  const double density[]={4.0,5.0,0.0};
  assert(Near(network.GetReactionRate(0,density,100.0,1.0),2.0));
  assert(Near(network.GetReactionRate(2,density,10.0,1.0),400.0));

  // Evaluate the first reaction in isolation to verify stoichiometric source
  // and loss accounting without contributions from the other channels.
  cReactionNetwork decayOnly(terms,3,reactions,1,3);
  double source[3],loss[3];
  decayOnly.EvaluateSourceLoss(density,100.0,1.0,source,loss);
  assert(Near(loss[0],2.0));
  assert(Near(source[1],2.0));
  assert(Near(source[2],4.0));

  const double sinkDensity[]={0.0,0.0,2.0};
  const cReaction sinkReaction={"ionization",0,1,1,0,
      RateLaw::HeliocentricInverseSquare,
      1.0e-6,0.0,0.0,1.0,0.0};
  cReactionNetwork sinkOnly(terms+8,1,&sinkReaction,1,3);
  assert(sinkOnly.Validate(&error));
  sinkOnly.EvaluateSourceLoss(sinkDensity,100.0,1.0,source,loss);
  assert(Near(loss[2],2.0e-6));
  assert(Near(source[0]+source[1]+source[2],0.0));

  const cStoichiometricTerm badTerms[]={{3,1.0},{1,1.0}};
  const cReaction badReaction={"bad",0,1,1,1,RateLaw::Constant,
      1.0,0.0,0.0,1.0,0.0};
  cReactionNetwork invalid(badTerms,2,&badReaction,1,3);
  assert(!invalid.Validate(&error));

  return 0;
}
