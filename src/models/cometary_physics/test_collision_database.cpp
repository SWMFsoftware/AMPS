#include "CollisionDatabase.h"

#include <cassert>
#include <cmath>
#include <string>

namespace {
bool Near(double a,double b) {
  return std::fabs(a-b)<=1.0e-12*(1.0+std::fabs(a)+std::fabs(b));
}
}

int main() {
  using namespace CometaryPhysics;

  // Exercise both law types, unordered lookup, the published temperature
  // interval contract, missing-pair reporting, and duplicate rejection.
  const cCollision entries[]={
      {0,0,CollisionLaw::TemperaturePowerLaw,
          1.66e-19,300.0,-0.6,50.0,300.0,"Crifo (1989)"},
      {0,1,CollisionLaw::Constant,
          3.2e-19,1.0,0.0,0.0,0.0,"Rubin et al. (2011), Table 3"}};
  cCollisionDatabase database(entries,2,3);
  std::string error;
  assert(database.Validate(&error));
  assert(database.Find(1,0)==entries+1);
  bool found=false;
  assert(Near(database.GetCrossSection(0,0,300.0,&found),1.66e-19));
  assert(found);
  assert(Near(database.GetCrossSection(0,0,75.0),
      1.66e-19*std::pow(75.0/300.0,-0.6)));
  bool temperatureInRange=true;
  assert(database.GetCrossSection(0,0,49.0,&found,&temperatureInRange)==0.0);
  assert(found && !temperatureInRange);
  assert(Near(database.GetCrossSection(1,0,0.0),3.2e-19));
  assert(database.GetCrossSection(1,2,100.0,&found)==0.0);
  assert(!found);

  const cCollision duplicate[]={entries[1],
      {1,0,CollisionLaw::Constant,3.2e-19,1.0,0.0,0.0,0.0,"duplicate"}};
  cCollisionDatabase invalid(duplicate,2,3);
  assert(!invalid.Validate(&error));
  return 0;
}
