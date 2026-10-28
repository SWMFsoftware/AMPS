// Empty fallback used only when compiling directly from the source tree.
// AMPS preprocessing replaces this file in the build tree from the selected
// PhysicsDatabase input record.
#ifndef _COMET_PHYSICS_GENERATED_H_
#define _COMET_PHYSICS_GENERATED_H_

namespace Comet { namespace Physics { namespace Generated {
static const char DatabaseName[]="NONE";
static const std::size_t CollisionCount=0;
static const CometaryPhysics::cCollision Collisions[]={
  // The sentinel only gives standard C++ a nonempty backing array.
  // CollisionCount=0 guarantees that it is never returned as physical data.
  {0,0,CometaryPhysics::CollisionLaw::Constant,
      1.0,1.0,0.0,0.0,0.0,"sentinel"}
};
}}} // namespace Comet::Physics::Generated

#endif
