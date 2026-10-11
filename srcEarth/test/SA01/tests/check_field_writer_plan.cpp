// Regression for aliasing between block-interior and neighboring ghost slots.
// No physical field is manufactured here: different numeric tags stand in for
// the distinct rounded coordinates each logical slot would use. Without a unique
// writer, reversing completion order changes the stored value (the old defect).
#include "../../../3d/FieldInitializationPlan.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

int main() {
  using Earth::Mode3D::UniqueFieldWriterPlan;
  int interiorBuffer=0,ghostOnlyBuffer=0,secondInterior=0;
  std::vector<int*> buffers={&interiorBuffer,&ghostOnlyBuffer,&interiorBuffer,
                            &secondInterior,&interiorBuffer,&ghostOnlyBuffer};
  UniqueFieldWriterPlan plan;
  for (std::size_t slot=0;slot<buffers.size();++slot)
    plan.Add(buffers[slot],slot,slot==2 || slot==3);
  plan.Add(nullptr,6,false);
  auto slots=plan.TakeSlots();
  assert((slots==std::vector<std::size_t>{2,1,3}));

  // The rejected legacy policy has schedule-dependent results on the same aliases.
  for (std::size_t slot=0;slot<buffers.size();++slot) *buffers[slot]=int(slot);
  const int oldForward=interiorBuffer;
  for (std::size_t slot=buffers.size();slot-- > 0;) *buffers[slot]=int(slot);
  assert(oldForward!=interiorBuffer);

  // Each selected physical buffer has exactly one writer; even adversarial
  // completion order leaves the owner-interior and ghost-only choices unchanged.
  for (auto slot:slots) *buffers[slot]=int(slot);
  assert(interiorBuffer==2 && ghostOnlyBuffer==1 && secondInterior==3);
  std::reverse(slots.begin(),slots.end());
  for (auto slot:slots) *buffers[slot]=int(slot);
  assert(interiorBuffer==2 && ghostOnlyBuffer==1 && secondInterior==3);

  // Distinct interior cells sharing storage indicate a mesh ownership defect,
  // not an alias the initializer may silently resolve through last-writer wins.
  UniqueFieldWriterPlan invalid;
  invalid.Add(&interiorBuffer,0,true);
  bool rejected=false;
  try { invalid.Add(&interiorBuffer,1,true); }
  catch (const std::runtime_error&) { rejected=true; }
  assert(rejected);
  std::cout << "Unique field writer alias/order/ownership regression: PASS\n";
}
