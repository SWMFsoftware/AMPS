// Deterministic ownership for standalone field initialization of aliased AMR cells.
//
// Neighboring blocks may refer to the SAME center-node buffer through an interior
// slot and one or more ghost slots. Splitting logical slots among workers does not
// split the physical storage: concurrent writes race, even if the resulting fields
// differ only through coordinate-rounding. This helper chooses one writer before
// workers start. It has no AMPS/MPI dependency so the aliasing contract is unit-testable.
#ifndef EARTH_MODE3D_FIELD_INITIALIZATION_PLAN_H
#define EARTH_MODE3D_FIELD_INITIALIZATION_PLAN_H

#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Earth { namespace Mode3D {

class UniqueFieldWriterPlan {
  struct Entry {
    std::size_t outputIndex;
    bool interior;
  };
  std::unordered_map<const void*,Entry> writers_;
  std::vector<std::size_t> slots_;

public:
  // Visit slots in deterministic tree/index order, never in pointer/hash order.
  // Prefer the owning interior cell's coordinates over an aliasing ghost's
  // coordinates. For ghost-only storage, keep the first deterministic slot.
  // Null slots are legitimate outside the mesh and do not create work items.
  void Add(const void* storage,std::size_t slot,bool interior) {
    if (!storage) return;
    const auto inserted=writers_.emplace(storage,Entry{slots_.size(),interior});
    if (inserted.second) {
      slots_.push_back(slot);
      return;
    }
    Entry& existing=inserted.first->second;
    if (interior && existing.interior && slots_[existing.outputIndex]!=slot)
      throw std::runtime_error("field initialization: distinct interior cells alias one buffer");
    if (interior && !existing.interior) {
      slots_[existing.outputIndex]=slot;
      existing.interior=true;
    }
  }

  // Release only the immutable work sequence to the workers. The temporary
  // pointer lookup can be destroyed before field evaluation starts; iteration
  // order, allocation addresses and worker scheduling cannot choose the writer.
  std::vector<std::size_t> TakeSlots() { return std::move(slots_); }
};

}}
#endif
