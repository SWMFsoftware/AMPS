#ifndef AMPS_CORE_TEST_STATE_H
#define AMPS_CORE_TEST_STATE_H

// Process-global state is part of PIC's production execution model.  Native
// module tests must therefore replace it transactionally: a case may fail by
// throwing from any assertion, but the next case must still observe exactly
// the state that existed before the fixture.  These small RAII guards keep the
// ownership/deallocation rules next to the snapshot logic and are shared by
// the particle-buffer, mover, collision, and photochemistry fixtures.

#include "pic.h"
#include "PhotolyticReactions.h"
#include "src/test/native/test_harness.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace AMPS {
namespace CoreTests {

// The anonymous T2 record extensions are intentionally owned by the layout
// translation unit.  Read-only accessors let split/merge tests state and
// verify the opaque-byte clone/survivor policy without duplicating offsets or
// making application extension semantics part of production PIC.
int T2AnonymousFieldCount();
long int T2AnonymousFieldOffset(int field);
int T2AnonymousFieldLength(int field);

namespace Fixture {

class ParticleBufferScope {
 public:
  explicit ParticleBufferScope(long int capacity)
      : data_(PIC::ParticleBuffer::ParticleDataBuffer),
        max_(PIC::ParticleBuffer::MaxNPart),
        count_(PIC::ParticleBuffer::NAllPart),
        first_(PIC::ParticleBuffer::FirstPBufferParticle),
        number_(PIC::ParticleBuffer::ParticleNumberTable),
        offset_(PIC::ParticleBuffer::ParticleOffsetTable),
        population_(PIC::ParticleBuffer::ParticlePopulationTable),
        thread_count_(PIC::ParticleBuffer::Thread::NTotalThreads),
        thread_available_(PIC::ParticleBuffer::Thread::AvailableParticleListLength),
        thread_first_(PIC::ParticleBuffer::Thread::FirstPBufferParticle) {
    PIC::ParticleBuffer::ParticleDataBuffer = NULL;
    PIC::ParticleBuffer::MaxNPart = 0;
    PIC::ParticleBuffer::NAllPart = 0;
    PIC::ParticleBuffer::FirstPBufferParticle = -1;
    PIC::ParticleBuffer::ParticleNumberTable = NULL;
    PIC::ParticleBuffer::ParticleOffsetTable = NULL;
    PIC::ParticleBuffer::ParticlePopulationTable = NULL;
    PIC::ParticleBuffer::Thread::NTotalThreads = 0;
    PIC::ParticleBuffer::Thread::AvailableParticleListLength = NULL;
    PIC::ParticleBuffer::Thread::FirstPBufferParticle = NULL;
    PIC::ParticleBuffer::Init(capacity);
  }

  ~ParticleBufferScope() {
    if (PIC::ParticleBuffer::ParticleNumberTable != NULL)
      amps_free_managed(PIC::ParticleBuffer::ParticleNumberTable);
    if (PIC::ParticleBuffer::ParticlePopulationTable != NULL)
      amps_free_managed(PIC::ParticleBuffer::ParticlePopulationTable);
    if (PIC::ParticleBuffer::ParticleOffsetTable != NULL)
      amps_free_managed(PIC::ParticleBuffer::ParticleOffsetTable);
    if (PIC::ParticleBuffer::ParticleDataBuffer != NULL) {
#if defined(__linux__) && _CUDA_MODE_ == _OFF_
      if (_ALIGN_STATE_VECTORS_BASE_ == 1) {
        amps_free_managed(PIC::ParticleBuffer::ParticleDataBuffer);
      } else {
        // ParticleBuffer::Init uses aligned_alloc on this host path.  The
        // matching release is free(), not the delete[] in amps_free_managed.
        std::free(PIC::ParticleBuffer::ParticleDataBuffer);
        PIC::ParticleBuffer::ParticleDataBuffer = NULL;
      }
#else
      amps_free_managed(PIC::ParticleBuffer::ParticleDataBuffer);
#endif
    }
    if (PIC::ParticleBuffer::Thread::AvailableParticleListLength != NULL) {
#if _CUDA_MODE_ == _ON_
      amps_free_managed(PIC::ParticleBuffer::Thread::AvailableParticleListLength);
#else
      delete[] PIC::ParticleBuffer::Thread::AvailableParticleListLength;
      PIC::ParticleBuffer::Thread::AvailableParticleListLength = NULL;
#endif
    }
    if (PIC::ParticleBuffer::Thread::FirstPBufferParticle != NULL) {
#if _CUDA_MODE_ == _ON_
      amps_free_managed(PIC::ParticleBuffer::Thread::FirstPBufferParticle);
#else
      delete[] PIC::ParticleBuffer::Thread::FirstPBufferParticle;
      PIC::ParticleBuffer::Thread::FirstPBufferParticle = NULL;
#endif
    }

    PIC::ParticleBuffer::ParticleDataBuffer = data_;
    PIC::ParticleBuffer::MaxNPart = max_;
    PIC::ParticleBuffer::NAllPart = count_;
    PIC::ParticleBuffer::FirstPBufferParticle = first_;
    PIC::ParticleBuffer::ParticleNumberTable = number_;
    PIC::ParticleBuffer::ParticleOffsetTable = offset_;
    PIC::ParticleBuffer::ParticlePopulationTable = population_;
    PIC::ParticleBuffer::Thread::NTotalThreads = thread_count_;
    PIC::ParticleBuffer::Thread::AvailableParticleListLength = thread_available_;
    PIC::ParticleBuffer::Thread::FirstPBufferParticle = thread_first_;
  }

 private:
  ParticleBufferScope(const ParticleBufferScope&);
  ParticleBufferScope& operator=(const ParticleBufferScope&);

  PIC::ParticleBuffer::byte* data_;
  long int max_, count_, first_;
  int *number_, *offset_;
  PIC::ParticleBuffer::cParticleTable* population_;
  int thread_count_;
  long int *thread_available_, *thread_first_;
};

class EmptyListHeadScope {
 public:
  explicit EmptyListHeadScope(long int& head) : head_(head), saved_(head) {
    // These fixtures intentionally own an otherwise empty list.  Detecting a
    // predecessor here is evidence of contamination and is safer than hiding
    // it by overwriting another case's particles.
    AMPS_EXPECT_EQ(saved_, -1);
  }
  ~EmptyListHeadScope() { head_ = saved_; }

 private:
  EmptyListHeadScope(const EmptyListHeadScope&);
  EmptyListHeadScope& operator=(const EmptyListHeadScope&);
  long int& head_;
  long int saved_;
};

class PhotolysisLifetimeScope {
 public:
  PhotolysisLifetimeScope()
      : table_(PIC::ChemicalReactions::PhotolyticReactions::ConstantTotalLifeTime) {
    AMPS_EXPECT_TRUE(table_ != NULL);
    values_.assign(table_, table_ + PIC::nTotalSpecies);
  }
  ~PhotolysisLifetimeScope() {
    if (table_ != NULL && values_.size() ==
                              static_cast<std::size_t>(PIC::nTotalSpecies))
      std::copy(values_.begin(), values_.end(), table_);
  }

 private:
  PhotolysisLifetimeScope(const PhotolysisLifetimeScope&);
  PhotolysisLifetimeScope& operator=(const PhotolysisLifetimeScope&);
  double* table_;
  std::vector<double> values_;
};

// NTC changes a cell's collecting-sample bytes, while its setup changes the
// block-local weights and time steps.  Preserve both as opaque/native state:
// copying the complete associated-data region also covers future collision
// sampling fields without embedding their offsets in the fixture.
class MeshCellStateScope {
 public:
  MeshCellStateScope(PIC::Mesh::cDataBlockAMR* block,
                     PIC::Mesh::cDataCenterNode* cell)
      : block_(block), cell_(cell) {
    AMPS_EXPECT_TRUE(block_ != NULL && cell_ != NULL);
    weights_.resize(PIC::nTotalSpecies);
    time_steps_.resize(PIC::nTotalSpecies);
    for (int species = 0; species < PIC::nTotalSpecies; ++species) {
      weights_[species] = block_->GetLocalParticleWeight(species);
      time_steps_[species] = block_->GetLocalTimeStep(species);
    }
    const int length = cell_->AssociatedDataLength();
    AMPS_EXPECT_TRUE(length >= 0);
    associated_data_.resize(static_cast<std::size_t>(length));
    if (length != 0) {
      AMPS_EXPECT_TRUE(cell_->GetAssociatedDataBufferPointer() != NULL);
      std::memcpy(&associated_data_[0], cell_->GetAssociatedDataBufferPointer(),
                  associated_data_.size());
    }
  }

  ~MeshCellStateScope() {
    // Restore storage, not setter semantics.  In species-global modes the
    // public setters are minimum reducers, and the legitimate pre-test state
    // can be the -1 "not initialized" sentinel.  Feeding that sentinel back
    // through a physical-value setter is both inexact and, after T2 domain
    // validation, intentionally rejected.  Direct restoration is confined to
    // this transactional test fixture and recreates the byte/global state that
    // the constructor observed.
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == \
    _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
    double* local_weights = reinterpret_cast<double*>(
        block_->GetAssociatedDataBufferPointer() +
        PIC::Mesh::cDataBlockAMR_static_data::LocalParticleWeightOffset);
    std::copy(weights_.begin(), weights_.end(), local_weights);
#elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == \
      _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    std::copy(weights_.begin(), weights_.end(),
              PIC::ParticleWeightTimeStep::GlobalParticleWeight);
#elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SINGLE_GLOBAL_PARTICLE_WEIGHT_
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[0] = weights_[0];
#else
#error Unsupported particle-weight mode in native test state fixture
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
    double* local_time_steps = reinterpret_cast<double*>(
        block_->GetAssociatedDataBufferPointer() +
        PIC::Mesh::cDataBlockAMR_static_data::LocalTimeStepOffset);
    std::copy(time_steps_.begin(), time_steps_.end(), local_time_steps);
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    std::copy(time_steps_.begin(), time_steps_.end(),
              PIC::ParticleWeightTimeStep::GlobalTimeStep);
#elif _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
    PIC::ParticleWeightTimeStep::GlobalTimeStep[0] = time_steps_[0];
#else
#error Unsupported time-step mode in native test state fixture
#endif
    if (!associated_data_.empty())
      std::memcpy(cell_->GetAssociatedDataBufferPointer(), &associated_data_[0],
                  associated_data_.size());
  }

 private:
  MeshCellStateScope(const MeshCellStateScope&);
  MeshCellStateScope& operator=(const MeshCellStateScope&);
  PIC::Mesh::cDataBlockAMR* block_;
  PIC::Mesh::cDataCenterNode* cell_;
  std::vector<double> weights_, time_steps_;
  std::vector<unsigned char> associated_data_;
};

template <typename T>
class ValueScope {
 public:
  explicit ValueScope(T& value) : value_(value), saved_(value) {}
  ~ValueScope() { value_ = saved_; }

 private:
  ValueScope(const ValueScope&);
  ValueScope& operator=(const ValueScope&);
  T& value_;
  T saved_;
};

}  // namespace Fixture
}  // namespace CoreTests
}  // namespace AMPS

#endif
