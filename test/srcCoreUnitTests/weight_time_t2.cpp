#include "pic.h"
#include "src/test/native/test_harness.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

// Statistical weight is the dimensionless number of physical particles
// represented by one model particle.  Time step is measured in seconds.  Both
// are process-global production state in the configured species-global modes,
// so each case restores every compiled species and remains order independent.
class GlobalWeightTimeScope {
 public:
  GlobalWeightTimeScope()
      : weights_(PIC::nTotalSpecies),
        time_steps_(PIC::nTotalSpecies),
        time_step_initialized_(
            PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized) {
    for (int species = 0; species < PIC::nTotalSpecies; ++species) {
      weights_[species] =
          PIC::ParticleWeightTimeStep::GlobalParticleWeight[species];
      time_steps_[species] =
          PIC::ParticleWeightTimeStep::GlobalTimeStep[species];

      // Negative one is the production "not initialized" sentinel.  Starting
      // from it is necessary because the global block setters retain the
      // smallest value they have seen rather than unconditionally assigning.
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[species] = -1.0;
      PIC::ParticleWeightTimeStep::GlobalTimeStep[species] = -1.0;
    }
  }

  ~GlobalWeightTimeScope() {
    for (int species = 0; species < PIC::nTotalSpecies; ++species) {
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[species] =
          weights_[species];
      PIC::ParticleWeightTimeStep::GlobalTimeStep[species] =
          time_steps_[species];
    }
    PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized =
        time_step_initialized_;
  }

 private:
  GlobalWeightTimeScope(const GlobalWeightTimeScope&);
  GlobalWeightTimeScope& operator=(const GlobalWeightTimeScope&);

  std::vector<double> weights_;
  std::vector<double> time_steps_;
  bool time_step_initialized_;
};

// This is production cDataBlockAMR storage without an AMR tree.  Global modes
// do not dereference associatedDataPointer, while local modes store one double
// per species at offsets reserved by PIC::Mesh::initCellSamplingDataBuffer().
// Supplying the complete associated-data span makes the same cases safe and
// meaningful when the test application is later configured for local modes.
class BlockStorage {
 public:
  BlockStorage() : bytes_(BlockBytes(), 0U) {
    block_.SetAssociatedDataBufferPointer(
        reinterpret_cast<char*>(&bytes_[0]));
  }

  PIC::Mesh::cDataBlockAMR* get() { return &block_; }

 private:
  BlockStorage(const BlockStorage&);
  BlockStorage& operator=(const BlockStorage&);

  std::size_t BlockBytes() {
    const int length = block_.AssociatedDataLength();
    AMPS_EXPECT_TRUE(length >= 0);
    // vector::operator[] requires a real element even when a global-mode block
    // reserves no associated bytes.
    return std::max<std::size_t>(1U, static_cast<std::size_t>(length));
  }

  PIC::Mesh::cDataBlockAMR block_;
  std::vector<unsigned char> bytes_;
};

#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == \
        _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_ && \
    _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
class LocalRootScope {
 public:
  explicit LocalRootScope(PIC::Mesh::cDataBlockAMR* block)
      : saved_(NULL) {
    AMPS_EXPECT_TRUE(PIC::Mesh::mesh != NULL);
    saved_ = PIC::Mesh::mesh->rootTree;
    root_.block = block;
    PIC::Mesh::mesh->rootTree = &root_;
  }
  ~LocalRootScope() { PIC::Mesh::mesh->rootTree = saved_; }

 private:
  LocalRootScope(const LocalRootScope&);
  LocalRootScope& operator=(const LocalRootScope&);
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> root_;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* saved_;
};
#endif

void RequireTwoSpecies() {
  AMPS_EXPECT_TRUE(PIC::nTotalSpecies >= 2);
}

void RequireThreeSpecies() {
  AMPS_EXPECT_TRUE(PIC::nTotalSpecies >= 3);
}

template <class Exception, class Callable>
void ExpectException(Callable callable) {
  bool rejected = false;
  try {
    callable();
  } catch (const Exception&) {
    rejected = true;
  }
  AMPS_EXPECT_TRUE(rejected);
}

}  // namespace

namespace AMPS {
namespace CoreTests {

// Referenced by the native test executable so archive linkage retains this
// translation unit and executes its static AMPS_TEST_CASE registrations.
void LinkWeightTimeT2Tests() {}

}  // namespace CoreTests
}  // namespace AMPS

// One source file is compiled under the supported global/local generated
// profiles.  Profile-qualified IDs keep the global registry unambiguous while
// allowing the independent runner to prove which compile-time branch was
// executed.  The established species-global IDs remain unchanged for
// backward-compatible T2 evidence and aliases.
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == \
        _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_ && \
    _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
#define AMPS_T2_WEIGHT_ROUNDTRIP_ID "PIC-WEIGHT-LOCAL-ROUNDTRIP-001"
#define AMPS_T2_TIMESTEP_ROUNDTRIP_ID "PIC-TIMESTEP-LOCAL-ROUNDTRIP-001"
#define AMPS_T2_INVALID_VALUE_ID \
  "PIC-WEIGHT-TIMESTEP-LOCAL-INVALID-VALUE-001"
#define AMPS_T2_INVALID_SPECIES_ID \
  "PIC-WEIGHT-TIMESTEP-LOCAL-INVALID-SPECIES-001"
#elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == \
          _SINGLE_GLOBAL_PARTICLE_WEIGHT_ && \
      _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
#define AMPS_T2_WEIGHT_ROUNDTRIP_ID \
  "PIC-WEIGHT-SINGLE-GLOBAL-ROUNDTRIP-001"
#define AMPS_T2_TIMESTEP_ROUNDTRIP_ID \
  "PIC-TIMESTEP-SINGLE-GLOBAL-ROUNDTRIP-001"
#define AMPS_T2_INVALID_VALUE_ID \
  "PIC-WEIGHT-TIMESTEP-SINGLE-GLOBAL-INVALID-VALUE-001"
#define AMPS_T2_INVALID_SPECIES_ID \
  "PIC-WEIGHT-TIMESTEP-SINGLE-GLOBAL-INVALID-SPECIES-001"
#else
#define AMPS_T2_WEIGHT_ROUNDTRIP_ID "PIC-WEIGHT-MODE-ROUNDTRIP-001"
#define AMPS_T2_TIMESTEP_ROUNDTRIP_ID "PIC-TIMESTEP-MODE-ROUNDTRIP-001"
#define AMPS_T2_INVALID_VALUE_ID "PIC-WEIGHT-TIMESTEP-INVALID-VALUE-001"
#define AMPS_T2_INVALID_SPECIES_ID "PIC-WEIGHT-TIMESTEP-INVALID-SPECIES-001"
#endif

AMPS_TEST_CASE(AMPS_T2_WEIGHT_ROUNDTRIP_ID,
               "particle weights round-trip with configured block/global scope") {
  RequireTwoSpecies();
  GlobalWeightTimeScope state;
  BlockStorage first;
  BlockStorage second;
  const int species0 = 0;
  const int species1 = 1;

#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == \
    _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  // In local mode, (block,species) is the ownership key.  Values in physical
  // particles/model-particle must remain independent across both dimensions.
  first.get()->SetLocalParticleWeight(64.0, species0);
  first.get()->SetLocalParticleWeight(125.0, species1);
  second.get()->SetLocalParticleWeight(32.0, species0);
  second.get()->SetLocalParticleWeight(250.0, species1);
  AMPS_EXPECT_NEAR(first.get()->GetLocalParticleWeight(species0), 64.0, 0.0);
  AMPS_EXPECT_NEAR(first.get()->GetLocalParticleWeight(species1), 125.0, 0.0);
  AMPS_EXPECT_NEAR(second.get()->GetLocalParticleWeight(species0), 32.0, 0.0);
  AMPS_EXPECT_NEAR(second.get()->GetLocalParticleWeight(species1), 250.0, 0.0);
#elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == \
      _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  // SetGlobalParticleWeight and the block setter are minimum reducers.  All
  // nodes see the same per-species value; a larger candidate cannot replace a
  // previously established global weight, while a smaller candidate must.
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(species0, 64.0, NULL);
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(species1, 125.0, NULL);
  AMPS_EXPECT_NEAR(first.get()->GetLocalParticleWeight(species0), 64.0, 0.0);
  AMPS_EXPECT_NEAR(second.get()->GetLocalParticleWeight(species1), 125.0, 0.0);
  second.get()->SetLocalParticleWeight(96.0, species0);
  AMPS_EXPECT_NEAR(first.get()->GetLocalParticleWeight(species0), 64.0, 0.0);
  second.get()->SetLocalParticleWeight(32.0, species0);
  AMPS_EXPECT_NEAR(first.get()->GetLocalParticleWeight(species0), 32.0, 0.0);
  AMPS_EXPECT_NEAR(first.get()->GetLocalParticleWeight(species1), 125.0, 0.0);
#elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SINGLE_GLOBAL_PARTICLE_WEIGHT_
  // Single-global mode aliases every species and block to global slot zero.
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(species0, 64.0, NULL);
  AMPS_EXPECT_NEAR(first.get()->GetLocalParticleWeight(species0), 64.0, 0.0);
  AMPS_EXPECT_NEAR(second.get()->GetLocalParticleWeight(species1), 64.0, 0.0);
  second.get()->SetLocalParticleWeight(32.0, species1);
  AMPS_EXPECT_NEAR(first.get()->GetLocalParticleWeight(species0), 32.0, 0.0);
  AMPS_EXPECT_NEAR(first.get()->GetLocalParticleWeight(species1), 32.0, 0.0);
#else
#error Unsupported AMPS particle-weight mode in T2 native tests
#endif
}

AMPS_TEST_CASE(AMPS_T2_TIMESTEP_ROUNDTRIP_ID,
               "time steps round-trip with configured block/global scope") {
  RequireTwoSpecies();
  GlobalWeightTimeScope state;
  BlockStorage first;
  BlockStorage second;
  const int species0 = 0;
  const int species1 = 1;

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  // Local time steps [s] are stored independently for every (block,species).
  first.get()->SetLocalTimeStep(0.50, species0);
  first.get()->SetLocalTimeStep(0.25, species1);
  second.get()->SetLocalTimeStep(0.125, species0);
  second.get()->SetLocalTimeStep(0.0625, species1);
  AMPS_EXPECT_NEAR(first.get()->GetLocalTimeStep(species0), 0.50, 0.0);
  AMPS_EXPECT_NEAR(first.get()->GetLocalTimeStep(species1), 0.25, 0.0);
  AMPS_EXPECT_NEAR(second.get()->GetLocalTimeStep(species0), 0.125, 0.0);
  AMPS_EXPECT_NEAR(second.get()->GetLocalTimeStep(species1), 0.0625, 0.0);
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  // The node setter reduces into one minimum stable step [s] per species.
  first.get()->SetLocalTimeStep(0.50, species0);
  first.get()->SetLocalTimeStep(0.25, species1);
  AMPS_EXPECT_NEAR(PIC::ParticleWeightTimeStep::GetGlobalTimeStep(species0),
                   0.50, 0.0);
  AMPS_EXPECT_NEAR(second.get()->GetLocalTimeStep(species1), 0.25, 0.0);
  second.get()->SetLocalTimeStep(0.75, species0);
  AMPS_EXPECT_NEAR(PIC::ParticleWeightTimeStep::GetGlobalTimeStep(species0),
                   0.50, 0.0);
  second.get()->SetLocalTimeStep(0.125, species0);
  AMPS_EXPECT_NEAR(first.get()->GetLocalTimeStep(species0), 0.125, 0.0);
  AMPS_EXPECT_NEAR(first.get()->GetLocalTimeStep(species1), 0.25, 0.0);
#elif _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
  // Single-global mode aliases every species and block to time-step slot zero.
  first.get()->SetLocalTimeStep(0.50, species1);
  AMPS_EXPECT_NEAR(PIC::ParticleWeightTimeStep::GetGlobalTimeStep(species0),
                   0.50, 0.0);
  AMPS_EXPECT_NEAR(second.get()->GetLocalTimeStep(species1), 0.50, 0.0);
  second.get()->SetLocalTimeStep(0.125, species0);
  AMPS_EXPECT_NEAR(first.get()->GetLocalTimeStep(species0), 0.125, 0.0);
  AMPS_EXPECT_NEAR(first.get()->GetLocalTimeStep(species1), 0.125, 0.0);
#else
#error Unsupported AMPS time-step mode in T2 native tests
#endif
}

#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == \
        _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_ && \
    _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
AMPS_TEST_CASE("PIC-WEIGHT-TIMESTEP-SPECIES-COPY-001",
               "species-global distribution copies scale only the target") {
  RequireThreeSpecies();
  GlobalWeightTimeScope state;
  BlockStorage block;
  const int source = 0;
  const int target = 1;
  const int control = 2;

  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(source, 40.0, NULL);
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(target, 7.0, NULL);
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(control, 99.0, NULL);
  block.get()->SetLocalTimeStep(0.20, source);
  block.get()->SetLocalTimeStep(0.03, target);
  block.get()->SetLocalTimeStep(0.75, control);

  // These production APIs define a distribution copy as target=coefficient
  // times source.  Source and unrelated species are invariants of the copy.
  PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(
      target, source, 2.5);
  PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(
      target, source, 0.5);
  AMPS_EXPECT_NEAR(
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[source], 40.0, 0.0);
  AMPS_EXPECT_NEAR(
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[target], 100.0, 0.0);
  AMPS_EXPECT_NEAR(
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[control], 99.0, 0.0);
  AMPS_EXPECT_NEAR(PIC::ParticleWeightTimeStep::GetGlobalTimeStep(source),
                   0.20, 0.0);
  AMPS_EXPECT_NEAR(PIC::ParticleWeightTimeStep::GetGlobalTimeStep(target),
                   0.10, 0.0);
  AMPS_EXPECT_NEAR(PIC::ParticleWeightTimeStep::GetGlobalTimeStep(control),
                   0.75, 0.0);
}
#endif

#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == \
        _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_ && \
    _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
AMPS_TEST_CASE("PIC-WEIGHT-TIMESTEP-LOCAL-RESCALE-001",
               "local distribution rescaling preserves weight over time step") {
  RequireThreeSpecies();
  GlobalWeightTimeScope state;
  BlockStorage block;
  LocalRootScope root(block.get());
  const int source = 0;
  const int target = 1;
  const int control = 2;
  block.get()->SetLocalParticleWeight(40.0, source);
  block.get()->SetLocalParticleWeight(7.0, target);
  block.get()->SetLocalParticleWeight(99.0, control);
  block.get()->SetLocalTimeStep(0.20, source);
  block.get()->SetLocalTimeStep(0.03, target);
  block.get()->SetLocalTimeStep(0.75, control);

  PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(
      target, source, 2.5);
  PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(
      target, source, 2.5);
  AMPS_EXPECT_NEAR(block.get()->GetLocalParticleWeight(source), 40.0, 0.0);
  AMPS_EXPECT_NEAR(block.get()->GetLocalParticleWeight(target), 100.0, 0.0);
  AMPS_EXPECT_NEAR(block.get()->GetLocalParticleWeight(control), 99.0, 0.0);
  AMPS_EXPECT_NEAR(block.get()->GetLocalTimeStep(source), 0.20, 0.0);
  AMPS_EXPECT_NEAR(block.get()->GetLocalTimeStep(target), 0.50, 0.0);
  AMPS_EXPECT_NEAR(block.get()->GetLocalTimeStep(control), 0.75, 0.0);
  AMPS_EXPECT_NEAR(
      block.get()->GetLocalParticleWeight(source) /
          block.get()->GetLocalTimeStep(source),
      block.get()->GetLocalParticleWeight(target) /
          block.get()->GetLocalTimeStep(target),
      0.0);
}
#endif

AMPS_TEST_CASE(AMPS_T2_INVALID_VALUE_ID,
               "weight/time APIs reject nonpositive and nonfinite scales") {
  GlobalWeightTimeScope state;
  BlockStorage block;
  const double invalid[4] = {
      0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::infinity()};

  for (int index = 0; index < 4; ++index) {
    const double value = invalid[index];
    ExpectException<std::invalid_argument>([&block, value]() {
      block.get()->SetLocalParticleWeight(value, 0);
    });
    ExpectException<std::invalid_argument>([&block, value]() {
      block.get()->SetLocalTimeStep(value, 0);
    });
    ExpectException<std::invalid_argument>([value]() {
      PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(0, value, NULL);
    });
    ExpectException<std::invalid_argument>([value]() {
      PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(
          0, 0, value);
    });
    ExpectException<std::invalid_argument>([value]() {
      PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(0, 0, value);
    });
  }

#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == \
    _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  PIC::ParticleWeightTimeStep::GlobalParticleWeight[0] = -1.0;
  ExpectException<std::invalid_argument>([]() {
    PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(
        0, 0, 1.0);
  });
  PIC::ParticleWeightTimeStep::GlobalParticleWeight[0] =
      std::numeric_limits<double>::max();
  ExpectException<std::invalid_argument>([]() {
    PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(
        0, 0, 2.0);
  });
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  PIC::ParticleWeightTimeStep::GlobalTimeStep[0] = -1.0;
  ExpectException<std::invalid_argument>([]() {
    PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(0, 0, 1.0);
  });
  PIC::ParticleWeightTimeStep::GlobalTimeStep[0] =
      std::numeric_limits<double>::max();
  ExpectException<std::invalid_argument>([]() {
    PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(0, 0, 2.0);
  });
#endif
}

AMPS_TEST_CASE(AMPS_T2_INVALID_SPECIES_ID,
               "weight/time APIs reject species indices before storage access") {
  GlobalWeightTimeScope state;
  BlockStorage block;
  const int invalid[2] = {-1, PIC::nTotalSpecies};
  for (int index = 0; index < 2; ++index) {
    const int species = invalid[index];
    ExpectException<std::out_of_range>([&block, species]() {
      block.get()->SetLocalParticleWeight(1.0, species);
    });
    ExpectException<std::out_of_range>([&block, species]() {
      (void)block.get()->GetLocalParticleWeight(species);
    });
    ExpectException<std::out_of_range>([&block, species]() {
      block.get()->SetLocalTimeStep(1.0, species);
    });
    ExpectException<std::out_of_range>([&block, species]() {
      (void)block.get()->GetLocalTimeStep(species);
    });
    ExpectException<std::out_of_range>([species]() {
      PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(species, 1.0, NULL);
    });
    ExpectException<std::out_of_range>([species]() {
      (void)PIC::ParticleWeightTimeStep::GetGlobalTimeStep(species);
    });
    ExpectException<std::out_of_range>([species]() {
      PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(
          species, 0, 1.0);
    });
    ExpectException<std::out_of_range>([species]() {
      PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(
          0, species, 1.0);
    });
    ExpectException<std::out_of_range>([species]() {
      PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(
          species, 0, 1.0);
    });
    ExpectException<std::out_of_range>([species]() {
      PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(
          0, species, 1.0);
    });
  }
}
