#include "pic.h"
#include "src/test/native/test_harness.h"
#include "core_test_state.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <vector>

namespace AMPS {
namespace CoreTests {

// This anchor follows the native-test archive contract.  T2 integration will
// add the call from main.cpp only in the profiles that include this object.
void LinkSplitMergeT2Tests() {}

}  // namespace CoreTests
}  // namespace AMPS

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_

namespace {

namespace PB = PIC::ParticleBuffer;
using AMPS::CoreTests::Fixture::ParticleBufferScope;

// Split/merge uses the process RNG to orient merged velocity pairs.  Preserve
// the complete calling-thread state so reverse/shuffled test order cannot
// perturb later stochastic physics cases.
class RandomStateScope {
 public:
  explicit RandomStateScope(unsigned long seed)
      : seed_(RandomNumberGenerator::rndLastSeed),
        generator_(RandomNumberGenerator::gen),
        spare_available_(RandomNumberGenerator::NormalSpareAvailable),
        spare_(RandomNumberGenerator::NormalSpare) {
    RandomNumberGenerator::rndLastSeed = seed;
    RandomNumberGenerator::gen.seed(static_cast<std::mt19937::result_type>(seed));
    RandomNumberGenerator::NormalSpareAvailable = false;
    RandomNumberGenerator::NormalSpare = 0.0;
  }

  ~RandomStateScope() {
    RandomNumberGenerator::rndLastSeed = seed_;
    RandomNumberGenerator::gen = generator_;
    RandomNumberGenerator::NormalSpareAvailable = spare_available_;
    RandomNumberGenerator::NormalSpare = spare_;
  }

 private:
  unsigned long seed_;
  std::mt19937 generator_;
  bool spare_available_;
  double spare_;
};

struct RecordState {
  long int pointer;
  int species;
  double weight;
  std::array<double, 3> x;
  std::array<double, 3> v;
  double rotational_energy;
  std::vector<double> vibrational_energy;
  double internal_energy;
};

struct Moments {
  int count;
  std::vector<int> species_count;
  double physical_weight;
  std::array<double, 3> weighted_position;
  std::array<double, 3> momentum;
  double translational_energy;
  double rotational_energy;
  std::vector<double> vibrational_energy;
  double internal_energy;
  std::array<double, 3> position_scale;
  std::array<double, 3> momentum_scale;
  double translational_scale;
  double rotational_scale;
  std::vector<double> vibrational_scale;
  double internal_scale;
  std::vector<RecordState> records;
};

bool IsValidSpecies(int species) {
  return species >= 0 && species < PIC::nTotalSpecies;
}

bool HasRotationalEnergy(int species) {
  return IsValidSpecies(species) &&
         PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_ >= 0 &&
         PIC::IDF::nTotalRotationalModes[species] > 0;
}

int VibrationalModeCount(int species) {
  if (!IsValidSpecies(species) ||
      PIC::IDF::LB::_VIBRATIONAL_ENERGY_OFFSET_ < 0)
    return 0;
  return std::max(PIC::IDF::nTotalVibtationalModes[species], 0);
}

double GetInternalEnergy(int species, PB::byte* data) {
  double result = 0.0;
  if (HasRotationalEnergy(species)) result += PIC::IDF::LB::GetRotE(data);
  for (int mode = 0; mode < VibrationalModeCount(species); ++mode)
    result += PIC::IDF::LB::GetVibE(mode, data);
  return result;
}

void SetInternalEnergy(int species, PB::byte* data, double base) {
  if (HasRotationalEnergy(species)) PIC::IDF::LB::SetRotE(base, data);
  const int modes = VibrationalModeCount(species);
  for (int mode = 0; mode < modes; ++mode)
    PIC::IDF::LB::SetVibE(base * (0.25 + 0.125 * mode), mode, data);
}

void SetOpaqueExtensionTag(PB::byte* data, int ordinal) {
  const PB::byte tag = static_cast<PB::byte>(0x20 + ordinal);
  for (int field = 0; field < AMPS::CoreTests::T2AnonymousFieldCount(); ++field)
    std::memset(data + AMPS::CoreTests::T2AnonymousFieldOffset(field), tag,
                static_cast<std::size_t>(
                    AMPS::CoreTests::T2AnonymousFieldLength(field)));
}

std::multiset<int> CaptureOpaqueExtensionTags(long int head) {
  std::multiset<int> tags;
  const int count = AMPS::CoreTests::T2AnonymousFieldCount();
  if (count == 0) return tags;
  for (long int pointer = head; pointer != -1; pointer = PB::GetNext(pointer)) {
    PB::byte* data = PB::GetParticleDataPointer(pointer);
    const int tag = static_cast<int>(
        data[AMPS::CoreTests::T2AnonymousFieldOffset(0)]);
    for (int field = 0; field < count; ++field) {
      const long int offset = AMPS::CoreTests::T2AnonymousFieldOffset(field);
      const int length = AMPS::CoreTests::T2AnonymousFieldLength(field);
      for (int byte = 0; byte < length; ++byte)
        AMPS_EXPECT_EQ(static_cast<int>(data[offset + byte]), tag);
    }
    tags.insert(tag);
  }
  return tags;
}

void ValidateList(long int head) {
  std::set<long int> visited;
  long int previous = -1;
  long int pointer = head;
  while (pointer != -1) {
    AMPS_EXPECT_TRUE(pointer >= 0);
    AMPS_EXPECT_TRUE(pointer < PB::GetMaxNPart());
    AMPS_EXPECT_TRUE(PB::IsParticleAllocated(pointer));
    AMPS_EXPECT_TRUE(visited.insert(pointer).second);
    AMPS_EXPECT_EQ(PB::GetPrev(pointer), previous);
    previous = pointer;
    pointer = PB::GetNext(pointer);
    AMPS_EXPECT_TRUE(visited.size() <=
                     static_cast<std::size_t>(PB::GetMaxNPart()));
  }
  AMPS_EXPECT_EQ(static_cast<long int>(visited.size()), PB::GetAllPartNum());
}

Moments Capture(long int head) {
  ValidateList(head);
  Moments result{};
  result.species_count.assign(PIC::nTotalSpecies, 0);
  for (long int pointer = head; pointer != -1; pointer = PB::GetNext(pointer)) {
    RecordState record{};
    record.pointer = pointer;
    record.species = PB::GetI(pointer);
    record.weight = PB::GetIndividualStatWeightCorrection(pointer);
    PB::GetX(record.x.data(), pointer);
    PB::GetV(record.v.data(), pointer);
    PB::byte* data = PB::GetParticleDataPointer(pointer);
    record.rotational_energy = HasRotationalEnergy(record.species) ?
        PIC::IDF::LB::GetRotE(data) : 0.0;
    record.vibrational_energy.resize(VibrationalModeCount(record.species));
    for (std::size_t mode = 0; mode < record.vibrational_energy.size(); ++mode)
      record.vibrational_energy[mode] =
          PIC::IDF::LB::GetVibE(static_cast<int>(mode), data);
    record.internal_energy = GetInternalEnergy(record.species, data);

    AMPS_EXPECT_TRUE(record.species >= 0 &&
                     record.species < PIC::nTotalSpecies);
    AMPS_EXPECT_TRUE(std::isfinite(record.weight) && record.weight > 0.0);
    ++result.count;
    ++result.species_count[record.species];
    result.physical_weight += record.weight;

    const double mass = PIC::MolecularData::GetMass(record.species);
    double speed2 = 0.0;
    for (int d = 0; d < 3; ++d) {
      AMPS_EXPECT_TRUE(std::isfinite(record.x[d]));
      AMPS_EXPECT_TRUE(std::isfinite(record.v[d]));
      const double weighted_x = record.weight * record.x[d];
      const double weighted_p = record.weight * mass * record.v[d];
      result.weighted_position[d] += weighted_x;
      result.momentum[d] += weighted_p;
      result.position_scale[d] += std::fabs(weighted_x);
      result.momentum_scale[d] += std::fabs(weighted_p);
      speed2 += record.v[d] * record.v[d];
    }
    const double translational = 0.5 * record.weight * mass * speed2;
    const double rotational = record.weight * record.rotational_energy;
    const double internal = record.weight * record.internal_energy;
    result.translational_energy += translational;
    result.rotational_energy += rotational;
    result.internal_energy += internal;
    result.translational_scale += std::fabs(translational);
    result.rotational_scale += std::fabs(rotational);
    result.internal_scale += std::fabs(internal);
    if (result.vibrational_energy.size() < record.vibrational_energy.size()) {
      result.vibrational_energy.resize(record.vibrational_energy.size(), 0.0);
      result.vibrational_scale.resize(record.vibrational_energy.size(), 0.0);
    }
    for (std::size_t mode = 0; mode < record.vibrational_energy.size(); ++mode) {
      const double vibrational =
          record.weight * record.vibrational_energy[mode];
      result.vibrational_energy[mode] += vibrational;
      result.vibrational_scale[mode] += std::fabs(vibrational);
    }
    result.records.push_back(record);
  }
  return result;
}

bool NearConserved(double before, double after, double scale,
                   double tolerance = 2.0e-11) {
  return std::isfinite(before) && std::isfinite(after) &&
         std::fabs(after-before) <=
             tolerance * std::max(scale, 1.0e-300);
}

double WeightedCentralVelocityMoment(const Moments& state) {
  long double total_weight = 0.0L;
  std::array<long double, 3> mean{};
  for (const RecordState& record : state.records) {
    const long double weight = record.weight;
    total_weight += weight;
    for (int d = 0; d < 3; ++d) mean[d] += weight * record.v[d];
  }
  for (int d = 0; d < 3; ++d) mean[d] /= total_weight;

  // Evaluate the independent oracle in the center-of-mass frame.  This avoids
  // the large-bulk cancellation that the production merge used to suffer and
  // makes loss of a cold population's finite velocity spread observable even
  // when total kinetic energy is dominated by uniform bulk motion.
  long double result = 0.0L;
  for (const RecordState& record : state.records) {
    const long double weight = record.weight;
    for (int d = 0; d < 3; ++d) {
      const long double delta =
          static_cast<long double>(record.v[d]) - mean[d];
      result += weight * delta * delta;
    }
  }
  return static_cast<double>(result);
}

void ExpectCentralVelocityMomentConservation(const Moments& before,
                                             const Moments& after,
                                             double tolerance = 2.0e-9) {
  const double initial = WeightedCentralVelocityMoment(before);
  const double final = WeightedCentralVelocityMoment(after);
  AMPS_EXPECT_TRUE(NearConserved(initial, final, std::fabs(initial),
                                 tolerance));
}

void ExpectCoreConservation(const Moments& before, const Moments& after,
                            int expected_count) {
  AMPS_EXPECT_EQ(after.count, expected_count);
  AMPS_EXPECT_EQ(after.species_count.size(), before.species_count.size());
  for (std::size_t species = 0; species < before.species_count.size(); ++species)
    AMPS_EXPECT_EQ(after.species_count[species],
                   before.species_count[species] +
                       (species == 0 ? expected_count-before.count : 0));
  AMPS_EXPECT_TRUE(NearConserved(before.physical_weight,
                                 after.physical_weight,
                                 std::fabs(before.physical_weight)));
  for (int d = 0; d < 3; ++d)
    AMPS_EXPECT_TRUE(NearConserved(before.momentum[d], after.momentum[d],
                                   before.momentum_scale[d]));
  AMPS_EXPECT_TRUE(NearConserved(before.translational_energy,
                                 after.translational_energy,
                                 before.translational_scale));
}

void ExpectExtendedConservation(const Moments& before, const Moments& after) {
  bool position_conserved = true;
  for (int d = 0; d < 3; ++d)
    position_conserved = position_conserved &&
        NearConserved(before.weighted_position[d], after.weighted_position[d],
                      before.position_scale[d]);
  bool reservoirs_conserved =
      NearConserved(before.rotational_energy, after.rotational_energy,
                    before.rotational_scale) &&
      before.vibrational_energy.size() == after.vibrational_energy.size();
  for (std::size_t mode = 0;
       reservoirs_conserved && mode < before.vibrational_energy.size(); ++mode)
    reservoirs_conserved = NearConserved(
        before.vibrational_energy[mode], after.vibrational_energy[mode],
        before.vibrational_scale[mode]);
  const bool internal_conserved =
      NearConserved(before.internal_energy, after.internal_energy,
                    before.internal_scale);
  if (!position_conserved || !reservoirs_conserved || !internal_conserved) {
    std::ostringstream message;
    message << "resampling did not conserve extended physical moments: "
            << "weighted_position=" << position_conserved
            << " individual_internal_reservoirs=" << reservoirs_conserved
            << " weighted_internal_energy=" << internal_conserved;
    throw AMPS::Test::Failure(message.str());
  }
}

long int AddParticle(long int& head, int species, int ordinal, double weight,
                     int velocity_group, bool heterogeneous_state) {
  const long int pointer = PB::GetNewParticle(head);
  PB::SetI(species, pointer);

  // Adjacent members of each pair/triple occupy the same velocity bin but are
  // not identical. This exercises the nonzero relative-energy algebra while
  // retaining a deterministic viable tuple at the initial bin resolution.
  const double parity = 0.08 * (ordinal % 3);
  double velocity[3] = {300.0 + 12.0 * velocity_group + parity,
                        25.0 + 0.4 * velocity_group - 0.5 * parity,
                        -7.0 + 0.2 * velocity_group + 0.25 * parity};
  double position[3] = {1.25, -2.5, 0.75};
  if (heterogeneous_state) {
    position[0] += 0.11 * ordinal;
    position[1] -= 0.07 * ordinal;
    position[2] += 0.03 * ordinal;
  }
  PB::SetV(velocity, pointer);
  PB::SetX(position, pointer);
  PB::SetIndividualStatWeightCorrection(weight, pointer);
  SetInternalEnergy(species, PB::GetParticleDataPointer(pointer),
                    heterogeneous_state ? 1.0e-21 * (ordinal + 1) : 2.5e-21);
  SetOpaqueExtensionTag(PB::GetParticleDataPointer(pointer), ordinal);
  return pointer;
}

void Populate(long int& head, int species, int count, bool variable_weight,
              int tuple_size, bool heterogeneous_state,
              std::vector<long int>* pointers = NULL) {
  for (int ordinal = 0; ordinal < count; ++ordinal) {
    const double weight = variable_weight ? 0.55 + 0.13 * ordinal : 1.0;
    const long int pointer = AddParticle(
        head, species, ordinal, weight, ordinal / tuple_size,
        heterogeneous_state);
    if (pointers != NULL) pointers->push_back(pointer);
  }
}

struct RawList {
  long int head;
  std::vector<long int> order;
  std::vector<std::vector<PB::byte>> bytes;
};

RawList CaptureRaw(long int head) {
  RawList result;
  result.head = head;
  const std::size_t length =
      static_cast<std::size_t>(PB::GetParticleDataLength());
  for (long int pointer = head; pointer != -1; pointer = PB::GetNext(pointer)) {
    result.order.push_back(pointer);
    result.bytes.push_back(std::vector<PB::byte>(length));
    std::memcpy(result.bytes.back().data(), PB::GetParticleDataPointer(pointer),
                length);
  }
  return result;
}

void ExpectRawUnchanged(const RawList& before, long int head) {
  const RawList after = CaptureRaw(head);
  AMPS_EXPECT_EQ(after.head, before.head);
  AMPS_EXPECT_EQ(after.order.size(), before.order.size());
  AMPS_EXPECT_EQ(after.bytes.size(), before.bytes.size());
  for (std::size_t i = 0; i < before.order.size(); ++i) {
    AMPS_EXPECT_EQ(after.order[i], before.order[i]);
    AMPS_EXPECT_EQ(after.bytes[i].size(), before.bytes[i].size());
    AMPS_EXPECT_TRUE(std::memcmp(after.bytes[i].data(), before.bytes[i].data(),
                                 before.bytes[i].size()) == 0);
  }
  ValidateList(head);
}

void RunNominal(bool variable_weight) {
  {
    ParticleBufferScope buffer(256);
    RandomStateScope random(1729);
    long int head = -1;
    Populate(head, 0, 6, variable_weight, 2, false);
    const Moments before = Capture(head);
    PIC::ParticleSplitting::SplitParticleList(0, head, 10);
    const Moments after = Capture(head);
    ExpectCoreConservation(before, after, 10);
    ExpectCentralVelocityMomentConservation(before, after);
    ExpectExtendedConservation(before, after);
  }
  {
    ParticleBufferScope buffer(256);
    RandomStateScope random(2718);
    long int head = -1;
    Populate(head, 0, 9, variable_weight, 3, false);
    const Moments before = Capture(head);
    PIC::ParticleSplitting::MergeParticleList(0, head, 5);
    const Moments after = Capture(head);
    ExpectCoreConservation(before, after, 5);
    ExpectCentralVelocityMomentConservation(before, after);
    ExpectExtendedConservation(before, after);
  }
}

void ExpectWeightedMeanState(const Moments& before,const Moments& after) {
  std::array<double,3> meanPosition{};
  for (int d=0;d<3;d++)
    meanPosition[d]=before.weighted_position[d]/before.physical_weight;
  const double meanRotationalEnergy=
      before.rotational_energy/before.physical_weight;
  std::vector<double> meanVibrationalEnergy(before.vibrational_energy.size());
  for (std::size_t mode=0;mode<meanVibrationalEnergy.size();++mode)
    meanVibrationalEnergy[mode]=
        before.vibrational_energy[mode]/before.physical_weight;
  const double meanInternalEnergy=
      before.internal_energy/before.physical_weight;

  // The generic Cartesian resampler reconstructs intensive core fields at
  // their pre-resampling weight average.  This makes every output independent
  // of which input record survived allocation/list mechanics and conserves
  // weighted position and each enabled LB energy reservoir exactly.
  for (const RecordState& record : after.records) {
    for (int d=0;d<3;d++)
      AMPS_EXPECT_NEAR(record.x[d],meanPosition[d],
                       2.0e-14*std::max(std::fabs(meanPosition[d]),1.0));
    AMPS_EXPECT_NEAR(
        record.rotational_energy,meanRotationalEnergy,
        2.0e-14*std::max(std::fabs(meanRotationalEnergy),1.0e-300));
    AMPS_EXPECT_EQ(record.vibrational_energy.size(),
                   meanVibrationalEnergy.size());
    for (std::size_t mode=0;mode<meanVibrationalEnergy.size();++mode)
      AMPS_EXPECT_NEAR(
          record.vibrational_energy[mode],meanVibrationalEnergy[mode],
          2.0e-14*std::max(std::fabs(meanVibrationalEnergy[mode]),1.0e-300));
    AMPS_EXPECT_NEAR(record.internal_energy,meanInternalEnergy,
                     2.0e-14*std::max(std::fabs(meanInternalEnergy),1.0e-300));
  }
}

}  // namespace

#if _STATE_VECTOR_MODE_ == _STATE_VECTOR_MODE_ALIGNED_
#define AMPS_T2_SPLIT_CASE(PACKED_ID, ALIGNED_ID) ALIGNED_ID
#else
#define AMPS_T2_SPLIT_CASE(PACKED_ID, ALIGNED_ID) PACKED_ID
#endif

AMPS_TEST_CASE(AMPS_T2_SPLIT_CASE(
                   "PIC-SPLIT-MERGE-EQUAL-WEIGHT-001",
                   "PIC-SPLIT-MERGE-ALIGNED-EQUAL-WEIGHT-001"),
               "equal-correction split/merge conserve physical moments") {
  RunNominal(false);
}

AMPS_TEST_CASE(AMPS_T2_SPLIT_CASE(
                   "PIC-SPLIT-MERGE-VARIABLE-WEIGHT-001",
                   "PIC-SPLIT-MERGE-ALIGNED-VARIABLE-WEIGHT-001"),
               "variable-correction split/merge conserve physical moments") {
  RunNominal(true);
}

AMPS_TEST_CASE(AMPS_T2_SPLIT_CASE(
                   "PIC-SPLIT-MERGE-BOUNDARY-001",
                   "PIC-SPLIT-MERGE-ALIGNED-BOUNDARY-001"),
               "empty, insufficient, unreachable, and unrelated populations are no-ops") {
  {
    ParticleBufferScope buffer(64);
    long int head = -1;
    PIC::ParticleSplitting::MergeParticleList(0, head, 1);
    PIC::ParticleSplitting::SplitParticleList(0, head, 3);
    AMPS_EXPECT_EQ(head, -1);
    ValidateList(head);
  }
  {
    ParticleBufferScope buffer(64);
    long int head = -1;
    AddParticle(head, 0, 0, 1.0, 0, true);
    const RawList before = CaptureRaw(head);
    PIC::ParticleSplitting::MergeParticleList(0, head, 1);
    PIC::ParticleSplitting::SplitParticleList(0, head, 3);
    ExpectRawUnchanged(before, head);
  }
  {
    ParticleBufferScope buffer(64);
    long int head = -1;
    const long int first = AddParticle(head, 0, 0, 1.0, 0, true);
    const long int second = AddParticle(head, 0, 1, 1.0, 0, true);
    double opposite[3] = {-300.0, -25.0, 7.0};
    PB::SetV(opposite, second);
    double positive[3] = {300.0, 25.0, -7.0};
    PB::SetV(positive, first);
    const RawList before = CaptureRaw(head);
    PIC::ParticleSplitting::MergeParticleList(0, head, 1);
    ExpectRawUnchanged(before, head);
    PIC::ParticleSplitting::SplitParticleList(0, head, 3);
    ExpectRawUnchanged(before, head);
  }
  {
    ParticleBufferScope buffer(64);
    long int head = -1;
    Populate(head, 1, 6, true, 2, true);
    const RawList before = CaptureRaw(head);
    PIC::ParticleSplitting::MergeParticleList(0, head, 2);
    PIC::ParticleSplitting::SplitParticleList(0, head, 12);
    ExpectRawUnchanged(before, head);
  }
}

AMPS_TEST_CASE(AMPS_T2_SPLIT_CASE(
                   "PIC-SPLIT-EXTENDED-STATE-001",
                   "PIC-SPLIT-ALIGNED-EXTENDED-STATE-001"),
               "split conserves and reconstructs weighted position/internal energy") {
  ParticleBufferScope buffer(64);
  long int head = -1;
  std::vector<long int> pointers;
  Populate(head, 0, 2, true, 2, true, &pointers);
  // With only two samples the production velocity range is defined by the
  // pair itself; any nonzero separation places them at opposite ends of the
  // initial bins.  Identical velocities form the intended viable split tuple
  // while heterogeneous position/internal/opaque state still exercises the
  // reconstruction policy independently of translational spread.
  double shared_velocity[3];
  PB::GetV(shared_velocity, pointers[0]);
  PB::SetV(shared_velocity, pointers[1]);
  const Moments before = Capture(head);
  const std::multiset<int> opaque_before = CaptureOpaqueExtensionTags(head);
  PIC::ParticleSplitting::SplitParticleList(0, head, 3);
  const Moments after = Capture(head);
  ExpectCoreConservation(before, after, 3);
  ExpectWeightedMeanState(before, after);
  ExpectExtendedConservation(before, after);
  if (!opaque_before.empty()) {
    AMPS_EXPECT_TRUE(opaque_before ==
                     (std::multiset<int>{0x20, 0x21}));
    const std::multiset<int> opaque_after = CaptureOpaqueExtensionTags(head);
    // The new particle clones the first member of the production-selected
    // tuple.  Intrusive-list order may select either original, so the policy
    // invariant is one complete duplicate of an existing opaque record, not
    // a particular ordinal chosen by the test.
    AMPS_EXPECT_TRUE(
        opaque_after == (std::multiset<int>{0x20, 0x20, 0x21}) ||
        opaque_after == (std::multiset<int>{0x20, 0x21, 0x21}));
  }
}

AMPS_TEST_CASE(AMPS_T2_SPLIT_CASE(
                   "PIC-MERGE-EXTENDED-STATE-001",
                   "PIC-MERGE-ALIGNED-EXTENDED-STATE-001"),
               "merge conserves cold spread, position, and internal reservoirs") {
  ParticleBufferScope buffer(64);
  RandomStateScope random(314159);
  long int head = -1;
  std::vector<long int> pointers;
  Populate(head, 0, 3, true, 3, true, &pointers);
  for (long int pointer : pointers) {
    double velocity[3];
    PB::GetV(velocity, pointer);
    velocity[0] += 1.0e9;
    velocity[1] -= 2.0e9;
    velocity[2] += 0.5e9;
    PB::SetV(velocity, pointer);
  }
  const Moments before = Capture(head);
  const std::multiset<int> opaque_before = CaptureOpaqueExtensionTags(head);
  PIC::ParticleSplitting::MergeParticleList(0, head, 2);
  const Moments after = Capture(head);
  ExpectCoreConservation(before, after, 2);
  ExpectCentralVelocityMomentConservation(before, after, 2.0e-6);
  ExpectWeightedMeanState(before, after);
  ExpectExtendedConservation(before, after);
  if (!opaque_before.empty()) {
    AMPS_EXPECT_TRUE(opaque_before ==
                     (std::multiset<int>{0x20, 0x21, 0x22}));
    AMPS_EXPECT_TRUE(CaptureOpaqueExtensionTags(head) ==
                     (std::multiset<int>{0x20, 0x21}));
  }
}

AMPS_TEST_CASE(AMPS_T2_SPLIT_CASE(
                   "PIC-WEIGHT-INDIVIDUAL-INVALID-001",
                   "PIC-WEIGHT-ALIGNED-INDIVIDUAL-INVALID-001"),
               "individual weight correction rejects nonphysical values") {
  ParticleBufferScope buffer(4);
  const long int particle = PB::GetNewParticle(true);
  PB::byte* data = PB::GetParticleDataPointer(particle);
  const double invalid[] = {
      0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::infinity()};
  for (std::size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); ++i) {
    bool pointer_rejected = false;
    bool record_rejected = false;
    try {
      PB::SetIndividualStatWeightCorrection(invalid[i], particle);
    } catch (const std::invalid_argument&) {
      pointer_rejected = true;
    }
    try {
      PB::SetIndividualStatWeightCorrection(invalid[i], data);
    } catch (const std::invalid_argument&) {
      record_rejected = true;
    }
    AMPS_EXPECT_TRUE(pointer_rejected);
    AMPS_EXPECT_TRUE(record_rejected);
  }

  const double invalid_stored = std::numeric_limits<double>::quiet_NaN();
  std::memcpy(data + _PIC_PARTICLE_DATA__WEIGHT_CORRECTION_OFFSET_,
              &invalid_stored, sizeof(invalid_stored));
  bool getter_rejected = false;
  try {
    (void)PB::GetIndividualStatWeightCorrection(data);
  } catch (const std::domain_error&) {
    getter_rejected = true;
  }
  AMPS_EXPECT_TRUE(getter_rejected);
}

#undef AMPS_T2_SPLIT_CASE

#endif  // individual particle weights enabled
