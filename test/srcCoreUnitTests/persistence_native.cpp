#include "pic.h"
#include "src/test/native/test_harness.h"
#include "core_test_mesh.h"
#include "core_test_state.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace AMPS {
namespace CoreTests {

// Static archives discard translation units with no referenced symbol.  The
// main program calls this anchor so the native restart case is registered in
// every compatible core-unit-test build.
void LinkNativePersistenceTests() {}

}  // namespace CoreTests
}  // namespace AMPS

namespace {

class RuntimeStateScope {
 public:
  RuntimeStateScope()
      : time_step_initialized_(
            PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized),
        random_seed_(RandomNumberGenerator::rndLastSeed),
        last_sample_length_(PIC::LastSampleLength),
        collecting_sample_counter_(PIC::CollectingSampleCounter),
        data_output_file_number_(PIC::DataOutputFileNumber) {
    AMPS_EXPECT_TRUE(PIC::ParticleWeightTimeStep::GlobalParticleWeight != NULL);
    weights_.assign(PIC::ParticleWeightTimeStep::GlobalParticleWeight,
                  PIC::ParticleWeightTimeStep::GlobalParticleWeight +
                      PIC::nTotalSpecies);
    time_steps_.assign(PIC::ParticleWeightTimeStep::GlobalTimeStep,
                       PIC::ParticleWeightTimeStep::GlobalTimeStep +
                           PIC::nTotalSpecies);
  }
  ~RuntimeStateScope() {
    std::copy(weights_.begin(), weights_.end(),
              PIC::ParticleWeightTimeStep::GlobalParticleWeight);
    std::copy(time_steps_.begin(), time_steps_.end(),
              PIC::ParticleWeightTimeStep::GlobalTimeStep);
    PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized =
        time_step_initialized_;
    RandomNumberGenerator::rndLastSeed = random_seed_;
    PIC::LastSampleLength = last_sample_length_;
    PIC::CollectingSampleCounter = collecting_sample_counter_;
    PIC::DataOutputFileNumber = data_output_file_number_;
  }

 private:
  std::vector<double> weights_, time_steps_;
  bool time_step_initialized_;
  unsigned long int random_seed_;
  long int last_sample_length_, collecting_sample_counter_;
  long int data_output_file_number_;
};

void SetOpaqueBytes(PIC::ParticleBuffer::byte* record,
                    std::vector<PIC::ParticleBuffer::byte>& expected) {
  expected.clear();
  for (int field = 0; field < AMPS::CoreTests::T2AnonymousFieldCount(); ++field) {
    const long int offset = AMPS::CoreTests::T2AnonymousFieldOffset(field);
    const int length = AMPS::CoreTests::T2AnonymousFieldLength(field);
    for (int index = 0; index < length; ++index) {
      const PIC::ParticleBuffer::byte value =
          static_cast<PIC::ParticleBuffer::byte>(0x31 + field + index);
      record[offset + index] = value;
      expected.push_back(value);
    }
  }
}

void CheckOpaqueBytes(const PIC::ParticleBuffer::byte* record,
                      const std::vector<PIC::ParticleBuffer::byte>& expected) {
  std::size_t cursor = 0;
  for (int field = 0; field < AMPS::CoreTests::T2AnonymousFieldCount(); ++field) {
    const long int offset = AMPS::CoreTests::T2AnonymousFieldOffset(field);
    const int length = AMPS::CoreTests::T2AnonymousFieldLength(field);
    for (int index = 0; index < length; ++index) {
      AMPS_EXPECT_TRUE(cursor < expected.size());
      AMPS_EXPECT_EQ(record[offset + index], expected[cursor++]);
    }
  }
  AMPS_EXPECT_EQ(cursor, expected.size());
}

}  // namespace

AMPS_TEST_CASE("PIC-RESTART-NATIVE-ROUNDTRIP-001",
               "native AMR particle restart preserves particle and generated-layout state") {
  // The T9 native fixture deliberately owns a one-rank contract.  Save/read
  // are collective production operations; a later T10 case will exercise
  // changed MPI ownership rather than letting rank zero run a collective
  // while other ranks are absent from this process-local harness.
  AMPS_EXPECT_EQ(PIC::nTotalThreads, 1);

  const double position[3] = {0.61, 0.57, 0.53};
  const double velocity[3] = {-120.0, 45.0, 7.5};
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
      AMPS::CoreTests::MoverFixtureNode(position);
  int i = -1, j = -1, k = -1;
  AMPS_EXPECT_TRUE(PIC::Mesh::mesh->FindCellIndex(
                       const_cast<double*>(position), i, j, k, node, false) >= 0);
  long int& head = node->block->FirstCellParticleTable[
      i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];

  AMPS::CoreTests::Fixture::ParticleBufferScope particles(128);
  AMPS::CoreTests::Fixture::EmptyListHeadScope empty_cell(head);
  RuntimeStateScope runtime_state;
  const double expected_weight = 8.75e19;
  const double expected_time_step = 0.03125;
  const unsigned long int expected_random_seed = 975318642ul;
  PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_] =
      expected_weight;
  PIC::ParticleWeightTimeStep::GlobalTimeStep[_H2O_SPEC_] =
      expected_time_step;
  PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized = true;
  RandomNumberGenerator::rndLastSeed = expected_random_seed;

  // Sampling restart is a separate AMPSPF1 kind.  Verify its real native
  // payload before the particle restart: completed sample bytes survive, the
  // two public counters return exactly, and recovery intentionally clears the
  // collecting counter so a partial pre-crash accumulation is not counted
  // twice.  The marker lives in the completed half of the double buffer and
  // is therefore not erased by Read()'s collecting-buffer flush.
  PIC::Mesh::cDataCenterNode* center =
      node->block->GetCenterNode(_getCenterNodeLocalNumber(i, j, k));
  AMPS_EXPECT_TRUE(center != NULL);
  const int sample_offset = PIC::Mesh::completedCellSampleDataPointerOffset;
  AMPS_EXPECT_TRUE(sample_offset >= 0 &&
      sample_offset + static_cast<int>(sizeof(double)) <=
      PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength);
  double* completed_marker = reinterpret_cast<double*>(
      center->GetAssociatedDataBufferPointer() + sample_offset);
  const double saved_marker = *completed_marker;
  const double expected_marker = 123.25;
  *completed_marker = expected_marker;
  PIC::LastSampleLength = 73;
  PIC::DataOutputFileNumber = 11;
  PIC::CollectingSampleCounter = 29;
  const std::string sampling_file = "native-sampling-restart.amps";
  std::remove(sampling_file.c_str());
  PIC::Restart::SamplingData::Save(sampling_file.c_str());
  *completed_marker = -1.0;
  PIC::LastSampleLength = -1;
  PIC::DataOutputFileNumber = -1;
  PIC::CollectingSampleCounter = -1;
  PIC::Restart::SamplingData::Read(sampling_file.c_str());
  AMPS_EXPECT_NEAR(*completed_marker, expected_marker, 0.0);
  AMPS_EXPECT_EQ(PIC::LastSampleLength, 73);
  AMPS_EXPECT_EQ(PIC::DataOutputFileNumber, 11);
  AMPS_EXPECT_EQ(PIC::CollectingSampleCounter, 0);
  AMPS_EXPECT_EQ(std::remove(sampling_file.c_str()), 0);

  const long int original = PIC::ParticleBuffer::GetNewParticle(head, true);
  PIC::ParticleBuffer::SetI(_H2O_SPEC_, original);
  PIC::ParticleBuffer::SetX(const_cast<double*>(position), original);
  PIC::ParticleBuffer::SetV(const_cast<double*>(velocity), original);
  PIC::ParticleBuffer::byte* original_record =
      PIC::ParticleBuffer::GetParticleDataPointer(original);

  const double rotational_energy = 2.75e-20;
  if (PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_ >= 0 &&
      PIC::IDF::nTotalRotationalModes[_H2O_SPEC_] > 0)
    PIC::IDF::LB::SetRotE(rotational_energy, original_record);
  std::vector<double> vibrational_energy;
  if (PIC::IDF::LB::_VIBRATIONAL_ENERGY_OFFSET_ >= 0) {
    const int modes = std::max(PIC::IDF::nTotalVibtationalModes[_H2O_SPEC_], 0);
    vibrational_energy.resize(modes);
    for (int mode = 0; mode < modes; ++mode) {
      vibrational_energy[mode] = (mode + 1) * 1.25e-21;
      PIC::IDF::LB::SetVibE(vibrational_energy[mode], mode, original_record);
    }
  }
  std::vector<PIC::ParticleBuffer::byte> opaque;
  SetOpaqueBytes(original_record, opaque);

  const std::string file = "native-particle-restart.amps";
  std::remove(file.c_str());
  PIC::Restart::SaveParticleData(file.c_str());
  AMPS_EXPECT_EQ(PIC::Restart::GetRestartFileParticleNumber(file.c_str()), 1);

  // Remove the in-memory source and perturb the global weight so successful
  // verification cannot accidentally read the pre-save state.  The restart
  // reader must allocate a new record and restore the serialized table.
  PIC::ParticleBuffer::DeleteParticle(original, head);
  AMPS_EXPECT_EQ(head, -1);
  PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_] = 1.0;
  PIC::ParticleWeightTimeStep::GlobalTimeStep[_H2O_SPEC_] = 99.0;
  PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized = false;
  RandomNumberGenerator::rndLastSeed = 1;
  PIC::Restart::ReadParticleData(file.c_str());

  AMPS_EXPECT_TRUE(head >= 0);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetNext(head), -1);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetI(head), _H2O_SPEC_);
  const double* restored_position = PIC::ParticleBuffer::GetX(head);
  const double* restored_velocity = PIC::ParticleBuffer::GetV(head);
  for (int dimension = 0; dimension < 3; ++dimension) {
    AMPS_EXPECT_NEAR(restored_position[dimension], position[dimension], 0.0);
    AMPS_EXPECT_NEAR(restored_velocity[dimension], velocity[dimension], 0.0);
  }
  AMPS_EXPECT_NEAR(
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_],
      expected_weight, 0.0);
  AMPS_EXPECT_NEAR(PIC::ParticleWeightTimeStep::GlobalTimeStep[_H2O_SPEC_],
                   expected_time_step, 0.0);
  AMPS_EXPECT_TRUE(PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized);
  AMPS_EXPECT_EQ(RandomNumberGenerator::rndLastSeed, expected_random_seed);

  PIC::ParticleBuffer::byte* restored_record =
      PIC::ParticleBuffer::GetParticleDataPointer(head);
  if (PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_ >= 0 &&
      PIC::IDF::nTotalRotationalModes[_H2O_SPEC_] > 0)
    AMPS_EXPECT_NEAR(PIC::IDF::LB::GetRotE(restored_record),
                     rotational_energy, 0.0);
  for (std::size_t mode = 0; mode < vibrational_energy.size(); ++mode)
    AMPS_EXPECT_NEAR(PIC::IDF::LB::GetVibE(static_cast<int>(mode),
                                            restored_record),
                     vibrational_energy[mode], 0.0);
  CheckOpaqueBytes(restored_record, opaque);

  const long int restored = head;
  PIC::ParticleBuffer::DeleteParticle(restored, head);
  AMPS_EXPECT_EQ(head, -1);
  AMPS_EXPECT_EQ(std::remove(file.c_str()), 0);
  *completed_marker = saved_marker;
}
