#include "pic.h"
#include "src/test/native/test_harness.h"
#include "core_test_state.h"

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

using AMPS::CoreTests::Fixture::ParticleBufferScope;

// A single source file is compiled in four T2 layouts.  Profile-qualified
// evidence IDs prevent the global runner from conflating a PASS in the packed
// fixed-weight ABI with the materially different aligned or individual-weight
// ABIs.
#if _STATE_VECTOR_MODE_ == _STATE_VECTOR_MODE_ALIGNED_
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
#define AMPS_LAYOUT_CASE_ID(SUFFIX) \
  "PIC-PBUFFER-LAYOUT-ALIGNED-INDIVIDUAL-" SUFFIX
#else
#define AMPS_LAYOUT_CASE_ID(SUFFIX) \
  "PIC-PBUFFER-LAYOUT-ALIGNED-FIXED-" SUFFIX
#endif
#else
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
#define AMPS_LAYOUT_CASE_ID(SUFFIX) \
  "PIC-PBUFFER-LAYOUT-PACKED-INDIVIDUAL-" SUFFIX
#else
#define AMPS_LAYOUT_CASE_ID(SUFFIX) \
  "PIC-PBUFFER-LAYOUT-PACKED-FIXED-" SUFFIX
#endif
#endif

namespace {

struct ExtensionField {
  long int offset;
  int length;
  std::size_t alignment;
};

ExtensionField extension_fields[] = {
    {-1, 3, 1}, {-1, 5, 2}, {-1, 7, 4}, {-1, 9, 8}, {-1, 17, 16}};
bool extension_fields_registered = false;

void RegisterExtensionFields() {
  if (extension_fields_registered) return;
#if _STATE_VECTOR_MODE_ == _STATE_VECTOR_MODE_ALIGNED_
  // The aligned profiles are the dedicated anonymous-extension fixtures.
  // Packed profiles intentionally retain the standard production layout so
  // both contracts are tested independently.
  for (std::size_t i = 0;
       i < sizeof(extension_fields) / sizeof(extension_fields[0]); ++i) {
    PIC::ParticleBuffer::RequestDataStorage(extension_fields[i].offset,
                                             extension_fields[i].length,
                                             extension_fields[i].alignment);
  }
#endif
  extension_fields_registered = true;
}

long int ExpectedBasicDataLength() {
#if _STATE_VECTOR_MODE_ == _STATE_VECTOR_MODE_ALIGNED_
  const long int species = 2 * static_cast<long int>(sizeof(long int));
  const long int velocity = species + sizeof(long int);
  const long int correction = velocity + 3 * sizeof(double);
  const long int position =
      correction + _PIC_PARTICLE_DATA__WEIGHT_CORRECTION_LENGTH_;
  return position + DIM * sizeof(double);
#else
  const long int species = 2 * static_cast<long int>(sizeof(long int));
  const long int velocity = species + sizeof(unsigned char);
  const long int position = velocity + 3 * sizeof(double);
  const long int correction = position + DIM * sizeof(double);
  return correction + _PIC_PARTICLE_DATA__WEIGHT_CORRECTION_LENGTH_;
#endif
}

long int ExpectedCompileTimeFullLength() {
  return ExpectedBasicDataLength() +
      _PIC_PARTICLE_DATA__DUST_GRAIN_MASS_LENGTH_ +
      _PIC_PARTICLE_DATA__DUST_GRAIN_RADIUS_LENGTH_ +
      _PIC_PARTICLE_DATA__DUST_GRAIN_CHARGE_LENGTH_ +
      _PIC_PARTICLE_DATA__MAGNETIC_MOMENT_LENGTH_ +
      _PIC_PARTICLE_DATA__INJECTION_FACE_LENGTH_ +
      _PIC_PARTICLE_DATA__WEIGHT_OVER_TIME_STEP_LENGTH_ +
      _PIC_PARTICLE_DATA__FIELD_LINE_ID_LENGTH_ +
      _PIC_PARTICLE_DATA__FIELD_LINE_COORD_LENGTH_ +
      _PIC_PARTICLE_DATA__V_PARALLEL_LENGTH_ +
      _PIC_PARTICLE_DATA__V_NORMAL_LENGTH_;
}

template <class Exception, class Function>
bool Throws(Function function) {
  try {
    function();
  } catch (const Exception&) {
    return true;
  }
  return false;
}

void ExpectExtensionFieldsZero(long int particle) {
#if _STATE_VECTOR_MODE_ == _STATE_VECTOR_MODE_ALIGNED_
  PIC::ParticleBuffer::byte* record =
      PIC::ParticleBuffer::GetParticleDataPointer(particle);
  for (std::size_t field = 0;
       field < sizeof(extension_fields) / sizeof(extension_fields[0]); ++field)
    for (int byte = 0; byte < extension_fields[field].length; ++byte)
      AMPS_EXPECT_EQ(record[extension_fields[field].offset + byte], 0);
#else
  (void)particle;
#endif
}

}  // namespace

namespace AMPS {
namespace CoreTests {

// main.cpp must call this anchor before PIC::Init_BeforeParser().  Registration
// is deliberately explicit rather than a static initializer: production model
// fields and test fields then share the documented pre-buffer lifecycle, while
// their relative order is deterministic across linkers and build machines.
void LinkParticleLayoutT2Tests() { RegisterExtensionFields(); }

int T2AnonymousFieldCount() {
#if _STATE_VECTOR_MODE_ == _STATE_VECTOR_MODE_ALIGNED_
  return static_cast<int>(sizeof(extension_fields) / sizeof(extension_fields[0]));
#else
  return 0;
#endif
}

long int T2AnonymousFieldOffset(int field) {
  if (field < 0 || field >= T2AnonymousFieldCount())
    throw std::out_of_range("anonymous particle field index");
  return extension_fields[field].offset;
}

int T2AnonymousFieldLength(int field) {
  if (field < 0 || field >= T2AnonymousFieldCount())
    throw std::out_of_range("anonymous particle field index");
  return extension_fields[field].length;
}

}  // namespace CoreTests
}  // namespace AMPS

AMPS_TEST_CASE(AMPS_LAYOUT_CASE_ID("OFFSETS-001"),
               "particle record offsets, extension alignment, and stride") {
  AMPS_EXPECT_EQ(_PIC_PARTICLE_DATA__BASIC_DATA_LENGTH_,
                 ExpectedBasicDataLength());
  // This equality independently detects an individual-weight block being
  // counted both in BASIC_DATA_LENGTH and again in optional-field offsets.
  AMPS_EXPECT_EQ(_PIC_PARTICLE_DATA__FULL_DATA_LENGTH_,
                 ExpectedCompileTimeFullLength());

  long int previous_end = _PIC_PARTICLE_DATA__FULL_DATA_LENGTH_;
#if _STATE_VECTOR_MODE_ == _STATE_VECTOR_MODE_ALIGNED_
  for (std::size_t field = 0;
       field < sizeof(extension_fields) / sizeof(extension_fields[0]); ++field) {
    AMPS_EXPECT_TRUE(extension_fields[field].offset >= previous_end);
    AMPS_EXPECT_EQ(extension_fields[field].offset %
                       static_cast<long int>(extension_fields[field].alignment),
                   0);
    previous_end = extension_fields[field].offset + extension_fields[field].length;
  }
#endif

  // The configured core-test application enables an IDF model whose accessors
  // dereference doubles directly.  This checks a real production extension in
  // addition to the deliberately awkward anonymous field sizes above.
  AMPS_EXPECT_TRUE(PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_ >= previous_end);
  AMPS_EXPECT_EQ(PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_ %
                     static_cast<int>(alignof(double)),
                 0);
  AMPS_EXPECT_EQ(PIC::IDF::LB::_VIBRATIONAL_ENERGY_OFFSET_,
                 PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_ +
                     static_cast<int>(sizeof(double)));

  ParticleBufferScope scope(4);
  AMPS_EXPECT_TRUE(PIC::ParticleBuffer::ParticleDataLength >= previous_end);
#if _STATE_VECTOR_MODE_ == _STATE_VECTOR_MODE_ALIGNED_
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::ParticleDataLength % 16, 0);
  for (long int particle = 0; particle < 2; ++particle) {
    PIC::ParticleBuffer::byte* record =
        PIC::ParticleBuffer::GetParticleDataPointer(particle);
    for (std::size_t field = 0;
         field < sizeof(extension_fields) / sizeof(extension_fields[0]); ++field) {
      const std::uintptr_t address = reinterpret_cast<std::uintptr_t>(
          record + extension_fields[field].offset);
      AMPS_EXPECT_EQ(address % extension_fields[field].alignment,
                     static_cast<std::uintptr_t>(0));
    }
  }
#else
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::ParticleDataLength %
                     _ALIGN_STATE_VECTORS_BASE_, 0);
#endif
}

AMPS_TEST_CASE(AMPS_LAYOUT_CASE_ID("CLONE-001"),
               "clone copies every enabled payload byte and preserves links") {
  ParticleBufferScope scope(8);
  long int head = -1;
  const long int source = PIC::ParticleBuffer::GetNewParticle(head, true);
  const long int copy = PIC::ParticleBuffer::GetNewParticle(head, true);
  PIC::ParticleBuffer::byte* source_data =
      PIC::ParticleBuffer::GetParticleDataPointer(source);
  PIC::ParticleBuffer::byte* copy_data =
      PIC::ParticleBuffer::GetParticleDataPointer(copy);

  // Fill the complete physical payload, including every compile-time and
  // runtime optional field.  Intrusive links are excluded because clone is
  // required to retain the target record's list membership.
  for (long int i = _PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_;
       i < PIC::ParticleBuffer::ParticleDataLength; ++i)
    source_data[i] = static_cast<PIC::ParticleBuffer::byte>(
        1 + (37 * i) % 251);
  PIC::ParticleBuffer::SetI(0, source_data);
  PIC::ParticleBuffer::SetInitFlag(true, source_data);
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(2.75, source_data);
#endif

  const long int saved_prev = PIC::ParticleBuffer::GetPrev(copy_data);
  const long int saved_next = PIC::ParticleBuffer::GetNext(copy_data);
  PIC::ParticleBuffer::CloneParticle(copy, source);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetPrev(copy_data), saved_prev);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetNext(copy_data), saved_next);
  AMPS_EXPECT_EQ(std::memcmp(
                     copy_data + _PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_,
                     source_data + _PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_,
                     static_cast<std::size_t>(PIC::ParticleBuffer::ParticleDataLength -
                         _PIC_PARTICLE_DATA__SPECIES_ID_OFFSET_)),
                 0);
}

AMPS_TEST_CASE(AMPS_LAYOUT_CASE_ID("REUSE-001"),
               "deleted records clear initialization and stale payload state") {
  ParticleBufferScope scope(4);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
  PIC::ParticleBuffer::byte* record =
      PIC::ParticleBuffer::GetParticleDataPointer(particle);
  PIC::ParticleBuffer::SetI(0, record);
  PIC::ParticleBuffer::SetInitFlag(true, record);
#if _STATE_VECTOR_MODE_ == _STATE_VECTOR_MODE_ALIGNED_
  for (std::size_t field = 0;
       field < sizeof(extension_fields) / sizeof(extension_fields[0]); ++field)
    std::memset(record + extension_fields[field].offset,
                static_cast<int>(0x41 + field),
                static_cast<std::size_t>(extension_fields[field].length));
#endif
  double position[3] = {9.0, -8.0, 7.0};
  double velocity[3] = {-6.0, 5.0, -4.0};
  PIC::ParticleBuffer::SetX(position, particle);
  PIC::ParticleBuffer::SetV(velocity, particle);
  if (PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_ >= 0)
    PIC::IDF::LB::SetRotE(7.5e-21, record);
  if (PIC::IDF::LB::_VIBRATIONAL_ENERGY_OFFSET_ >= 0 &&
      PIC::IDF::nTotalVibtationalModes[0] > 0)
    PIC::IDF::LB::SetVibE(3.25e-21, 0, record);

  PIC::ParticleBuffer::DeleteParticle(particle);
  AMPS_EXPECT_FALSE(PIC::ParticleBuffer::TestInitFlag(record));
  const long int reused = PIC::ParticleBuffer::GetNewParticle(true);
  AMPS_EXPECT_EQ(reused, particle);
  AMPS_EXPECT_FALSE(PIC::ParticleBuffer::TestInitFlag(record));
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetPrev(record), -1);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetNext(record), -1);
  const double* reused_position = PIC::ParticleBuffer::GetX(reused);
  const double* reused_velocity = PIC::ParticleBuffer::GetV(reused);
  for (int dimension = 0; dimension < 3; ++dimension) {
    AMPS_EXPECT_NEAR(reused_position[dimension], 0.0, 0.0);
    AMPS_EXPECT_NEAR(reused_velocity[dimension], 0.0, 0.0);
  }
  if (PIC::IDF::LB::_ROTATIONAL_ENERGY_OFFSET_ >= 0)
    AMPS_EXPECT_NEAR(PIC::IDF::LB::GetRotE(record), 0.0, 0.0);
  if (PIC::IDF::LB::_VIBRATIONAL_ENERGY_OFFSET_ >= 0 &&
      PIC::IDF::nTotalVibtationalModes[0] > 0)
    AMPS_EXPECT_NEAR(PIC::IDF::LB::GetVibE(0, record), 0.0, 0.0);
  ExpectExtensionFieldsZero(reused);
}

AMPS_TEST_CASE(AMPS_LAYOUT_CASE_ID("LIST-001"),
               "middle deletion and free-record reuse preserve list integrity") {
  ParticleBufferScope scope(8);
  long int head = -1;
  const long int first = PIC::ParticleBuffer::GetNewParticle(head, true);
  const long int middle = PIC::ParticleBuffer::GetNewParticle(head, true);
  const long int last = PIC::ParticleBuffer::GetNewParticle(head, true);
  AMPS_EXPECT_EQ(head, last);
  PIC::ParticleBuffer::DeleteParticle(middle, head);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetNext(last), first);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetPrev(first), last);

  const long int reused = PIC::ParticleBuffer::GetNewParticle(head, true);
  AMPS_EXPECT_EQ(reused, middle);
  std::vector<long int> visited;
  for (long int current = head; current >= 0;
       current = PIC::ParticleBuffer::GetNext(current)) {
    AMPS_EXPECT_TRUE(visited.size() < 3);
    for (std::size_t i = 0; i < visited.size(); ++i)
      AMPS_EXPECT_TRUE(visited[i] != current);
    if (!visited.empty())
      AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetPrev(current), visited.back());
    else
      AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetPrev(current), -1);
    AMPS_EXPECT_TRUE(PIC::ParticleBuffer::IsParticleAllocated(current));
    visited.push_back(current);
  }
  AMPS_EXPECT_EQ(visited.size(), static_cast<std::size_t>(3));
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 3);
}

AMPS_TEST_CASE(AMPS_LAYOUT_CASE_ID("CAPACITY-001"),
               "empty and exactly-full buffer states remain consistent") {
  ParticleBufferScope scope(4);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 0);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::FirstPBufferParticle, 0);
  std::vector<long int> particles;
  for (int i = 0; i < 4; ++i)
    particles.push_back(PIC::ParticleBuffer::GetNewParticle(true));
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 4);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::FirstPBufferParticle, -1);
  for (std::size_t i = 0; i < particles.size(); ++i)
    AMPS_EXPECT_TRUE(PIC::ParticleBuffer::IsParticleAllocated(particles[i]));
  for (std::size_t i = 0; i < particles.size(); ++i)
    PIC::ParticleBuffer::DeleteParticle(particles[i]);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 0);
  AMPS_EXPECT_TRUE(PIC::ParticleBuffer::FirstPBufferParticle >= 0);
}

AMPS_TEST_CASE(AMPS_LAYOUT_CASE_ID("INVALID-001"),
               "invalid extension and record access fail without memory access") {
  long int offset = -73;
  AMPS_EXPECT_TRUE(Throws<std::invalid_argument>([&]() {
    PIC::ParticleBuffer::RequestDataStorage(offset, 0, 1);
  }));
  AMPS_EXPECT_EQ(offset, -73);
  AMPS_EXPECT_TRUE(Throws<std::invalid_argument>([&]() {
    PIC::ParticleBuffer::RequestDataStorage(offset, 8, 3);
  }));
  AMPS_EXPECT_EQ(offset, -73);

  ParticleBufferScope scope(4);
  AMPS_EXPECT_TRUE(Throws<std::logic_error>([&]() {
    PIC::ParticleBuffer::RequestDataStorage(offset, 8, 8);
  }));
  AMPS_EXPECT_EQ(offset, -73);
  AMPS_EXPECT_TRUE(Throws<std::out_of_range>([]() {
    PIC::ParticleBuffer::GetParticleDataPointer(-1);
  }));
  AMPS_EXPECT_TRUE(Throws<std::out_of_range>([]() {
    PIC::ParticleBuffer::GetParticleDataPointer(
        PIC::ParticleBuffer::GetMaxNPart());
  }));
}

AMPS_TEST_CASE(AMPS_LAYOUT_CASE_ID("EXHAUSTION-001"),
               "one allocation beyond capacity reaches the fatal guard") {
  // This is intentionally a death case.  The outer Python runner must launch
  // this ID in its own mpiexec process and accept it only when the process is
  // nonzero and emits the exact production diagnostic below.  Keeping the
  // abort out of the multi-case process prevents MPI_Abort from suppressing
  // receipts for unrelated tests.
  ParticleBufferScope scope(2);
  PIC::ParticleBuffer::GetNewParticle(true);
  PIC::ParticleBuffer::GetNewParticle(true);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 2);
  // The generic AMPS error code defaults to zero.  Give MPI_Abort an explicit
  // nonzero function code so the outer death-case runner can require a
  // nonzero process result rather than relying on launcher-specific treatment
  // of MPI_Abort(...,0).
  SetExitErrorCode(__LINE__, 97);
  PIC::ParticleBuffer::GetNewParticle(true);
  throw AMPS::Test::Failure(
      "full particle buffer unexpectedly allowed another allocation");
}

#undef AMPS_LAYOUT_CASE_ID
