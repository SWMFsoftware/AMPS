#include "pic.h"
#include "src/test/native/test_harness.h"
#include "core_test_state.h"

using AMPS::CoreTests::Fixture::ParticleBufferScope;

namespace AMPS {
namespace CoreTests {

// See main.cpp.  The function intentionally has no run-time work; its link
// reference makes this translation unit and its AutoRegister objects members
// of the final executable even though the application is linked as an archive.
void LinkParticleBufferTests() {}

}  // namespace CoreTests
}  // namespace AMPS

AMPS_TEST_CASE("PIC-PBUFFER-INIT-001", "particle-buffer allocation and deletion counters") {
  // Allocation with the initialization flag must produce an addressable
  // record and increment NAllPart exactly once.  Deletion must clear the
  // allocation marker and reverse that accounting without rebuilding state.
  ParticleBufferScope scope(64);
  AMPS_EXPECT_TRUE(PIC::ParticleBuffer::ParticleDataBuffer != NULL);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetMaxNPart(), 64);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 0);
  const long int particle = PIC::ParticleBuffer::GetNewParticle(true);
  AMPS_EXPECT_TRUE(particle >= 0);
  AMPS_EXPECT_TRUE(PIC::ParticleBuffer::IsParticleAllocated(particle));
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 1);
  PIC::ParticleBuffer::DeleteParticle(particle);
  AMPS_EXPECT_FALSE(PIC::ParticleBuffer::IsParticleAllocated(particle));
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetAllPartNum(), 0);
}

AMPS_TEST_CASE("PIC-PBUFFER-STATE-001", "particle species, position, velocity, and clone state") {
  // Coordinates are deliberately asymmetric/nonintegral so swapped fields,
  // dimensions, or byte offsets cannot accidentally satisfy the comparison.
  ParticleBufferScope scope(64);
  const long int source = PIC::ParticleBuffer::GetNewParticle(true);
  double position[3] = {1.25, -2.5, 7.0};
  double velocity[3] = {-4.0, 5.5, 6.25};
  PIC::ParticleBuffer::SetI(0, source);
  PIC::ParticleBuffer::SetX(position, source);
  PIC::ParticleBuffer::SetV(velocity, source);
  const long int copy = PIC::ParticleBuffer::GetNewParticle(true);
  PIC::ParticleBuffer::CloneParticle(copy, source);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetI(copy), 0);
  const double* copied_position = PIC::ParticleBuffer::GetX(copy);
  const double* copied_velocity = PIC::ParticleBuffer::GetV(copy);
  for (int dimension = 0; dimension < 3; ++dimension) {
    AMPS_EXPECT_NEAR(copied_position[dimension], position[dimension], 0.0);
    AMPS_EXPECT_NEAR(copied_velocity[dimension], velocity[dimension], 0.0);
  }
}

AMPS_TEST_CASE("PIC-PBUFFER-LIST-001", "intrusive particle-list insertion and exclusion") {
  // GetNewParticle(first,...) prepends to an intrusive list.  Removing its
  // head must make the prior head visible and preserve the terminal -1 link.
  ParticleBufferScope scope(64);
  long int first = -1;
  const long int first_inserted = PIC::ParticleBuffer::GetNewParticle(first, true);
  const long int second_inserted = PIC::ParticleBuffer::GetNewParticle(first, true);
  AMPS_EXPECT_EQ(first, second_inserted);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetNext(second_inserted), first_inserted);
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetNext(first_inserted), -1);
  PIC::ParticleBuffer::ExcludeParticleFromList(second_inserted, first);
  AMPS_EXPECT_EQ(first, first_inserted);
}

AMPS_TEST_CASE("PIC-PBUFFER-SIGNATURE-001", "canonical and list-inclusive signatures preserve full-width legacy bytes without changing records") {
  ParticleBufferScope scope(16);
  const long int source=PIC::ParticleBuffer::GetNewParticle(true);
  auto* record=PIC::ParticleBuffer::GetParticleDataPointer(source);
  const std::size_t length=PIC::ParticleBuffer::ParticleDataLength;
  // Pattern every compiled byte, including anonymous extensions and padding.
  // The record is not attached to a mesh list; links are deliberate byte data
  // for the inclusion/exclusion contract and restored by the fixture lifetime.
  for (std::size_t i=0;i<length;i++) record[i]=static_cast<unsigned char>(31*i+7);
  PIC::ParticleBuffer::SetPrev(-1,record);
  PIC::ParticleBuffer::SetNext(-1,record);
  PIC::ParticleBuffer::SetParticleAllocated(record);
  const std::vector<unsigned char> before(record,record+length);
  std::vector<unsigned char> canonical(length,0);
  // This unmodified clone API is the previous canonicalization algorithm:
  // zero destination bytes, copy all source bytes, preserve destination links.
  PIC::ParticleBuffer::CloneParticle(canonical.data(),record);
  auto reference=[](const unsigned char* bytes,std::size_t count,unsigned long state) {
    for (std::size_t i=0;i<count;i++) {
      unsigned long term=(((state>>24)^bytes[i])&255UL)<<24;
      for (int bit=0;bit<8;bit++) term=(term&0x80000000UL) ? (term<<1)^0x04c11db7UL:term<<1;
      state=(state<<8)^term;
    }
    return state;
  };
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetParticleSignature(source,false),
      reference(canonical.data(),length,0));
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetParticleSignature(source,true),
      reference(before.data(),length,0));
  CRC32 incremental; incremental.crc_accum=~0UL;
  AMPS_EXPECT_EQ(PIC::ParticleBuffer::GetParticleSignature(source,&incremental,false),
      reference(canonical.data(),length,~0UL));
  AMPS_EXPECT_TRUE(std::memcmp(record,before.data(),length)==0);
}
