#ifndef SRCMOON_SPECIES_QUALIFICATION_RANDOM_H
#define SRCMOON_SPECIES_QUALIFICATION_RANDOM_H
#include <cstdint>
#include <string>
#include <vector>
namespace Moon { namespace Species { namespace QualificationRandom {
// Semantic domains, never species slots or allocation/rank-dependent numbers.
enum class Purpose : std::uint64_t { SourceClock=1,Channel=2,ProductClock=3,Kinematics=4 };
struct Key {
  std::uint64_t campaignSeed=0,originalParent=0;
  // Matched timestep/MPI trials deliberately share a cohort ID; a distinct
  // physical cohort declares another ID before execution, never after scoring.
  std::string cohort;
  // NEW daughters append caller-declared event/daughter ordinals. Ownership
  // clones preserve the exact path, purpose and accepted draw counter.
  std::vector<std::uint64_t> genealogy;
};
// Stateless, transactional pure calls. They do not advance any clock/counter,
// access PIC/MPI, draw another event on retry or install a production provider.
// HMAC-SHA256 construction/encoding is normative in qualification/CONTRACT.md.
bool DigestHex(const Key&,Purpose,std::uint64_t drawIndex,std::string*,std::string* error);
// First52 bits -> exact midpoint (k+1/2)*2^-52, strictly0<u<1. The finite-grid
// law and omitted extreme optical-depth tail are explicit qualification scope.
bool OpenUniform(const Key&,Purpose,std::uint64_t drawIndex,double*,std::string* error);
} } }
#endif
