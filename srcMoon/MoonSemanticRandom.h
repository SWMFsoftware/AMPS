#ifndef SRCMOON_SEMANTIC_RANDOM_H
#define SRCMOON_SEMANTIC_RANDOM_H
#include "MoonReactionController.h"
#include "species/QualificationRandom.h"
#include <array>
#include <cstdint>
#include <string>
#include <vector>
namespace Moon { namespace Chemistry { namespace SemanticRandom {
namespace Q=Species::QualificationRandom;
struct Configuration {
  std::uint64_t seed=0;
  std::string cohort;
  std::size_t genealogyCapacity=0; // Explicit storage admission, NOT physics cap.
  // Stored in the opaque configuration fingerprint, never in the random key.
  // Changing log precision therefore preserves every parent/channel draw.
  Species::ReactionExecution::OpticalDepthPrecision clockPrecision=Species::ReactionExecution::OpticalDepthPrecision::Binary64Log;
};
struct State {
  std::uint64_t originalParent=0;
  std::vector<std::uint64_t> genealogy;
  std::array<std::uint64_t,4> acceptedDraws{{0,0,0,0}};
};
// Pure supplied-byte codec and key derivation, shared by actual PIC storage.
// No mutation on invalid shape/capacity/configuration. Bytes contain no pointer.
bool RecordSize(const Configuration&,std::size_t*,std::string* error);
bool Encode(const Configuration&,const State&,std::vector<unsigned char>*,std::string* error);
bool Decode(const Configuration&,const unsigned char*,std::size_t,State*,std::string* error);
bool Draw(const Configuration&,Q::Purpose,State*,double*,std::string* error);
bool ChildState(const Configuration&,const State&,std::uint64_t productOrdinal,
    bool productClockAlreadyDrawn,State*,std::string* error);
// Prospective child draw0 is PURE: no parent or future child counter changes.
bool ProspectiveProductClock(const Configuration&,const State&,std::uint64_t productOrdinal,
    double*,std::string* error);

// Opaque per-particle bridge: register before ParticleBuffer::Init; all actual
// PIC clones/MPI transfers preserve full genealogy/counters automatically.
bool RegisterStorage(const Configuration&,std::string* error);
bool Registered();
bool InitializeSource(long int,std::uint64_t originalGlobalParent,std::string* error);
bool Read(long int,State*,std::string* error);
bool ProspectiveProductClockForParticle(long int,std::uint64_t ordinal,double*,std::string* error);
bool DrawForParticle(Controller::DrawPurpose,long int,std::size_t,void*,double*,std::string*);
// Initially scoped to terminal single-ion Na channel. Its provider must use
// ProspectiveProductClock ordinal0 once; this initializer installs counter1.
// Multi-product callers use their explicit canonical ordinal and identical law.
bool InitializeStagedProduct(const Species::ReactionExecution::EventPlan&,
    std::size_t,long int,long int,void*,std::string*);
}}}
#endif
