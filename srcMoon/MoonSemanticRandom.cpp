#include "MoonSemanticRandom.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
namespace Moon { namespace Chemistry { namespace SemanticRandom {
namespace {
constexpr std::uint64_t Magic=UINT64_C(0x4d3353454d4b4559);
constexpr std::size_t HeaderWords=9; // magic,config-fingerprint,parent,depth,4 counters,capacity
long int offset=-1;std::unique_ptr<const Configuration> active;
bool Reject(std::string* e,const std::string& s){if(e)*e=s;return false;}
void Put(std::vector<unsigned char>& b,std::size_t i,std::uint64_t v){
  // Byte-level big endian format is deterministic and never relies on alignment.
  for(unsigned j=0;j<8;++j)b[8*i+j]=static_cast<unsigned char>(v>>(56-8*j));}
std::uint64_t Get(const unsigned char*b,std::size_t i){std::uint64_t v=0;for(unsigned j=0;j<8;++j)v=(v<<8)|b[8*i+j];return v;}
bool Fingerprint(const Configuration& c,std::uint64_t* out,std::string* e){
  Q::Key key;key.campaignSeed=c.seed;key.cohort=c.cohort;key.originalParent=0;
  // Capacity and seed/cohort are part of opaque record interpretation. This
  // integrity tag is not a substitute for the full retained SHA provenance.
  key.genealogy={static_cast<std::uint64_t>(c.genealogyCapacity)};
  // Domain-separate only the storage fingerprint. Every physical draw uses
  // the same seed, semantic key and counters regardless of clock arithmetic.
  if(c.clockPrecision!=Species::ReactionExecution::OpticalDepthPrecision::Binary64Log)key.genealogy.push_back(1);
  std::string hash;
  if(!Q::DigestHex(key,Q::Purpose::SourceClock,0,&hash,e))return false;
  std::uint64_t x=0;for(unsigned i=0;i<16;++i){const char h=hash[i];x=(x<<4)|static_cast<unsigned>(h<='9'?h-'0':h-'a'+10);}*out=x;return true;
}
bool CheckParticle(long int p,std::string* e){
  return (active&&offset>=0&&p>=0&&p<PIC::ParticleBuffer::GetMaxNPart()&&PIC::ParticleBuffer::IsParticleAllocated(p))||
    Reject(e,"semantic RNG needs registered storage and an actually allocated particle");}
bool Write(long int p,const State& state,std::string* e){
  if(!CheckParticle(p,e))return false;std::vector<unsigned char> bytes;
  if(!Encode(*active,state,&bytes,e))return false;
  std::memcpy(PIC::ParticleBuffer::GetParticleDataPointer(p)+offset,bytes.data(),bytes.size());return true;
}
bool Purpose(Controller::DrawPurpose p,Q::Purpose* q,std::string* e){
  switch(p){case Controller::DrawPurpose::SourceClock:*q=Q::Purpose::SourceClock;return true;
    case Controller::DrawPurpose::Channel:*q=Q::Purpose::Channel;return true;
    case Controller::DrawPurpose::ProductClock:*q=Q::Purpose::ProductClock;return true;
    case Controller::DrawPurpose::Kinematics:*q=Q::Purpose::Kinematics;return true;}
  return Reject(e,"unknown controller semantic RNG purpose");
}
}
bool RecordSize(const Configuration& c,std::size_t* n,std::string* e){
  if(!n||!e)return false;
  if(!c.genealogyCapacity||c.genealogyCapacity>(static_cast<std::size_t>(std::numeric_limits<int>::max())/8)-HeaderWords)
    return Reject(e,"semantic genealogy capacity must be explicit positive and fit PIC integer record size; increase capacity only within available memory");
  std::uint64_t tag=0;if(!Fingerprint(c,&tag,e))return false;
  *n=8*(HeaderWords+c.genealogyCapacity);e->clear();return true;
}
bool Encode(const Configuration& c,const State& s,std::vector<unsigned char>* out,std::string* e){
  if(!out||!e)return false;std::size_t n=0;if(!RecordSize(c,&n,e))return false;
  if(s.genealogy.size()>c.genealogyCapacity)return Reject(e,"full semantic genealogy exceeds declared opaque record capacity. No reaction was published. Increase genealogy capacity, rebuild/restart from the declared initial state; never truncate/fold the genealogy or discard a branch");
  std::uint64_t tag=0;if(!Fingerprint(c,&tag,e))return false;std::vector<unsigned char> b(n,0);
  Put(b,0,Magic);Put(b,1,tag);Put(b,2,s.originalParent);Put(b,3,s.genealogy.size());
  for(unsigned j=0;j<4;++j)Put(b,4+j,s.acceptedDraws[j]);Put(b,8,c.genealogyCapacity);
  for(std::size_t j=0;j<s.genealogy.size();++j)Put(b,HeaderWords+j,s.genealogy[j]);
  out->swap(b);e->clear();return true;
}
bool Decode(const Configuration& c,const unsigned char* b,std::size_t length,State* out,std::string* e){
  if(!out||!e)return false;std::size_t n=0;if(!RecordSize(c,&n,e))return false;std::uint64_t tag=0;
  if(!Fingerprint(c,&tag,e))return false;
  if(!b||length!=n||Get(b,0)!=Magic||Get(b,1)!=tag||Get(b,8)!=c.genealogyCapacity||Get(b,3)>c.genealogyCapacity)
    return Reject(e,"semantic particle record magic/configuration/capacity/size is incompatible; ownership clones must preserve complete opaque bytes");
  State s;s.originalParent=Get(b,2);const auto depth=static_cast<std::size_t>(Get(b,3));
  for(unsigned j=0;j<4;++j)s.acceptedDraws[j]=Get(b,4+j);
  for(std::size_t j=0;j<depth;++j)s.genealogy.push_back(Get(b,HeaderWords+j));
  for(std::size_t j=depth;j<c.genealogyCapacity;++j)if(Get(b,HeaderWords+j)!=0)
    return Reject(e,"semantic unused genealogy bytes are not canonical zero");
  *out=std::move(s);e->clear();return true;
}
bool Draw(const Configuration& c,Q::Purpose purpose,State* state,double* out,std::string* e){
  if(!state||!out||!e)return false;
  const auto raw=static_cast<std::uint64_t>(purpose);if(raw<1||raw>4)return Reject(e,"unknown semantic RNG purpose");
  std::size_t bytes=0;if(!RecordSize(c,&bytes,e))return false;
  if(state->genealogy.size()>c.genealogyCapacity)return Reject(e,"semantic key exceeds configured genealogy capacity");
  const auto index=static_cast<std::size_t>(raw-1);
  if(state->acceptedDraws[index]==std::numeric_limits<std::uint64_t>::max())
    return Reject(e,"semantic draw counter exhausted uint64 domain; restart with a predeclared finite cohort/stream plan rather than reusing a draw");
  Q::Key key;key.campaignSeed=c.seed;key.cohort=c.cohort;key.originalParent=state->originalParent;key.genealogy=state->genealogy;
  double value=0;if(!Q::OpenUniform(key,purpose,state->acceptedDraws[index],&value,e))return false;
  ++state->acceptedDraws[index];*out=value;e->clear();return true;
}
bool ChildState(const Configuration& c,const State& parent,std::uint64_t ordinal,bool used,State* out,std::string* e){
  if(!out||!e)return false;std::size_t n=0;if(!RecordSize(c,&n,e))return false;
  if(parent.genealogy.size()>=c.genealogyCapacity)return Reject(e,"reaction cannot append the complete daughter genealogy within declared storage. No product/ledger may be published. Increase genealogyCapacity and restart; never censor the branch, shorten/hash-fold the key or redraw the event");
  State child;child.originalParent=parent.originalParent;child.genealogy=parent.genealogy;child.genealogy.push_back(ordinal);
  if(used)child.acceptedDraws[2]=1;*out=std::move(child);e->clear();return true;
}
bool ProspectiveProductClock(const Configuration& c,const State& p,std::uint64_t ordinal,double* out,std::string* e){
  if(!out||!e)return false;State child;if(!ChildState(c,p,ordinal,false,&child,e))return false;
  return Draw(c,Q::Purpose::ProductClock,&child,out,e);
}
bool RegisterStorage(const Configuration& c,std::string* e){
  if(!e)return false;if(active||offset>=0||PIC::ParticleBuffer::ParticleDataBuffer)
    return Reject(e,"register semantic RNG storage exactly once before allocating the PIC particle buffer");
  std::size_t n=0;if(!RecordSize(c,&n,e))return false;
  auto copy=std::unique_ptr<const Configuration>(new Configuration(c));
  PIC::ParticleBuffer::RequestDataStorage(offset,static_cast<int>(n),1);active=std::move(copy);e->clear();return true;
}
bool Registered(){return static_cast<bool>(active);}
bool InitializeSource(long int p,std::uint64_t id,std::string* e){State s;s.originalParent=id;return Write(p,s,e);}
bool Read(long int p,State* out,std::string* e){
  if(!out||!e||!CheckParticle(p,e))return false;std::size_t n=0;if(!RecordSize(*active,&n,e))return false;
  return Decode(*active,PIC::ParticleBuffer::GetParticleDataPointer(p)+offset,n,out,e);
}
bool ProspectiveProductClockForParticle(long int ptr,std::uint64_t ordinal,double* out,std::string* e){
  State parent;if(!Read(ptr,&parent,e))return false;
  return ProspectiveProductClock(*active,parent,ordinal,out,e);
}
bool DrawForParticle(Controller::DrawPurpose p,long int ptr,std::size_t,void*,double* out,std::string* e){
  if(!out||!e)return false;
  if(p==Controller::DrawPurpose::ProductClock)
    return Reject(e,"daughter clocks require ProspectiveProductClockForParticle with explicit canonical product ordinal; parent-key ProductClock is forbidden");
  Q::Purpose q;if(!Purpose(p,&q,e))return false;State state;
  if(!Read(ptr,&state,e))return false;double value=0;if(!Draw(*active,q,&state,&value,e)||!Write(ptr,state,e))return false;
  *out=value;e->clear();return true;
}
bool InitializeStagedProduct(const Species::ReactionExecution::EventPlan& plan,std::size_t ordinal,long int p,long int child,void*,std::string* e){
  State parent,next;if(!Read(p,&parent,e)||!ChildState(*active,parent,ordinal,true,&next,e))return false;
  double u=0;Species::ReactionExecution::OpticalDepthClock expected;
  if(ordinal>=plan.products().size()||!ProspectiveProductClock(*active,parent,ordinal,&u,e)||
      !Species::ReactionExecution::OpticalDepthClock::FromOpenUniform(u,&expected,e,active->clockPrecision))return false;
  const auto& actual=plan.products()[ordinal].clock;
  if(!actual.fresh()||actual.consumed()||!actual.SamePendingDepth(expected))
    return Reject(e,"staged daughter's supplied optical clock was not constructed from its declared full semantic key/draw0; publication refused");
  return Write(child,next,e);
}
}}}
