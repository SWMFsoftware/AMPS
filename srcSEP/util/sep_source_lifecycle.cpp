#include "sep_source_lifecycle.h"

#include <algorithm>
#include <limits>
#include <utility>

// One monotone event index belongs to each physical (kind,line,species)
// channel. Source calls reserve events collectively before any owner-rank
// return; empty events still advance their channel. Macroparticle identities
// and purpose-tagged draws exclude sampled momentum, MPI rank and wall time.
// These value/I/O helpers do not identify moving/remapped lines on their own:
// the native bridge supplies and revalidates its exact fixed-geometry namespace.
// FILE serialization uses explicit little-endian uint64 words and staged reads;
// the enclosing full checkpoint owns atomic publication and checksum validation.
namespace SEP { namespace SourceLifecycle {
namespace {
Transport::Status Invalid(const char* text) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument, text);
}
bool Less(const Channel& a, const Channel& b) {
  if (a.kind != b.kind) return a.kind < b.kind;
  if (a.line != b.line) return a.line < b.line;
  return a.species < b.species;
}
bool Valid(const Channel& channel) {
  const auto kind=static_cast<std::uint64_t>(channel.kind);
  return kind >= 1 && kind <= 7;
}
void Append(std::uint64_t word, std::uint64_t* hash) {
  // Explicit byte order makes keys independent of host integer endianness.
  for (unsigned byte=0;byte<8;++byte) {
    *hash ^= (word>>(8*byte))&UINT64_C(255);
    *hash *= UINT64_C(1099511628211);
  }
}
std::uint64_t Hash(const Event& event,std::uint64_t macro,std::uint64_t purpose) {
  std::uint64_t hash=UINT64_C(1469598103934665603);
  Append(UINT64_C(0x534550534f555243),&hash);
  Append(1,&hash);Append(event.campaign,&hash);
  Append(static_cast<std::uint64_t>(event.channel.kind),&hash);
  Append(event.channel.line,&hash);Append(event.channel.species,&hash);
  Append(event.index,&hash);Append(macro,&hash);Append(purpose,&hash);
  return hash == 0 ? UINT64_C(1) : hash;
}
bool Put(FILE* file,std::uint64_t value) {
  unsigned char bytes[8];
  for(unsigned i=0;i<8;++i)bytes[i]=static_cast<unsigned char>(value>>(8*i));
  return std::fwrite(bytes,1,8,file)==8;
}
bool Get(FILE* file,std::uint64_t* value) {
  unsigned char bytes[8];
  if(std::fread(bytes,1,8,file)!=8)return false;
  *value=0;
  for(unsigned i=0;i<8;++i)*value|=static_cast<std::uint64_t>(bytes[i])<<(8*i);
  return true;
}
Transport::Status Validate(const State& state) {
  if(state.schema!=1||state.fixedLineLayout.empty())
    return Invalid("incompatible source scheduler schema or missing fixed-line namespace");
  for(std::size_t i=0;i<state.counters.size();++i) {
    if(!Valid(state.counters[i].channel)||
       (i&&!Less(state.counters[i-1].channel,state.counters[i].channel)))
      return Invalid("source scheduler channels must be valid, unique, and canonically ordered");
  }
  return Transport::Status::Ok();
}
}

Transport::Status Scheduler::Initialize(std::uint64_t campaign,const std::string& layout) {
  if(layout.empty())return Invalid("source scheduler needs an explicit fixed-line namespace");
  if(initialized_) {
    if(state_.campaign!=campaign||state_.fixedLineLayout!=layout)
      return Invalid("source scheduler campaign or physical line layout changed; stable remapping is unsupported");
    return Transport::Status::Ok();
  }
  state_.campaign=campaign;state_.fixedLineLayout=layout;initialized_=true;
  return Transport::Status::Ok();
}
Transport::Status Scheduler::Begin(const Channel& channel,Event* event) {
  // Keep canonical channel ordering independent of source call enumeration.
  // Exhaustion is an explicit error before wrapping an existing event history;
  // rebuilding event zero would reuse physical births and random realizations.
  if(!initialized_||!event||!Valid(channel))return Invalid("uninitialized or invalid source event");
  auto at=std::lower_bound(state_.counters.begin(),state_.counters.end(),channel,
      [](const Counter& counter,const Channel& key){return Less(counter.channel,key);});
  if(at==state_.counters.end()||Less(channel,at->channel)) {
    Counter counter;counter.channel=channel;
    at=state_.counters.insert(at,counter);
  }
  if(at->next==std::numeric_limits<std::uint64_t>::max())
    return Transport::Status::Error(Transport::StatusCode::StepUnderflow,
        "source event identity counter exhausted; refusing wraparound");
  Event staged;staged.channel=channel;staged.campaign=state_.campaign;staged.index=at->next;
  ++at->next;*event=staged;
  return Transport::Status::Ok();
}
State Scheduler::Capture() const {return state_;}
void Scheduler::Swap(Scheduler* other) {
  std::swap(initialized_,other->initialized_);std::swap(state_,other->state_);
}
Transport::Status Scheduler::Restore(const State& state,std::uint64_t campaign,
                                      const std::string& layout) {
  // Compatibility means the same campaign and physical line namespace, not
  // merely an empty current population. Native restore stages a private owner
  // before swapping it into the collective checkpoint commit.
  const auto status=Validate(state);if(!status.ok())return status;
  if(state.campaign!=campaign||state.fixedLineLayout!=layout)
    return Invalid("source checkpoint campaign or fixed physical line namespace mismatch");
  state_=state;initialized_=true;return Transport::Status::Ok();
}
std::uint64_t ParticleIdentity(const Event& event,std::uint64_t macro) {
  return Hash(event,macro,UINT64_C(0x42495254484944));
}
Transport::KeyedRandomStream Random(const Event& event,std::uint64_t macro,
                                    std::uint64_t purpose) {
  return Transport::KeyedRandomStream(event.campaign,ParticleIdentity(event,macro),
      Hash(event,macro,purpose),event.index);
}
Transport::Status WriteState(FILE* file,const State& state) {
  const auto status=Validate(state);if(!status.ok())return status;
  if(!file||!Put(file,UINT64_C(0x5345505352433031))||!Put(file,state.schema)||
      !Put(file,state.campaign)||!Put(file,state.fixedLineLayout.size())||
      std::fwrite(state.fixedLineLayout.data(),1,state.fixedLineLayout.size(),file)!=state.fixedLineLayout.size()||
      !Put(file,state.counters.size()))return Invalid("cannot write source scheduler checkpoint");
  for(const auto& counter:state.counters)
    if(!Put(file,static_cast<std::uint64_t>(counter.channel.kind))||
       !Put(file,counter.channel.line)||!Put(file,counter.channel.species)||!Put(file,counter.next))
      return Invalid("cannot write source scheduler counters");
  return Transport::Status::Ok();
}
Transport::Status ReadState(FILE* file,State* state) {
  State staged;std::uint64_t magic=0,length=0,count=0;
  if(!file||!state||!Get(file,&magic)||magic!=UINT64_C(0x5345505352433031)||
      !Get(file,&staged.schema)||staged.schema!=1||!Get(file,&staged.campaign)||!Get(file,&length))
    return Invalid("incompatible or truncated source scheduler checkpoint header");
  // Check the remaining file extent before allocating untrusted lengths. This
  // format is embedded only in seekable, validated corrected-model files.
  const long begin=std::ftell(file);
  if(begin<0||std::fseek(file,0,SEEK_END)!=0)return Invalid("source checkpoint must be seekable");
  const long end=std::ftell(file);
  if(end<begin||std::fseek(file,begin,SEEK_SET)!=0||length>static_cast<std::uint64_t>(end-begin))
    return Invalid("source checkpoint line namespace exceeds payload");
  staged.fixedLineLayout.resize(static_cast<std::size_t>(length));
  if((length&&std::fread(&staged.fixedLineLayout[0],1,length,file)!=length)||!Get(file,&count))
    return Invalid("truncated source checkpoint line namespace");
  const long countersBegin=std::ftell(file);
  if(countersBegin<0||end<countersBegin||count>static_cast<std::uint64_t>(end-countersBegin)/32)
    return Invalid("source checkpoint counter table exceeds payload");
  staged.counters.resize(static_cast<std::size_t>(count));
  for(auto& counter:staged.counters) {
    std::uint64_t kind=0;
    if(!Get(file,&kind)||!Get(file,&counter.channel.line)||!Get(file,&counter.channel.species)||!Get(file,&counter.next))
      return Invalid("truncated source checkpoint counter table");
    counter.channel.kind=static_cast<Kind>(kind);
  }
  const auto status=Validate(staged);if(!status.ok())return status;
  *state=std::move(staged);return Transport::Status::Ok();
}
} }  // namespace SEP::SourceLifecycle
