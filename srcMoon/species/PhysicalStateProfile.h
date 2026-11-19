#ifndef SRCMOON_SPECIES_PHYSICAL_STATE_PROFILE_H
#define SRCMOON_SPECIES_PHYSICAL_STATE_PROFILE_H

#include "ReducedKinematics.h"
#include <map>
#include <string>
#include <vector>

namespace Moon { namespace Species { namespace PhysicalStateProfile {
namespace K=ReducedKinematics;

// Source authentication is independent of source authority/use permission.
// These are explicit original research or separately permitted source paths;
// this loader neither imports them into a managed database nor grants use.
struct SourceFile {std::string id,path,sha256,role;};
struct PinnedCsv {std::string path,sha256;};
struct Files {
  std::string schema="srcMoon-M3-state-profile-v1";
  PinnedCsv isotopes,constants,states;
  std::vector<SourceFile> sources;
};
struct Constant {
  K::Real value=0;std::string unit,sourceIdentity,sourceRecord;
};
struct StateMetadata {
  std::string id,label,kind,missingReason;
  std::vector<ReactionExecution::NucleusCount> nuclei;
  std::int64_t charge=0;
  bool available=false,physicalRestMassAvailable=false;
};

// A loaded profile is an authenticated declaration, NOT a qualified reaction
// family, physical-error budget, permitted managed import, or I21--I24 PASS.
// All constructors keep physical activation outside this component. Unknown
// states/energies/rest masses never become zero or the additive reference mass.
class Registry {
 public:
  bool Load(const Files&,std::string* error);
  bool GetState(const std::string& id,K::ResolvedState* output,std::string* error)const;
  bool GetPhysicalRestMass(const std::string& id,K::MassValue* output,std::string* error)const;
  bool GetConstant(const std::string& id,const std::string& expectedUnit,
      Constant* output,std::string* error)const;
  bool loaded()const{return loaded_;}
  const Files& files()const{return files_;}
  const std::vector<StateMetadata>& states()const{return metadata_;}
  const std::vector<K::IsotopeMass>& isotopes()const{return isotopes_;}
  const std::string& identity()const{return identity_;}
 private:
  bool loaded_=false;Files files_;std::string identity_;
  std::vector<K::IsotopeMass> isotopes_;
  std::vector<StateMetadata> metadata_;
  std::map<std::string,K::ResolvedState> resolved_;
  std::map<std::string,K::MassValue> physical_;
  std::map<std::string,Constant> constants_;
};
}}}
#endif
