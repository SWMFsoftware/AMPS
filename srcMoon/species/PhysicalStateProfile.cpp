#include "PhysicalStateProfile.h"
#include "Photochemistry.h" // SAME shared production source/CSV authentication
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <utility>

namespace Moon { namespace Species { namespace PhysicalStateProfile {
namespace {
namespace D=ReactionData;
bool Fail(const std::string& message,std::string* error){if(error)*error=message;return false;}
bool Identifier(const std::string& text) {
  return !text.empty()&&std::all_of(text.begin(),text.end(),[](unsigned char c){
    return (c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='_'||c=='.'||c==':'||c=='-'||c=='+';});
}
bool Number(const std::string& text,K::Real* output,std::string* error) {
  if(text.empty())return Fail("missing numerical profile field; no zero substitution",error);
  std::istringstream stream(text);stream.imbue(std::locale::classic());
  K::Real value=0;stream>>std::noskipws>>value;
  if(!stream||stream.peek()!=std::char_traits<char>::eof()||!std::isfinite(value))
    return Fail("invalid/unrepresentable finite profile number: "+text,error);
  *output=value;return true;
}
bool Integer(const std::string& text,std::int64_t* output,bool positive,std::string* error) {
  if(text.empty()||text=="-0"||text.front()=='+'||
      (text.size()>1&&text.front()=='0')||
      (text.size()>2&&text[0]=='-'&&text[1]=='0'))
    return Fail("noncanonical profile integer: "+text,error);
  std::size_t begin=text.front()=='-'?1:0;
  if(begin==text.size()||!std::all_of(text.begin()+begin,text.end(),[](unsigned char c){return c>='0'&&c<='9';}))
    return Fail("invalid profile integer: "+text,error);
  std::istringstream stream(text);stream.imbue(std::locale::classic());std::int64_t value=0;
  stream>>std::noskipws>>value;
  if(!stream||stream.peek()!=std::char_traits<char>::eof()||(positive&&value<=0))
    return Fail("profile integer outside its valid range: "+text,error);
  *output=value;return true;
}
bool Read(const PinnedCsv& file,const D::Row& header,D::Csv* output,std::string* error) {
  std::string bytes,digest;
  if(!D::DigestSyntax(file.sha256)||!D::ReadFileBytes(file.path,0,&bytes,error)||
      !D::BytesSha256(bytes,&digest,error))return Fail("profile CSV requires explicit path/hash: "+file.path,error);
  if(digest!=file.sha256)return Fail("profile CSV SHA-256 changed: "+file.path,error);
  // Authenticate the SAME byte snapshot that the production CSV reader parses.
  return D::ParseCsv(bytes,header,output,error);
}
bool Nuclei(const std::string& text,const std::set<std::string>& isotopes,
    std::vector<ReactionExecution::NucleusCount>* output,std::string* error) {
  if(text.empty())return true; // Only an explicitly classified free electron may lack nuclei.
  std::set<std::string> seen;std::size_t begin=0;
  do {
    const auto end=text.find(';',begin);const auto term=text.substr(begin,end==std::string::npos?end:end-begin);
    const auto colon=term.find(':');std::int64_t count=0;
    if(colon==std::string::npos||term.find(':',colon+1)!=std::string::npos)
      return Fail("isotope stoichiometry requires explicit isotope:count pairs",error);
    const auto id=term.substr(0,colon);
    if(!isotopes.count(id)||!seen.insert(id).second||!Integer(term.substr(colon+1),&count,true,error))
      return Fail("unknown/duplicate isotope or invalid count in state: "+term,error);
    output->push_back({id,count});if(end==std::string::npos)break;begin=end+1;
  }while(begin<=text.size());
  return true;
}
using Sources=std::map<std::string,SourceFile>;
bool Source(const Sources& sources,const std::string& id,const std::string& record,
    const std::set<std::string>& allowed,std::string* identity,std::string* role,std::string* error) {
  const auto found=sources.find(id);
  if(found==sources.end()||record.empty()||!allowed.count(found->second.role))
    return Fail("missing source/record or incompatible source role: "+id,error);
  *identity=id+"@sha256:"+found->second.sha256+"#"+record;
  if(role)*role=found->second.role;
  return true;
}
bool Energy(const Sources& sources,const std::string& value,const std::string& source,
    const std::string& record,bool excitation,const StateMetadata& metadata,
    K::EnergyValue* output,std::string* error) {
  K::EnergyValue energy;std::string role;
  if(!Number(value,&energy.joules,error)||
      !Source(sources,source,record,{"energy","reference-definition"},&energy.sourceIdentity,&role,error)||
      (excitation&&energy.joules<0))return Fail("invalid explicit state energy/source: "+metadata.id,error);
  if(role=="reference-definition") {
    const bool electron=metadata.kind=="free-electron"&&metadata.nuclei.empty()&&metadata.charge==-1;
    const bool neutralAtom=metadata.nuclei.size()==1&&metadata.nuclei[0].count==1&&metadata.charge==0;
    if(energy.joules!=0||(!excitation&&!electron&&!neutralAtom)||
        (excitation&&metadata.kind!="ground"&&!electron))
      return Fail("reference definition supplies only approved atom/electron or ground-state exact zeros: "+metadata.id,error);
  }
  *output=std::move(energy);return true;
}
}

bool Registry::Load(const Files& files,std::string* error) {
  try {
    if(files.schema!="srcMoon-M3-state-profile-v1")return Fail("unsupported physical-state profile schema",error);
    Registry candidate;candidate.files_=files;Sources sources;
    for(const auto& source:files.sources) {
      if(!Identifier(source.id)||!D::DigestSyntax(source.sha256)||
          !std::set<std::string>{"isotope-mass","constants","energy","reference-definition"}.count(source.role)||
          !sources.emplace(source.id,source).second)
        return Fail("invalid/duplicate explicit physical-state source: "+source.id,error);
      if(!D::VerifyFileSha256(source.path,source.sha256,error))return false;
    }
    D::Csv masses,constants,states;
    if(!Read(files.isotopes,{"isotope","neutral_mass_kg","source_id","source_record"},&masses,error)||
        !Read(files.constants,{"name","value","unit","source_id","source_record"},&constants,error)||
        !Read(files.states,{"state_id","state_label","state_kind","nuclei","charge","formation_J",
        "formation_source_id","formation_source_record","excitation_J","excitation_source_id",
        "excitation_source_record","physical_rest_mass_kg","physical_mass_source_id",
        "physical_mass_source_record","status","missing_reason"},&states,error))return false;
    std::set<std::string> isotopeIds;
    for(const auto& row:masses) {
      K::IsotopeMass isotope;isotope.identity=row[0];
      if(!Identifier(row[0])||!isotopeIds.insert(row[0]).second||
          !Number(row[1],&isotope.neutralAtom.kg,error)||isotope.neutralAtom.kg<=0||
          !Source(sources,row[2],row[3],{"isotope-mass"},&isotope.neutralAtom.sourceIdentity,nullptr,error))
        return Fail("invalid source-bound neutral isotope mass: "+row[0],error);
      candidate.isotopes_.push_back(std::move(isotope));
    }
    if(isotopeIds.empty())return Fail("physical-state profile has no isotope inventory",error);
    for(const auto& row:constants) {
      Constant value;value.unit=row[2];value.sourceRecord=row[4];
      if(row[0].empty()||row[2].empty()||!Number(row[1],&value.value,error)||value.value<=0||
          !Source(sources,row[3],row[4],{"constants"},&value.sourceIdentity,nullptr,error)||
          !candidate.constants_.emplace(row[0],value).second)
        return Fail("invalid/duplicate source-bound positive constant: "+row[0],error);
    }
    const auto electron=candidate.constants_.find("electron mass");
    if(electron==candidate.constants_.end()||electron->second.unit!="kg")
      return Fail("explicit source electron mass in kg is required",error);
    const K::MassValue electronMass{electron->second.value,electron->second.sourceIdentity};
    std::set<std::string> stateIds;
    for(const auto& row:states) {
      StateMetadata metadata;metadata.id=row[0];metadata.label=row[1];metadata.kind=row[2];metadata.missingReason=row[15];
      if(!Identifier(metadata.id)||metadata.label.empty()||!stateIds.insert(metadata.id).second||
          !std::set<std::string>{"ground","excited","unspecified","free-electron"}.count(metadata.kind)||
          !Integer(row[4],&metadata.charge,false,error)||
          !Nuclei(row[3],isotopeIds,&metadata.nuclei,error))
        return Fail("invalid/duplicate explicit state identity: "+row[0],error);
      const bool electronState=metadata.kind=="free-electron";
      if((electronState&&(!metadata.nuclei.empty()||metadata.charge!=-1))||
          (!electronState&&metadata.nuclei.empty()))
        return Fail("state nuclei/charge contradict explicit free-electron classification: "+metadata.id,error);
      if(row[14]=="unavailable") {
        if(metadata.missingReason.empty()||std::any_of(row.begin()+5,row.begin()+14,[](const std::string& field){return !field.empty();}))
          return Fail("unavailable state must preserve a reason and no invented energy/mass fields: "+metadata.id,error);
        candidate.metadata_.push_back(std::move(metadata));continue;
      }
      if(row[14]!="available"||!metadata.missingReason.empty()||metadata.kind=="unspecified")
        return Fail("unsupported state status or unresolved available state: "+metadata.id,error);
      K::StateDefinition definition;definition.identity=metadata.id;definition.state=metadata.label;
      definition.nuclei=metadata.nuclei;definition.elementaryCharge=metadata.charge;
      if(!Energy(sources,row[5],row[6],row[7],false,metadata,&definition.formation,error)||
          !Energy(sources,row[8],row[9],row[10],true,metadata,&definition.excitation,error))return false;
      K::ResolvedState resolved;
      if(!K::ResolveState(definition,candidate.isotopes_,electronMass,&resolved,error))return false;
      candidate.resolved_.emplace(metadata.id,std::move(resolved));metadata.available=true;
      if(!row[11].empty()||!row[12].empty()||!row[13].empty()) {
        K::MassValue physical;
        if(!Number(row[11],&physical.kg,error)||physical.kg<=0||
            !Source(sources,row[12],row[13],{"isotope-mass","constants","energy"},&physical.sourceIdentity,nullptr,error))
          return Fail("invalid optional physical rest mass; no reference-mass fallback: "+metadata.id,error);
        candidate.physical_.emplace(metadata.id,std::move(physical));metadata.physicalRestMassAvailable=true;
      }
      candidate.metadata_.push_back(std::move(metadata));
    }
    if(candidate.metadata_.empty())return Fail("physical-state profile has no state declarations",error);
    // Identity is independent of site paths and binds exact profile/source
    // bytes. It is a reproducibility identity, never an approval signature.
    std::string identity=files.schema+"\n"+files.isotopes.sha256+"\n"+files.constants.sha256+"\n"+files.states.sha256+"\n";
    for(const auto& source:sources)identity+=source.first+"\n"+source.second.sha256+"\n"+source.second.role+"\n";
    if(!D::BytesSha256(identity,&candidate.identity_,error))return false;
    candidate.loaded_=true;*this=std::move(candidate);if(error)error->clear();return true;
  }catch(const std::exception& ex){return Fail(std::string("physical-state profile allocation/parser failure: ")+ex.what(),error);}
}

bool Registry::GetState(const std::string& id,K::ResolvedState* output,std::string* error)const {
  if(!output)return Fail("state output storage is required",error);
  if(!loaded_)return Fail("physical-state profile unavailable: not loaded",error);
  const auto found=resolved_.find(id);
  if(found==resolved_.end()) {
    const auto metadata=std::find_if(metadata_.begin(),metadata_.end(),[&](const StateMetadata& row){return row.id==id;});
    return Fail("physical state unavailable: "+id+(metadata==metadata_.end()?" (not declared)":": "+metadata->missingReason),error);
  }
  try{K::ResolvedState candidate=found->second;*output=std::move(candidate);if(error)error->clear();return true;}
  catch(const std::exception& ex){return Fail(std::string("state copy allocation failed: ")+ex.what(),error);}
}
bool Registry::GetPhysicalRestMass(const std::string& id,K::MassValue* output,std::string* error)const {
  if(!output)return Fail("physical-rest-mass output storage is required",error);
  const auto found=physical_.find(id);
  if(!loaded_||found==physical_.end())return Fail("physical rest mass unavailable; additive reference mass is not a substitute: "+id,error);
  try{K::MassValue candidate=found->second;*output=std::move(candidate);if(error)error->clear();return true;}
  catch(const std::exception& ex){return Fail(std::string("physical mass copy allocation failed: ")+ex.what(),error);}
}
bool Registry::GetConstant(const std::string& id,const std::string& expectedUnit,Constant* output,std::string* error)const {
  if(!output)return Fail("constant output storage is required",error);
  const auto found=constants_.find(id);
  if(!loaded_||found==constants_.end()||expectedUnit.empty()||found->second.unit!=expectedUnit)
    return Fail("constant unavailable or explicit native unit mismatch: "+id,error);
  try{Constant candidate=found->second;*output=std::move(candidate);if(error)error->clear();return true;}
  catch(const std::exception& ex){return Fail(std::string("constant copy allocation failed: ")+ex.what(),error);}
}
}}}
