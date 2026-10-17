#include "DamageResponseLibrary.hh"

#include "ContentHash.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace shield_sim {
namespace {

namespace fs=std::filesystem;

std::string ReadBytes(const std::string& path){
  std::ifstream input(path,std::ios::binary);
  if(!input) throw JsonError("cannot open damage-response table '"+path+"'");
  std::ostringstream bytes;
  bytes<<input.rdbuf();
  if(!input.good() && !input.eof())
    throw JsonError("failed while reading damage-response table '"+path+"'");
  return bytes.str();
}

std::string RequiredString(const Json::Object& object,const std::string& key,
                           const std::string& context){
  const auto found=object.find(key);
  if(found==object.end()) throw JsonError(context+" requires '"+key+"'");
  const std::string value=found->second.AsString(context+"."+key);
  if(value.empty()) throw JsonError(context+"."+key+" must not be empty");
  return value;
}

double RequiredPositive(const Json::Object& object,const std::string& key,
                        const std::string& context){
  const auto found=object.find(key);
  if(found==object.end()) throw JsonError(context+" requires '"+key+"'");
  const double value=found->second.AsNumber(context+"."+key);
  if(!std::isfinite(value) || value<=0.0)
    throw JsonError(context+"."+key+" must be finite and > 0");
  return value;
}

double ParseNumber(const std::string& token,const std::string& context){
  std::size_t used=0;
  double result=0.0;
  try { result=std::stod(token,&used); }
  catch(...) { throw JsonError(context+" is not a number"); }
  if(used!=token.size() || !std::isfinite(result))
    throw JsonError(context+" must be one finite number");
  return result;
}

std::vector<std::string> SplitCsv(const std::string& line){
  std::vector<std::string> fields;
  std::size_t begin=0;
  for(;;){
    const std::size_t comma=line.find(',',begin);
    fields.push_back(line.substr(begin,comma==std::string::npos
        ? std::string::npos : comma-begin));
    if(comma==std::string::npos) break;
    begin=comma+1;
  }
  return fields;
}

std::vector<std::vector<std::string>> ParseCsv(const std::string& bytes,
                                               const std::string& path){
  std::istringstream input(bytes);
  std::vector<std::vector<std::string>> rows;
  std::string line;
  std::size_t lineNumber=0;
  while(std::getline(input,line)){
    ++lineNumber;
    if(!line.empty() && line.back()=='\r') line.pop_back();
    if(line.empty()) throw JsonError(path+":"+std::to_string(lineNumber)+
                                    ": blank rows are not permitted");
    rows.push_back(SplitCsv(line));
  }
  if(rows.size()<2) throw JsonError(path+" must contain a header and data rows");
  return rows;
}

double IntegrateLogLogSegment(double e0,double e1,double s0,double s1,
                              double lower,double upper){
  // S(E)=S0*(E/E0)^m is exactly the interpolation contract. Integrating that
  // power law avoids assigning an arbitrary bin-center NIEL to a wide source
  // bin. The logarithmic branch is the analytic m=-1 limit.
  const double m=std::log(s1/s0)/std::log(e1/e0);
  if(std::abs(m+1.0)<1.0e-12) return s0*e0*std::log(upper/lower);
  return s0*e0*(std::pow(upper/e0,m+1.0)-std::pow(lower/e0,m+1.0))/(m+1.0);
}

} // namespace

std::string DamageResponseStatusName(DamageResponseStatus status){
  switch(status){
    case DamageResponseStatus::Available: return "available";
    case DamageResponseStatus::UnsupportedSpecies: return "unsupported_species";
    case DamageResponseStatus::OutOfDomain: return "out_of_domain";
  }
  return "out_of_domain";
}

DamageResponseLibrary DamageResponseLibrary::Load(const std::string& manifestPath){
  DamageResponseLibrary library;
  library.manifestPath_=fs::absolute(manifestPath).lexically_normal().string();
  library.canonicalManifest_=Json::ParseFile(library.manifestPath_);
  const auto& root=library.canonicalManifest_.AsObject("damage_response_library");
  RejectUnknownFields(root,{"schema_version","library_id","library_version","target",
      "out_of_domain_policy","charged_particle_interpolation","neutron_integration",
      "n_eq","tables"},"damage_response_library");
  if(RequiredString(root,"schema_version","damage_response_library")!=
     "shieldSim.damage_response_library.v1")
    throw JsonError("damage_response_library.schema_version is unsupported");
  library.libraryId_=RequiredString(root,"library_id","damage_response_library");
  library.libraryVersion_=RequiredString(root,"library_version","damage_response_library");
  if(RequiredString(root,"out_of_domain_policy","damage_response_library")!=
     "not_assessable_no_extrapolation")
    throw JsonError("damage_response_library requires not_assessable_no_extrapolation");
  if(RequiredString(root,"charged_particle_interpolation","damage_response_library")!=
     "log_log_between_positive_nodes")
    throw JsonError("damage_response_library charged-particle interpolation must be log-log");
  if(RequiredString(root,"neutron_integration","damage_response_library")!=
     "ASTM_groupwise_explicit_edges")
    throw JsonError("damage_response_library neutron integration must be ASTM groupwise");

  const auto& target=root.at("target").AsObject("damage_response_library.target");
  RejectUnknownFields(target,{"material_id","description","atomic_weight_g_mol"},
                      "damage_response_library.target");
  library.targetMaterialId_=RequiredString(target,"material_id","damage_response_library.target");
  (void)RequiredString(target,"description","damage_response_library.target");
  const double atomicWeight=RequiredPositive(target,"atomic_weight_g_mol",
      "damage_response_library.target");

  const auto& neq=root.at("n_eq").AsObject("damage_response_library.n_eq");
  RejectUnknownFields(neq,{"convention","reference_damage_function_MeV_mb",
      "reference_niel_MeV_cm2_g","avogadro_constant_mol_inverse","supported_species"},
      "damage_response_library.n_eq");
  if(RequiredString(neq,"convention","damage_response_library.n_eq")!=
     "ASTM_E722_19_silicon_neutron_contribution_only")
    throw JsonError("damage_response_library n_eq convention is unsupported");
  if(std::abs(RequiredPositive(neq,"reference_damage_function_MeV_mb",
      "damage_response_library.n_eq")-95.0)>1.0e-12)
    throw JsonError("damage_response_library E722 reference must be 95 MeV-mb");
  library.referenceNeutronNiel_=RequiredPositive(neq,"reference_niel_MeV_cm2_g",
      "damage_response_library.n_eq");
  const double avogadro=RequiredPositive(neq,"avogadro_constant_mol_inverse",
      "damage_response_library.n_eq");
  // E722's 95 MeV-mb reference is a microscopic damage function. Convert it
  // independently to mass NIEL: 95 MeV * 1e-27 cm2 * N_A / M_Si. Reject a
  // manifest that states an inconsistent denominator, since such a mismatch
  // would scale every reported neutron-equivalent fluence.
  const double derivedReference=95.0e-27*avogadro/atomicWeight;
  if(std::abs(library.referenceNeutronNiel_-derivedReference)>
     1.0e-12*derivedReference)
    throw JsonError("damage_response_library n_eq reference_niel_MeV_cm2_g is "
                    "inconsistent with 95 MeV-mb, N_A, and the target atomic weight");
  const auto& neqSpecies=neq.at("supported_species").AsArray(
      "damage_response_library.n_eq.supported_species");
  if(neqSpecies.size()!=1 || neqSpecies.front().AsString()!="neutron")
    throw JsonError("damage_response_library n_eq supports neutron contributions only");

  const fs::path directory=fs::path(library.manifestPath_).parent_path();
  const auto& tables=root.at("tables").AsArray("damage_response_library.tables");
  if(tables.empty()) throw JsonError("damage_response_library.tables must not be empty");
  for(std::size_t index=0;index<tables.size();++index){
    const std::string context="damage_response_library.tables["+std::to_string(index)+"]";
    const auto& object=tables[index].AsObject(context);
    RejectUnknownFields(object,{"species","representation","energy_coordinate","energy_unit",
        "response_unit","path","sha256","domain_min_MeV","domain_max_MeV","source"},context);
    Table table;
    table.species=RequiredString(object,"species",context);
    table.representation=RequiredString(object,"representation",context);
    if(table.species!="proton" && table.species!="alpha" && table.species!="neutron")
      throw JsonError(context+".species is outside the Stage-5 support matrix");
    const std::string requiredRepresentation=table.species=="neutron"
        ? "group_constant" : "positive_point_nodes";
    if(table.representation!=requiredRepresentation)
      throw JsonError(context+" has an inconsistent species/representation pair");
    if(RequiredString(object,"energy_coordinate",context)!="total_kinetic_energy" ||
       RequiredString(object,"energy_unit",context)!="MeV" ||
       RequiredString(object,"response_unit",context)!="MeV cm2/g")
      throw JsonError(context+" requires total kinetic energy [MeV] and NIEL [MeV cm2/g]");
    table.domainMinMeV=RequiredPositive(object,"domain_min_MeV",context);
    table.domainMaxMeV=RequiredPositive(object,"domain_max_MeV",context);
    if(table.domainMaxMeV<=table.domainMinMeV)
      throw JsonError(context+" has an invalid domain");
    table.contentHash=RequiredString(object,"sha256",context);
    table.path=(directory/RequiredString(object,"path",context)).lexically_normal().string();
    const auto& source=object.at("source").AsObject(context+".source");
    std::vector<std::string> sourceFields={"name","version","retrieved_utc","url","license_status",
        "raw_response_sha256","parameters"};
    if(table.species=="neutron"){
      sourceFields.push_back("edge_source"); sourceFields.push_back("edge_source_url");
      sourceFields.push_back("edge_source_sha256");
    }
    RejectUnknownFields(source,sourceFields,context+".source");
    for(const char* field:{"name","version","retrieved_utc","url","license_status",
                           "raw_response_sha256","parameters"})
      (void)RequiredString(source,field,context+".source");
    if(table.species=="neutron")
      for(const char* field:{"edge_source","edge_source_url","edge_source_sha256"})
        (void)RequiredString(source,field,context+".source");
    const std::string bytes=ReadBytes(table.path);
    if(Sha256(bytes)!=table.contentHash)
      throw JsonError(context+" SHA-256 does not match table '"+table.path+"'");
    const auto rows=ParseCsv(bytes,table.path);
    if(table.species=="neutron"){
      const std::vector<std::string> header={"lower_MeV","upper_MeV","representative_MeV",
          "niel_MeV_cm2_g","hardness_factor_to_1MeV_neutron"};
      if(rows.front()!=header) throw JsonError(table.path+" has an unexpected neutron header");
      for(std::size_t row=1;row<rows.size();++row){
        if(rows[row].size()!=5) throw JsonError(table.path+":"+std::to_string(row+1)+
                                               ": expected five fields");
        const std::string line=table.path+":"+std::to_string(row+1);
        const double lower=ParseNumber(rows[row][0],line);
        const double upper=ParseNumber(rows[row][1],line);
        const double representative=ParseNumber(rows[row][2],line);
        const double niel=ParseNumber(rows[row][3],line);
        const double hardness=ParseNumber(rows[row][4],line);
        if(lower<=0.0 || upper<=lower || representative<=lower || representative>=upper ||
           niel<=0.0 || hardness<=0.0)
          throw JsonError(line+" contains an invalid E722 group");
        if(!table.upperMeV.empty() && lower!=table.upperMeV.back())
          throw JsonError(line+" does not continue the preceding E722 group exactly");
        if(std::abs(hardness-niel/library.referenceNeutronNiel_)>
           5.0e-12*std::max(1.0,std::abs(hardness)))
          throw JsonError(line+" hardness factor is inconsistent with the 95 MeV-mb reference");
        table.lowerMeV.push_back(lower); table.upperMeV.push_back(upper);
        table.energyMeV.push_back(representative);
        table.nielMeVCm2PerG.push_back(niel); table.neutronHardness.push_back(hardness);
      }
      if(table.lowerMeV.front()!=table.domainMinMeV ||
         table.upperMeV.back()!=table.domainMaxMeV)
        throw JsonError(context+" domain does not equal its explicit neutron group edges");
    } else {
      const std::vector<std::string> header={"energy_MeV","niel_MeV_cm2_g"};
      if(rows.front()!=header) throw JsonError(table.path+" has an unexpected charged-particle header");
      for(std::size_t row=1;row<rows.size();++row){
        if(rows[row].size()!=2) throw JsonError(table.path+":"+std::to_string(row+1)+
                                               ": expected two fields");
        const std::string line=table.path+":"+std::to_string(row+1);
        const double energy=ParseNumber(rows[row][0],line);
        const double niel=ParseNumber(rows[row][1],line);
        if(energy<=0.0 || niel<=0.0 ||
           (!table.energyMeV.empty() && energy<=table.energyMeV.back()))
          throw JsonError(line+" violates positive, strictly increasing log-log nodes");
        table.energyMeV.push_back(energy); table.nielMeVCm2PerG.push_back(niel);
      }
      if(table.energyMeV.size()<2 || table.energyMeV.front()!=table.domainMinMeV ||
         table.energyMeV.back()!=table.domainMaxMeV)
        throw JsonError(context+" domain does not equal its charged-particle endpoint nodes");
    }
    if(!library.tables_.emplace(table.species,std::move(table)).second)
      throw JsonError("damage_response_library contains duplicate species tables");
  }
  const std::set<std::string> expected={"proton","alpha","neutron"};
  std::set<std::string> actual;
  for(const auto& item:library.tables_) actual.insert(item.first);
  if(actual!=expected)
    throw JsonError("damage_response_library must contain exactly proton, alpha, and neutron tables");
  library.manifestHash_=Sha256(library.canonicalManifest_.DumpCanonical());
  return library;
}

const DamageResponseLibrary::Table* DamageResponseLibrary::FindTable(
    const std::string& species) const{
  const auto found=tables_.find(species);
  return found==tables_.end()?nullptr:&found->second;
}

bool DamageResponseLibrary::SupportsSpecies(const std::string& species) const{
  return FindTable(species)!=nullptr;
}

DamageResponseValue DamageResponseLibrary::Evaluate(const std::string& species,
                                                     double energyMeV) const{
  DamageResponseValue result;
  const Table* table=FindTable(species);
  if(!table){
    result.status=DamageResponseStatus::UnsupportedSpecies;
    result.reason="no approved natural-silicon NIEL table for species '"+species+"'";
    return result;
  }
  if(!std::isfinite(energyMeV) || energyMeV<table->domainMinMeV ||
     energyMeV>table->domainMaxMeV){
    result.status=DamageResponseStatus::OutOfDomain;
    result.reason="energy is outside the pinned "+species+" NIEL domain; extrapolation is forbidden";
    return result;
  }
  result.status=DamageResponseStatus::Available;
  if(species=="neutron"){
    auto upper=std::upper_bound(table->lowerMeV.begin(),table->lowerMeV.end(),energyMeV);
    std::size_t index=upper==table->lowerMeV.begin()?0:
        static_cast<std::size_t>(upper-table->lowerMeV.begin()-1);
    if(index>=table->upperMeV.size()) index=table->upperMeV.size()-1;
    if(energyMeV>table->upperMeV[index] ||
       (energyMeV==table->upperMeV[index] && index+1<table->upperMeV.size())){
      result.status=DamageResponseStatus::OutOfDomain;
      result.reason="energy falls in a gap between pinned E722 groups";
      return result;
    }
    result.intervalIndex=index;
    result.nielMeVCm2PerG=table->nielMeVCm2PerG[index];
    result.neutronHardnessFactor=table->neutronHardness[index];
    return result;
  }
  auto upper=std::lower_bound(table->energyMeV.begin(),table->energyMeV.end(),energyMeV);
  if(upper!=table->energyMeV.end() && *upper==energyMeV){
    result.intervalIndex=static_cast<std::size_t>(upper-table->energyMeV.begin());
    result.nielMeVCm2PerG=table->nielMeVCm2PerG[result.intervalIndex];
    return result;
  }
  const std::size_t high=static_cast<std::size_t>(upper-table->energyMeV.begin());
  const std::size_t low=high-1;
  const double fraction=std::log(energyMeV/table->energyMeV[low])/
      std::log(table->energyMeV[high]/table->energyMeV[low]);
  result.intervalIndex=low;
  result.nielMeVCm2PerG=std::exp(std::log(table->nielMeVCm2PerG[low])+fraction*
      std::log(table->nielMeVCm2PerG[high]/table->nielMeVCm2PerG[low]));
  return result;
}

DamageFoldResult DamageResponseLibrary::FoldPiecewiseConstantDifferential(
    const std::string& species,const std::vector<DifferentialFluenceBin>& bins) const{
  DamageFoldResult result;
  const Table* table=FindTable(species);
  if(!table){
    result.status=DamageResponseStatus::UnsupportedSpecies;
    result.reason="no approved natural-silicon NIEL table for species '"+species+"'";
    return result;
  }
  for(std::size_t binIndex=0;binIndex<bins.size();++binIndex){
    const auto& bin=bins[binIndex];
    if(!std::isfinite(bin.lowerMeV) || !std::isfinite(bin.upperMeV) ||
       !std::isfinite(bin.differentialFluencePerCm2PerMeV) ||
       bin.lowerMeV<=0.0 || bin.upperMeV<=bin.lowerMeV ||
       bin.differentialFluencePerCm2PerMeV<0.0)
      throw std::invalid_argument("invalid explicit differential-fluence bin at index "+
                                  std::to_string(binIndex));
    if(bin.differentialFluencePerCm2PerMeV==0.0) continue;
    result.fluencePerCm2+=bin.differentialFluencePerCm2PerMeV*(bin.upperMeV-bin.lowerMeV);
    if(bin.lowerMeV<table->domainMinMeV || bin.upperMeV>table->domainMaxMeV){
      result.status=DamageResponseStatus::OutOfDomain;
      result.reason="nonzero source bin extends outside the pinned "+species+
          " NIEL domain; the bin was not clipped or extrapolated";
      return result;
    }
    double responseIntegral=0.0;
    if(species=="neutron"){
      for(std::size_t group=0;group<table->lowerMeV.size();++group){
        const double lower=std::max(bin.lowerMeV,table->lowerMeV[group]);
        const double upper=std::min(bin.upperMeV,table->upperMeV[group]);
        if(upper>lower)
          responseIntegral+=(upper-lower)*table->nielMeVCm2PerG[group];
      }
    } else {
      for(std::size_t node=0;node+1<table->energyMeV.size();++node){
        const double lower=std::max(bin.lowerMeV,table->energyMeV[node]);
        const double upper=std::min(bin.upperMeV,table->energyMeV[node+1]);
        if(upper>lower)
          responseIntegral+=IntegrateLogLogSegment(table->energyMeV[node],
              table->energyMeV[node+1],table->nielMeVCm2PerG[node],
              table->nielMeVCm2PerG[node+1],lower,upper);
      }
    }
    result.dddMeVPerG+=bin.differentialFluencePerCm2PerMeV*responseIntegral;
  }
  if(species=="neutron")
    result.neutronEquivalentFluencePerCm2=result.dddMeVPerG/referenceNeutronNiel_;
  result.reason="all nonzero fluence lies inside the pinned response domain";
  return result;
}

} // namespace shield_sim
