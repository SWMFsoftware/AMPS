#include "SourceCatalog.h"
#include "../src/general/input/amps_input.h"
#include "../src/general/input/amps_input_digest.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <locale>
#include <set>
#include <sstream>

namespace Enceladus { namespace ForwardModel {
// This module resolves data ownership, not plume physics. The configured
// root is supplied by the PIC adapter; all catalog paths below it are local
// immutable assets in IAU_ENCELADUS and SI units. Source-law, chemistry and
// emission validators remain in ForwardModel.cpp. No fitting procedure or
// observational nuisance background is interpreted here.
namespace {
namespace In=::AMPS::Input;
namespace fs=std::filesystem;
using ::AMPS::InputDigest::Sha256Hex;
std::string Unquote(std::string x) {
  if(x.size()>1&&((x.front()=='"'&&x.back()=='"')||(x.front()=='\''&&x.back()=='\'')))return x.substr(1,x.size()-2);
  return x;
}
Status Read(const std::string& path,std::string* out) {
  std::ifstream f(path,std::ios::binary);
  if(!f)return Status::Error("cannot read asset: "+path);
  std::ostringstream bytes;bytes<<f.rdbuf();
  if(f.bad())return Status::Error("asset read failed: "+path);
  *out=bytes.str();return Status::Ok();
}
bool Inside(const fs::path& root,const fs::path& path) {
  auto p=path.begin();
  for(auto r=root.begin();r!=root.end();++r,++p)if(p==path.end()||*p!=*r)return false;
  return true;
}
Status RelativeAsset(const fs::path& root,const std::string& text,fs::path* result) {
  // Canonical component comparison (rather than a string prefix) prevents
  // both ../ traversal and a symlink into a neighbouring data directory.
  // Explicit errors are retained until the application's MPI status gate.
  const fs::path relative=Unquote(text);
  if(relative.empty()||relative.is_absolute())return Status::Error("asset path must be nonempty and relative: "+text);
  for(const auto& part:relative)if(part=="..")return Status::Error("asset path cannot traverse '..': "+text);
  std::error_code error;
  const auto p=fs::canonical(root/relative,error);
  if(error)return Status::Error("missing asset: "+(root/relative).string());
  if(!Inside(root,p))return Status::Error("asset symlink escapes publication root: "+text);
  *result=p;return Status::Ok();
}
bool HasInclude(const std::string& bytes) {
  std::istringstream f(bytes);std::string line;
  while(std::getline(f,line)) {
    const auto p=line.find_first_not_of(" \t\r");
    if(p!=std::string::npos) {
      std::string prefix=line.substr(p,8);
      std::transform(prefix.begin(),prefix.end(),prefix.begin(),[](unsigned char c){return char(std::tolower(c));});
      if(prefix=="#include")return true;
    }
  }
  return false;
}
Status SingleApplication(const RuntimeSelection& selection,In::Document* document,const In::Block** app) {
  auto parsed=In::ParseFile(selection.inputFile,document);
  if(!parsed.ok())return Status::Error(parsed.message);
  auto status=In::RequireNoAssignments(document->root);
  if(!status.ok())return Status::Error(status.message);
  if(document->root.blocks.size()!=1||document->root.blocks.front().kind!="application"||document->root.blocks.front().name!="enceladus")
    return Status::Error("runtime input requires exactly one application: enceladus block");
  *app=&document->root.blocks.front();return Status::Ok();
}
void Set(In::Block* block,const std::string& key,const std::string& value) {
  for(auto& a:block->assignments)if(a.key==key){a.value=value;return;}
  In::Assignment a;a.key=key;a.value=value;a.where=block->begin;block->assignments.push_back(a);
}
void Remove(In::Block* block,const std::string& key) {
  block->assignments.erase(std::remove_if(block->assignments.begin(),block->assignments.end(),[&](const In::Assignment& a){return a.key==key;}),block->assignments.end());
}
Status Require(const In::Block& b,const std::vector<std::string>& keys) {
  for(const auto& k:keys)if(!b.Find(k))return Status::Error(b.begin.Prefix()+"missing catalog key "+k);
  return Status::Ok();
}
Status OwnedRuntimeKeys(const In::Block& block) {
  // Static ownership is checked even in unselected cases, while physical
  // values, species capabilities and external assets are selected-case only.
  const std::set<std::string> keys={"schema_version","source_set","simulated_species","injected_species",
    "chemistry_species","selected_reactions","reaction_products","sample_length_steps","total_iterations",
    "random_seed","particle_splitting_merging","solar_activity","plasma_background_mode","volume_field_output",
    "reproducible_parallel_run",
    "flyby_density_output","flyby_mesh_refinement","flyby_refinement_target_cell_size_m","flyby_refinement_radius_m",
    "flyby_refinement_time_min_s","flyby_refinement_time_max_s","split_position_shift_cell_diagonal_fraction",
    "split_velocity_shift_thermal_speed_fraction","min_particles_per_cell","max_particles_per_cell",
    "nonuniform_split_position_shift","plasma_cache_directory","ices_installation_directory","ices_case",
    "chemistry_weight_reference_number_density_m-3","chemistry_weight_reference_relative_speed_m_s-1",
    "mesh.domain_min_m","mesh.domain_max_m","mesh.background_cell_size_m","mesh.surface_cell_size_north_m",
    "mesh.surface_cell_size_south_m","mesh.south_polar_half_angle_deg","mesh.south_polar_cell_size_m"};
  for(const auto& a:block.assignments) {
    bool valid=keys.count(a.key);
    for(const auto& prefix:{"total_source_rate.","primary_macroparticles_per_step.","characteristic_speed.","courant.","plasma_ion_fraction."})
      if(a.key.rfind(prefix,0)==0)valid=true;
    if(!valid)return Status::Error(a.where.Prefix()+"unowned case-runtime key '"+a.key+"'; sources/observation paths belong to the source catalog");
  }
  for(const auto& child:block.blocks) {
    if(child.kind!="case"&&child.kind!="refinement")return Status::Error(child.begin.Prefix()+"case input accepts only case/refinement blocks; observations belong to the catalog");
    if(child.kind=="case") {
      if(child.Find("schema_version")||child.Find("source_set"))return Status::Error(child.begin.Prefix()+"case cannot override schema_version/source_set");
      for(const auto& nested:child.blocks)if(nested.kind=="case")return Status::Error(nested.begin.Prefix()+"nested cases are forbidden");
      auto s=OwnedRuntimeKeys(child);if(!s.success)return s;
    }
  }
  return Status::Ok();
}
std::string Json(const std::string& x) {
  std::ostringstream o;o<<'"';
  for(unsigned char c:x) {
    if(c=='"'||c=='\\')o<<'\\'<<char(c);
    else if(c=='\n')o<<"\\n";
    else if(c=='\r')o<<"\\r";
    else if(c=='\t')o<<"\\t";
    else if(c<32)o<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;
    else o<<char(c);
  }
  return o.str()+'"';
}
}

Status ResolveRun(const RuntimeSelection& selection,const std::string& configuredRoot,ResolvedRun* output) try {
  // Resolve into a private value. Caller-visible configuration is assigned
  // only after every selected source/observation byte and physical control
  // has passed validation. Unselected cases do not acquire external assets.
  if(!output)return Status::Error("resolved-run output is null");
  ResolvedRun run;In::Document runtime;const In::Block* app=nullptr;
  auto status=SingleApplication(selection,&runtime,&app);if(!status.success)return status;
  In::Block effective;
  auto resolved=In::ResolveCase(*app,selection.caseId,&effective);
  if(!resolved.ok())return Status::Error(resolved.message);
  run.caseId=selection.caseId;
  std::string identity=In::ToCanonicalText(effective);
  // Bind actual input/include bytes too. Logical relative paths are portable
  // across rank-local mounts; absolute diagnostic locations are not physics.
  for(const auto& item:runtime.sourceBytes) {
    const auto logical=fs::path(item.first).lexically_relative(fs::path(selection.inputFile).parent_path()).generic_string();
    identity+="\ninput:"+logical+":"+Sha256Hex(item.second);
    // Retain raw root/include files as well as the resolved canonical deck.
    // This binds provenance to what was parsed, including comments and CRLF;
    // native readers never need to reopen these mutable upstream paths.
    run.assetBytes[item.first]=item.second;
    run.assetHashes[item.first]=Sha256Hex(item.second);
    run.assetLogicalPaths[item.first]="runtime_input/"+logical;
  }
  auto capture=[&](const std::string& path,const std::string& logical,const std::string& expected) {
    std::string bytes;auto s=Read(path,&bytes);if(!s.success)return s;
    const auto digest=Sha256Hex(bytes);
    if(!expected.empty()&&digest!=Unquote(expected))return Status::Error("asset SHA-256 mismatch: "+logical);
    run.assetBytes[path]=bytes;run.assetHashes[path]=digest;run.assetLogicalPaths[path]=logical;
    return Status::Ok();
  };
  if(selection.caseId.empty()) {
    if(!effective.Find("schema_version")||Unquote(effective.Find("schema_version")->value)!="2")
      return Status::Error("case-free legacy input requires schema_version=2");
    status=ParseConfiguration(selection.inputFile,&run.configuration);if(!status.success)return status;
    status=ParseSourceDistribution(run.configuration.sourceDistributionFile,&run.distribution);if(!status.success)return status;
    status=capture(run.configuration.sourceDistributionFile,"source_distribution","");if(!status.success)return status;
    run.sourceBytes=run.assetBytes.at(run.configuration.sourceDistributionFile);
    run.canonicalInput=In::ToCanonicalText(effective);
  } else {
    status=ParseCaseIdentity(selection.caseId,&run.measurement,&run.basis);if(!status.success)return status;
    status=OwnedRuntimeKeys(*app);if(!status.success)return status;
    if(!app->Find("schema_version")||Unquote(app->Find("schema_version")->value)!="3")return Status::Error("case mode requires application schema_version=3");
    if(!app->Find("source_set"))return Status::Error("case mode requires source_set");
    run.sourceSet=Unquote(app->Find("source_set")->value);
    std::error_code error;
    const auto root=fs::canonical(configuredRoot,error);
    if(configuredRoot.empty()||error||!fs::is_directory(root))return Status::Error("missing configured model-data root; use Config.pl -model-data-path=PATH");
    fs::path publication;status=RelativeAsset(root,run.sourceSet,&publication);if(!status.success)return status;
    run.dataRoot=root.string();
    if(!fs::is_directory(publication))return Status::Error("source_set must name a directory");
    fs::path catalogFile;status=RelativeAsset(publication,"enceladus_cases.in",&catalogFile);if(!status.success)return status;
    const auto catalogPath=catalogFile.string();
    std::string catalogBytes;status=Read(catalogPath,&catalogBytes);if(!status.success)return status;
    if(HasInclude(catalogBytes))return Status::Error("published case catalog must be self-contained without #include");
    In::Document catalog;
    auto p=In::ParseText(catalogBytes,catalogPath,publication.string(),&catalog);
    if(!p.ok())return Status::Error(p.message);
    if(!catalog.root.assignments.empty()||catalog.root.blocks.size()!=1||catalog.root.blocks[0].kind!="case_catalog"||catalog.root.blocks[0].name!="enceladus")
      return Status::Error("catalog requires exactly case_catalog: enceladus");
    const auto& c=catalog.root.blocks[0];
    status=Require(c,{"schema_version","publication_id","frame","reference_radius_m"});if(!status.success)return status;
    const std::set<std::string> rootKeys={"schema_version","publication_id","frame","reference_radius_m"};
    for(const auto& a:c.assignments)if(!rootKeys.count(a.key))return Status::Error(a.where.Prefix()+"unknown catalog key "+a.key);
    if(Unquote(c.Find("schema_version")->value)!="1"||Unquote(c.Find("frame")->value)!="IAU_ENCELADUS")return Status::Error("catalog schema/frame mismatch");
    const In::Block* selected=nullptr;
    for(const auto& b:c.blocks) {
      if(b.kind!="case")return Status::Error(b.begin.Prefix()+"catalog accepts only case children");
      status=ParseCaseIdentity(b.name,nullptr,nullptr);if(!status.success)return status;
      if(b.name==selection.caseId)selected=&b;
    }
    if(!selected)return Status::Error("no catalog case published for "+selection.caseId);
    In::Block chosen=*selected;
    status=Require(chosen,{"measurement_id","instrument","source_basis","source_distribution_file","source_sha256","available_species","provenance_file","provenance_sha256"});
    if(!status.success)return status;
    const std::set<std::string> caseKeys={"measurement_id","instrument","source_basis","source_distribution_file","source_sha256","available_species","provenance_file","provenance_sha256","flyby_trajectory_file","flyby_trajectory_sha256"};
    for(const auto& a:chosen.assignments)if(!caseKeys.count(a.key)&&a.key.rfind("unavailable_source.",0)!=0)
      return Status::Error(a.where.Prefix()+"unknown catalog case key "+a.key);
    if(Unquote(chosen.Find("measurement_id")->value)!=run.measurement||Unquote(chosen.Find("source_basis")->value)!=run.basis)
      return Status::Error("catalog measurement/basis does not match case ID");
    const auto instrument=Unquote(chosen.Find("instrument")->value);
    if(instrument!="INMS"&&instrument!="UVIS")return Status::Error("catalog instrument must be INMS or UVIS");
    run.publication=Unquote(c.Find("publication_id")->value);
    run.assetBytes[catalogPath]=catalogBytes;run.assetHashes[catalogPath]=Sha256Hex(catalogBytes);run.assetLogicalPaths[catalogPath]="enceladus_cases.in";
    auto asset=[&](In::Block* owner,const std::string& key,const std::string& hashKey) {
      auto s=Require(*owner,{key,hashKey});if(!s.success)return s;
      const auto logical=Unquote(owner->Find(key)->value);fs::path path;
      s=RelativeAsset(publication,logical,&path);if(!s.success)return s;
      if(!fs::is_regular_file(path))return Status::Error("asset must be a regular file: "+logical);
      s=capture(path.string(),logical,owner->Find(hashKey)->value);if(!s.success)return s;
      Set(owner,key,path.string());Remove(owner,hashKey);return Status::Ok();
    };
    // The fingerprint binds logical catalog data before paths are absolutized.
    identity+="\ncatalog:"+In::ToCanonicalText(c);
    status=asset(&chosen,"source_distribution_file","source_sha256");if(!status.success)return status;
    const auto source=Unquote(chosen.Find("source_distribution_file")->value);
    run.sourceBytes=run.assetBytes.at(source);
    if(HasInclude(run.sourceBytes))return Status::Error("published source distribution must not use #include");
    status=ParseSourceDistributionText(run.sourceBytes,source,&run.distribution);if(!status.success)return status;
    if(run.distribution.sourceCase!=run.basis||run.distribution.frame!=Unquote(c.Find("frame")->value))return Status::Error("source basis/frame disagrees with catalog");
    double reference=0;
    try{std::size_t end=0;reference=std::stod(c.Find("reference_radius_m")->value,&end);if(end!=c.Find("reference_radius_m")->value.size())return Status::Error("invalid catalog radius");}
    catch(...){return Status::Error("invalid catalog radius");}
    if(reference!=run.distribution.stripeReferenceRadiusM)return Status::Error("catalog/source reference radius mismatch");
    // available_species is redundant intentionally: it diagnoses malformed
    // publication metadata rather than silently ignoring extra strengths.
    std::set<std::string> species;
    std::istringstream symbols(Unquote(chosen.Find("available_species")->value));std::string symbol;
    while(std::getline(symbols,symbol,',')) {
      const auto first=symbol.find_first_not_of(" \t"),last=symbol.find_last_not_of(" \t");
      if(first==std::string::npos||!species.insert(symbol.substr(first,last-first+1)).second)return Status::Error("invalid/duplicate available_species");
    }
    std::set<std::string> actual;
    for(const auto& x:run.distribution.sources.front().strength)actual.insert(x.first);
    if(actual!=species)return Status::Error("catalog available_species disagrees with source strengths");
    status=asset(&chosen,"provenance_file","provenance_sha256");if(!status.success)return status;
    Set(&effective,"source_distribution_file",source);Remove(&effective,"source_set");
    if(effective.Find("flyby_density_output")&&Unquote(effective.Find("flyby_density_output")->value)=="true") {
      status=asset(&chosen,"flyby_trajectory_file","flyby_trajectory_sha256");if(!status.success)return status;
      Set(&effective,"flyby_id",run.measurement);Set(&effective,"flyby_trajectory_file",chosen.Find("flyby_trajectory_file")->value);
    }
    for(auto observation:chosen.blocks) {
      if(observation.kind!="observation")return Status::Error(observation.begin.Prefix()+"catalog case accepts only observation blocks");
      if(!observation.Find("event_id")||Unquote(observation.Find("event_id")->value)!=run.measurement||
         !observation.Find("instrument")||Unquote(observation.Find("instrument")->value)!=instrument)
        return Status::Error("catalog observation measurement/instrument mismatch");
      status=asset(&observation,"data_file","data_sha256");if(!status.success)return status;
      if(observation.Find("observer_geometry_file")) {
        status=asset(&observation,"observer_geometry_file","observer_geometry_sha256");if(!status.success)return status;
      }
      // Unsupported species requests are errors later; unrelated comparison
      // products of a valid multispecies bundle are omitted deliberately.
      if(observation.Find("comparison_species")) {
        const auto s=Unquote(observation.Find("comparison_species")->value);
        const auto injected=effective.Find("simulated_species");
        if(injected) {
          std::istringstream list(Unquote(injected->value));std::string item;bool active=false;
          while(std::getline(list,item,',')) {item=In::detail::Trim(item);if(item==s)active=true;}
          if(s!="none"&&!active)continue;
        }
      }
      effective.blocks.push_back(observation);
    }
    if(chosen.blocks.empty())return Status::Error("published measurement case has no observations");
    status=ParseResolvedConfiguration(effective,selection.inputFile,&run.configuration);if(!status.success)return status;
    run.canonicalInput=In::ToCanonicalText(effective);
  }
  // Legacy paths are already resolved by their existing owner. Snapshot all
  // selected tables so native observation readers use precisely these bytes.
  auto legacyAsset=[&](const std::string& path,const std::string& logical) {
    return path.empty()||run.assetBytes.count(path)?Status::Ok():capture(path,logical,"");
  };
  status=legacyAsset(run.configuration.flybyTrajectoryFile,"trajectory");if(!status.success)return status;
  for(const auto& o:run.configuration.observations) {
    status=legacyAsset(o.dataFile,"observation/"+o.id);if(!status.success)return status;
    status=legacyAsset(o.observerGeometryFile,"observer/"+o.id);if(!status.success)return status;
  }
  status=ValidateConfiguration(run.configuration,run.distribution);if(!status.success)return status;
  std::map<std::string,std::string> logicalHashes;
  for(const auto& h:run.assetHashes)logicalHashes[run.assetLogicalPaths.at(h.first)]=h.second;
  for(const auto& h:logicalHashes)identity+="\nasset:"+h.first+":"+h.second;
  run.fingerprint=Sha256Hex(identity);*output=std::move(run);return Status::Ok();
} catch(const fs::filesystem_error& error) {
  // Filesystem accessibility may differ on MPI ranks. Return a typed input
  // failure so every rank still reaches the caller's collective reduction;
  // an uncaught exception here could strand the other ranks in MPI.
  return Status::Error(std::string("case filesystem access failed: ")+error.what());
}

Status ListRuntimeCases(const RuntimeSelection& selection,const std::string& root,std::vector<CaseInventory>* output) {
  In::Document input;const In::Block* app=nullptr;auto s=SingleApplication(selection,&input,&app);if(!s.success)return s;
  output->clear();
  for(const auto& c:app->blocks)if(c.kind=="case") {
    RuntimeSelection selected=selection;selected.caseId=c.name;ResolvedRun run;
    s=ResolveRun(selected,root,&run);
    const auto state=s.success?"READY":s.message.find("no catalog case")!=std::string::npos?"MISSING_PUBLICATION":
      s.message.find("missing")!=std::string::npos||s.message.find("cannot read")!=std::string::npos?"MISSING_DATA":
      s.message.find("no strengths")!=std::string::npos?"DEFERRED_SOURCE":"INVALID";
    output->push_back({c.name,state,s.success?DescribeRun(run):s.message});
  }
  if(output->empty())return Status::Error("input contains no named cases");
  return Status::Ok();
}
std::string DescribeRun(const ResolvedRun& run) {
  std::ostringstream out;out.imbue(std::locale::classic());out<<std::setprecision(17);
  out<<"case="<<run.caseId<<" publication="<<run.publication<<" source_basis="<<run.distribution.sourceCase
     <<" source_mode="<<Name(run.distribution.mode)<<" source_count="<<run.distribution.sources.size()
     <<" observations="<<run.configuration.observations.size()<<" fingerprint="<<run.fingerprint<<"\n";
  for(const auto& symbol:run.configuration.injectedSpecies) {
    double q=0;for(const auto& source:run.distribution.sources) {const auto a=source.strength.find(symbol);if(a!=source.strength.end())q+=a->second;}
    if(run.distribution.mode==SourceMode::Shape)q=run.configuration.totalSourceRatePerS.at(symbol);
    out<<"species="<<symbol<<" Q_s-1="<<q<<" particles_per_step="<<run.configuration.primaryMacroparticlesPerStep.at(symbol)<<"\n";
  }
  if(run.configuration.mesh.explicitConfiguration) {
    const auto& m=run.configuration.mesh;
    out<<"domain_min_m="<<m.minimumM.x<<","<<m.minimumM.y<<","<<m.minimumM.z
       <<" domain_max_m="<<m.maximumM.x<<","<<m.maximumM.y<<","<<m.maximumM.z
       <<" background_cell_size_m="<<m.backgroundCellSizeM<<"\n";
  }
  out<<"global_dt_s=pending_realized_AMR statistical_weights=pending_realized_AMR\n";
  return out.str();
}
Status WriteRunSnapshot(const ResolvedRun& run,const std::string& directory,const NumericalPlan* plan) {
  // Called by rank zero only. Assets are the already verified in-memory
  // bytes, so this receipt cannot silently capture a later file revision.
  // Final dt/W are absent until the actual AMR minimum edge is known;
  // an atomic manifest rename exposes each complete phase to report readers.
  const fs::path root=directory,inputs=root/"case_inputs";std::error_code ec;
  fs::create_directories(inputs,ec);if(ec)return Status::Error("cannot create case snapshot directory: "+ec.message());
  auto write=[&](const fs::path& p,const std::string& bytes) {
    std::ofstream f(p,std::ios::binary);f<<bytes;return bool(f);
  };
  if(!write(inputs/"resolved.in",run.canonicalInput)||!write(inputs/"sources.in",run.sourceBytes))return Status::Error("cannot write source/input snapshots");
  std::ostringstream o;o.imbue(std::locale::classic());o<<std::setprecision(17);
  o<<"{\n\"schema_version\":1,\"case_id\":"<<Json(run.caseId)<<",\"publication\":"<<Json(run.publication)
   <<",\"measurement_id\":"<<Json(run.measurement)<<",\"source_basis\":"<<Json(run.distribution.sourceCase)
   <<",\"source_mode\":"<<Json(Name(run.distribution.mode))<<",\"frame\":"<<Json(run.distribution.frame)
   <<",\"reference_radius_m\":"<<run.distribution.stripeReferenceRadiusM
   <<",\"fingerprint\":"<<Json(run.fingerprint)<<",\"data_root\":"<<Json(run.dataRoot)
   <<",\"source_set\":"<<Json(run.sourceSet)<<",\"phase\":"<<Json(run.completedIterations?"complete":plan?"mesh_ready":"input_ready")
   <<",\"completed_iterations\":"<<run.completedIterations<<",\"mpi_ranks\":"<<run.mpiRanks
   <<",\"executable_sha256\":"<<Json(run.executableSha256)<<",\"assets\":[";
  std::size_t index=0;
  for(const auto& asset:run.assetBytes) {
    const auto name="asset_"+std::to_string(index)+".dat";
    if(!write(inputs/name,asset.second))return Status::Error("cannot snapshot "+asset.first);
    if(index++)o<<",";
    o<<"{\"path\":"<<Json(asset.first)<<",\"logical_path\":"<<Json(run.assetLogicalPaths.at(asset.first))
     <<",\"sha256\":"<<Json(run.assetHashes.at(asset.first))<<",\"snapshot\":"<<Json("case_inputs/"+name)<<"}";
  }
  o<<"]";
  if(plan) {
    o<<",\"global_dt_s\":"<<plan->globalTimeStepS<<",\"minimum_cell_size_m\":"<<plan->minimumAllocatedCellSizeM
     <<",\"limiting_species\":"<<Json(plan->limitingSpecies)<<",\"species\":{";bool first=true;
    for(const auto& q:plan->physicalSourceRatePerS) {
      if(!first)o<<",";first=false;
      o<<Json(q.first)<<":{\"Q_s-1\":"<<q.second<<",\"weight\":"<<plan->statisticalWeight.at(q.first)
       <<",\"particles_per_step\":"<<plan->expectedMacroparticlesPerStep.at(q.first)<<"}";
    }
    o<<"}";
  }
  o<<"}\n";
  if(!write(root/"run_manifest.json.tmp",o.str()))return Status::Error("cannot write run manifest");
  fs::rename(root/"run_manifest.json.tmp",root/"run_manifest.json",ec);
  return ec?Status::Error("cannot publish run manifest: "+ec.message()):Status::Ok();
}
} }
