// Explicit conditional physical Na cohort. Normal amps_init/SPICE/mesh and
// actual production services are used. Current missing Na+ state data =>77
// BEFORE expensive initialization, never a synthetic state or zero rate.
#include "MoonSodiumCohort.h"
#include "MoonInput.h"
#include "MoonTransportClock.h"
#include "Moon.h"
#include "MoonSurfaceEnvironment.h"
#include "zero_force_recipe.h"
#include "collective_admission.h"
#include "src/general/input/amps_input.h"
#include <filesystem>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>
void amps_init();
#define MOON_COHORT_STRING_INNER(...) #__VA_ARGS__
#define MOON_COHORT_STRING(...) MOON_COHORT_STRING_INNER(__VA_ARGS__)
const char* NativeAccelerationMacro=MOON_COHORT_STRING(_PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(a,s,p,x,v,n));
namespace S=Moon::Chemistry::SodiumCohort;
namespace E=Moon::Species::PhysicalReactionEvents;
namespace Q=Moon::Chemistry::SemanticRandom;
namespace D=Moon::Species::ReactionData;
namespace fs=std::filesystem;
namespace {
using Node=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;
struct Input {
 S::Request source;Moon::Runtime::Configuration application;
 std::string applicationPath,caseKind,cohort,identity;
 std::vector<std::string> documents;
 std::uint64_t parents=0,steps=0;double duration=0,horizon=0,representedNumber=0;
 std::array<double,3> position{},velocity{};
 Moon::Species::ReactionExecution::ClosureBound numberBound;
 Moon::Species::ReducedKinematics::Bound rateBound;
 // Explicit numerical weights for the approved 1:1/1:8/8:1 reference matrix.
 // Production mode preserves normal initialized weights; represented numbers
 // are always independently checked against the actual W*correction product.
 std::string referenceWeightMode="production";
 double parentReferenceWeight=0,ionReferenceWeight=0;
};
std::string Json(const std::string&s){std::ostringstream o;o<<'"';for(unsigned char c:s){
 if(c=='"'||c=='\\')o<<'\\'<<c;else if(c=='\n')o<<"\\n";else if(c=='\r')o<<"\\r";else if(c=='\t')o<<"\\t";
 else if(c<32)o<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<static_cast<unsigned>(c)<<std::dec;else o<<c;}o<<'"';return o.str();}
struct Owner {
 const AMPS::Input::Block& block;std::set<std::string> used;
 std::string Text(const std::string&key){const auto*a=block.Find(key);if(!a)throw std::runtime_error(block.begin.Prefix()+"missing mandatory key "+key);used.insert(key);return a->value;}
 std::string Path(const std::string&key){auto s=Text(key);if(s.size()>=2&&s.front()=='"'&&s.back()=='"')s=s.substr(1,s.size()-2);
 const auto*a=block.Find(key);fs::path p(s);if(p.is_relative())p=fs::path(a->where.file).parent_path()/p;return fs::absolute(p).lexically_normal().string();}
 long double Number(const std::string&key){auto s=Text(key);std::size_t i=0;long double v=0;try{v=std::stold(s,&i);}catch(...){throw std::runtime_error("invalid numerical value: "+key);}
 if(i!=s.size()||!std::isfinite(v))throw std::runtime_error("nonfinite/partial numerical value: "+key);return v;}
 double Double(const std::string&key){const long double v=Number(key);const double x=static_cast<double>(v);if(!std::isfinite(x)||(v!=0&&x==0))throw std::runtime_error("unrepresentable binary64 input: "+key);return x;}
 std::uint64_t Uint(const std::string&key){const auto s=Text(key);if(s.empty()||(s.size()>1&&s.front()=='0'))throw std::runtime_error("noncanonical unsigned integer: "+key);
 for(char c:s)if(c<'0'||c>'9')throw std::runtime_error("unsigned integer required: "+key);
 std::size_t i=0;auto v=std::stoull(s,&i);if(i!=s.size())throw std::runtime_error("invalid uint64: "+key);return v;}
 void Complete(){for(const auto&a:block.assignments)if(!used.count(a.key))throw std::runtime_error(a.where.Prefix()+"unknown owned key "+a.key);}
};
Input Parse(const std::string& path,bool sourceOnly){
 AMPS::Input::Document doc;auto st=AMPS::Input::ParseFile(path,&doc);if(!st.ok())throw std::runtime_error(st.message);
 if(!doc.root.assignments.empty()||doc.root.blocks.size()!=1||doc.root.blocks[0].kind!="sodium_qualification"||doc.root.blocks[0].name.empty())throw std::runtime_error("one named sodium_qualification block required");
 const auto& b=doc.root.blocks[0];
 const std::set<std::string> known={"application_input","campaign_seed","case_kind","channel_id","constants_file","constants_sha256","crosswalk_file","crosswalk_sha256","earth_occultation","electron_maximum_relative_acceleration_error","electron_maximum_relative_inertial_mass_error","electron_maximum_relative_relativistic_kinetic_error","electron_state_id","energy_absolute_bound","energy_reference_identity","energy_relative_bound","exposure_absolute_depth_tolerance","exposure_maximum_geometry_intervals","exposure_maximum_segments","exposure_relative_depth_tolerance","exposure_root_time_tolerance_s","genealogy_capacity","hazard_maximum_root_refinements","inertial_frame","inertial_origin","ion_compiled_species","ion_maximum_relative_acceleration_error","ion_maximum_relative_inertial_mass_error","ion_maximum_relative_relativistic_kinetic_error","ion_state_id","isotopes_file","isotopes_sha256","mass_absolute_bound","mass_relative_bound","maximum_classical_speed_ms","model_identity","number_absolute_bound","number_relative_bound","momentum_absolute_bound","momentum_relative_bound","original_parents","parent_compiled_species","parent_maximum_relative_acceleration_error","parent_maximum_relative_inertial_mass_error","parent_maximum_relative_relativistic_kinetic_error","parent_state_id","physical_horizon_s","physical_step_s","physical_steps","planner_energy_absolute_bound","planner_energy_relative_bound","planner_mass_absolute_bound","planner_mass_relative_bound","planner_momentum_absolute_bound","planner_momentum_relative_bound","profile_metadata_file","profile_metadata_sha256","profile_native_manifest_file","profile_native_manifest_sha256","represented_number_per_parent","rotation_absolute_bound","rotation_relative_bound","rate_absolute_bound","rate_relative_bound","schema","solar_endpoint","source_position_lso_m_0","source_position_lso_m_1","source_position_lso_m_2","source_reference_distance_identity","source_reference_distance_m","source_velocity_lso_ms_0","source_velocity_lso_ms_1","source_velocity_lso_ms_2","states_file","states_sha256","threshold_absolute_bound","threshold_relative_bound","transport_policy_evidence_file","transport_policy_evidence_sha256","transport_policy_identity","unit_vector_tolerance"};
 for(const auto&a:b.assignments)if(!known.count(a.key)&&a.key!="photon_energy_policy"&&a.key!="reference_weight_mode"&&a.key!="parent_reference_weight"&&a.key!="ion_reference_weight")throw std::runtime_error(a.where.Prefix()+"unknown sodium_qualification key "+a.key);
 Owner o{b,{}};Input c;c.cohort=b.name;c.documents=doc.files;
 if(o.Text("schema")!="srcMoon-native-sodium-cohort-v1")throw std::runtime_error("unsupported native sodium cohort schema");
 c.applicationPath=o.Path("application_input");std::string error;
 if(!Moon::Runtime::ParseApplicationInput(c.applicationPath,&c.application,&error))throw std::runtime_error(error);
 auto&r=c.source;r.profileMetadataPath=o.Path("profile_metadata_file");r.profileMetadataSha256=o.Text("profile_metadata_sha256");
 r.profileNativeManifestPath=o.Path("profile_native_manifest_file");r.profileNativeManifestSha256=o.Text("profile_native_manifest_sha256");
 r.profile.isotopes={o.Path("isotopes_file"),o.Text("isotopes_sha256")};r.profile.constants={o.Path("constants_file"),o.Text("constants_sha256")};r.profile.states={o.Path("states_file"),o.Text("states_sha256")};
 for(const auto&child:b.blocks){if(child.kind!="source"||child.name.empty()||!child.blocks.empty())throw std::runtime_error(child.begin.Prefix()+"only named source children allowed");
  Owner item{child,{}};r.profile.sources.push_back({child.name,item.Path("path"),item.Text("sha256"),item.Text("role")});item.Complete();}
 if(r.profile.sources.empty())throw std::runtime_error("physical state profile needs its explicitly pinned source files");
 r.parentStateId=o.Text("parent_state_id");r.ionStateId=o.Text("ion_state_id");r.electronStateId=o.Text("electron_state_id");
 if(sourceOnly)return c; // Source availability phase: no numeric/physical default chosen.
 const auto parent=o.Uint("parent_compiled_species"),ion=o.Uint("ion_compiled_species");
 if(parent>std::numeric_limits<int>::max()||ion>std::numeric_limits<int>::max())throw std::runtime_error("compiled species index outside int range");
 r.parentCompiledSpecies=static_cast<int>(parent);r.ionCompiledSpecies=static_cast<int>(ion);
 r.crosswalkPath=o.Path("crosswalk_file");r.crosswalkSha256=o.Text("crosswalk_sha256");
 r.random.seed=o.Uint("campaign_seed");r.random.cohort=b.name;const auto capacity=o.Uint("genealogy_capacity");
 if(capacity>std::numeric_limits<std::size_t>::max())throw std::runtime_error("genealogy capacity outside native size domain");r.random.genealogyCapacity=static_cast<std::size_t>(capacity);
 c.parents=o.Uint("original_parents");c.steps=o.Uint("physical_steps");c.duration=o.Double("physical_step_s");c.horizon=o.Double("physical_horizon_s");c.representedNumber=o.Double("represented_number_per_parent");
 c.numberBound.absolute=o.Double("number_absolute_bound");c.numberBound.relative=o.Double("number_relative_bound");
 c.rateBound.absolute=o.Number("rate_absolute_bound");c.rateBound.relative=o.Number("rate_relative_bound");
 if(!c.parents||!c.steps||!(c.duration>0)||!(c.horizon>0)||!(c.representedNumber>0)||static_cast<long double>(c.steps)*c.duration!=c.horizon)throw std::runtime_error("positive complete finite cohort and exact common-step horizon required; declare representable equal steps, never truncate last physical interval");
 for(int k=0;k<3;++k){c.position[k]=o.Double("source_position_lso_m_"+std::to_string(k));c.velocity[k]=o.Double("source_velocity_lso_ms_"+std::to_string(k));}
 c.caseKind=o.Text("case_kind");if(c.caseKind!="homogeneous-fixed-position"&&c.caseKind!="trajectory-integrated-hazard")throw std::runtime_error("unknown explicit cohort case_kind");
 if(c.caseKind=="homogeneous-fixed-position")for(double v:c.velocity)if(v!=0)throw std::runtime_error("homogeneous fixed-position parent requires explicitly zero initial velocity");
 r.inertialOrigin=o.Text("inertial_origin");r.inertialFrame=o.Text("inertial_frame");r.sourceReferenceDistanceM=o.Double("source_reference_distance_m");r.sourceReferenceDistanceIdentity=o.Text("source_reference_distance_identity");
 const auto occult=o.Text("earth_occultation");if(occult!="on"&&occult!="off")throw std::runtime_error("earth_occultation must explicitly be on/off");r.earthOccultation=occult=="on";
 r.exposure.rootTimeToleranceS=o.Double("exposure_root_time_tolerance_s");r.exposure.absoluteDepthTolerance=o.Number("exposure_absolute_depth_tolerance");r.exposure.relativeDepthTolerance=o.Number("exposure_relative_depth_tolerance");
 r.exposure.maximumSegments=o.Uint("exposure_maximum_segments");r.exposure.maximumGeometryIntervals=o.Uint("exposure_maximum_geometry_intervals");r.maximumRootRefinements=o.Uint("hazard_maximum_root_refinements");
 auto&d=r.controls;d.identity=o.Text("model_identity");d.channelId=o.Text("channel_id");d.energyReferenceIdentity=o.Text("energy_reference_identity");
 if(b.Find("photon_energy_policy")){
  const auto energyPolicy=o.Text("photon_energy_policy");
  if(energyPolicy=="source-effective-photon")d.photonEnergyPolicy=Moon::Species::ReducedKinematics::PhotonEnergyPolicy::SourceEffectivePhoton;
  else if(energyPolicy!="source-equals-state")throw std::runtime_error("photon_energy_policy must be source-equals-state or source-effective-photon");
 }
 if(b.Find("reference_weight_mode"))c.referenceWeightMode=o.Text("reference_weight_mode");
 if(c.referenceWeightMode=="explicit"){
  c.parentReferenceWeight=o.Double("parent_reference_weight");c.ionReferenceWeight=o.Double("ion_reference_weight");
  if(!(c.parentReferenceWeight>0)||!(c.ionReferenceWeight>0))throw std::runtime_error("explicit reference weights must be finite positive represented-number weights");
 }else if(c.referenceWeightMode!="production"||b.Find("parent_reference_weight")||b.Find("ion_reference_weight")){
  throw std::runtime_error("reference_weight_mode must be production or explicit; only explicit mode admits caller-supplied reference weights");
 }
 const auto endpoint=o.Text("solar_endpoint");if(endpoint!="quiet"&&endpoint!="active")throw std::runtime_error("separate quiet/active endpoint required");d.solar=endpoint=="quiet"?Moon::Species::Photochemistry::SolarCase::Quiet:Moon::Species::Photochemistry::SolarCase::Active;
 for(auto pair:{std::make_pair("parent",&d.parent),std::make_pair("ion",&d.ion),std::make_pair("electron",&d.electron)}){
  const std::string prefix=pair.first;auto& limits=pair.second->modelErrorLimits;
  limits.maximumRelativeInertialMassError=o.Number(prefix+"_maximum_relative_inertial_mass_error");limits.maximumRelativeAccelerationError=o.Number(prefix+"_maximum_relative_acceleration_error");limits.maximumRelativeRelativisticKineticError=o.Number(prefix+"_maximum_relative_relativistic_kinetic_error");}
 d.numerical.maximumClassicalSpeedMS=o.Number("maximum_classical_speed_ms");d.numerical.unitVectorTolerance=o.Number("unit_vector_tolerance");
 for(auto pair:{std::make_pair("mass",&d.numerical.referenceMassKg),std::make_pair("momentum",&d.numerical.momentumKgMS),std::make_pair("energy",&d.numerical.energyJ),std::make_pair("threshold",&d.numerical.thresholdJ),std::make_pair("rotation",&d.rotation)}){
  pair.second->absolute=o.Number(std::string(pair.first)+"_absolute_bound");pair.second->relative=o.Number(std::string(pair.first)+"_relative_bound");}
 auto& policy=d.plannerAuthorization;policy.identity=o.Text("transport_policy_identity");policy.evidenceIdentity=o.Text("transport_policy_evidence_sha256");
 const auto policyPath=o.Path("transport_policy_evidence_file");if(!D::VerifyFileSha256(policyPath,policy.evidenceIdentity,&error))throw std::runtime_error(error);
 // This is admission to the already approved EXPERIMENTAL reduced model, not
 // a formal physical validation flag. Actual source states/mass/Q/errors are
 // independently refused by the production composer before any startup.
 policy.qualified=true;policy.eventFrame=Moon::Frames::SolarOrbital;policy.maximumClassicalSpeedMS=static_cast<double>(d.numerical.maximumClassicalSpeedMS);
 for(auto pair:{std::make_pair("planner_mass",&policy.bounds.modelMassKg),std::make_pair("planner_momentum",&policy.bounds.momentumKgMS),std::make_pair("planner_energy",&policy.bounds.energyJ)}){
  pair.second->absolute=o.Double(std::string(pair.first)+"_absolute_bound");pair.second->relative=o.Double(std::string(pair.first)+"_relative_bound");}
 o.Complete();
 if(c.application.particleTracingMover=="legacy"||c.application.tracingClockMode!="common-time"||c.application.simulationEnvironmentMode!="frozen"||c.application.m2EnvironmentEnabled)throw std::runtime_error("conditional physical cohort requires explicit first/second-order common-time frozen geometry with M2 disabled");
 c.identity=fs::absolute(path).string();return c;
}
void InputReceipt(const Input& c,const fs::path& out,bool sourceOnly,const Moon::Species::TrajectoryExposure::Field* geometry=nullptr){
 if(PIC::ThisThread!=0)return;
 std::ofstream f(out/"input_receipt.json");f<<std::setprecision(21);
 f<<"{\n  \"schema\": \"srcMoon-native-sodium-input-receipt-v1\",\n  \"profile\": {\"path\": "<<Json(c.source.profileMetadataPath)<<", \"sha256\": "<<Json(c.source.profileMetadataSha256)<<"},\n  \"mapping\": {\"parent\": "<<Json(c.source.parentStateId)<<", \"ion\": "<<Json(c.source.ionStateId)<<", \"electron\": "<<Json(c.source.electronStateId);
 if(!sourceOnly)f<<", \"event_frame\": \"LSO\", \"profile_sha256\": "<<Json(c.source.profileMetadataSha256)<<", \"ion_compiled_species\": "<<c.source.ionCompiledSpecies<<", \"ion_catalog_species\": 1, \"energy_reference_identity\":"<<Json(c.source.controls.energyReferenceIdentity)<<", \"photon_energy_policy\": "<<Json(c.source.controls.photonEnergyPolicy==Moon::Species::ReducedKinematics::PhotonEnergyPolicy::SourceEqualsState?"source-equals-state":"source-effective-photon")<<", \"source_crosswalk\":{\"path\":"<<Json(c.source.crosswalkPath)<<",\"sha256\":"<<Json(c.source.crosswalkSha256)<<"}";
 f<<"},\n  \"phase\": "<<Json(sourceOnly?"source-availability-only":"explicit-experimental-physical-input");
 if(!sourceOnly){f<<",\n  \"bounds\": {";bool first=true;for(const auto& item:{std::make_pair("mass",c.source.controls.plannerAuthorization.bounds.modelMassKg),std::make_pair("momentum",c.source.controls.plannerAuthorization.bounds.momentumKgMS),std::make_pair("energy",c.source.controls.plannerAuthorization.bounds.energyJ)}){
  if(!first)f<<',';first=false;f<<Json(item.first)<<": {\"absolute\": \""<<item.second.absolute<<"\", \"relative\": \""<<item.second.relative<<"\"}";}
 f<<",\"rate_s\":{\"absolute\":\""<<c.rateBound.absolute<<"\",\"relative\":\""<<c.rateBound.relative<<"\"},\"number\":{\"absolute\":\""<<c.numberBound.absolute<<"\",\"relative\":\""<<c.numberBound.relative<<"\"},\"rotation\":{\"absolute\":\""<<c.source.controls.rotation.absolute<<"\",\"relative\":\""<<c.source.controls.rotation.relative<<"\"}";
 f<<"},\n  \"original_parents\": "<<c.parents<<",\n  \"physical_steps\": "<<c.steps<<",\n  \"physical_step_s\": \""<<c.duration<<"\",\n  \"physical_horizon_s\": \""<<c.horizon<<"\",\n  \"initial_represented_number\": \""<<c.representedNumber<<"\",\n  \"case_kind\": "<<Json(c.caseKind)<<",\n  \"cohort_id\": "<<Json(c.cohort)<<",\n  \"seed\": "<<c.source.random.seed<<",\n  \"source_reference_distance_m\": \""<<c.source.sourceReferenceDistanceM<<"\",\n  \"source_reference_distance_identity\": "<<Json(c.source.sourceReferenceDistanceIdentity);}
 if(!sourceOnly){
  f<<",\n  \"reference_weight_mode\": "<<Json(c.referenceWeightMode);
  if(c.referenceWeightMode=="explicit")f<<",\n  \"configured_parent_reference_weight\": \""<<c.parentReferenceWeight<<"\",\n  \"configured_ion_reference_weight\": \""<<c.ionReferenceWeight<<"\"";
 }
 if(geometry){
  const auto* transport=Moon::GetCommonTransportClock();if(!transport)throw std::runtime_error("actual native common transport clock unavailable after normal startup");
  f<<",\n  \"transport_origin_et_s\": \""<<transport->originEtS()<<"\"";
  f<<",\n  \"geometry_basis_et_s\": \""<<geometry->frozenGeometryEpochS<<"\",\n  \"sun_position_lso_m\": [\""<<geometry->sunPositionM[0]<<"\",\""<<geometry->sunPositionM[1]<<"\",\""<<geometry->sunPositionM[2]<<"\"],\n  \"source_position_lso_m\": [\""<<c.position[0]<<"\",\""<<c.position[1]<<"\",\""<<c.position[2]<<"\"],\n  \"actual_source_distance_m\": \""<<std::hypot(geometry->sunPositionM[0]-c.position[0],geometry->sunPositionM[1]-c.position[1],geometry->sunPositionM[2]-c.position[2])<<"\",\n  \"total_acceleration_macro\": "<<Json(NativeAccelerationMacro);
  f<<",\n  \"solar_endpoint\": "<<Json(c.source.controls.solar==Moon::Species::Photochemistry::SolarCase::Quiet?"quiet":"active")<<",\n  \"source_channels\": [";
  const auto* source=Moon::Runtime::GetPhotochemistryRegistry();std::vector<Moon::Species::Photochemistry::Channel> channels;std::string error;
  if(!source||!source->ChannelsForParent("Na",&channels,&error))throw std::runtime_error(error);
  bool channelFirst=true;for(const auto& channel:channels){if(!channelFirst)f<<',';channelFirst=false;
   f<<"{\"id\":"<<Json(channel.id)<<",\"quiet_rate_1au_s\":\""<<channel.rate1AuS[0]<<"\",\"active_rate_1au_s\":\""<<channel.rate1AuS[1]<<"\",\"threshold_ev\":\""<<channel.thresholdEnergyEv<<"\",\"quiet_mean_excess_ev\":\""<<channel.meanExcessEnergyEv[0]<<"\",\"active_mean_excess_ev\":\""<<channel.meanExcessEnergyEv[1]<<"\"}";}
  f<<"]";
  const auto endpoint=c.source.controls.solar==Moon::Species::Photochemistry::SolarCase::Quiet?0:1;
  Moon::Species::TrajectoryExposure::Arc point;point.positionM=c.position;Moon::Species::TrajectoryExposure::History history;
  if(channels.size()!=1||!Moon::Species::TrajectoryExposure::BuildPhotonHistory(point,c.duration,*geometry,{{0,channels.front().rate1AuS[endpoint]}},{},c.source.exposure,&history,&error)||history.segments.empty())throw std::runtime_error("initial homogeneous point geometry unavailable: "+error);
  const bool visible=history.segments.front().channels.front().beginRateS>0;
  f<<",\n  \"point_sun_visible\":"<<(visible?1:0)<<",\n  \"homogeneous_reference\": {\"sun_position_lso_m\":[\""<<geometry->sunPositionM[0]<<"\",\""<<geometry->sunPositionM[1]<<"\",\""<<geometry->sunPositionM[2]<<"\"],\"source_position_lso_m\":[\""<<c.position[0]<<"\",\""<<c.position[1]<<"\",\""<<c.position[2]<<"\"],\"point_sun_visible\":"<<(visible?1:0)<<",\"source_reference_distance_m\":\""<<geometry->referenceDistanceM<<"\",\"source_reference_distance_identity\":"<<Json(c.source.sourceReferenceDistanceIdentity)<<",\"solar_endpoint\":"<<Json(endpoint==0?"quiet":"active")<<",\"total_acceleration_macro\":"<<Json(NativeAccelerationMacro)<<",\"geometry_basis_epoch_et_s\":\""<<geometry->frozenGeometryEpochS<<"\"}";
 }
 f<<",\n  \"input_documents\": [";
 auto documents=c.documents;AMPS::Input::Document application;
 const auto parsed=AMPS::Input::ParseFile(c.applicationPath,&application);if(!parsed.ok())throw std::runtime_error(parsed.message);
 for(const auto& path:application.files)if(std::find(documents.begin(),documents.end(),path)==documents.end())documents.push_back(path);
 bool first=true;for(const auto& path:documents){std::string digest,error;if(!D::FileSha256(path,&digest,&error))throw std::runtime_error(error);if(!first)f<<',';first=false;f<<"{\"path\":"<<Json(path)<<",\"sha256\":"<<Json(digest)<<'}';}
 f<<"]\n}\n";
}
void Emit(const fs::path& out,const char* status,const std::string& reason,int code){
 if(PIC::ThisThread!=0)return;
 std::ofstream f(out/"result.json");f<<"{\n  \"schema\": \"srcMoon-native-sodium-cohort-result-v1\",\n  \"status\": "<<Json(status)<<",\n  \"native_returncode\": "<<code<<",\n  \"scope\": \"conditional physical cohort execution; no statistical or formal M3 qualification\",\n  \"formal_M3_status\": \"SKIPPED / NOT QUALIFIED\",\n  \"reason\": "<<Json(reason)<<"\n}\n";
}
void EachParticle(const std::function<void(long int,Node*)>& fn){
 for(unsigned i=0;i<PIC::DomainBlockDecomposition::nLocalBlocks;++i){auto* n=PIC::DomainBlockDecomposition::BlockTable[i];if(!n||!n->block)continue;
 for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)for(long int p=n->block->FirstCellParticleTable[cell];p!=-1;p=PIC::ParticleBuffer::GetNext(p))fn(p,n);}
}
int Execute(const Input& c,std::shared_ptr<S::Prepared> prepared,const fs::path& out){
 auto status=prepared->RegisterForNormalStartup();if(!status.ok())throw std::runtime_error(status.message);
 amps_init();PIC::DomainBlockDecomposition::UpdateBlockTable();InputReceipt(c,out,false,&prepared->geometry());
 if(PIC::nTotalThreads!=1&&PIC::nTotalThreads!=8)throw std::runtime_error("declared qualification layouts are1/8MPI withOMP1");
 if(c.caseKind=="homogeneous-fixed-position"){
  if(PIC::Mover::ForceCutSurface::Acceleration!=nullptr||!Moon::SodiumCohortInput::ExactZeroAccelerationMacro(NativeAccelerationMacro))
    throw std::runtime_error("homogeneous source must use the actual configured core zero-acceleration macro with no overriding callback; a marker flag alone is not proof");
#if !defined(_MOON_COHORT_ZERO_FORCE_REFERENCE_) || _MOON_COHORT_ZERO_FORCE_REFERENCE_ != 1
  throw std::runtime_error("homogeneous reference requires the explicitly reviewed zero-force compile recipe and effective macro proof; a normal Moon-force binary cannot serve this case");
#endif
 }
 // The explicit physical cohort dt replaces only this isolated reference
 // experiment's actual native dt table; it is never inferred from species.
 // In global modes the normal public SetLocalTimeStep API is a MINIMUM
 // reduction. It cannot increase the ordinary species flight/gyro step to
 // this declared zero-force reference interval. Reset the global reduction
 // accumulators to their normal uninitialized sentinel BEFORE calling that
 // same API, after amps_init and before the first common step or particle.
 // Without this reset Na and Na+ retain different normal startup minima and
 // BeginCommonTransportStep correctly refuses the mixed clocks. Local tables
 // are direct assignments and require no reset. No normal application mode,
 // common-time admission rule, core setter, or accepted clock is changed.
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
 for(int species=0;species<PIC::nTotalSpecies;++species)
   PIC::ParticleWeightTimeStep::GlobalTimeStep[species]=-1.0;
#elif _SIMULATION_TIME_STEP_MODE_ == _SINGLE_GLOBAL_TIME_STEP_
 PIC::ParticleWeightTimeStep::GlobalTimeStep[0]=-1.0;
#elif _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
 throw std::runtime_error("explicit cohort timestep requires a supported actual PIC timestep table mode");
#endif
 for(unsigned i=0;i<PIC::DomainBlockDecomposition::nLocalBlocks;++i){auto* n=PIC::DomainBlockDecomposition::BlockTable[i];if(!n||!n->block)continue;for(int species=0;species<PIC::nTotalSpecies;++species)n->block->SetLocalTimeStep(c.duration,species);}
 // Read the actual public values back; a requested input or marker is not
 // proof that the common physical interval is installed. The normal common
 // clock independently performs its collective equality check at Begin.
 for(unsigned i=0;i<PIC::DomainBlockDecomposition::nLocalBlocks;++i){auto* n=PIC::DomainBlockDecomposition::BlockTable[i];if(!n||!n->block)continue;
   for(int species=0;species<PIC::nTotalSpecies;++species)
     if(n->block->GetLocalTimeStep(species)!=c.duration)
       throw std::runtime_error("actual native timestep differs from the explicitly declared cohort interval before particle creation");
 }
 // Exercise the normal PIC weight tables through their public setters. This
 // is the declared numerical weight-ratio experiment, not a daughter dt factor
 // or a change to the normal application's weight policy. Source number and
 // all live product readbacks remain subject to the existing exact checks.
 if(c.referenceWeightMode=="explicit")for(unsigned index=0;index<PIC::DomainBlockDecomposition::nLocalBlocks;++index){
  auto* block=PIC::DomainBlockDecomposition::BlockTable[index];if(!block||!block->block)continue;
  block->block->SetLocalParticleWeight(c.parentReferenceWeight,c.source.parentCompiledSpecies);
  block->block->SetLocalParticleWeight(c.ionReferenceWeight,c.source.ionCompiledSpecies);
 }
 if(c.parents>static_cast<std::uint64_t>(std::numeric_limits<int>::max()))throw std::runtime_error("parent reporting exceeds MPI int-count resource domain; use a separately declared cluster cohort/report partition, never drop parents");
 // A declared exospheric source cannot begin inside/on the solid. Sphere
 // uses its actual compiled radius. LOLA is the production radially generated
 // closed mesh; query its SAME body-fixed BVH along the source radial ray,
 // rather than substituting a reference sphere or a test-only intersection.
 const double sourceRadius=std::hypot(c.position[0],c.position[1],c.position[2]);
 if(!(sourceRadius>0))throw std::runtime_error("source point is at/inside the lunar centre");
 if(c.application.surfaceGeometry==Moon::Runtime::SurfaceGeometry::Sphere){
  if(!(sourceRadius>_RADIUS_(_TARGET_)+c.application.tracingGeometryToleranceM))throw std::runtime_error("source point is inside/on the actual spherical lunar boundary; declare an exterior source point");
 }else{
#if defined(_NO_SPICE_CALLS_)
  throw std::runtime_error("LOLA source exterior classification requires furnished actual frame SPICE");
#else
  const auto* terrain=Moon::SurfaceEnvironment::ProductionTerrain();if(!terrain)throw std::runtime_error("actual LOLA BVH unavailable for source admission");
  SpiceDouble rotation[3][3];pxform_c(Moon::Frames::SolarOrbital,Moon::Frames::BodyFixed,prepared->geometry().frozenGeometryEpochS,rotation);
  if(failed_c()){reset_c();throw std::runtime_error("cannot transform explicit source radial query into actual LOLA frame");}
  std::array<double,3> direction{};for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)direction[i]+=rotation[i][j]*c.position[j]/sourceRadius;
  Moon::SurfaceEnvironment::RayHit hit;std::string error;
  if(!terrain->IntersectRay({{0,0,0}},direction,0,sourceRadius,std::numeric_limits<std::size_t>::max(),&hit,&error))throw std::runtime_error(error);
  if(!hit.hit||!(hit.distanceM+c.application.tracingGeometryToleranceM<sourceRadius))throw std::runtime_error("source point is inside/on the actual radially generated LOLA surface; no reference-sphere approximation or source relocation is applied");
#endif
 }
 auto* node=PIC::Mesh::mesh->findTreeNode(const_cast<double*>(c.position.data()));if(!node)throw std::runtime_error("explicit source point is outside actual normal mesh");
 if(node->Thread==PIC::ThisThread){
  const auto free=PIC::ParticleBuffer::GetMaxNPart()-PIC::ParticleBuffer::GetAllPartNum();
  if(free<=0||c.parents>=static_cast<std::uint64_t>(free))throw std::runtime_error("complete initial cohort plus at least one atomic product reservation does not fit actual particle buffer. Increase configured buffer capacity or reduce the predeclared parent count; no partial cohort is accepted");
  int i,j,k;PIC::Mesh::mesh->FindCellIndex(const_cast<double*>(c.position.data()),i,j,k,node,false);long&head=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
  for(std::uint64_t id=0;id<c.parents;++id){const auto p=PIC::ParticleBuffer::GetNewParticle();auto*b=PIC::ParticleBuffer::GetParticleDataPointer(p);
   PIC::ParticleBuffer::SetI(c.source.parentCompiledSpecies,b);double x[3],v[3];for(int axis=0;axis<3;++axis){x[axis]=c.position[axis];v[axis]=c.velocity[axis];}
   PIC::ParticleBuffer::SetX(x,b);PIC::ParticleBuffer::SetV(v,b);PIC::ParticleBuffer::SetIndividualStatWeightCorrection(c.representedNumber/node->block->GetLocalParticleWeight(c.source.parentCompiledSpecies),b);
   const double actualNumber=static_cast<double>(static_cast<long double>(node->block->GetLocalParticleWeight(c.source.parentCompiledSpecies))*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p));
   if(actualNumber!=c.representedNumber)throw std::runtime_error("requested represented parent number is not reproducible by actual native W*individual correction. Choose an explicitly representable cohort weight/reference-weight combination; no source weight is silently changed");
   PIC::Mover::ForceCutSurface::ResetParticleState(b);Moon::Chemistry::Controller::ResetForNewParticle(b);PIC::Mover::ForceCutSurface::SetParticleEpochS(p,Moon::GetCommonTransportClock()->originEtS());
   std::string error;if(!prepared->InitializeSourceKey(p,id,&error))throw std::runtime_error(error);
   PIC::ParticleBuffer::SetPrev(-1,b);PIC::ParticleBuffer::SetNext(head,b);if(head!=-1)PIC::ParticleBuffer::SetPrev(p,head);head=p;
  }
 }
 for(std::uint64_t step=0;step<c.steps;++step){std::string error;if(!Moon::BeginCommonTransportStep(Moon::GetCommonTransportClock()->originEtS(),&error))throw std::runtime_error(error);
  Moon::PrepareProductionParticleEnvironment();
  double moving=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
  if(!Moon::CommitCommonTransportStep(&error))throw std::runtime_error(error);
 }
 const std::size_t n=static_cast<std::size_t>(c.parents);std::vector<int> local(n,0),all(n,0);std::vector<double> times(n,0),allTimes(n,0);
 EachParticle([&](long p,Node*){Q::State key;std::string error;if(!Q::Read(p,&key,&error))throw std::runtime_error(error);
  if(key.originalParent>=c.parents)throw std::runtime_error("unexpected original parent identity");
  if(PIC::ParticleBuffer::GetI(p)==c.source.parentCompiledSpecies){if(!key.genealogy.empty())throw std::runtime_error("unexpected regenerated neutral genealogy");++local[key.originalParent];}
 });
 const auto& receipts=prepared->composedReceipts();if(receipts.size()!=prepared->ledger().events.size())throw std::runtime_error("composed source records are not paired one-for-one with accepted atomic ledger events");
 for(const auto& event:receipts){if(!event.key.genealogy.empty()||event.key.originalParent>=c.parents)throw std::runtime_error("invalid physical parent event identity");const auto id=static_cast<std::size_t>(event.key.originalParent);local[id]+=2;times[id]=static_cast<double>(event.relativeEventTimeS);}
 MPI_Reduce(local.data(),all.data(),static_cast<int>(n),MPI_INT,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);MPI_Reduce(times.data(),allTimes.data(),static_cast<int>(n),MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
 std::ofstream events(out/("events.rank"+std::to_string(PIC::ThisThread)+".csv"));events<<std::setprecision(21)<<"original_parent_id,relative_event_time_s,channel_id,parent_represented_number,bath_charge_elementary,physical_event_epoch_s\n";
 for(std::size_t i=0;i<receipts.size();++i){const auto&e=receipts[i];events<<e.key.originalParent<<','<<e.relativeEventTimeS<<','<<e.receipt.channelId<<','<<prepared->ledger().events[i].ledger.parentNumber<<','<<e.receipt.representedElectronNumberChange<<','<<e.receipt.frame.physicalEventEpochEtS<<'\n';}
 // Full actual source/parent/product/bath receipt; no Boolean conservation
 // flags substitute for an independent reconstruction from these SI moments.
 std::ofstream metrics(out/("event_metrics.rank"+std::to_string(PIC::ThisThread)+".csv"));metrics<<std::setprecision(21)<<"original_parent_id,key,value\n";
 std::ofstream hexadecimal(out/("event_metrics.rank"+std::to_string(PIC::ThisThread)+".hex.csv"));hexadecimal<<std::hexfloat<<"original_parent_id,key,value\n";
 for(std::size_t j=0;j<receipts.size();++j){const auto&e=receipts[j];const auto& receipt=e.receipt;const auto& ledger=prepared->ledger().events[j].ledger;
  auto scalar=[&](const std::string& key,long double value){metrics<<e.key.originalParent<<','<<key<<','<<value<<'\n';hexadecimal<<e.key.originalParent<<','<<key<<','<<value<<'\n';};
  auto vector=[&](const std::string& key,const auto& value){for(unsigned k=0;k<3;++k)scalar(key+std::to_string(k),value[k]);};
  const auto& initial=prepared->ledger().events[j].initialProducts;
  if(initial.size()!=1)throw std::runtime_error("accepted atomic event lacks exactly one ACTUAL initial PIC product snapshot");
  scalar("actual_ion_child_count",initial.size());
  const auto& live=initial.front();
  scalar("actual_parent_compiled_mass_kg",PIC::MolecularData::GetMass(c.source.parentCompiledSpecies));
  scalar("actual_ion_compiled_mass_kg",PIC::MolecularData::GetMass(live.compiledSpecies));
  scalar("actual_parent_charge_c",PIC::MolecularData::GetElectricCharge(c.source.parentCompiledSpecies));
  scalar("actual_ion_charge_c",PIC::MolecularData::GetElectricCharge(live.compiledSpecies));
  scalar("planned_interval_start_s",e.ionAtEvent.interval.startS());scalar("planned_interval_duration_s",e.ionAtEvent.interval.durationS());scalar("planned_interval_elapsed_s",e.ionAtEvent.interval.elapsedS());
  scalar("actual_ion_compiled_index",live.compiledSpecies);scalar("actual_ion_catalog_species",live.catalogSpecies);
  scalar("actual_ion_reference_weight",live.referenceWeight);scalar("actual_ion_individual_correction",live.individualCorrection);
  scalar("actual_ion_internal_j",live.internalEnergyJ);scalar("actual_ion_epoch_et_s",live.physicalEpochS);
  scalar("actual_ion_interval_start_s",live.interval.startS());scalar("actual_ion_interval_duration_s",live.interval.durationS());scalar("actual_ion_interval_elapsed_s",live.interval.elapsedS());
  vector("actual_ion_x",live.positionM);vector("actual_ion_v",live.velocityMS);
  vector("parent_event_x",e.parentAtEvent.positionM);vector("ion_event_x",e.ionAtEvent.positionM);
  scalar("parent_number",ledger.parentNumber);scalar("ion_number",e.ionAtEvent.representedNumber);
  scalar("parent_reference_weight",e.parentAtEvent.referenceWeight);scalar("parent_individual_correction",e.parentAtEvent.individualCorrection);
  scalar("ion_reference_weight",e.ionAtEvent.referenceWeight);scalar("ion_individual_correction",e.ionAtEvent.individualCorrection);
  scalar("parent_internal_j",e.parentAtEvent.internalEnergyJ);scalar("ion_internal_j",e.ionAtEvent.internalEnergyJ);
  vector("parent_event_v",e.parentAtEvent.velocityMS);vector("ion_event_v",e.ionAtEvent.velocityMS);
  vector("parent_inertial_v",receipt.parentInertialVelocityMS);vector("electron_event_v",receipt.outgoingElectron.eventVelocityMS);vector("electron_inertial_v",receipt.outgoingElectron.inertialVelocityMS);
  vector("photon_inertial_p",receipt.incomingPhoton.inertialMomentumKgMS);scalar("photon_inertial_j",receipt.incomingPhoton.inertialEnergyJ);
  vector("bath_event_p",receipt.bathEventMomentumChangeKgMS);scalar("bath_event_j",receipt.bathEventEnergyChangeJ);scalar("bath_mass_kg",receipt.bathModelMassChangeKg);
  scalar("electron_state_j",receipt.outgoingElectron.stateEnergyJ);scalar("source_threshold_j",receipt.nativeThresholdJ);scalar("source_mean_excess_j",receipt.matchedSourceMeanExcessJ);scalar("source_q_j",receipt.stateQJ);
  scalar("physical_threshold_j",receipt.physicalThresholdJ);scalar("physical_minus_source_threshold_j",receipt.physicalMinusSourceThresholdJ);
  scalar("source_rate_1au_s",receipt.sourceRate1AuS);scalar("source_rate_event_s",receipt.eventRateS);scalar("accepted_rate_s",e.acceptedInstantaneousRateS);
  scalar("relative_event_time_s",e.relativeEventTimeS);scalar("event_epoch_et_s",receipt.frame.physicalEventEpochEtS);scalar("geometry_basis_et_s",receipt.frame.geometryBasisEpochEtS);
  vector("velocity_offset_event",receipt.frame.velocityOffsetEventMS);for(unsigned k=0;k<9;++k)scalar("inertial_to_event_r"+std::to_string(k),receipt.frame.inertialToEventRotation[k]);
  scalar("ledger_parent_mass_kg",ledger.parentModelMassKg);scalar("ledger_ion_mass_kg",ledger.productModelMassKg);scalar("ledger_bath_mass_kg",ledger.bathModelMassChangeKg);
  vector("ledger_tracked_dp",ledger.trackedMomentumChangeKgMS);vector("ledger_bath_dp",ledger.bathMomentumChangeKgMS);
  scalar("ledger_tracked_de_j",ledger.trackedEnergyChangeJ);scalar("ledger_bath_de_j",ledger.bathEnergyChangeJ);
  scalar("tracked_charge_change_e",ledger.trackedChargeChangeElementary);scalar("bath_charge_change_e",ledger.bathChargeChangeElementary);
  scalar("exposure_error_bound",prepared->ledger().events[j].exposureApproximationErrorBound);scalar("clock_storage_error_bound",prepared->ledger().events[j].clockStorageErrorBound);
 }
 int valid=1;
 if(PIC::ThisThread==0){std::ofstream rows(out/"parent_records.csv");rows<<std::setprecision(17)<<"original_parent_id,status,relative_event_time_s,channel_id,initial_represented_number\n";
  for(std::size_t id=0;id<n;++id){if(all[id]!=1&&all[id]!=2){valid=0;break;}rows<<id<<','<<(all[id]==1?"SURVIVED":"REACTED")<<',';if(all[id]==2)rows<<allTimes[id];rows<<',';if(all[id]==2)rows<<c.source.controls.channelId;rows<<','<<c.representedNumber<<'\n';}
 }
 MPI_Bcast(&valid,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);if(!valid)throw std::runtime_error("initial parent missing/duplicated or physically censored: homogeneous binomial scoring is refused; retain boundary/event ledgers and choose a separately declared independent transport oracle");
 if(PIC::ThisThread==0)Emit(out,"PASS","actual declared finite cohort completed; statistics/numerical bias/budget/predecessor qualification remain separate",0);
 return 0;
}
}
int main(int argc,char**argv){
 std::string input,mode;fs::path out;
 for(int i=1;i<argc;++i){const std::string key=argv[i];if(i+1>=argc){std::cerr<<"each option needs a value\n";return 2;}const std::string value=argv[++i];if(key=="--input")input=value;else if(key=="--output-dir")out=value;else if(key=="--mode")mode=value;else return 2;}
 if(input.empty()||out.empty()||(mode!="preflight"&&mode!="run"&&mode!="source-preflight"))return 2;
 bool ownOutput=false;
 try{
  // Lightweight MPI setup permits one writer and coherent preflight outcomes;
  // it allocates no mesh/SPICE/particles. amps_init's guarded re-entry uses the
  // SAME normal communicator later. Semantic keys never use this rank.
  PIC::InitMPI();int ready=1;
  if(PIC::ThisThread==0){if(fs::exists(out))ready=0;else {fs::create_directories(out);ownOutput=true;}}
  MPI_Bcast(&ready,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  if(!ready)throw std::runtime_error("output directory already exists; preserve previous attempt");
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  auto c=Parse(input,mode=="source-preflight");InputReceipt(c,out,mode=="source-preflight");std::string error;
  if(mode=="source-preflight"){
    S::P::Registry profile;const auto status=S::Prepared::InspectProfile(c.source,&profile);
    const auto admission=Moon::SodiumCohortInput::CollectiveAdmission(MPI_GLOBAL_COMMUNICATOR,
        status.ok()?0:status.code==E::Code::Unavailable?1:2,status.message,"source-profile");
    const int code=admission.severity==1?77:admission.severity==2?2:0;
    if(PIC::ThisThread==0)Emit(out,code==0?"PASS":code==77?"SKIPPED":"ERROR",code==0?"actual source states/rest masses inspected; composer/numerics/transport NOT evaluated":admission.diagnostic,code);
    MPI_Finalize();return code;
  }
  if(!Moon::Runtime::InstallConfiguration(c.application,&error))throw std::runtime_error(error);
  const auto* source=Moon::Runtime::GetPhotochemistryRegistry();
  const auto sourceAdmission=Moon::SodiumCohortInput::CollectiveAdmission(MPI_GLOBAL_COMMUNICATOR,
      source?0:1,source?"":"explicit frozen normal PHIDRATES service unavailable","runtime-source-registry");
  if(sourceAdmission.severity){const int code=sourceAdmission.severity==1?77:2;Emit(out,code==77?"SKIPPED":"ERROR",sourceAdmission.diagnostic,code);MPI_Finalize();return code;}
  std::shared_ptr<S::Prepared> prepared;const auto status=S::Prepared::Create(*source,c.source,&prepared);
  const auto preparedAdmission=Moon::SodiumCohortInput::CollectiveAdmission(MPI_GLOBAL_COMMUNICATOR,
      status.ok()?0:status.code==E::Code::Unavailable?1:2,status.message,"physical-composer");
  if(preparedAdmission.severity){const int code=preparedAdmission.severity==1?77:2;Emit(out,code==77?"SKIPPED":"ERROR",preparedAdmission.diagnostic,code);MPI_Finalize();return code;}
  if(mode=="preflight"){Emit(out,"PASS","actual source/state/composer preparation evaluated; no simulation/statistical/formal qualification",0);MPI_Finalize();return 0;}
  const int code=Execute(c,prepared,out);MPI_Finalize();return code;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';if(ownOutput)Emit(out,"ERROR",e.what(),2);
  int initialized=0;MPI_Initialized(&initialized);if(initialized)MPI_Abort(MPI_GLOBAL_COMMUNICATOR,2);return 2;}
}
