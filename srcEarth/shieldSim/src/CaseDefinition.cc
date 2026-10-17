#include "CaseDefinition.hh"

#include "ContentHash.hh"
#include "PublishedSourceAdapters.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <set>
#include <sstream>
#include <tuple>

namespace shield_sim {
namespace {

namespace fs = std::filesystem;

std::string RequiredString(const Json::Object& object, const std::string& field,
                           const std::string& context) {
  const auto found = object.find(field);
  if (found == object.end()) throw JsonError(context + " is missing required field '" + field + "'");
  const std::string value = found->second.AsString(context + '.' + field);
  if (value.empty()) throw JsonError(context + '.' + field + " must not be empty");
  return value;
}

double FiniteNumber(const Json& value, const std::string& context) {
  const double result = value.AsNumber(context);
  if (!std::isfinite(result)) throw JsonError(context + " must be finite");
  return result;
}

double PositiveNumber(const Json& value, const std::string& context) {
  const double result = FiniteNumber(value, context);
  if (result <= 0.0) throw JsonError(context + " must be > 0");
  return result;
}

long PositiveInteger(const Json& value, const std::string& context) {
  const double number = PositiveNumber(value, context);
  if (std::floor(number) != number || number > static_cast<double>(std::numeric_limits<long>::max()))
    throw JsonError(context + " must be a positive integer in range");
  return static_cast<long>(number);
}

long NonNegativeInteger(const Json& value, const std::string& context) {
  const double number = FiniteNumber(value, context);
  if (number < 0.0 || std::floor(number) != number ||
      number > static_cast<double>(std::numeric_limits<long>::max()))
    throw JsonError(context + " must be a non-negative integer in range");
  return static_cast<long>(number);
}

double LengthMm(const Json& value, const std::string& context, bool allowZero = false) {
  const auto& object = value.AsObject(context);
  RejectUnknownFields(object, {"value","unit"}, context);
  const double number = FiniteNumber(value.At("value", context), context + ".value");
  if (number < 0.0 || (!allowZero && number == 0.0))
    throw JsonError(context + ".value must be " + std::string(allowZero ? ">= 0" : "> 0"));
  const std::string unit = RequiredString(object, "unit", context);
  if (unit == "mm") return number;
  if (unit == "cm") return number * 10.0;
  if (unit == "m") return number * 1000.0;
  throw JsonError(context + ".unit must be mm, cm, or m");
}

double PositionMm(const Json& value, const std::string& context) {
  const auto& object = value.AsObject(context);
  RejectUnknownFields(object, {"value","unit"}, context);
  const double number = FiniteNumber(value.At("value", context), context + ".value");
  const std::string unit = RequiredString(object, "unit", context);
  if (unit == "mm") return number;
  if (unit == "cm") return number * 10.0;
  if (unit == "m") return number * 1000.0;
  throw JsonError(context + ".unit must be mm, cm, or m");
}

double EnergyMeV(const Json& value, const std::string& context) {
  const auto& object = value.AsObject(context);
  RejectUnknownFields(object, {"value","unit"}, context);
  const double number = PositiveNumber(value.At("value", context), context + ".value");
  const std::string unit = RequiredString(object, "unit", context);
  if (unit == "MeV") return number;
  if (unit == "GeV") return number * 1000.0;
  throw JsonError(context + ".unit must be MeV or GeV");
}

Json Measure(double value, const std::string& unit);

double AngleRad(const Json& value, const std::string& context) {
  const auto& object=value.AsObject(context);
  RejectUnknownFields(object,{"value","unit"},context);
  const double number=FiniteNumber(value.At("value",context),context+".value");
  const std::string unit=RequiredString(object,"unit",context);
  if(unit=="rad") return number;
  if(unit=="deg") return number*3.14159265358979323846/180.0;
  throw JsonError(context+".unit must be rad or deg");
}

RigidTransform ParseTransform(const Json& value,const std::string& context){
  const auto& object=value.AsObject(context);
  RejectUnknownFields(object,{"translation","rotation_quaternion_local_to_world"},context);
  const auto& translation=object.at("translation").AsObject(context+".translation");
  RejectUnknownFields(translation,{"x","y","z"},context+".translation");
  const auto& rotation=object.at("rotation_quaternion_local_to_world")
      .AsObject(context+".rotation_quaternion_local_to_world");
  RejectUnknownFields(rotation,{"w","x","y","z"},context+".rotation_quaternion_local_to_world");
  RigidTransform result;
  result.translationMm={PositionMm(translation.at("x"),context+".translation.x"),
                        PositionMm(translation.at("y"),context+".translation.y"),
                        PositionMm(translation.at("z"),context+".translation.z")};
  try {
    result.rotation=Quaterniond{
      FiniteNumber(rotation.at("w"),context+".rotation.w"),
      FiniteNumber(rotation.at("x"),context+".rotation.x"),
      FiniteNumber(rotation.at("y"),context+".rotation.y"),
      FiniteNumber(rotation.at("z"),context+".rotation.z")}.Canonical();
  } catch(const std::invalid_argument& error){
    throw JsonError(context+": "+error.what());
  }
  return result;
}

Json TransformToJson(const RigidTransform& transform){
  const Quaterniond q=transform.rotation.Canonical();
  return Json(Json::Object{
    {"rotation_quaternion_local_to_world",Json(Json::Object{
      {"w",q.w},{"x",q.x},{"y",q.y},{"z",q.z}})},
    {"translation",Json(Json::Object{
      {"x",Measure(transform.translationMm.x,"mm")},
      {"y",Measure(transform.translationMm.y,"mm")},
      {"z",Measure(transform.translationMm.z,"mm")}})}
  });
}

// Exact corner containment is sufficient for an oriented box inside the
// axis-aligned world: a convex box is inside a convex world whenever all eight
// vertices are inside.  Stage-3 tests use this independently of Geant4's own
// overlap/navigation diagnostics so malformed inputs fail before transport.
bool BoxInsideWorld(const RigidTransform& transform,double sx,double sy,double sz,
                    double wx,double wy,double wz){
  constexpr double tolerance=1.0e-9;
  for(double x:{-sx/2.0,sx/2.0}) for(double y:{-sy/2.0,sy/2.0})
    for(double z:{-sz/2.0,sz/2.0}){
      const Vector3d p=transform.TransformPoint({x,y,z});
      if(std::abs(p.x)>wx/2.0+tolerance || std::abs(p.y)>wy/2.0+tolerance ||
         std::abs(p.z)>wz/2.0+tolerance) return false;
    }
  return true;
}

// Separating-axis theorem for two oriented boxes.  Touching faces are allowed
// (the shield and an absorber may be contiguous), while positive-volume
// overlap is rejected.  The small epsilon handles nearly parallel axes without
// converting a legitimate shared face into a false overlap.
bool BoxesOverlap(const RigidTransform& a,const Vector3d& aHalf,
                  const RigidTransform& b,const Vector3d& bHalf){
  const auto A=a.Axes(),B=b.Axes();
  const double ah[3]={aHalf.x,aHalf.y,aHalf.z};
  const double bh[3]={bHalf.x,bHalf.y,bHalf.z};
  double R[3][3],AbsR[3][3];
  for(int i=0;i<3;++i) for(int j=0;j<3;++j){
    R[i][j]=Dot(A[i],B[j]);
    AbsR[i][j]=std::abs(R[i][j])+1.0e-12;
  }
  const Vector3d delta=b.translationMm-a.translationMm;
  const double t[3]={Dot(delta,A[0]),Dot(delta,A[1]),Dot(delta,A[2])};
  for(int i=0;i<3;++i){
    const double rb=bh[0]*AbsR[i][0]+bh[1]*AbsR[i][1]+bh[2]*AbsR[i][2];
    if(std::abs(t[i])>=ah[i]+rb-1.0e-9) return false;
  }
  for(int j=0;j<3;++j){
    const double ta=std::abs(t[0]*R[0][j]+t[1]*R[1][j]+t[2]*R[2][j]);
    const double ra=ah[0]*AbsR[0][j]+ah[1]*AbsR[1][j]+ah[2]*AbsR[2][j];
    if(ta>=ra+bh[j]-1.0e-9) return false;
  }
  for(int i=0;i<3;++i) for(int j=0;j<3;++j){
    const double ra=ah[(i+1)%3]*AbsR[(i+2)%3][j]+ah[(i+2)%3]*AbsR[(i+1)%3][j];
    const double rb=bh[(j+1)%3]*AbsR[i][(j+2)%3]+bh[(j+2)%3]*AbsR[i][(j+1)%3];
    const double value=std::abs(t[(i+2)%3]*R[(i+1)%3][j]-t[(i+1)%3]*R[(i+2)%3][j]);
    if(value>=ra+rb-1.0e-9) return false;
  }
  return true;
}

double ArealDensity(const Json& value, const std::string& context) {
  const auto& object = value.AsObject(context);
  RejectUnknownFields(object, {"value","unit"}, context);
  const double number = PositiveNumber(value.At("value", context), context + ".value");
  if (RequiredString(object, "unit", context) != "g/cm2")
    throw JsonError(context + ".unit must be exactly 'g/cm2'");
  return number;
}

Json Measure(double value, const std::string& unit) {
  return Json(Json::Object{{"unit", Json(unit)}, {"value", Json(value)}});
}

std::vector<std::string> StringArray(const Json& value, const std::string& context) {
  std::vector<std::string> result;
  std::set<std::string> seen;
  for (const auto& item : value.AsArray(context)) {
    const std::string text = item.AsString(context + "[]");
    if (text.empty()) throw JsonError(context + " contains an empty string");
    if (!seen.insert(text).second) throw JsonError(context + " contains duplicate '" + text + "'");
    result.push_back(text);
  }
  return result;
}

std::string ResolvePath(const fs::path& directory, const std::string& path) {
  const fs::path input(path);
  return (input.is_absolute() ? input : directory / input).lexically_normal().string();
}

AccuracyRequirement ParseAccuracyRequirement(const Json& value, const std::string& context,
                                             const AccuracyRequirement& base, bool requireAll) {
  const auto& object = value.AsObject(context);
  RejectUnknownFields(object, requireAll
      ? std::vector<std::string>{"target","scope","coverage","systematic_budget"}
      : std::vector<std::string>{"target","scope","coverage","systematic_budget","quantity_overrides","sampling"}, context);
  AccuracyRequirement result = base;
  if (requireAll && (!object.count("target") || !object.count("scope") || !object.count("coverage")))
    throw JsonError(context + " must replace target, scope, and coverage together");
  if (object.count("target")) {
    result.target = PositiveNumber(object.at("target"), context + ".target");
    if (result.target >= 1.0) throw JsonError(context + ".target must be < 1");
  }
  if (object.count("scope")) {
    result.scope = object.at("scope").AsString(context + ".scope");
    if (result.scope != "combined" && result.scope != "statistical")
      throw JsonError(context + ".scope must be combined or statistical");
  }
  if (object.count("coverage")) {
    const auto& coverage = object.at("coverage").AsObject(context + ".coverage");
    RejectUnknownFields(coverage,
      {"standard_k","standard_probability","interval_probability","interval_method"}, context + ".coverage");
    result.standardK = PositiveNumber(coverage.at("standard_k"), context + ".coverage.standard_k");
    result.standardProbability = PositiveNumber(coverage.at("standard_probability"), context + ".coverage.standard_probability");
    result.intervalProbability = PositiveNumber(coverage.at("interval_probability"), context + ".coverage.interval_probability");
    result.intervalMethod = RequiredString(coverage, "interval_method", context + ".coverage");
    if (result.standardProbability >= 1.0 || result.intervalProbability >= 1.0)
      throw JsonError(context + ".coverage probabilities must be < 1");
  }
  if (object.count("systematic_budget")) {
    const auto& budget=object.at("systematic_budget").AsObject(context+".systematic_budget");
    RejectUnknownFields(budget,{"status","relative_standard_uncertainty","known_relative_floor",
                                "components","correlations","pair_policy"},
                        context+".systematic_budget");
    const std::string status=RequiredString(budget,"status",context+".systematic_budget");
    result.systematicComponents.clear();
    result.systematicCorrelations.clear();
    result.systematicPairPolicy.clear();
    if(status=="complete"){
      const bool aggregate=budget.count("relative_standard_uncertainty");
      const bool structured=budget.count("components") || budget.count("correlations") ||
                            budget.count("pair_policy");
      if(aggregate==structured)
        throw JsonError(context+".systematic_budget complete status requires exactly one of the legacy aggregate value or the structured components/correlations form");
      if(aggregate){
        result.systematicRelativeUncertainty=FiniteNumber(
            budget.at("relative_standard_uncertainty"),context+".systematic_budget.relative_standard_uncertainty");
        if(result.systematicRelativeUncertainty<0.0)
          throw JsonError(context+".systematic_budget.relative_standard_uncertainty must be >= 0");
      } else {
        if(!budget.count("components") || !budget.count("correlations") ||
           !budget.count("pair_policy"))
          throw JsonError(context+".systematic_budget structured form requires components, correlations, and pair_policy");
        result.systematicPairPolicy=RequiredString(budget,"pair_policy",context+".systematic_budget");
        if(result.systematicPairPolicy!="all_pairs_explicit")
          throw JsonError(context+".systematic_budget.pair_policy must be all_pairs_explicit");
        const auto& components=budget.at("components").AsArray(context+".systematic_budget.components");
        if(components.empty())
          throw JsonError(context+".systematic_budget.components must not be empty");
        const std::set<std::string> categories={"material_mass","source","numerical","physics_reference"};
        std::map<std::string,std::size_t> indices;
        std::set<std::string> representedCategories;
        for(std::size_t i=0;i<components.size();++i){
          const std::string itemContext=context+".systematic_budget.components["+
              std::to_string(i)+"]";
          const auto& item=components[i].AsObject(itemContext);
          RejectUnknownFields(item,{"id","category","relative_standard_uncertainty","provenance"},itemContext);
          SystematicUncertaintyComponent component;
          component.id=RequiredString(item,"id",itemContext);
          component.category=RequiredString(item,"category",itemContext);
          component.provenance=RequiredString(item,"provenance",itemContext);
          component.relativeStandardUncertainty=FiniteNumber(
              item.at("relative_standard_uncertainty"),itemContext+".relative_standard_uncertainty");
          if(component.relativeStandardUncertainty<0.0)
            throw JsonError(itemContext+".relative_standard_uncertainty must be >= 0");
          if(!categories.count(component.category))
            throw JsonError(itemContext+".category is not a Stage-4 uncertainty category");
          if(!indices.emplace(component.id,i).second)
            throw JsonError(itemContext+" duplicates component id '"+component.id+"'");
          representedCategories.insert(component.category);
          result.systematicComponents.push_back(std::move(component));
        }
        if(representedCategories!=categories)
          throw JsonError(context+".systematic_budget complete structured form must represent material_mass, source, numerical, and physics_reference categories");

        const auto& correlations=budget.at("correlations").AsArray(context+".systematic_budget.correlations");
        std::set<std::pair<std::size_t,std::size_t>> declaredPairs;
        std::vector<std::vector<double>> correlationMatrix(
            components.size(),std::vector<double>(components.size(),0.0));
        for(std::size_t i=0;i<components.size();++i) correlationMatrix[i][i]=1.0;
        for(std::size_t i=0;i<correlations.size();++i){
          const std::string itemContext=context+".systematic_budget.correlations["+
              std::to_string(i)+"]";
          const auto& item=correlations[i].AsObject(itemContext);
          RejectUnknownFields(item,{"component_i","component_j","correlation_coefficient","rationale"},itemContext);
          SystematicCorrelation correlation;
          correlation.firstComponent=RequiredString(item,"component_i",itemContext);
          correlation.secondComponent=RequiredString(item,"component_j",itemContext);
          correlation.rationale=RequiredString(item,"rationale",itemContext);
          correlation.coefficient=FiniteNumber(item.at("correlation_coefficient"),
                                                itemContext+".correlation_coefficient");
          if(correlation.coefficient < -1.0 || correlation.coefficient > 1.0)
            throw JsonError(itemContext+".correlation_coefficient must be in [-1,1]");
          if(!indices.count(correlation.firstComponent) || !indices.count(correlation.secondComponent))
            throw JsonError(itemContext+" refers to an unknown uncertainty component");
          std::size_t a=indices.at(correlation.firstComponent);
          std::size_t b=indices.at(correlation.secondComponent);
          if(a==b) throw JsonError(itemContext+" must name two different components");
          if(a>b) std::swap(a,b);
          if(!declaredPairs.emplace(a,b).second)
            throw JsonError(itemContext+" duplicates an unordered component pair");
          correlationMatrix[a][b]=correlationMatrix[b][a]=correlation.coefficient;
          result.systematicCorrelations.push_back(std::move(correlation));
        }
        const std::size_t expectedPairs=components.size()*(components.size()-1)/2;
        if(declaredPairs.size()!=expectedPairs)
          throw JsonError(context+".systematic_budget all_pairs_explicit requires exactly "+
                          std::to_string(expectedPairs)+" unique correlations");

        // Cholesky with a small numerical tolerance rejects a correlation
        // matrix that cannot represent any real covariance model. A zero
        // pivot is permitted only when the remaining row is also dependent.
        std::vector<std::vector<double>> lower(
            components.size(),std::vector<double>(components.size(),0.0));
        for(std::size_t i=0;i<components.size();++i){
          for(std::size_t j=0;j<=i;++j){
            double value=correlationMatrix[i][j];
            for(std::size_t k=0;k<j;++k) value-=lower[i][k]*lower[j][k];
            if(i==j){
              if(value < -1.0e-12)
                throw JsonError(context+".systematic_budget correlation matrix is not positive semidefinite");
              lower[i][j]=std::sqrt(std::max(0.0,value));
            } else if(lower[j][j]>1.0e-12) {
              lower[i][j]=value/lower[j][j];
            } else if(std::abs(value)>1.0e-10) {
              throw JsonError(context+".systematic_budget correlation matrix is not positive semidefinite");
            }
          }
        }
        double variance=0.0;
        for(std::size_t i=0;i<components.size();++i){
          const double ui=result.systematicComponents[i].relativeStandardUncertainty;
          variance+=ui*ui;
          for(std::size_t j=i+1;j<components.size();++j)
            variance+=2.0*correlationMatrix[i][j]*ui*
                result.systematicComponents[j].relativeStandardUncertainty;
        }
        if(variance < -1.0e-14)
          throw JsonError(context+".systematic_budget produces a negative combined variance");
        result.systematicRelativeUncertainty=std::sqrt(std::max(0.0,variance));
      }
      result.systematicBudgetComplete=true;
    } else if(status=="incomplete"){
      if(budget.count("relative_standard_uncertainty") || budget.count("components") ||
         budget.count("correlations") || budget.count("pair_policy"))
        throw JsonError(context+".systematic_budget incomplete status cannot claim a complete covariance model");
      result.systematicBudgetComplete=false;
    } else {
      throw JsonError(context+".systematic_budget.status must be complete or incomplete");
    }
    if(budget.count("known_relative_floor")){
      result.knownSystematicFloor=FiniteNumber(
          budget.at("known_relative_floor"),context+".systematic_budget.known_relative_floor");
      if(result.knownSystematicFloor<0.0)
        throw JsonError(context+".systematic_budget.known_relative_floor must be >= 0");
      result.hasKnownSystematicFloor=true;
    }
  }
  return result;
}

Json AccuracyToJson(const AccuracyRequirement& requirement) {
  Json::Object coverage;
  coverage["interval_method"] = requirement.intervalMethod;
  coverage["interval_probability"] = requirement.intervalProbability;
  coverage["standard_k"] = requirement.standardK;
  coverage["standard_probability"] = requirement.standardProbability;
  Json::Object budget;
  budget["status"]=requirement.systematicBudgetComplete ? "complete" : "incomplete";
  if(requirement.systematicBudgetComplete){
    budget["combined_relative_standard_uncertainty"]=requirement.systematicRelativeUncertainty;
    if(requirement.systematicComponents.empty()){
      // Legacy case-v1 representation is retained verbatim for reproducible
      // hashes. Stage-4 case-v2 inputs never take this branch.
      budget["relative_standard_uncertainty"]=requirement.systematicRelativeUncertainty;
    } else {
      Json::Array components;
      for(const auto& component:requirement.systematicComponents)
        components.emplace_back(Json::Object{{"category",component.category},{"id",component.id},
          {"provenance",component.provenance},
          {"relative_standard_uncertainty",component.relativeStandardUncertainty}});
      Json::Array correlations;
      for(const auto& correlation:requirement.systematicCorrelations)
        correlations.emplace_back(Json::Object{{"component_i",correlation.firstComponent},
          {"component_j",correlation.secondComponent},
          {"correlation_coefficient",correlation.coefficient},{"rationale",correlation.rationale}});
      budget["components"]=Json(std::move(components));
      budget["correlations"]=Json(std::move(correlations));
      budget["pair_policy"]=requirement.systematicPairPolicy;
    }
  }
  if(requirement.hasKnownSystematicFloor)
    budget["known_relative_floor"]=requirement.knownSystematicFloor;
  return Json(Json::Object{{"coverage", Json(std::move(coverage))},
                           {"scope", requirement.scope},
                           {"systematic_budget",Json(std::move(budget))},
                           {"target", requirement.target}});
}

} // namespace

double CaseDefinition::ShieldTotalThicknessMm() const {
  double result = 0.0;
  for (const auto& layer : shieldLayers) result += layer.thicknessMm;
  return result;
}

double CaseDefinition::ShieldTotalArealDensityGPerCm2() const {
  double result = 0.0;
  for (const auto& layer : shieldLayers) result += layer.arealDensityGPerCm2;
  return result;
}

const AccuracyRequirement& CaseDefinition::AccuracyFor(const std::string& quantity) const {
  const auto found = accuracy.quantityOverrides.find(quantity);
  return found == accuracy.quantityOverrides.end() ? accuracy.global : found->second;
}

std::shared_ptr<CaseDefinition> CaseDefinition::Load(const std::string& path) {
  auto result = std::shared_ptr<CaseDefinition>(new CaseDefinition());
  result->casePath = path;
  const fs::path caseDirectory = fs::absolute(fs::path(path)).parent_path();
  const Json input = Json::ParseFile(path);
  const auto& root = input.AsObject("case");
  RejectUnknownFields(root,
    {"schema_version","case_id","purpose","material_library","material_overrides","geometry",
     "source","quantities","damage_response","accuracy","transport","output"}, "case");
  result->schemaVersion=RequiredString(root,"schema_version","case");
  if(result->schemaVersion!="shieldSim.case.v1" &&
     result->schemaVersion!="shieldSim.case.v2")
    throw JsonError("case.schema_version must be 'shieldSim.case.v1' or 'shieldSim.case.v2'");
  const bool stage3Schema=result->schemaVersion=="shieldSim.case.v2";
  result->caseId = RequiredString(root, "case_id", "case");
  result->purpose = RequiredString(root, "purpose", "case");
  if (result->purpose != "spacecraft_anomaly_environment_effect_assessment" &&
      result->purpose != "preliminary_mission_environment_effect_assessment" &&
      result->purpose != "verification_test")
    throw JsonError("case.purpose is outside the approved model scope");

  const auto& lib = input.At("material_library", "case").AsObject("case.material_library");
  RejectUnknownFields(lib, {"path","library_id","library_version"}, "case.material_library");
  const std::string libraryPath = ResolvePath(caseDirectory, RequiredString(lib, "path", "case.material_library"));
  const Json::Array noOverrides;
  const Json::Array& overrides = root.count("material_overrides")
      ? root.at("material_overrides").AsArray("case.material_overrides") : noOverrides;
  result->materials = std::make_shared<MaterialLibrary>(MaterialLibrary::Load(libraryPath, overrides));
  if (result->materials->Id() != RequiredString(lib, "library_id", "case.material_library") ||
      result->materials->Version() != RequiredString(lib, "library_version", "case.material_library"))
    throw JsonError("case.material_library ID/version does not match the loaded file");

  const auto& geometry = input.At("geometry", "case").AsObject("case.geometry");
  RejectUnknownFields(geometry, {"coordinate_frame","world","source_plane","shield_stack","absorbers"}, "case.geometry");
  const auto& frame = input.At("geometry", "case").At("coordinate_frame", "case.geometry")
                            .AsObject("case.geometry.coordinate_frame");
  RejectUnknownFields(frame, stage3Schema
      ? std::vector<std::string>{"name","handedness","origin","local_surface_normal"}
      : std::vector<std::string>{"name","slab_normal","origin"},
      "case.geometry.coordinate_frame");
  result->coordinateFrame = RequiredString(frame, "name", "case.geometry.coordinate_frame");
  if(stage3Schema){
    if(RequiredString(frame,"handedness","case.geometry.coordinate_frame")!="right")
      throw JsonError("case-v2 geometry requires coordinate_frame.handedness='right'");
    if(RequiredString(frame,"local_surface_normal","case.geometry.coordinate_frame")!="+z")
      throw JsonError("case-v2 slab providers require local_surface_normal='+z'");
  } else if (RequiredString(frame, "slab_normal", "case.geometry.coordinate_frame") != "+z") {
    throw JsonError("Stage-1 slab geometry requires coordinate_frame.slab_normal='+z'");
  }
  if (RequiredString(frame, "origin", "case.geometry.coordinate_frame") != "world_center")
    throw JsonError("case geometry requires coordinate_frame.origin='world_center'");

  const auto& world = input.At("geometry", "case").At("world", "case.geometry").AsObject("case.geometry.world");
  RejectUnknownFields(world, {"size_x","size_y","size_z","center_x","center_y","center_z"}, "case.geometry.world");
  result->worldSizeXmm = LengthMm(world.at("size_x"), "case.geometry.world.size_x");
  result->worldSizeYmm = LengthMm(world.at("size_y"), "case.geometry.world.size_y");
  result->worldSizeZmm = LengthMm(world.at("size_z"), "case.geometry.world.size_z");
  result->worldCenterXmm = PositionMm(world.at("center_x"), "case.geometry.world.center_x");
  result->worldCenterYmm = PositionMm(world.at("center_y"), "case.geometry.world.center_y");
  result->worldCenterZmm = PositionMm(world.at("center_z"), "case.geometry.world.center_z");
  if (result->worldCenterXmm != 0.0 || result->worldCenterYmm != 0.0 || result->worldCenterZmm != 0.0)
    throw JsonError("versioned-case Geant4 world center must be explicitly [0,0,0] mm");

  const auto& sourcePlane = input.At("geometry", "case").At("source_plane", "case.geometry")
                                  .AsObject("case.geometry.source_plane");
  RejectUnknownFields(sourcePlane, stage3Schema
      ? std::vector<std::string>{"size_x","size_y","transform"}
      : std::vector<std::string>{"size_x","size_y","center_x","center_y","z"},
      "case.geometry.source_plane");
  result->sourceSizeXmm = LengthMm(sourcePlane.at("size_x"), "case.geometry.source_plane.size_x");
  result->sourceSizeYmm = LengthMm(sourcePlane.at("size_y"), "case.geometry.source_plane.size_y");
  if(stage3Schema){
    result->sourceTransform=ParseTransform(sourcePlane.at("transform"),"case.geometry.source_plane.transform");
    result->sourceCenterXmm=result->sourceTransform.translationMm.x;
    result->sourceCenterYmm=result->sourceTransform.translationMm.y;
    result->sourceZmm=result->sourceTransform.translationMm.z;
  } else {
    result->sourceCenterXmm = PositionMm(sourcePlane.at("center_x"), "case.geometry.source_plane.center_x");
    result->sourceCenterYmm = PositionMm(sourcePlane.at("center_y"), "case.geometry.source_plane.center_y");
    result->sourceZmm = PositionMm(sourcePlane.at("z"), "case.geometry.source_plane.z");
    result->sourceTransform.translationMm={result->sourceCenterXmm,result->sourceCenterYmm,result->sourceZmm};
  }

  const auto& stack = input.At("geometry", "case").At("shield_stack", "case.geometry")
                            .AsObject("case.geometry.shield_stack");
  RejectUnknownFields(stack, stage3Schema
      ? std::vector<std::string>{"size_x","size_y","transform","layers"}
      : std::vector<std::string>{"size_x","size_y","center_x","center_y","front_z","layers"},
      "case.geometry.shield_stack");
  result->shieldSizeXmm = LengthMm(stack.at("size_x"), "case.geometry.shield_stack.size_x");
  result->shieldSizeYmm = LengthMm(stack.at("size_y"), "case.geometry.shield_stack.size_y");
  if(stage3Schema){
    result->shieldTransform=ParseTransform(stack.at("transform"),"case.geometry.shield_stack.transform");
    result->shieldCenterXmm=result->shieldTransform.translationMm.x;
    result->shieldCenterYmm=result->shieldTransform.translationMm.y;
    result->shieldFrontZmm=result->shieldTransform.translationMm.z;
  } else {
    result->shieldCenterXmm = PositionMm(stack.at("center_x"), "case.geometry.shield_stack.center_x");
    result->shieldCenterYmm = PositionMm(stack.at("center_y"), "case.geometry.shield_stack.center_y");
    result->shieldFrontZmm = PositionMm(stack.at("front_z"), "case.geometry.shield_stack.front_z");
    result->shieldTransform.translationMm={result->shieldCenterXmm,result->shieldCenterYmm,result->shieldFrontZmm};
  }
  const auto& layers = stack.at("layers").AsArray("case.geometry.shield_stack.layers");
  if (layers.empty() || layers.size() > 5)
    throw JsonError("case.geometry.shield_stack.layers must contain one through five layers");
  std::set<std::string> layerIds;
  for (std::size_t i = 0; i < layers.size(); ++i) {
    const std::string context = "case.geometry.shield_stack.layers[" + std::to_string(i) + ']';
    const auto& layer = layers[i].AsObject(context);
    RejectUnknownFields(layer, {"layer_id","material","thickness"}, context);
    ShieldLayerDefinition parsed;
    parsed.id = RequiredString(layer, "layer_id", context);
    if (!layerIds.insert(parsed.id).second) throw JsonError(context + ": duplicate layer_id '" + parsed.id + "'");
    const auto& material = result->materials->Resolve(RequiredString(layer, "material", context), "shield");
    parsed.materialId = material.id;
    parsed.materialHash = material.contentHash;
    parsed.densityGPerCm3 = material.densityGPerCm3;
    const auto& thickness = layer.at("thickness").AsObject(context + ".thickness");
    RejectUnknownFields(thickness, {"type","value","unit"}, context + ".thickness");
    parsed.inputThicknessType = RequiredString(thickness, "type", context + ".thickness");
    Json thicknessMeasure(Json::Object{{"value", thickness.at("value")}, {"unit", thickness.at("unit")}});
    if (parsed.inputThicknessType == "geometric") {
      parsed.thicknessMm = LengthMm(thicknessMeasure, context + ".thickness");
      parsed.arealDensityGPerCm2 = material.densityGPerCm3 * parsed.thicknessMm / 10.0;
    } else if (parsed.inputThicknessType == "areal_density") {
      parsed.arealDensityGPerCm2 = ArealDensity(thicknessMeasure, context + ".thickness");
      parsed.thicknessMm = 10.0 * parsed.arealDensityGPerCm2 / material.densityGPerCm3;
    } else {
      throw JsonError(context + ".thickness.type must be geometric or areal_density");
    }
    if (!std::isfinite(parsed.thicknessMm) || parsed.thicknessMm <= 0.0 || parsed.thicknessMm > 100.0 + 1e-12)
      throw JsonError(context + " derived geometric thickness must be in (0, 100] mm");
    result->shieldLayers.push_back(parsed);
  }

  const auto& absorberValues = input.At("geometry", "case").At("absorbers", "case.geometry")
                                     .AsArray("case.geometry.absorbers");
  if (absorberValues.empty()) throw JsonError("case.geometry.absorbers must not be empty");
  std::set<std::string> absorberIds;
  for (std::size_t i = 0; i < absorberValues.size(); ++i) {
    const std::string context = "case.geometry.absorbers[" + std::to_string(i) + ']';
    const auto& object = absorberValues[i].AsObject(context);
    RejectUnknownFields(object, {"absorber_id","material","geometry"}, context);
    SlabAbsorberDefinition parsed;
    parsed.id = RequiredString(object, "absorber_id", context);
    if (!absorberIds.insert(parsed.id).second) throw JsonError(context + ": duplicate absorber_id '" + parsed.id + "'");
    const auto& material = result->materials->Resolve(RequiredString(object, "material", context), "absorber");
    parsed.materialId = material.id;
    parsed.materialHash = material.contentHash;
    parsed.densityGPerCm3 = material.densityGPerCm3;
    const auto& shape = object.at("geometry").AsObject(context + ".geometry");
    RejectUnknownFields(shape, stage3Schema
        ? std::vector<std::string>{"type","size_x","size_y","thickness","transform"}
        : std::vector<std::string>{"type","size_x","size_y","center_x","center_y","thickness","front_z"},
        context + ".geometry");
    if (RequiredString(shape, "type", context + ".geometry") != "slab")
      throw JsonError(context + ".geometry.type must be 'slab'; future providers retain a later validation gate");
    parsed.sizeXmm = LengthMm(shape.at("size_x"), context + ".geometry.size_x");
    parsed.sizeYmm = LengthMm(shape.at("size_y"), context + ".geometry.size_y");
    parsed.thicknessMm = LengthMm(shape.at("thickness"), context + ".geometry.thickness");
    if(stage3Schema){
      parsed.transform=ParseTransform(shape.at("transform"),context+".geometry.transform");
      parsed.centerXmm=parsed.transform.translationMm.x;
      parsed.centerYmm=parsed.transform.translationMm.y;
      parsed.frontZmm=parsed.transform.translationMm.z;
    } else {
      parsed.centerXmm = PositionMm(shape.at("center_x"), context + ".geometry.center_x");
      parsed.centerYmm = PositionMm(shape.at("center_y"), context + ".geometry.center_y");
      parsed.frontZmm = PositionMm(shape.at("front_z"), context + ".geometry.front_z");
      parsed.transform.translationMm={parsed.centerXmm,parsed.centerYmm,parsed.frontZmm};
    }
    Json::Object hashInput{{"absorber_id", parsed.id},
      {"transform",TransformToJson(parsed.transform)},
      {"material_hash", parsed.materialHash}, {"size_x_mm", parsed.sizeXmm},
      {"size_y_mm", parsed.sizeYmm}, {"thickness_mm", parsed.thicknessMm}, {"type", "slab"}};
    parsed.geometryHash = Sha256(Json(std::move(hashInput)).DumpCanonical());
    result->absorbers.push_back(parsed);
  }

  const auto& source = input.At("source", "case").AsObject("case.source");
  RejectUnknownFields(source, {"mode","spectrum","energy_limits"}, "case.source");
  result->sourceMode = RequiredString(source, "mode", "case.source");
  if (result->sourceMode != "beam" && result->sourceMode != "isotropic" &&
      !(stage3Schema && result->sourceMode=="directional_plane"))
    throw JsonError("case.source.mode must be beam, isotropic, or case-v2 directional_plane");
  const auto& spectrum = source.at("spectrum").AsObject("case.source.spectrum");
  result->sourceFormat = RequiredString(spectrum, "format", "case.source.spectrum");
  if (result->sourceFormat == "legacy_three_column") {
    RejectUnknownFields(spectrum, {"format","path","metadata_path"}, "case.source.spectrum");
    result->spectrum = SpectrumDataset::LoadLegacyThreeColumn(
        ResolvePath(caseDirectory, RequiredString(spectrum, "path", "case.source.spectrum")),
        ResolvePath(caseDirectory, RequiredString(spectrum, "metadata_path", "case.source.spectrum")));
    if (result->spectrum.angularModel != "declared_by_case" &&
        result->spectrum.angularModel != result->sourceMode)
      throw JsonError("case.source.mode conflicts with spectrum metadata.angular_model");
    result->sourceComponentId = "proton_and_alpha_shared_grid";
  } else if (result->sourceFormat == "ams02_hapi_component") {
    // This is deliberately a narrow validation adapter, not a generic license
    // to turn any ingested observation into a boundary condition.  The pinned
    // AMS product has one unambiguous transport species: its proton component.
    // Inclusive helium has no isotope mixture in the product and therefore
    // must not be silently represented as Geant4 alpha particles.
    RejectUnknownFields(spectrum,
        {"format","info_path","data_path","interpretation_path","component_id"},
        "case.source.spectrum");
    if (result->purpose != "verification_test")
      throw JsonError("ams02_hapi_component is admitted only for verification_test cases, not mission environments");
    const std::string componentId = RequiredString(
        spectrum, "component_id", "case.source.spectrum");
    SpectrumDataset published = AMS02HapiAdapter::Load(
        ResolvePath(caseDirectory, RequiredString(spectrum, "info_path", "case.source.spectrum")),
        ResolvePath(caseDirectory, RequiredString(spectrum, "data_path", "case.source.spectrum")),
        ResolvePath(caseDirectory, RequiredString(spectrum, "interpretation_path", "case.source.spectrum")));
    const auto selected = std::find_if(published.components.begin(), published.components.end(),
        [&](const SpectrumComponent& component) { return component.componentId == componentId; });
    if (selected == published.components.end())
      throw JsonError("case.source.spectrum.component_id is absent from the AMS dataset");
    if (selected->componentId != "proton" || selected->speciesId != "proton" ||
        selected->geant4Particle != "proton" || selected->atomicNumber != 1 ||
        selected->massNumber != 1)
      throw JsonError("Stage-2 AMS transport validation admits only the explicit A=1 proton component");
    if (selected->energyVariable != "kinetic_energy_per_nucleon" ||
        selected->energyUnit != "MeV/n")
      throw JsonError("AMS proton transport requires kinetic_energy_per_nucleon in MeV/n");
    if (selected->valueType != "differential_intensity" ||
        selected->valueUnit != "1/(cm2 sr s (MeV/n))")
      throw JsonError("AMS proton transport requires differential intensity in 1/(cm2 sr s (MeV/n))");
    if (selected->angularRepresentation != "isotropic_intensity_at_instrument" ||
        result->sourceMode != "isotropic")
      throw JsonError("AMS isotropic intensity validation requires case.source.mode='isotropic'");
    const std::size_t count = selected->value.size();
    if (count == 0 || selected->binLower.size() != count ||
        selected->binUpper.size() != count)
      throw JsonError("AMS proton component has inconsistent or empty bin arrays");
    if (!selected->valid.empty() && selected->valid.size() != count)
      throw JsonError("AMS proton component validity-mask size differs from its bin count");
    const SpectrumComponent transportedComponent = *selected;

    // For A=1 only, kinetic energy per nucleon is numerically equal to total
    // proton kinetic energy.  Preserve exact published edges and use a
    // piecewise-constant differential intensity inside each bin.  No analogous
    // conversion is made for mixed-isotope components.
    published.energyMeV.clear();
    published.binWidthMeV.clear();
    published.protonValue.clear();
    published.alphaValue.clear();
    published.transportBinLowerMeV.clear();
    published.transportBinUpperMeV.clear();
    for (std::size_t i = 0; i < count; ++i) {
      const double lower = selected->binLower[i];
      const double upper = selected->binUpper[i];
      const double value = selected->value[i];
      if (!std::isfinite(lower) || !std::isfinite(upper) || lower <= 0.0 || upper <= lower ||
          !std::isfinite(value) || value < 0.0)
        throw JsonError("AMS proton component contains an invalid transport bin");
      if (!selected->valid.empty() && !selected->valid[i])
        throw JsonError("AMS proton component contains a bin not approved by its validity mask");
      published.transportBinLowerMeV.push_back(lower);
      published.transportBinUpperMeV.push_back(upper);
      published.energyMeV.push_back(0.5 * (lower + upper));
      published.binWidthMeV.push_back(upper - lower);
      published.protonValue.push_back(value);
      published.alphaValue.push_back(0.0);
    }
    published.transportComponentId = componentId;
    // Keep only the transported component in the live representation.  The
    // canonical metadata still retains and hashes the complete source record.
    published.components = {transportedComponent};
    result->spectrum = std::move(published);
    result->sourceComponentId = componentId;
  } else if(result->sourceFormat=="monoenergetic_response_cell") {
    if(!stage3Schema)
      throw JsonError("monoenergetic_response_cell requires shieldSim.case.v2");
    if(result->purpose!="verification_test")
      throw JsonError("monoenergetic response cells are verification infrastructure, not mission environments");
    if(result->sourceMode!="directional_plane")
      throw JsonError("monoenergetic_response_cell requires source.mode='directional_plane'");
    RejectUnknownFields(spectrum,{"format","incident_id","particle","energy","direction"},
                        "case.source.spectrum");
    auto& incident=result->responseIncident;
    incident.enabled=true;
    incident.incidentId=RequiredString(spectrum,"incident_id","case.source.spectrum");
    const auto& particle=spectrum.at("particle").AsObject("case.source.spectrum.particle");
    RejectUnknownFields(particle,{"name","Z","A"},"case.source.spectrum.particle");
    incident.particleName=RequiredString(particle,"name","case.source.spectrum.particle");
    incident.atomicNumber=static_cast<int>(NonNegativeInteger(particle.at("Z"),"case.source.spectrum.particle.Z"));
    incident.massNumber=static_cast<int>(NonNegativeInteger(particle.at("A"),"case.source.spectrum.particle.A"));
    if(incident.massNumber<incident.atomicNumber)
      throw JsonError("case.source.spectrum.particle requires A >= Z");
    if((incident.particleName=="proton" && (incident.atomicNumber!=1 || incident.massNumber!=1)) ||
       (incident.particleName=="alpha" && (incident.atomicNumber!=2 || incident.massNumber!=4)) ||
       (incident.particleName=="neutron" && (incident.atomicNumber!=0 || incident.massNumber!=1)) ||
       (incident.particleName=="gamma" && (incident.atomicNumber!=0 || incident.massNumber!=0)) ||
       (incident.particleName!="proton" && incident.particleName!="alpha" &&
        incident.particleName!="neutron" && incident.particleName!="gamma" &&
        incident.particleName!="ion") ||
       (incident.particleName=="ion" && (incident.atomicNumber==0 || incident.massNumber==0)))
      throw JsonError("response particle name/Z/A must be proton/1/1, alpha/2/4, neutron/0/1, gamma/0/0, or ion with explicit positive Z/A");
    const auto& energy=spectrum.at("energy").AsObject("case.source.spectrum.energy");
    RejectUnknownFields(energy,{"coordinate","value","unit"},"case.source.spectrum.energy");
    incident.energyCoordinate=RequiredString(energy,"coordinate","case.source.spectrum.energy");
    incident.energyCoordinateValue=PositiveNumber(energy.at("value"),"case.source.spectrum.energy.value");
    incident.energyCoordinateUnit=RequiredString(energy,"unit","case.source.spectrum.energy");
    double multiplier=1.0;
    if(incident.energyCoordinate=="total_kinetic_energy"){
      if(incident.energyCoordinateUnit=="MeV") multiplier=1.0;
      else if(incident.energyCoordinateUnit=="GeV") multiplier=1000.0;
      else throw JsonError("total_kinetic_energy response coordinate requires MeV or GeV");
    } else if(incident.energyCoordinate=="kinetic_energy_per_nucleon"){
      if(incident.particleName=="gamma" || incident.particleName=="neutron")
        throw JsonError("gamma and neutron response cells require total_kinetic_energy");
      if(incident.energyCoordinateUnit=="MeV/n") multiplier=incident.massNumber;
      else if(incident.energyCoordinateUnit=="GeV/n") multiplier=1000.0*incident.massNumber;
      else throw JsonError("kinetic_energy_per_nucleon response coordinate requires MeV/n or GeV/n");
    } else {
      throw JsonError("response energy coordinate must be total_kinetic_energy or kinetic_energy_per_nucleon");
    }
    incident.totalKineticEnergyMeV=incident.energyCoordinateValue*multiplier;
    const auto& direction=spectrum.at("direction").AsObject("case.source.spectrum.direction");
    RejectUnknownFields(direction,{"mu","phi"},"case.source.spectrum.direction");
    incident.mu=FiniteNumber(direction.at("mu"),"case.source.spectrum.direction.mu");
    if(incident.mu<=0.0 || incident.mu>1.0)
      throw JsonError("case.source.spectrum.direction.mu must be in (0,1]");
    incident.phiRad=AngleRad(direction.at("phi"),"case.source.spectrum.direction.phi");
    incident.phiRad=std::fmod(incident.phiRad,2.0*3.14159265358979323846);
    if(incident.phiRad<0.0) incident.phiRad+=2.0*3.14159265358979323846;
    result->sourceComponentId=incident.incidentId;
    result->spectrum.id="stage3_monoenergetic_response_cell";
    result->spectrum.scienceStatus="verification_only_not_a_physical_environment";
    result->spectrum.valueType="one_primary_per_history";
    result->spectrum.valueUnit="1";
    result->spectrum.contentHash=Sha256(spectrum.at("particle").DumpCanonical()+
        spectrum.at("energy").DumpCanonical()+spectrum.at("direction").DumpCanonical());
    result->protonMinMeV=result->protonMaxMeV=incident.particleName=="proton"
        ? incident.totalKineticEnergyMeV : 1.0;
    result->alphaMinMeV=result->alphaMaxMeV=incident.particleName=="alpha"
        ? incident.totalKineticEnergyMeV : 1.0;
  } else if(result->sourceFormat=="discrete_response_nodes") {
    if(!stage3Schema || result->purpose!="verification_test" ||
       result->sourceMode!="directional_plane")
      throw JsonError("discrete_response_nodes requires a case-v2 verification_test directional_plane source");
    RejectUnknownFields(spectrum,{"format","nodes"},"case.source.spectrum");
    const auto& nodes=spectrum.at("nodes").AsArray("case.source.spectrum.nodes");
    if(nodes.empty()) throw JsonError("case.source.spectrum.nodes must not be empty");
    std::set<std::string> ids;
    std::set<std::tuple<std::string,int,int,double,double,double>> coordinates;
    double weightSum=0.0;
    for(std::size_t i=0;i<nodes.size();++i){
      const std::string context="case.source.spectrum.nodes["+std::to_string(i)+"]";
      const auto& node=nodes[i].AsObject(context);
      RejectUnknownFields(node,{"incident_id","particle","energy","direction","weight"},context);
      WeightedResponseIncident weighted;
      auto& incident=weighted.incident;
      incident.enabled=true;
      incident.incidentId=RequiredString(node,"incident_id",context);
      if(!ids.insert(incident.incidentId).second)
        throw JsonError(context+" has a duplicate incident_id");
      weighted.weight=PositiveNumber(node.at("weight"),context+".weight");
      weightSum+=weighted.weight;
      const auto& particle=node.at("particle").AsObject(context+".particle");
      RejectUnknownFields(particle,{"name","Z","A"},context+".particle");
      incident.particleName=RequiredString(particle,"name",context+".particle");
      incident.atomicNumber=static_cast<int>(NonNegativeInteger(particle.at("Z"),context+".particle.Z"));
      incident.massNumber=static_cast<int>(NonNegativeInteger(particle.at("A"),context+".particle.A"));
      if(incident.massNumber<incident.atomicNumber ||
         (incident.particleName=="proton" && (incident.atomicNumber!=1 || incident.massNumber!=1)) ||
         (incident.particleName=="alpha" && (incident.atomicNumber!=2 || incident.massNumber!=4)) ||
         (incident.particleName=="neutron" && (incident.atomicNumber!=0 || incident.massNumber!=1)) ||
         (incident.particleName=="gamma" && (incident.atomicNumber!=0 || incident.massNumber!=0)) ||
         (incident.particleName!="proton" && incident.particleName!="alpha" &&
          incident.particleName!="neutron" && incident.particleName!="gamma" &&
          incident.particleName!="ion") ||
         (incident.particleName=="ion" && (incident.atomicNumber==0 || incident.massNumber==0)))
        throw JsonError(context+" has an inconsistent particle name/Z/A");
      const auto& energy=node.at("energy").AsObject(context+".energy");
      RejectUnknownFields(energy,{"coordinate","value","unit"},context+".energy");
      incident.energyCoordinate=RequiredString(energy,"coordinate",context+".energy");
      incident.energyCoordinateValue=PositiveNumber(energy.at("value"),context+".energy.value");
      incident.energyCoordinateUnit=RequiredString(energy,"unit",context+".energy");
      if(incident.energyCoordinate=="total_kinetic_energy"){
        if(incident.energyCoordinateUnit=="MeV") incident.totalKineticEnergyMeV=incident.energyCoordinateValue;
        else if(incident.energyCoordinateUnit=="GeV") incident.totalKineticEnergyMeV=1000.0*incident.energyCoordinateValue;
        else throw JsonError(context+" total energy requires MeV or GeV");
      } else if(incident.energyCoordinate=="kinetic_energy_per_nucleon"){
        if(incident.particleName=="gamma" || incident.particleName=="neutron")
          throw JsonError(context+" gamma and neutron nodes require total_kinetic_energy");
        if(incident.energyCoordinateUnit=="MeV/n")
          incident.totalKineticEnergyMeV=incident.massNumber*incident.energyCoordinateValue;
        else if(incident.energyCoordinateUnit=="GeV/n")
          incident.totalKineticEnergyMeV=1000.0*incident.massNumber*incident.energyCoordinateValue;
        else throw JsonError(context+" per-nucleon energy requires MeV/n or GeV/n");
      } else throw JsonError(context+" has unsupported energy coordinate");
      const auto& direction=node.at("direction").AsObject(context+".direction");
      RejectUnknownFields(direction,{"mu","phi"},context+".direction");
      incident.mu=FiniteNumber(direction.at("mu"),context+".direction.mu");
      if(incident.mu<=0.0 || incident.mu>1.0)
        throw JsonError(context+".direction.mu must be in (0,1]");
      incident.phiRad=AngleRad(direction.at("phi"),context+".direction.phi");
      incident.phiRad=std::fmod(incident.phiRad,2.0*3.14159265358979323846);
      if(incident.phiRad<0.0) incident.phiRad+=2.0*3.14159265358979323846;
      if(!coordinates.emplace(incident.particleName,incident.atomicNumber,incident.massNumber,
                              incident.totalKineticEnergyMeV,incident.mu,incident.phiRad).second)
        throw JsonError(context+" duplicates a species/energy/direction coordinate");
      result->directIncidents.push_back(std::move(weighted));
    }
    if(std::abs(weightSum-1.0)>1.0e-12)
      throw JsonError("case.source.spectrum.nodes weights must sum to one within 1e-12");
    result->sourceComponentId="declared_discrete_node_mixture";
    result->spectrum.id="stage3_discrete_response_nodes";
    result->spectrum.scienceStatus="verification_only_not_a_physical_environment";
    result->spectrum.valueType="probability_mass_per_node";
    result->spectrum.valueUnit="1";
    result->spectrum.contentHash=Sha256(spectrum.at("nodes").DumpCanonical());
    result->protonMinMeV=result->alphaMinMeV=std::numeric_limits<double>::infinity();
    result->protonMaxMeV=result->alphaMaxMeV=0.0;
    for(const auto& weighted:result->directIncidents){
      const auto& incident=weighted.incident;
      if(incident.particleName=="proton"){
        result->protonMinMeV=std::min(result->protonMinMeV,incident.totalKineticEnergyMeV);
        result->protonMaxMeV=std::max(result->protonMaxMeV,incident.totalKineticEnergyMeV);
      } else if(incident.particleName=="alpha"){
        result->alphaMinMeV=std::min(result->alphaMinMeV,incident.totalKineticEnergyMeV);
        result->alphaMaxMeV=std::max(result->alphaMaxMeV,incident.totalKineticEnergyMeV);
      }
    }
    if(!std::isfinite(result->protonMinMeV))
      result->protonMinMeV=result->protonMaxMeV=1.0;
    if(!std::isfinite(result->alphaMinMeV))
      result->alphaMinMeV=result->alphaMaxMeV=1.0;
  } else {
    throw JsonError("case.source.spectrum.format must be legacy_three_column, ams02_hapi_component, monoenergetic_response_cell, or discrete_response_nodes");
  }
  if(!result->responseIncident.enabled && result->directIncidents.empty()){
    if(!source.count("energy_limits"))
      throw JsonError("case.source.energy_limits is required for spectrum-sampled sources");
    const auto& limits = source.at("energy_limits").AsObject("case.source.energy_limits");
    RejectUnknownFields(limits, {"proton","alpha"}, "case.source.energy_limits");
    const auto& proton = limits.at("proton").AsObject("case.source.energy_limits.proton");
    const auto& alpha = limits.at("alpha").AsObject("case.source.energy_limits.alpha");
    RejectUnknownFields(proton, {"min","max"}, "case.source.energy_limits.proton");
    RejectUnknownFields(alpha, {"min","max"}, "case.source.energy_limits.alpha");
    result->protonMinMeV = EnergyMeV(proton.at("min"), "case.source.energy_limits.proton.min");
    result->protonMaxMeV = EnergyMeV(proton.at("max"), "case.source.energy_limits.proton.max");
    result->alphaMinMeV = EnergyMeV(alpha.at("min"), "case.source.energy_limits.alpha.min");
    result->alphaMaxMeV = EnergyMeV(alpha.at("max"), "case.source.energy_limits.alpha.max");
    if (result->protonMinMeV > result->protonMaxMeV || result->alphaMinMeV > result->alphaMaxMeV)
      throw JsonError("case.source energy min must be <= max; ranges are never reordered");
  } else if(source.count("energy_limits")) {
    throw JsonError("Stage-3 declared response sources must not also provide energy_limits");
  }
  if (!result->spectrum.transportBinLowerMeV.empty()) {
    bool overlaps = false;
    for (std::size_t i = 0; i < result->spectrum.transportBinLowerMeV.size(); ++i)
      overlaps = overlaps || (std::max(result->protonMinMeV,
                                       result->spectrum.transportBinLowerMeV[i]) <
                              std::min(result->protonMaxMeV,
                                       result->spectrum.transportBinUpperMeV[i]));
    if (!overlaps)
      throw JsonError("case.source proton energy limits do not overlap the selected published component");
  }

  result->quantities = StringArray(input.At("quantities", "case"), "case.quantities");
  if (result->quantities.empty()) throw JsonError("case.quantities must not be empty");
  for (const auto& quantity : result->quantities) {
    const bool allowed=stage3Schema
        ? quantity=="transport_raw" || quantity=="TID" || quantity=="DDD" ||
          quantity=="n_eq"
        : quantity=="TID" || quantity=="H100/10";
    if (!allowed)
      throw JsonError("quantity '" + quantity + "' is not enabled in this case schema; case-v2 permits transport_raw, Stage-4 finite-slab TID, and Stage-5 silicon DDD/neutron-only n_eq while LET, H100/10, crew-health, and device/organ claims retain later gates");
  }

  const bool requestsDdd=std::find(result->quantities.begin(),result->quantities.end(),"DDD")!=
      result->quantities.end();
  const bool requestsNeq=std::find(result->quantities.begin(),result->quantities.end(),"n_eq")!=
      result->quantities.end();
  const bool requestsDamage=requestsDdd || requestsNeq;
  if(requestsDamage){
    if(!stage3Schema)
      throw JsonError("DDD and n_eq require shieldSim.case.v2");
    if(!root.count("damage_response"))
      throw JsonError("Stage-5 DDD/n_eq requires case.damage_response");
    const auto& damage=root.at("damage_response").AsObject("case.damage_response");
    RejectUnknownFields(damage,{"path","library_id","library_version"},"case.damage_response");
    const std::string damagePath=ResolvePath(caseDirectory,
        RequiredString(damage,"path","case.damage_response"));
    result->damageResponse=std::make_shared<DamageResponseLibrary>(
        DamageResponseLibrary::Load(damagePath));
    if(result->damageResponse->LibraryId()!=
         RequiredString(damage,"library_id","case.damage_response") ||
       result->damageResponse->LibraryVersion()!=
         RequiredString(damage,"library_version","case.damage_response"))
      throw JsonError("case.damage_response ID/version does not match the loaded manifest");
    for(const auto& absorber:result->absorbers)
      if(absorber.materialId!=result->damageResponse->TargetMaterialId())
        throw JsonError("Stage-5 DDD/n_eq absorber '"+absorber.id+"' uses material_id '"+
            absorber.materialId+"'; the approved response target is '"+
            result->damageResponse->TargetMaterialId()+"'");
    const auto validateDamageIncident=[&](const ResponseIncidentDefinition& incident,
                                          const std::string& context){
      if(!result->damageResponse->SupportsSpecies(incident.particleName))
        throw JsonError(context+" particle '"+incident.particleName+
            "' is outside the Stage-5 proton/alpha/neutron incident support matrix");
    };
    if(result->responseIncident.enabled)
      validateDamageIncident(result->responseIncident,"case.source.spectrum");
    for(std::size_t i=0;i<result->directIncidents.size();++i)
      validateDamageIncident(result->directIncidents[i].incident,
          "case.source.spectrum.nodes["+std::to_string(i)+"]");
  } else if(root.count("damage_response")) {
    throw JsonError("case.damage_response is permitted only when DDD or n_eq is requested");
  }

  // Stage 4 validates absorbed dose only for proton and He-4 irradiation.
  // Gamma response nodes remain useful for the independent XCOM attenuation
  // check, and explicit heavy-ion nodes remain part of the Stage-3 raw-score
  // infrastructure, but neither may silently acquire a validated TID claim.
  // Spectrum-sampled sources already contain only proton and alpha columns.
  const bool requestsTid=std::find(result->quantities.begin(),result->quantities.end(),"TID")!=
      result->quantities.end();
  const auto tidSpeciesIsValidated=[](const ResponseIncidentDefinition& incident){
    return incident.particleName=="proton" || incident.particleName=="alpha";
  };
  if(stage3Schema && requestsTid && result->responseIncident.enabled &&
     !tidSpeciesIsValidated(result->responseIncident))
    throw JsonError("case-v2 TID is validated only for proton and alpha response incidents; gamma and explicit-ion nodes are transport_raw verification inputs");
  if(stage3Schema && requestsTid){
    for(const auto& weighted:result->directIncidents)
      if(!tidSpeciesIsValidated(weighted.incident))
        throw JsonError("case-v2 TID discrete_response_nodes may contain only proton and alpha incidents; gamma and explicit-ion nodes are transport_raw verification inputs");

    // The Stage-4 material reference matrix covers these three absorber
    // identities. User overrides may change their declared physical
    // properties (and are carried into the material hash and uncertainty
    // record), but a different material identity needs a later validation
    // gate. Shield layers are intentionally not restricted by this absorber-
    // response boundary.
    const std::set<std::string> validatedAbsorberMaterials={
      "Al_natural_G4_11_4","Si_natural_G4_11_4","Water_G4_11_4"};
    for(const auto& absorber:result->absorbers)
      if(!validatedAbsorberMaterials.count(absorber.materialId))
        throw JsonError("case-v2 TID absorber '"+absorber.id+"' uses material_id '"+
            absorber.materialId+"'; Stage 4 validates only Al_natural_G4_11_4, Si_natural_G4_11_4, and Water_G4_11_4 absorbers");
  }

  const auto& accuracy = input.At("accuracy", "case").AsObject("case.accuracy");
  RejectUnknownFields(accuracy, {"target","scope","coverage","systematic_budget",
                                 "quantity_overrides","sampling"}, "case.accuracy");
  result->accuracy.targetWasDefaulted = !accuracy.count("target");
  result->accuracy.global = ParseAccuracyRequirement(input.At("accuracy", "case"), "case.accuracy",
                                                     AccuracyRequirement{}, false);
  if (accuracy.count("quantity_overrides")) {
    const auto& overridesObject = accuracy.at("quantity_overrides").AsObject("case.accuracy.quantity_overrides");
    for (const auto& item : overridesObject) {
      if (std::find(result->quantities.begin(), result->quantities.end(), item.first) == result->quantities.end())
        throw JsonError("case.accuracy.quantity_overrides names unrequested quantity '" + item.first + "'");
      result->accuracy.quantityOverrides.emplace(item.first,
          ParseAccuracyRequirement(item.second, "case.accuracy.quantity_overrides." + item.first,
                                   result->accuracy.global, true));
    }
  }
  if(stage3Schema && std::find(result->quantities.begin(),result->quantities.end(),"TID")!=
       result->quantities.end()){
    const auto checkTidBudget=[&](const AccuracyRequirement& requirement,
                                  const std::string& context){
      if(requirement.systematicBudgetComplete && requirement.systematicComponents.empty())
        throw JsonError(context+" case-v2 TID rejects a legacy aggregate uncertainty; declare the four Stage-4 component categories and all pair correlations");
      if(requirement.scope=="combined" && requirement.systematicBudgetComplete &&
         requirement.intervalMethod!="student_t_mc_plus_normal_systematics_quadrature")
        throw JsonError(context+" combined Stage-4 TID requires coverage.interval_method='student_t_mc_plus_normal_systematics_quadrature'");
      if(requirement.scope=="combined" && requirement.systematicBudgetComplete &&
         (std::abs(requirement.standardK-1.0)>1.0e-12 ||
          std::abs(requirement.intervalProbability-0.95)>1.0e-12))
        throw JsonError(context+" combined Stage-4 TID requires standard_k=1 and interval_probability=0.95");
    };
    checkTidBudget(result->accuracy.global,"case.accuracy");
    const auto override=result->accuracy.quantityOverrides.find("TID");
    if(override!=result->accuracy.quantityOverrides.end())
      checkTidBudget(override->second,"case.accuracy.quantity_overrides.TID");
  }
  if(stage3Schema && requestsDamage){
    const auto checkDamageBudget=[&](const AccuracyRequirement& requirement,
                                     const std::string& context){
      // The reference tables publish no complete uncertainty model that this
      // code can legitimately invent. An incomplete budget is accepted but
      // will make combined accuracy not_assessable in the result. A claimed
      // complete budget must instead carry a case-owned reference term and its
      // provenance through the existing structured covariance contract.
      if(!requirement.systematicBudgetComplete) return;
      if(requirement.systematicComponents.empty())
        throw JsonError(context+" Stage-5 DDD/n_eq rejects a legacy aggregate uncertainty; a case-supplied physics_reference component with provenance is required");
      const bool hasPhysicsReference=std::any_of(requirement.systematicComponents.begin(),
          requirement.systematicComponents.end(),[](const auto& component){
            return component.category=="physics_reference" && !component.provenance.empty();
          });
      if(!hasPhysicsReference)
        throw JsonError(context+" Stage-5 DDD/n_eq requires a physics_reference uncertainty component with provenance");
      if(requirement.scope=="combined" &&
         requirement.intervalMethod!="student_t_mc_plus_normal_systematics_quadrature")
        throw JsonError(context+" combined Stage-5 DDD/n_eq requires coverage.interval_method='student_t_mc_plus_normal_systematics_quadrature'");
      if(requirement.scope=="combined" &&
         (std::abs(requirement.standardK-1.0)>1.0e-12 ||
          std::abs(requirement.intervalProbability-0.95)>1.0e-12))
        throw JsonError(context+" combined Stage-5 DDD/n_eq requires standard_k=1 and interval_probability=0.95");
    };
    for(const std::string quantity:{"DDD","n_eq"}){
      if(std::find(result->quantities.begin(),result->quantities.end(),quantity)==
         result->quantities.end()) continue;
      const auto override=result->accuracy.quantityOverrides.find(quantity);
      checkDamageBudget(override==result->accuracy.quantityOverrides.end()
          ? result->accuracy.global : override->second,
          override==result->accuracy.quantityOverrides.end()
          ? "case.accuracy" : "case.accuracy.quantity_overrides."+quantity);
    }
  }
  const auto& sampling = accuracy.at("sampling").AsObject("case.accuracy.sampling");
  RejectUnknownFields(sampling, {"min_histories","max_histories","batch_size"}, "case.accuracy.sampling");
  result->accuracy.minHistories = PositiveInteger(sampling.at("min_histories"), "case.accuracy.sampling.min_histories");
  result->accuracy.maxHistories = PositiveInteger(sampling.at("max_histories"), "case.accuracy.sampling.max_histories");
  result->accuracy.batchSize = PositiveInteger(sampling.at("batch_size"), "case.accuracy.sampling.batch_size");
  if (result->accuracy.minHistories > result->accuracy.maxHistories)
    throw JsonError("case.accuracy.sampling requires min_histories <= max_histories");
  if (result->accuracy.batchSize > result->accuracy.maxHistories)
    throw JsonError("case.accuracy.sampling.batch_size must be <= max_histories");

  const auto& transport = input.At("transport", "case").AsObject("case.transport");
  RejectUnknownFields(transport, {"physics_list","ion_physics","em_physics","production_cut","max_step","random_seed"}, "case.transport");
  result->physicsList = RequiredString(transport, "physics_list", "case.transport");
  const std::set<std::string> lists = {"FTFP_BERT","FTFP_BERT_HP","Shielding","QGSP_BIC_HP"};
  if (!lists.count(result->physicsList)) throw JsonError("case.transport.physics_list is unsupported");
  if(stage3Schema && !transport.count("ion_physics"))
    throw JsonError("case-v2 transport requires explicit ion_physics");
  if(transport.count("ion_physics"))
    result->ionPhysics=RequiredString(transport,"ion_physics","case.transport");
  if(result->ionPhysics!="reference_list_default" && result->ionPhysics!="G4IonQMDPhysics")
    throw JsonError("case.transport.ion_physics must be reference_list_default or G4IonQMDPhysics");
  if(stage3Schema && !transport.count("em_physics"))
    throw JsonError("case-v2 transport requires explicit em_physics");
  if(transport.count("em_physics"))
    result->emPhysics=RequiredString(transport,"em_physics","case.transport");
  if(result->emPhysics!="reference_list_default" &&
     result->emPhysics!="G4EmStandardPhysics_option4")
    throw JsonError("case.transport.em_physics must be reference_list_default or G4EmStandardPhysics_option4");
  result->productionCutMm = LengthMm(transport.at("production_cut"), "case.transport.production_cut");
  result->maxStepMm = LengthMm(transport.at("max_step"), "case.transport.max_step");
  result->randomSeed = PositiveInteger(transport.at("random_seed"), "case.transport.random_seed");
  const auto& output = input.At("output", "case").AsObject("case.output");
  RejectUnknownFields(output, {"prefix"}, "case.output");
  result->outputPrefix = RequiredString(output, "prefix", "case.output");

  // Validate every explicitly placed object against the world and against the
  // ordered upstream-to-downstream layout.  Absorbers are checked separately
  // because each is transported in an independent world, not as a stack.
  if(!stage3Schema){
    const double wx0 = result->worldCenterXmm - result->worldSizeXmm / 2.0;
    const double wx1 = result->worldCenterXmm + result->worldSizeXmm / 2.0;
    const double wy0 = result->worldCenterYmm - result->worldSizeYmm / 2.0;
    const double wy1 = result->worldCenterYmm + result->worldSizeYmm / 2.0;
    const double wz0 = result->worldCenterZmm - result->worldSizeZmm / 2.0;
    const double wz1 = result->worldCenterZmm + result->worldSizeZmm / 2.0;
    auto transverseInside = [&](double cx,double cy,double sx, double sy) {
      return cx - sx / 2.0 >= wx0 && cx + sx / 2.0 <= wx1 &&
             cy - sy / 2.0 >= wy0 && cy + sy / 2.0 <= wy1;
    };
    if (!transverseInside(result->sourceCenterXmm,result->sourceCenterYmm,
                          result->sourceSizeXmm, result->sourceSizeYmm) ||
        result->sourceZmm <= wz0 || result->sourceZmm >= wz1)
      throw JsonError("case.geometry.source_plane must lie strictly inside the world");
    if (!transverseInside(result->shieldCenterXmm,result->shieldCenterYmm,
                          result->shieldSizeXmm, result->shieldSizeYmm))
      throw JsonError("case.geometry.shield_stack transverse bounds exceed the world");
    const double shieldRear = result->shieldFrontZmm + result->ShieldTotalThicknessMm();
    if (result->shieldFrontZmm <= wz0 || shieldRear >= wz1)
      throw JsonError("case.geometry.shield_stack must lie strictly inside the world");
    if (result->sourceZmm >= result->shieldFrontZmm)
      throw JsonError("case.geometry.source_plane.z must be upstream of shield_stack.front_z");
    for (const auto& absorber : result->absorbers) {
      if (!transverseInside(absorber.centerXmm,absorber.centerYmm,
                            absorber.sizeXmm, absorber.sizeYmm) || absorber.frontZmm < shieldRear ||
          absorber.frontZmm + absorber.thicknessMm >= wz1)
        throw JsonError("absorber '" + absorber.id + "' is outside the world or overlaps the shield");
    }
  } else {
    if(!BoxInsideWorld(result->sourceTransform,result->sourceSizeXmm,result->sourceSizeYmm,0.0,
                       result->worldSizeXmm,result->worldSizeYmm,result->worldSizeZmm))
      throw JsonError("case.geometry.source_plane transformed corners exceed the world");
    RigidTransform shieldBox=result->shieldTransform;
    shieldBox.translationMm=result->shieldTransform.TransformPoint(
        {0.0,0.0,result->ShieldTotalThicknessMm()/2.0});
    if(!BoxInsideWorld(shieldBox,result->shieldSizeXmm,result->shieldSizeYmm,
                       result->ShieldTotalThicknessMm(),result->worldSizeXmm,
                       result->worldSizeYmm,result->worldSizeZmm))
      throw JsonError("case.geometry.shield_stack transformed bounds exceed the world");
    const Vector3d shieldNormal=result->shieldTransform.Axes()[2];
    if(Dot(result->sourceTransform.translationMm-result->shieldTransform.translationMm,
           shieldNormal)>=-1.0e-9)
      throw JsonError("case.geometry.source_plane origin must be upstream of the shield reference surface");
    for(const auto& absorber:result->absorbers){
      RigidTransform absorberBox=absorber.transform;
      absorberBox.translationMm=absorber.transform.TransformPoint(
          {0.0,0.0,absorber.thicknessMm/2.0});
      if(!BoxInsideWorld(absorberBox,absorber.sizeXmm,absorber.sizeYmm,absorber.thicknessMm,
                         result->worldSizeXmm,result->worldSizeYmm,result->worldSizeZmm))
        throw JsonError("absorber '"+absorber.id+"' transformed bounds exceed the world");
      if(Dot(absorber.transform.translationMm-result->shieldTransform.translationMm,
             shieldNormal)<result->ShieldTotalThicknessMm()-1.0e-9)
        throw JsonError("absorber '"+absorber.id+"' reference surface is not downstream of the shield");
      if(BoxesOverlap(shieldBox,{result->shieldSizeXmm/2.0,result->shieldSizeYmm/2.0,
                                 result->ShieldTotalThicknessMm()/2.0},
                      absorberBox,{absorber.sizeXmm/2.0,absorber.sizeYmm/2.0,
                                   absorber.thicknessMm/2.0}))
        throw JsonError("absorber '"+absorber.id+"' overlaps the shield");
    }
  }

  // Serialize all resolved values, including defaults and both forms of layer
  // thickness.  The archived representation is therefore sufficient to
  // reproduce what was transported without consulting parser defaults.
  Json::Array layerJson;
  for (const auto& layer : result->shieldLayers) {
    layerJson.emplace_back(Json::Object{{"areal_density", Measure(layer.arealDensityGPerCm2, "g/cm2")},
      {"input_thickness_type", layer.inputThicknessType}, {"layer_id", layer.id},
      {"material_hash", layer.materialHash}, {"material_id", layer.materialId},
      {"thickness", Measure(layer.thicknessMm, "mm")}});
  }
  Json::Array absorberJson;
  // Absorbers are independent alternative cases, so their array order is not
  // physical. Canonicalize by stable ID; shield-layer order above deliberately
  // remains untouched because that order is physical.
  std::vector<SlabAbsorberDefinition> canonicalAbsorbers=result->absorbers;
  std::sort(canonicalAbsorbers.begin(),canonicalAbsorbers.end(),
            [](const auto& a,const auto& b){ return a.id<b.id; });
  for (const auto& absorber : canonicalAbsorbers) {
    absorberJson.emplace_back(Json::Object{{"absorber_id", absorber.id},
      {"center_x",Measure(absorber.centerXmm,"mm")}, {"center_y",Measure(absorber.centerYmm,"mm")},
      {"front_z", Measure(absorber.frontZmm, "mm")},
      {"geometry_hash", absorber.geometryHash}, {"geometry_type", "slab"}, {"material_hash", absorber.materialHash},
      {"material_id", absorber.materialId}, {"density_g_cm3",absorber.densityGPerCm3},
      {"size_x", Measure(absorber.sizeXmm, "mm")},
      {"size_y", Measure(absorber.sizeYmm, "mm")}, {"thickness", Measure(absorber.thicknessMm, "mm")},
      {"transform",TransformToJson(absorber.transform)}});
  }
  Json::Object geometryJson{
    {"absorbers", Json(std::move(absorberJson))}, {"coordinate_frame", result->coordinateFrame},
    {"shield_front_z", Measure(result->shieldFrontZmm, "mm")}, {"shield_layers", Json(std::move(layerJson))},
    {"shield_center_x",Measure(result->shieldCenterXmm,"mm")},
    {"shield_center_y",Measure(result->shieldCenterYmm,"mm")},
    {"shield_transform",TransformToJson(result->shieldTransform)},
    {"shield_size_x", Measure(result->shieldSizeXmm, "mm")}, {"shield_size_y", Measure(result->shieldSizeYmm, "mm")},
    {"source_size_x", Measure(result->sourceSizeXmm, "mm")}, {"source_size_y", Measure(result->sourceSizeYmm, "mm")},
    {"source_center_x",Measure(result->sourceCenterXmm,"mm")},
    {"source_center_y",Measure(result->sourceCenterYmm,"mm")},
    {"source_transform",TransformToJson(result->sourceTransform)},
    {"source_z", Measure(result->sourceZmm, "mm")}, {"world_center_x", Measure(result->worldCenterXmm, "mm")},
    {"world_center_y", Measure(result->worldCenterYmm, "mm")}, {"world_center_z", Measure(result->worldCenterZmm, "mm")},
    {"world_size_x", Measure(result->worldSizeXmm, "mm")}, {"world_size_y", Measure(result->worldSizeYmm, "mm")},
    {"world_size_z", Measure(result->worldSizeZmm, "mm")}};
  result->geometryHash = Sha256(Json(geometryJson).DumpCanonical());

  Json::Object quantityOverrides;
  for (const auto& item : result->accuracy.quantityOverrides)
    quantityOverrides[item.first] = AccuracyToJson(item.second);
  Json::Object accuracyJson = AccuracyToJson(result->accuracy.global).AsObject();
  accuracyJson["quantity_overrides"] = Json(std::move(quantityOverrides));
  accuracyJson["sampling"] = Json(Json::Object{{"batch_size", static_cast<double>(result->accuracy.batchSize)},
    {"max_histories", static_cast<double>(result->accuracy.maxHistories)},
    {"min_histories", static_cast<double>(result->accuracy.minHistories)}});
  accuracyJson["target_source"] = result->accuracy.targetWasDefaulted ? "default_0.10" : "case_input";
  Json::Array quantityJson;
  for (const auto& quantity : result->quantities) quantityJson.emplace_back(quantity);
  Json::Object sourceJson{{"alpha_max", Measure(result->alphaMaxMeV, "MeV")},
      {"alpha_min", Measure(result->alphaMinMeV, "MeV")}, {"mode", result->sourceMode},
      {"format",result->sourceFormat}, {"transport_component_id",result->sourceComponentId},
      {"metadata",result->spectrum.canonicalMetadata}, {"metadata_path",result->spectrum.metadataPath},
      {"proton_max", Measure(result->protonMaxMeV, "MeV")}, {"proton_min", Measure(result->protonMinMeV, "MeV")},
      {"science_status",result->spectrum.scienceStatus}, {"spectrum_hash", result->spectrum.contentHash},
      {"spectrum_id", result->spectrum.id}, {"spectrum_path",result->spectrum.spectrumPath},
      {"value_type",result->spectrum.valueType}, {"value_unit",result->spectrum.valueUnit}};
  if(result->responseIncident.enabled){
    sourceJson["response_incident"]=Json(Json::Object{
      {"A",result->responseIncident.massNumber},
      {"Z",result->responseIncident.atomicNumber},
      {"energy_coordinate",result->responseIncident.energyCoordinate},
      {"energy_coordinate_unit",result->responseIncident.energyCoordinateUnit},
      {"energy_coordinate_value",result->responseIncident.energyCoordinateValue},
      {"incident_id",result->responseIncident.incidentId},
      {"mu",result->responseIncident.mu},
      {"particle",result->responseIncident.particleName},
      {"phi_rad",result->responseIncident.phiRad},
      {"total_kinetic_energy_MeV",result->responseIncident.totalKineticEnergyMeV}});
  }
  if(!result->directIncidents.empty()){
    Json::Array nodes;
    for(const auto& weighted:result->directIncidents){
      const auto& incident=weighted.incident;
      nodes.emplace_back(Json::Object{{"A",incident.massNumber},{"Z",incident.atomicNumber},
        {"energy_coordinate",incident.energyCoordinate},
        {"energy_coordinate_unit",incident.energyCoordinateUnit},
        {"energy_coordinate_value",incident.energyCoordinateValue},
        {"incident_id",incident.incidentId},{"mu",incident.mu},
        {"particle",incident.particleName},{"phi_rad",incident.phiRad},
        {"total_kinetic_energy_MeV",incident.totalKineticEnergyMeV},
        {"weight",weighted.weight}});
    }
    sourceJson["direct_incident_nodes"]=Json(std::move(nodes));
  }
  Json::Object resolved{
    {"input_schema_version",result->schemaVersion},
    {"accuracy", Json(std::move(accuracyJson))}, {"case_id", result->caseId},
    {"geometry", Json(std::move(geometryJson))}, {"geometry_hash", result->geometryHash},
    {"material_library", Json(Json::Object{{"canonical_records",result->materials->CanonicalDocument()},
      {"content_hash",result->materials->ContentHash()}, {"library_id",result->materials->Id()},
      {"library_version",result->materials->Version()}, {"source_path",result->materials->SourcePath()}})},
    {"purpose", result->purpose},
    {"quantities", Json(std::move(quantityJson))},
    {"schema_version",stage3Schema ? "shieldSim.resolved_case.v2" : "shieldSim.resolved_case.v1"},
    {"source", Json(std::move(sourceJson))},
    {"transport", Json(Json::Object{{"max_step", Measure(result->maxStepMm, "mm")},
      {"em_physics",result->emPhysics},{"ion_physics",result->ionPhysics},{"physics_list", result->physicsList},
      {"production_cut", Measure(result->productionCutMm, "mm")},
      {"random_seed", static_cast<double>(result->randomSeed)}})}};
  if(result->damageResponse){
    resolved["damage_response"]=Json(Json::Object{
      {"canonical_manifest",result->damageResponse->CanonicalManifest()},
      {"library_id",result->damageResponse->LibraryId()},
      {"library_version",result->damageResponse->LibraryVersion()},
      {"manifest_hash",result->damageResponse->ManifestHash()},
      {"manifest_path",result->damageResponse->ManifestPath()},
      {"target_material_id",result->damageResponse->TargetMaterialId()}});
  }
  result->resolvedJson = Json(std::move(resolved));
  result->caseHash = Sha256(result->resolvedJson.DumpCanonical());
  result->resolvedJson.AsObject()["case_hash"] = result->caseHash;
  return result;
}

} // namespace shield_sim
