#include "CaseDefinition.hh"

#include "ContentHash.hh"
#include "PublishedSourceAdapters.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <set>
#include <sstream>

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
    RejectUnknownFields(budget,{"status","relative_standard_uncertainty","known_relative_floor"},
                        context+".systematic_budget");
    const std::string status=RequiredString(budget,"status",context+".systematic_budget");
    if(status=="complete"){
      if(!budget.count("relative_standard_uncertainty"))
        throw JsonError(context+".systematic_budget complete status requires relative_standard_uncertainty");
      result.systematicRelativeUncertainty=FiniteNumber(
          budget.at("relative_standard_uncertainty"),context+".systematic_budget.relative_standard_uncertainty");
      if(result.systematicRelativeUncertainty<0.0)
        throw JsonError(context+".systematic_budget.relative_standard_uncertainty must be >= 0");
      result.systematicBudgetComplete=true;
    } else if(status=="incomplete"){
      if(budget.count("relative_standard_uncertainty"))
        throw JsonError(context+".systematic_budget cannot claim a combined uncertainty while incomplete");
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
  if(requirement.systematicBudgetComplete)
    budget["relative_standard_uncertainty"]=requirement.systematicRelativeUncertainty;
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
     "source","quantities","accuracy","transport","output"}, "case");
  if (RequiredString(root, "schema_version", "case") != "shieldSim.case.v1")
    throw JsonError("case.schema_version must be 'shieldSim.case.v1'");
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
  RejectUnknownFields(frame, {"name","slab_normal","origin"}, "case.geometry.coordinate_frame");
  result->coordinateFrame = RequiredString(frame, "name", "case.geometry.coordinate_frame");
  if (RequiredString(frame, "slab_normal", "case.geometry.coordinate_frame") != "+z")
    throw JsonError("Stage-1 slab geometry requires coordinate_frame.slab_normal='+z'");
  if (RequiredString(frame, "origin", "case.geometry.coordinate_frame") != "world_center")
    throw JsonError("Stage-1 slab geometry requires coordinate_frame.origin='world_center'");

  const auto& world = input.At("geometry", "case").At("world", "case.geometry").AsObject("case.geometry.world");
  RejectUnknownFields(world, {"size_x","size_y","size_z","center_x","center_y","center_z"}, "case.geometry.world");
  result->worldSizeXmm = LengthMm(world.at("size_x"), "case.geometry.world.size_x");
  result->worldSizeYmm = LengthMm(world.at("size_y"), "case.geometry.world.size_y");
  result->worldSizeZmm = LengthMm(world.at("size_z"), "case.geometry.world.size_z");
  result->worldCenterXmm = PositionMm(world.at("center_x"), "case.geometry.world.center_x");
  result->worldCenterYmm = PositionMm(world.at("center_y"), "case.geometry.world.center_y");
  result->worldCenterZmm = PositionMm(world.at("center_z"), "case.geometry.world.center_z");
  if (result->worldCenterXmm != 0.0 || result->worldCenterYmm != 0.0 || result->worldCenterZmm != 0.0)
    throw JsonError("Stage-1 Geant4 world center must be explicitly [0,0,0] mm");

  const auto& sourcePlane = input.At("geometry", "case").At("source_plane", "case.geometry")
                                  .AsObject("case.geometry.source_plane");
  RejectUnknownFields(sourcePlane, {"size_x","size_y","center_x","center_y","z"}, "case.geometry.source_plane");
  result->sourceSizeXmm = LengthMm(sourcePlane.at("size_x"), "case.geometry.source_plane.size_x");
  result->sourceSizeYmm = LengthMm(sourcePlane.at("size_y"), "case.geometry.source_plane.size_y");
  result->sourceCenterXmm = PositionMm(sourcePlane.at("center_x"), "case.geometry.source_plane.center_x");
  result->sourceCenterYmm = PositionMm(sourcePlane.at("center_y"), "case.geometry.source_plane.center_y");
  result->sourceZmm = PositionMm(sourcePlane.at("z"), "case.geometry.source_plane.z");

  const auto& stack = input.At("geometry", "case").At("shield_stack", "case.geometry")
                            .AsObject("case.geometry.shield_stack");
  RejectUnknownFields(stack, {"size_x","size_y","center_x","center_y","front_z","layers"}, "case.geometry.shield_stack");
  result->shieldSizeXmm = LengthMm(stack.at("size_x"), "case.geometry.shield_stack.size_x");
  result->shieldSizeYmm = LengthMm(stack.at("size_y"), "case.geometry.shield_stack.size_y");
  result->shieldCenterXmm = PositionMm(stack.at("center_x"), "case.geometry.shield_stack.center_x");
  result->shieldCenterYmm = PositionMm(stack.at("center_y"), "case.geometry.shield_stack.center_y");
  result->shieldFrontZmm = PositionMm(stack.at("front_z"), "case.geometry.shield_stack.front_z");
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
    const auto& shape = object.at("geometry").AsObject(context + ".geometry");
    RejectUnknownFields(shape, {"type","size_x","size_y","center_x","center_y","thickness","front_z"}, context + ".geometry");
    if (RequiredString(shape, "type", context + ".geometry") != "slab")
      throw JsonError(context + ".geometry.type must be 'slab' in Stage 1");
    parsed.sizeXmm = LengthMm(shape.at("size_x"), context + ".geometry.size_x");
    parsed.sizeYmm = LengthMm(shape.at("size_y"), context + ".geometry.size_y");
    parsed.centerXmm = PositionMm(shape.at("center_x"), context + ".geometry.center_x");
    parsed.centerYmm = PositionMm(shape.at("center_y"), context + ".geometry.center_y");
    parsed.thicknessMm = LengthMm(shape.at("thickness"), context + ".geometry.thickness");
    parsed.frontZmm = PositionMm(shape.at("front_z"), context + ".geometry.front_z");
    Json::Object hashInput{{"absorber_id", parsed.id}, {"center_x_mm",parsed.centerXmm},
      {"center_y_mm",parsed.centerYmm}, {"front_z_mm", parsed.frontZmm},
      {"material_hash", parsed.materialHash}, {"size_x_mm", parsed.sizeXmm},
      {"size_y_mm", parsed.sizeYmm}, {"thickness_mm", parsed.thicknessMm}, {"type", "slab"}};
    parsed.geometryHash = Sha256(Json(std::move(hashInput)).DumpCanonical());
    result->absorbers.push_back(parsed);
  }

  const auto& source = input.At("source", "case").AsObject("case.source");
  RejectUnknownFields(source, {"mode","spectrum","energy_limits"}, "case.source");
  result->sourceMode = RequiredString(source, "mode", "case.source");
  if (result->sourceMode != "beam" && result->sourceMode != "isotropic")
    throw JsonError("case.source.mode must be beam or isotropic");
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
  } else {
    throw JsonError("case.source.spectrum.format must be legacy_three_column or ams02_hapi_component");
  }
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
    if (quantity != "TID" && quantity != "H100/10")
      throw JsonError("quantity '" + quantity + "' is not enabled in Stage 1; DDD, n_eq, and LET require their later validation stages");
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
  RejectUnknownFields(transport, {"physics_list","production_cut","max_step","random_seed"}, "case.transport");
  result->physicsList = RequiredString(transport, "physics_list", "case.transport");
  const std::set<std::string> lists = {"FTFP_BERT","FTFP_BERT_HP","Shielding","QGSP_BIC_HP"};
  if (!lists.count(result->physicsList)) throw JsonError("case.transport.physics_list is unsupported");
  result->productionCutMm = LengthMm(transport.at("production_cut"), "case.transport.production_cut");
  result->maxStepMm = LengthMm(transport.at("max_step"), "case.transport.max_step");
  result->randomSeed = PositiveInteger(transport.at("random_seed"), "case.transport.random_seed");
  const auto& output = input.At("output", "case").AsObject("case.output");
  RejectUnknownFields(output, {"prefix"}, "case.output");
  result->outputPrefix = RequiredString(output, "prefix", "case.output");

  // Validate every explicitly placed object against the world and against the
  // ordered upstream-to-downstream layout.  Absorbers are checked separately
  // because each is transported in an independent world, not as a stack.
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
      {"material_id", absorber.materialId}, {"size_x", Measure(absorber.sizeXmm, "mm")},
      {"size_y", Measure(absorber.sizeYmm, "mm")}, {"thickness", Measure(absorber.thicknessMm, "mm")}});
  }
  Json::Object geometryJson{
    {"absorbers", Json(std::move(absorberJson))}, {"coordinate_frame", result->coordinateFrame},
    {"shield_front_z", Measure(result->shieldFrontZmm, "mm")}, {"shield_layers", Json(std::move(layerJson))},
    {"shield_center_x",Measure(result->shieldCenterXmm,"mm")},
    {"shield_center_y",Measure(result->shieldCenterYmm,"mm")},
    {"shield_size_x", Measure(result->shieldSizeXmm, "mm")}, {"shield_size_y", Measure(result->shieldSizeYmm, "mm")},
    {"source_size_x", Measure(result->sourceSizeXmm, "mm")}, {"source_size_y", Measure(result->sourceSizeYmm, "mm")},
    {"source_center_x",Measure(result->sourceCenterXmm,"mm")},
    {"source_center_y",Measure(result->sourceCenterYmm,"mm")},
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
  Json::Object resolved{
    {"accuracy", Json(std::move(accuracyJson))}, {"case_id", result->caseId},
    {"geometry", Json(std::move(geometryJson))}, {"geometry_hash", result->geometryHash},
    {"material_library", Json(Json::Object{{"canonical_records",result->materials->CanonicalDocument()},
      {"content_hash",result->materials->ContentHash()}, {"library_id",result->materials->Id()},
      {"library_version",result->materials->Version()}, {"source_path",result->materials->SourcePath()}})},
    {"purpose", result->purpose},
    {"quantities", Json(std::move(quantityJson))}, {"schema_version", "shieldSim.resolved_case.v1"},
    {"source", Json(Json::Object{{"alpha_max", Measure(result->alphaMaxMeV, "MeV")},
      {"alpha_min", Measure(result->alphaMinMeV, "MeV")}, {"mode", result->sourceMode},
      {"format",result->sourceFormat}, {"transport_component_id",result->sourceComponentId},
      {"metadata",result->spectrum.canonicalMetadata}, {"metadata_path",result->spectrum.metadataPath},
      {"proton_max", Measure(result->protonMaxMeV, "MeV")}, {"proton_min", Measure(result->protonMinMeV, "MeV")},
      {"science_status",result->spectrum.scienceStatus}, {"spectrum_hash", result->spectrum.contentHash},
      {"spectrum_id", result->spectrum.id}, {"spectrum_path",result->spectrum.spectrumPath},
      {"value_type",result->spectrum.valueType}, {"value_unit",result->spectrum.valueUnit}})},
    {"transport", Json(Json::Object{{"max_step", Measure(result->maxStepMm, "mm")},
      {"physics_list", result->physicsList}, {"production_cut", Measure(result->productionCutMm, "mm")},
      {"random_seed", static_cast<double>(result->randomSeed)}})}};
  result->resolvedJson = Json(std::move(resolved));
  result->caseHash = Sha256(result->resolvedJson.DumpCanonical());
  result->resolvedJson.AsObject()["case_hash"] = result->caseHash;
  return result;
}

} // namespace shield_sim
