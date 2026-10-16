#include "MaterialLibrary.hh"

#include "ContentHash.hh"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <sstream>

namespace shield_sim {
namespace {

constexpr double kFractionTolerance = 1.0e-9;

const Json::Object& Object(const Json& value, const std::string& context) {
  return value.AsObject(context);
}

std::string RequiredString(const Json::Object& object, const std::string& field,
                           const std::string& context) {
  const auto found = object.find(field);
  if (found == object.end()) throw JsonError(context + " is missing required field '" + field + "'");
  const std::string value = found->second.AsString(context + '.' + field);
  if (value.empty()) throw JsonError(context + '.' + field + " must not be empty");
  return value;
}

double RequiredFinitePositive(const Json::Object& object, const std::string& field,
                              const std::string& context) {
  const auto found = object.find(field);
  if (found == object.end()) throw JsonError(context + " is missing required field '" + field + "'");
  const double value = found->second.AsNumber(context + '.' + field);
  if (!std::isfinite(value) || value <= 0.0)
    throw JsonError(context + '.' + field + " must be finite and > 0");
  return value;
}

std::vector<std::string> StringArray(const Json& value, const std::string& context) {
  std::vector<std::string> result;
  std::set<std::string> unique;
  for (const auto& item : value.AsArray(context)) {
    const std::string text = item.AsString(context + "[]");
    if (text.empty()) throw JsonError(context + " contains an empty string");
    if (!unique.insert(text).second) throw JsonError(context + " contains duplicate '" + text + "'");
    result.push_back(text);
  }
  return result;
}

MaterialStatus ParseStatus(const std::string& value, const std::string& context) {
  if (value == "reference") return MaterialStatus::Reference;
  if (value == "trade_study") return MaterialStatus::TradeStudy;
  if (value == "template") return MaterialStatus::Template;
  if (value == "project_specific") return MaterialStatus::ProjectSpecific;
  throw JsonError(context + " has unsupported material status '" + value + "'");
}

void ValidateUncertainty(const Json& value,const std::string& context){
  const auto& object=value.AsObject(context);
  RejectUnknownFields(object,{"status","type","value","unit"},context);
  if(object.count("status")){
    if(RequiredString(object,"status",context)!="not_provided" || object.size()!=1)
      throw JsonError(context+" status form must be exactly {'status':'not_provided'}");
    return;
  }
  const std::string type=RequiredString(object,"type",context);
  if(type!="standard" && type!="range" && type!="expanded")
    throw JsonError(context+".type must be standard, range, or expanded");
  const double uncertainty=object.at("value").AsNumber(context+".value");
  if(!std::isfinite(uncertainty) || uncertainty<0.0)
    throw JsonError(context+".value must be finite and >= 0");
  RequiredString(object,"unit",context);
}

// Atomic weights are used only to deterministically convert an atom-count
// recipe to elemental mass fractions.  The values are conventional standard
// atomic weights (or representative mass numbers for elements without a
// standard weight); transport itself uses the selected Geant4 element data.
const std::map<std::string, double>& AtomicWeights() {
  static const std::map<std::string, double> values = {
    {"H",1.008},{"He",4.002602},{"Li",6.94},{"Be",9.0121831},{"B",10.81},{"C",12.011},
    {"N",14.007},{"O",15.999},{"F",18.998403163},{"Ne",20.1797},{"Na",22.98976928},{"Mg",24.305},
    {"Al",26.9815385},{"Si",28.085},{"P",30.973761998},{"S",32.06},{"Cl",35.45},{"Ar",39.948},
    {"K",39.0983},{"Ca",40.078},{"Sc",44.955908},{"Ti",47.867},{"V",50.9415},{"Cr",51.9961},
    {"Mn",54.938044},{"Fe",55.845},{"Co",58.933194},{"Ni",58.6934},{"Cu",63.546},{"Zn",65.38},
    {"Ga",69.723},{"Ge",72.630},{"As",74.921595},{"Se",78.971},{"Br",79.904},{"Kr",83.798},
    {"Rb",85.4678},{"Sr",87.62},{"Y",88.90584},{"Zr",91.224},{"Nb",92.90637},{"Mo",95.95},
    {"Tc",98},{"Ru",101.07},{"Rh",102.90550},{"Pd",106.42},{"Ag",107.8682},{"Cd",112.414},
    {"In",114.818},{"Sn",118.710},{"Sb",121.760},{"Te",127.60},{"I",126.90447},{"Xe",131.293},
    {"Cs",132.90545196},{"Ba",137.327},{"La",138.90547},{"Ce",140.116},{"Pr",140.90766},{"Nd",144.242},
    {"Pm",145},{"Sm",150.36},{"Eu",151.964},{"Gd",157.25},{"Tb",158.92535},{"Dy",162.500},
    {"Ho",164.93033},{"Er",167.259},{"Tm",168.93422},{"Yb",173.045},{"Lu",174.9668},{"Hf",178.49},
    {"Ta",180.94788},{"W",183.84},{"Re",186.207},{"Os",190.23},{"Ir",192.217},{"Pt",195.084},
    {"Au",196.966569},{"Hg",200.592},{"Tl",204.38},{"Pb",207.2},{"Bi",208.98040},{"Po",209},
    {"At",210},{"Rn",222},{"Fr",223},{"Ra",226},{"Ac",227},{"Th",232.0377},{"Pa",231.03588},
    {"U",238.02891},{"Np",237},{"Pu",244},{"Am",243},{"Cm",247},{"Bk",247},{"Cf",251},{"Es",252},
    {"Fm",257},{"Md",258},{"No",259},{"Lr",266},{"Rf",267},{"Db",268},{"Sg",269},{"Bh",270},
    {"Hs",277},{"Mt",278},{"Ds",281},{"Rg",282},{"Cn",285},{"Nh",286},{"Fl",289},{"Mc",290},
    {"Lv",293},{"Ts",294},{"Og",294}
  };
  return values;
}

bool HasRole(const MaterialRecord& record, const std::string& role) {
  return std::find(record.roles.begin(), record.roles.end(), role) != record.roles.end();
}

std::pair<std::string,int> ParseIsotopeId(const std::string& id,const std::string& context){
  const std::size_t dash=id.rfind('-');
  if(dash==std::string::npos || dash==0 || dash+1>=id.size())
    throw JsonError(context+" isotope must use '<element>-<mass_number>', for example B-10");
  const std::string symbol=id.substr(0,dash);
  if(!AtomicWeights().count(symbol)) throw JsonError(context+" has unknown isotope element '"+symbol+"'");
  int mass=0;
  try {
    std::size_t used=0;
    mass=std::stoi(id.substr(dash+1),&used);
    if(used!=id.size()-dash-1) throw std::invalid_argument("trailing");
  } catch(...) { throw JsonError(context+" has invalid isotope mass number in '"+id+"'"); }
  if(mass<=0) throw JsonError(context+" isotope mass number must be > 0");
  return {symbol,mass};
}

void NormalizeAndCheck(std::map<std::string, double>& values, const std::string& context) {
  double sum = 0.0;
  for (const auto& item : values) {
    if (!std::isfinite(item.second) || item.second < 0.0)
      throw JsonError(context + " contains a negative or non-finite fraction");
    sum += item.second;
  }
  if (!std::isfinite(sum) || sum <= 0.0) throw JsonError(context + " has zero total composition");
  for (auto& item : values) item.second /= sum;
}

} // namespace

std::string MaterialStatusName(MaterialStatus status) {
  switch (status) {
    case MaterialStatus::Reference: return "reference";
    case MaterialStatus::TradeStudy: return "trade_study";
    case MaterialStatus::Template: return "template";
    case MaterialStatus::ProjectSpecific: return "project_specific";
  }
  return "unknown";
}

MaterialRecord MaterialLibrary::ParseRecord(const Json& value, const std::string& context) {
  const auto& object = Object(value, context);
  RejectUnknownFields(object,
    {"material_id","revision","aliases","label","status","roles","density","state",
     "composition","geant4_name","mean_excitation_energy","thermal_scattering_law",
     "damage_response_links","provenance","applicability_and_limitations","extends"}, context);

  MaterialRecord record;
  record.id = RequiredString(object, "material_id", context);
  record.revision = RequiredString(object, "revision", context);
  record.label = RequiredString(object, "label", context);
  record.status = ParseStatus(RequiredString(object, "status", context), context + ".status");
  record.roles = StringArray(value.At("roles", context), context + ".roles");
  if (record.roles.empty()) throw JsonError(context + ".roles must not be empty");
  for (const auto& role : record.roles)
    if (role != "shield" && role != "absorber")
      throw JsonError(context + ".roles contains unsupported role '" + role + "'");
  if (object.count("aliases")) record.aliases = StringArray(object.at("aliases"), context + ".aliases");
  if (object.count("extends")) record.parent = object.at("extends").AsString(context + ".extends");
  if (object.count("geant4_name")) record.geant4Name = object.at("geant4_name").AsString(context + ".geant4_name");
  if (object.count("mean_excitation_energy")) {
    const auto& excitation=object.at("mean_excitation_energy").AsObject(context+".mean_excitation_energy");
    RejectUnknownFields(excitation,{"value","unit"},context+".mean_excitation_energy");
    record.meanExcitationEnergyEV=RequiredFinitePositive(excitation,"value",context+".mean_excitation_energy");
    if(RequiredString(excitation,"unit",context+".mean_excitation_energy")!="eV")
      throw JsonError(context+".mean_excitation_energy.unit must be exactly 'eV'");
    record.hasMeanExcitationEnergy=true;
  }
  if(object.count("thermal_scattering_law")){
    const auto& thermal=object.at("thermal_scattering_law").AsObject(context+".thermal_scattering_law");
    RejectUnknownFields(thermal,{"status"},context+".thermal_scattering_law");
    if(RequiredString(thermal,"status",context+".thermal_scattering_law")!="none")
      throw JsonError(context+".thermal_scattering_law associations are not enabled by the Stage-1 transport factory");
  }
  if(object.count("damage_response_links")){
    const auto& links=object.at("damage_response_links").AsArray(context+".damage_response_links");
    for(std::size_t i=0;i<links.size();++i){
      const std::string lc=context+".damage_response_links["+std::to_string(i)+']';
      const auto& link=links[i].AsObject(lc);
      RejectUnknownFields(link,{"quantity","dataset_id","version","target_material_hash",
                                "incident_species","energy_domain","interpolation","out_of_range_policy"},lc);
      RequiredString(link,"quantity",lc); RequiredString(link,"dataset_id",lc);
      RequiredString(link,"version",lc); RequiredString(link,"target_material_hash",lc);
      RequiredString(link,"interpolation",lc); RequiredString(link,"out_of_range_policy",lc);
      if(!link.count("incident_species") || link.at("incident_species").AsArray(lc+".incident_species").empty())
        throw JsonError(lc+".incident_species must not be empty");
      if(!link.count("energy_domain") || !link.at("energy_domain").IsObject())
        throw JsonError(lc+".energy_domain must be an object");
    }
  }
  record.applicability = RequiredString(object, "applicability_and_limitations", context);

  if (!object.count("provenance")) throw JsonError(context + " is missing required field 'provenance'");
  record.provenance = object.at("provenance");
  const auto& provenance = record.provenance.AsArray(context + ".provenance");
  if (provenance.empty()) throw JsonError(context + ".provenance must not be empty");
  for (std::size_t i = 0; i < provenance.size(); ++i) {
    const auto& p = provenance[i].AsObject(context + ".provenance[]");
    RejectUnknownFields(p, {"properties","source","version","retrieval_date","transformation","license"},
                        context + ".provenance[]");
    RequiredString(p, "source", context + ".provenance[]");
    RequiredString(p, "version", context + ".provenance[]");
    if (!p.count("properties") || StringArray(p.at("properties"), context + ".provenance[].properties").empty())
      throw JsonError(context + ".provenance[] must name at least one property");
  }

  // Templates deliberately omit unresolved density/composition.  Every other
  // status is transport-usable and therefore must be fully explicit.
  if (record.status == MaterialStatus::Template) {
    if (object.count("density") || object.count("composition"))
      throw JsonError(context + ": template records must not masquerade as complete transport materials");
    record.canonicalRecord = value;
    record.contentHash = Sha256(value.DumpCanonical());
    return record;
  }

  const auto& density = value.At("density", context).AsObject(context + ".density");
  RejectUnknownFields(density, {"value","unit","uncertainty","reference_condition"}, context + ".density");
  record.densityGPerCm3 = RequiredFinitePositive(density, "value", context + ".density");
  if (RequiredString(density, "unit", context + ".density") != "g/cm3")
    throw JsonError(context + ".density.unit must be exactly 'g/cm3'");
  if (!density.count("uncertainty"))
    throw JsonError(context + ".density must explicitly contain uncertainty metadata");
  ValidateUncertainty(density.at("uncertainty"),context+".density.uncertainty");

  const auto& state = value.At("state", context).AsObject(context + ".state");
  RejectUnknownFields(state, {"phase","temperature","pressure"}, context + ".state");
  record.phase = RequiredString(state, "phase", context + ".state");
  if (record.phase != "solid" && record.phase != "liquid" && record.phase != "gas")
    throw JsonError(context + ".state.phase must be solid, liquid, or gas");
  const auto& temperature = value.At("state", context).At("temperature", context + ".state")
                                  .AsObject(context + ".state.temperature");
  RejectUnknownFields(temperature, {"value","unit"}, context + ".state.temperature");
  record.temperatureK = RequiredFinitePositive(temperature, "value", context + ".state.temperature");
  if (RequiredString(temperature, "unit", context + ".state.temperature") != "K")
    throw JsonError(context + ".state.temperature.unit must be exactly 'K'");
  if(state.count("pressure")){
    const auto& pressure=state.at("pressure").AsObject(context+".state.pressure");
    RejectUnknownFields(pressure,{"value","unit"},context+".state.pressure");
    record.pressurePa=RequiredFinitePositive(pressure,"value",context+".state.pressure");
    const std::string pressureUnit=RequiredString(pressure,"unit",context+".state.pressure");
    if(pressureUnit=="kPa") record.pressurePa*=1000.0;
    else if(pressureUnit!="Pa")
      throw JsonError(context+".state.pressure.unit must be Pa or kPa");
    record.hasPressure=true;
  }
  if(record.phase=="gas" && !record.hasPressure)
    throw JsonError(context+".state.pressure is required for a gas");

  const auto& composition = value.At("composition", context).AsObject(context + ".composition");
  RejectUnknownFields(composition, {"basis","components","residual_policy","isotopes"}, context + ".composition");
  record.compositionBasis = RequiredString(composition, "basis", context + ".composition");
  if (record.compositionBasis != "elemental_mass_fraction" &&
      record.compositionBasis != "atom_count" &&
      record.compositionBasis != "constituent_mass_fraction" &&
      record.compositionBasis != "constituent_volume_fraction")
    throw JsonError(context + ".composition.basis is unsupported");
  record.residualPolicy = RequiredString(composition, "residual_policy", context + ".composition");
  if (record.residualPolicy != "exact")
    throw JsonError(context + ".composition.residual_policy must be 'exact' in Stage 1");
  const std::string isotopeMode=composition.count("isotopes")
      ? composition.at("isotopes").AsString(context + ".composition.isotopes") : "";
  const bool constituentBasis=record.compositionBasis=="constituent_mass_fraction" ||
                              record.compositionBasis=="constituent_volume_fraction";
  if ((!constituentBasis && isotopeMode!="natural_geant4" && isotopeMode!="explicit") ||
      (constituentBasis && !isotopeMode.empty()))
    throw JsonError(context + ".composition.isotopes must be 'natural_geant4' or 'explicit' for elemental recipes and omitted for constituent recipes");

  const auto& components = composition.at("components").AsArray(context + ".composition.components");
  if (components.empty()) throw JsonError(context + ".composition.components must not be empty");
  std::set<std::string> seen;
  double sum = 0.0;
  for (std::size_t i = 0; i < components.size(); ++i) {
    const std::string cc = context + ".composition.components[" + std::to_string(i) + ']';
    const auto& component = components[i].AsObject(cc);
    RejectUnknownFields(component, {"element","isotope","material_id","value","uncertainty"}, cc);
    int kinds = static_cast<int>(component.count("element")) +
                static_cast<int>(component.count("isotope")) +
                static_cast<int>(component.count("material_id"));
    if (kinds != 1) throw JsonError(cc + " must contain exactly one of element, isotope, or material_id");
    MaterialComponent parsed;
    if (component.count("element")) {
      parsed.kind = MaterialComponent::Kind::Element;
      parsed.id = RequiredString(component, "element", cc);
      if (!AtomicWeights().count(parsed.id)) throw JsonError(cc + " has unknown element '" + parsed.id + "'");
    } else if (component.count("isotope")) {
      parsed.kind = MaterialComponent::Kind::Isotope;
      parsed.id = RequiredString(component, "isotope", cc);
      (void)ParseIsotopeId(parsed.id,cc);
      if(isotopeMode!="explicit")
        throw JsonError(cc+" uses an explicit isotope but composition.isotopes is not 'explicit'");
    } else {
      parsed.kind = MaterialComponent::Kind::Material;
      parsed.id = RequiredString(component, "material_id", cc);
    }
    if (!seen.insert(std::to_string(static_cast<int>(parsed.kind)) + ':' + parsed.id).second)
      throw JsonError(cc + " duplicates component '" + parsed.id + "'");
    parsed.value = RequiredFinitePositive(component, "value", cc);
    if(component.count("uncertainty"))
      ValidateUncertainty(component.at("uncertainty"),cc+".uncertainty");
    sum += parsed.value;
    record.components.push_back(parsed);
  }

  const bool fractionBasis = record.compositionBasis != "atom_count";
  if (fractionBasis && std::abs(sum - 1.0) > kFractionTolerance)
    throw JsonError(context + ".composition fractions must sum to 1 within 1e-9; got " + std::to_string(sum));
  for (const auto& component : record.components) {
    const bool wantsMaterial = record.compositionBasis == "constituent_mass_fraction" ||
                               record.compositionBasis == "constituent_volume_fraction";
    if (wantsMaterial != (component.kind == MaterialComponent::Kind::Material))
      throw JsonError(context + ".composition component kind does not match basis '" + record.compositionBasis + "'");
  }

  record.canonicalRecord = value;
  record.contentHash = Sha256(value.DumpCanonical());
  return record;
}

Json MaterialLibrary::MergeOverride(const MaterialRecord& parent, const Json& overrideValue,
                                    const std::string& context) {
  const auto& overrideObject = overrideValue.AsObject(context);
  RejectUnknownFields(overrideObject,
    {"material_id","revision","extends","aliases","label","status","roles","density","state",
     "composition","geant4_name","mean_excitation_energy","thermal_scattering_law",
     "damage_response_links","provenance","applicability_and_limitations"}, context);
  RequiredString(overrideObject, "material_id", context);
  RequiredString(overrideObject, "revision", context);
  const std::string extends = RequiredString(overrideObject, "extends", context);
  if (extends != parent.id + '@' + parent.revision)
    throw JsonError(context + ".extends must pin exact parent '" + parent.id + '@' + parent.revision + "'");
  if (!overrideObject.count("status"))
    throw JsonError(context + " must explicitly set status");
  const std::string status = overrideObject.at("status").AsString(context + ".status");
  if (status != "trade_study" && status != "project_specific")
    throw JsonError(context + ".status must be trade_study or project_specific");
  if (!overrideObject.count("provenance"))
    throw JsonError(context + " must provide provenance for the override");

  Json::Object merged = parent.canonicalRecord.AsObject(context + " parent");
  // Aliases are never inherited: a case-local record must be selected through
  // its new stable ID, avoiding an alias that silently changes meaning.
  merged.erase("aliases");
  // Composition changes invalidate inherited response/interaction metadata.
  // The override may reintroduce them only by supplying explicit compatible
  // fields, which are included in its new content hash.
  if(overrideObject.count("composition")){
    if(!overrideObject.count("mean_excitation_energy")) merged.erase("mean_excitation_energy");
    if(!overrideObject.count("thermal_scattering_law")) merged.erase("thermal_scattering_law");
    if(!overrideObject.count("damage_response_links")) merged.erase("damage_response_links");
  }
  for (const auto& field : overrideObject) merged[field.first] = field.second;
  return Json(std::move(merged));
}

void MaterialLibrary::Insert(MaterialRecord record, const std::string& context) {
  if (records_.count(record.id)) throw JsonError(context + ": duplicate material_id '" + record.id + "'");
  if (aliases_.count(record.id)) throw JsonError(context + ": material_id collides with alias '" + record.id + "'");
  const std::string id = record.id;
  for (const auto& alias : record.aliases) {
    if (records_.count(alias) || aliases_.count(alias))
      throw JsonError(context + ": duplicate/colliding material alias '" + alias + "'");
    aliases_[alias] = id;
  }
  records_.emplace(id, std::move(record));
}

MaterialLibrary MaterialLibrary::Load(const std::string& path, const Json::Array& caseOverrides) {
  MaterialLibrary library;
  library.sourcePath_ = path;
  Json document = Json::ParseFile(path);
  const auto& root = document.AsObject("material library");
  RejectUnknownFields(root, {"schema_version","library_id","library_version","materials"}, "material library");
  if (RequiredString(root, "schema_version", "material library") != "shieldSim.material_library.v1")
    throw JsonError("material library schema_version must be 'shieldSim.material_library.v1'");
  library.id_ = RequiredString(root, "library_id", "material library");
  library.version_ = RequiredString(root, "library_version", "material library");
  const auto& records = document.At("materials", "material library").AsArray("material library.materials");
  if (records.empty()) throw JsonError("material library.materials must not be empty");
  for (std::size_t i = 0; i < records.size(); ++i) {
    const std::string context = "material library.materials[" + std::to_string(i) + ']';
    library.Insert(ParseRecord(records[i], context), context);
  }

  // Verify all constituent recipes immediately so a bad unused material
  // cannot hide in a library that otherwise appears valid.
  for (const auto& item : library.records_) {
    if (item.second.status != MaterialStatus::Template)
      (void)library.ElementalMassFractions(item.second);
  }

  for (std::size_t i = 0; i < caseOverrides.size(); ++i) {
    const std::string context = "case.material_overrides[" + std::to_string(i) + ']';
    const auto& object = caseOverrides[i].AsObject(context);
    if(object.count("extends")){
      const std::string extends = RequiredString(object, "extends", context);
      const std::size_t at = extends.rfind('@');
      if (at == std::string::npos) throw JsonError(context + ".extends must be '<material_id>@<revision>'");
      const std::string parentId = extends.substr(0, at);
      const auto found = library.records_.find(parentId);
      if (found == library.records_.end()) throw JsonError(context + " extends unknown material '" + parentId + "'");
      Json merged = MergeOverride(found->second, caseOverrides[i], context);
      library.Insert(ParseRecord(merged, context + " resolved"), context);
    } else {
      MaterialRecord created=ParseRecord(caseOverrides[i],context);
      if(created.status!=MaterialStatus::TradeStudy && created.status!=MaterialStatus::ProjectSpecific)
        throw JsonError(context+" new case material must have trade_study or project_specific status");
      library.Insert(std::move(created),context);
    }
  }

  // The library hash includes both the shared document and case-local resolved
  // records.  This ensures a density/composition override invalidates cached
  // geometry and response data even if the original file is unchanged.
  Json::Array resolved;
  for (const auto& item : library.records_) resolved.push_back(item.second.canonicalRecord);
  Json::Object canonical;
  canonical["schema_version"] = "shieldSim.resolved_material_library.v1";
  canonical["library_id"] = library.id_;
  canonical["library_version"] = library.version_;
  canonical["materials"] = Json(std::move(resolved));
  library.canonicalDocument_ = Json(std::move(canonical));
  library.contentHash_ = Sha256(library.canonicalDocument_.DumpCanonical());
  return library;
}

const MaterialRecord& MaterialLibrary::ById(const std::string& id) const {
  const auto found = records_.find(id);
  if (found == records_.end()) throw JsonError("unknown material_id '" + id + "'");
  return found->second;
}

const MaterialRecord& MaterialLibrary::Resolve(const std::string& idOrAlias,
                                               const std::string& requiredRole) const {
  std::string id = idOrAlias;
  const auto alias = aliases_.find(idOrAlias);
  if (alias != aliases_.end()) id = alias->second;
  const MaterialRecord& record = ById(id);
  if (record.status == MaterialStatus::Template)
    throw JsonError("material '" + idOrAlias + "' resolves to disabled template '" + record.id + "'");
  if (!HasRole(record, requiredRole))
    throw JsonError("material '" + record.id + "' is not approved for role '" + requiredRole + "'");
  return record;
}

std::map<std::string, double> MaterialLibrary::Flatten(const MaterialRecord& record,
                                                       std::vector<std::string>& stack) const {
  if (std::find(stack.begin(), stack.end(), record.id) != stack.end()) {
    std::ostringstream path;
    for (const auto& id : stack) path << id << " -> ";
    throw JsonError("constituent material cycle: " + path.str() + record.id);
  }
  stack.push_back(record.id);
  std::map<std::string, double> result;

  if (record.compositionBasis == "elemental_mass_fraction") {
    for (const auto& component : record.components)
      result[(component.kind==MaterialComponent::Kind::Isotope ? "isotope:" : "")+component.id] += component.value;
  } else if (record.compositionBasis == "atom_count") {
    for (const auto& component : record.components){
      if(component.kind==MaterialComponent::Kind::Isotope){
        const auto isotope=ParseIsotopeId(component.id,"material '"+record.id+"'");
        result["isotope:"+component.id] += component.value*isotope.second;
      } else {
        result[component.id] += component.value * AtomicWeights().at(component.id);
      }
    }
  } else {
    double constituentTotal = 0.0;
    std::vector<double> weights;
    for (const auto& component : record.components) {
      const MaterialRecord& child = ById(component.id);
      if (child.status == MaterialStatus::Template)
        throw JsonError("material '" + record.id + "' uses unresolved template constituent '" + child.id + "'");
      double weight = component.value;
      if (record.compositionBasis == "constituent_volume_fraction")
        weight *= child.densityGPerCm3;
      weights.push_back(weight);
      constituentTotal += weight;
    }
    if (constituentTotal <= 0.0) throw JsonError("material '" + record.id + "' has zero constituent mass");
    for (std::size_t i = 0; i < record.components.size(); ++i) {
      const auto childFractions = Flatten(ById(record.components[i].id), stack);
      for (const auto& element : childFractions)
        result[element.first] += (weights[i] / constituentTotal) * element.second;
    }
  }

  stack.pop_back();
  NormalizeAndCheck(result, "material '" + record.id + "'");
  return result;
}

std::map<std::string, double> MaterialLibrary::ElementalMassFractions(const MaterialRecord& record) const {
  std::vector<std::string> stack;
  return Flatten(record, stack);
}

std::vector<std::string> MaterialLibrary::MaterialIds() const {
  std::vector<std::string> ids;
  for (const auto& item : records_) ids.push_back(item.first);
  return ids;
}

} // namespace shield_sim
