#include "MaterialCatalog.hh"

/* ============================================================================
 * MaterialCatalog.cc
 *
 * Material implementation and reference notes
 * -------------------------------------------
 * This file is the single place where shieldSim maps user-facing material names
 * to Geant4 materials.  The logic is deliberately centralized because material
 * names are used by several otherwise independent parts of the code:
 *
 *   - CLI parsing and --help/--list-materials output,
 *   - detector construction,
 *   - dose-vs-thickness sweep mode,
 *   - configuration printouts and Tecplot metadata.
 *
 * Material categories
 * -------------------
 * 1. Geant4/NIST aliases
 *      Entries such as Al, Cu, W, Ta, Kapton, and Water map directly to G4_*
 *      material names.  Their density/composition is therefore the value in the
 *      Geant4 material database used by the installed Geant4 version.
 *
 * 2. ShieldSim custom engineering approximations
 *      Entries such as HDPE, BPE, Kevlar, CFRP, SiCComposite, LunarRegolith,
 *      and MarsRegolith are built explicitly below.  They are intended for
 *      preliminary shielding trade studies.  The comments above each builder
 *      state whether the composition is defined by atom count or by mass
 *      fraction and give the reference assumptions.
 *
 * References used for the default catalog
 * ---------------------------------------
 * The help text prints compact reference notes.  The longer URLs are listed
 * here so the material-property provenance stays with the source code.
 *
 * [G4-MAT]  Geant4 Book for Application Developers, Appendix: Geant4 Material
 *           Database, material names and predefined G4_* materials:
 *           https://geant4.web.cern.ch/documentation/dev/bfad_html/
 *           ForApplicationDevelopers/Appendix/materialNames.html
 *
 * [NIST-XCOM] NIST XCOM / X-Ray Mass Attenuation Coefficients, material
 *           constants and elemental densities used in photon-interaction
 *           evaluations:
 *           https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
 *
 * [NIST-KAPTON] NIST STAR material composition for Kapton polyimide film,
 *           density 1.420 g/cm3 and H/C/N/O mass fractions:
 *           https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=179
 *
 * [G4-PE]   Geant4/NIST material database includes G4_POLYETHYLENE with a
 *           polyethylene formula and density near 0.94 g/cm3.  ShieldSim uses
 *           0.95 g/cm3 for the HDPE approximation unless changed below.
 *
 * [FSRI-HDPE] Fire Safety Research Institute material database, measured HDPE
 *           sheet density near 971 kg/m3.  This is used as a sanity check for
 *           the HDPE-density scale, not as a unique mission-material standard:
 *           https://materials.fsri.org/materialdetail/high-density-polyethylene-hdpe
 *
 * [DUPONT-KEVLAR] DuPont Kevlar Aramid Fiber Technical Guide, typical Kevlar
 *           fiber density 1.44 g/cm3.  The polymer repeat unit is represented
 *           here as C14H10N2O2 for transport-element composition.
 *
 * [SIC]     PubChem / engineering handbooks give SiC density about 3.21 g/cm3.
 *           ShieldSim uses stoichiometric SiC for the ceramic filler phase.
 *
 * [LUNAR-SOURCEBOOK] Lunar Sourcebook chapters on lunar regolith composition
 *           and physical properties; lunar soils vary substantially with site,
 *           maturity, compaction, and grain-size distribution:
 *           https://www.lpi.usra.edu/publications/books/lunar_sourcebook/
 *
 * [MARS-APXS] Mars Pathfinder / MER / Curiosity APXS literature gives in-situ
 *           bulk chemistry of Martian soils and rocks.  ShieldSim's default is
 *           only a rough basaltic/regolith elemental mixture.
 *
 * How to add a new material
 * -------------------------
 * See the long comment in include/MaterialCatalog.hh.  In short:
 *   - add an entry in ShieldMaterialCatalog();
 *   - if it is a G4_* material, no builder is needed;
 *   - if it is a new custom material, implement BuildNewMaterial() and register
 *     it in BuildCustomMaterial();
 *   - keep comments, density units, fraction conventions, and references with
 *     the builder so future users can audit the model.
 * ========================================================================== */

#include <G4Element.hh>
#include <G4Exception.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4SystemOfUnits.hh>

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

namespace {

std::string NormalizeToken(const std::string& s){
  // Convert user input to a forgiving lookup token.  The user may type any of
  // these and get the same material:
  //   "HDPE", "high-density polyethylene", "High Density Polyethylene".
  // Only letters and digits are kept; spaces, slashes, hyphens, and underscores
  // are ignored.  The catalog can therefore support user-friendly names without
  // duplicating many exact spellings in the parser.
  std::string out;
  out.reserve(s.size());
  for(unsigned char ch : s){
    if(std::isalnum(ch)) out.push_back(static_cast<char>(std::tolower(ch)));
  }
  return out;
}

bool StartsWithG4(const std::string& s){
  // Geant4/NIST material names have the conventional prefix G4_.  We pass these
  // through directly so users can select materials that are not in the shorter
  // shieldSim catalog, e.g. G4_Pb, G4_STAINLESS-STEEL, or G4_AIR.
  return s.size()>=3 && s[0]=='G' && s[1]=='4' && s[2]=='_';
}

G4Element* El(const char* symbol){
  // Convenience wrapper for the Geant4/NIST element builder.  Elements are
  // cached internally by Geant4, so repeated calls are cheap and safe.
  return G4NistManager::Instance()->FindOrBuildElement(symbol);
}

G4Material* Existing(const char* name){
  // Return a previously-created material without raising a Geant4 exception if
  // it is absent.  This is essential in sweep mode because the geometry can be
  // rebuilt many times within one executable run.
  return G4Material::GetMaterial(name,false);
}

G4Material* BuildHDPE(){
  const char* name="ShieldSim_HDPE";
  if(auto* m=Existing(name)) return m;

  // High-density polyethylene approximation.
  //   Composition convention: atom counts in the repeat unit, CH2.
  //   Density: 0.95 g/cm3, representative of polyethylene/HDPE shielding
  //            stock and close to the Geant4 G4_POLYETHYLENE and measured HDPE
  //            density scale.  Replace with supplier data when available.
  auto* m = new G4Material(name,0.95*g/cm3,2);
  m->AddElement(El("C"),1);
  m->AddElement(El("H"),2);
  return m;
}

G4Material* BuildBPE(){
  const char* name="ShieldSim_BPE_5pctB";
  if(auto* m=Existing(name)) return m;

  // Borated polyethylene approximation.
  //   Composition convention: mass fractions.
  //   Default: 95 wt% of the HDPE approximation above + 5 wt% natural boron.
  //   Density: engineering placeholder 0.98 g/cm3.  Commercial borated-PE
  //            density depends on boron loading, polymer grade, and processing.
  //            Use measured supplier density for mission-quality analysis.
  auto* hdpe = BuildHDPE();
  auto* m = new G4Material(name,0.98*g/cm3,2);
  m->AddMaterial(hdpe,0.95);
  m->AddElement(El("B"),0.05);
  return m;
}

G4Material* BuildKevlar(){
  const char* name="ShieldSim_Kevlar";
  if(auto* m=Existing(name)) return m;

  // Para-aramid / Kevlar approximation.
  //   Composition convention: atom counts in the polymer repeat unit
  //                           C14 H10 N2 O2.
  //   Density: 1.44 g/cm3, typical Kevlar fiber density from technical guides.
  //   Note: woven Kevlar cloth, resin-impregnated laminates, or aramid composites
  //         can have different effective densities and should be defined as
  //         custom project-specific mixtures if needed.
  auto* m = new G4Material(name,1.44*g/cm3,4);
  m->AddElement(El("C"),14);
  m->AddElement(El("H"),10);
  m->AddElement(El("N"), 2);
  m->AddElement(El("O"), 2);
  return m;
}

G4Material* BuildEpoxy(){
  const char* name="ShieldSim_Epoxy";
  if(auto* m=Existing(name)) return m;

  // Generic epoxy-resin approximation used only as a matrix phase for simple
  // composite examples below.
  //   Composition convention: atom counts in a representative epoxy-like repeat
  //                           unit C21 H25 O5.
  //   Density: 1.20 g/cm3 engineering placeholder.
  //   Replace this with the actual resin formulation if the composite material
  //   is important to a quantitative result.
  auto* m = new G4Material(name,1.20*g/cm3,3);
  m->AddElement(El("C"),21);
  m->AddElement(El("H"),25);
  m->AddElement(El("O"), 5);
  return m;
}

G4Material* BuildCFRP(){
  const char* name="ShieldSim_CFRP";
  if(auto* m=Existing(name)) return m;

  // Carbon fiber / epoxy composite approximation.
  //   Composition convention: mass fractions.
  //   Default: 60 wt% elemental carbon fiber + 40 wt% generic epoxy matrix.
  //   Density: 1.60 g/cm3 engineering placeholder.
  //   This is a first-order shielding trade-study material.  Real CFRP depends
  //   on fiber volume fraction, layup, void fraction, resin, cure, and additives.
  auto* epoxy = BuildEpoxy();
  auto* m = new G4Material(name,1.60*g/cm3,2);
  m->AddElement(El("C"),0.60);
  m->AddMaterial(epoxy,0.40);
  return m;
}

G4Material* BuildSiC(){
  const char* name="ShieldSim_SiC";
  if(auto* m=Existing(name)) return m;

  // Silicon carbide filler/ceramic phase.
  //   Composition convention: atom counts, SiC.
  //   Density: 3.21 g/cm3 for dense SiC.  This helper is not exposed directly in
  //            the CLI; it is used by BuildSiCCompositePlastic().
  auto* m = new G4Material(name,3.21*g/cm3,2);
  m->AddElement(El("Si"),1);
  m->AddElement(El("C"), 1);
  return m;
}

G4Material* BuildSiCCompositePlastic(){
  const char* name="ShieldSim_SiCCompositePlastic";
  if(auto* m=Existing(name)) return m;

  // SiC-filled plastic composite approximation.
  //   Composition convention: mass fractions.
  //   Default: 70 wt% dense SiC + 30 wt% generic epoxy-like polymer matrix.
  //   Density: 2.20 g/cm3 engineering placeholder.
  //   Change both the SiC loading and density for a specific plastic composite.
  auto* sic   = BuildSiC();
  auto* epoxy = BuildEpoxy();
  auto* m = new G4Material(name,2.20*g/cm3,2);
  m->AddMaterial(sic,  0.70);
  m->AddMaterial(epoxy,0.30);
  return m;
}

G4Material* BuildLunarRegolith(){
  const char* name="ShieldSim_LunarRegolith";
  if(auto* m=Existing(name)) return m;

  // Approximate bulk lunar-regolith elemental composition.
  //   Composition convention: mass fractions; values sum to 1.0.
  //   Density: 1.60 g/cm3 representative loose/bulk density for transport trade
  //            studies.  Actual lunar regolith density varies strongly with
  //            depth, compaction, porosity, mare/highland provenance, glass
  //            content, agglutinates, and grain size.
  //   Replace this definition with a sample-specific or simulant-specific
  //   composition whenever the shielding result depends on regolith chemistry.
  auto* m = new G4Material(name,1.60*g/cm3,9);
  m->AddElement(El("O" ),0.430);
  m->AddElement(El("Si"),0.210);
  m->AddElement(El("Al"),0.100);
  m->AddElement(El("Ca"),0.080);
  m->AddElement(El("Mg"),0.060);
  m->AddElement(El("Fe"),0.090);
  m->AddElement(El("Ti"),0.020);
  m->AddElement(El("Na"),0.005);
  m->AddElement(El("K" ),0.005);
  return m;
}

G4Material* BuildMarsRegolith(){
  const char* name="ShieldSim_MarsRegolith";
  if(auto* m=Existing(name)) return m;

  // Approximate bulk Martian-regolith/soil elemental composition.
  //   Composition convention: mass fractions; values sum to 1.0.
  //   Density: 1.50 g/cm3 loose-material placeholder.  Real Martian soils and
  //            rocks vary by site, grain-size distribution, cementation,
  //            hydration/sulfate/chloride content, and compaction state.
  //   Replace this definition with APXS-derived site-specific chemistry or a
  //   measured simulant formulation for quantitative studies.
  auto* m = new G4Material(name,1.50*g/cm3,9);
  m->AddElement(El("O" ),0.440);
  m->AddElement(El("Si"),0.210);
  m->AddElement(El("Fe"),0.140);
  m->AddElement(El("Mg"),0.060);
  m->AddElement(El("Al"),0.050);
  m->AddElement(El("Ca"),0.050);
  m->AddElement(El("S" ),0.030);
  m->AddElement(El("Na"),0.010);
  m->AddElement(El("K" ),0.010);
  return m;
}

G4Material* BuildCustomMaterial(const std::string& canonical){
  // Dispatch table for ShieldSim_* materials.  Every custom material listed in
  // ShieldMaterialCatalog() must appear here.  If it does not, the material will
  // be advertised in --help but will fail at geometry construction.
  if(canonical=="ShieldSim_HDPE")                 return BuildHDPE();
  if(canonical=="ShieldSim_BPE_5pctB")            return BuildBPE();
  if(canonical=="ShieldSim_Kevlar")               return BuildKevlar();
  if(canonical=="ShieldSim_CFRP")                 return BuildCFRP();
  if(canonical=="ShieldSim_SiCCompositePlastic")  return BuildSiCCompositePlastic();
  if(canonical=="ShieldSim_LunarRegolith")        return BuildLunarRegolith();
  if(canonical=="ShieldSim_MarsRegolith")         return BuildMarsRegolith();
  return nullptr;
}

std::map<std::string,std::string> MakeAliasMap(){
  // Build a normalized alias -> canonical-name map from the catalog itself.
  // This keeps the parser data-driven: adding an alias to a catalog entry is
  // enough to make the CLI accept it.
  std::map<std::string,std::string> out;
  for(const auto& e : ShieldMaterialCatalog()){
    out[NormalizeToken(e.key)] = e.canonicalName;
    out[NormalizeToken(e.displayName)] = e.canonicalName;
    out[NormalizeToken(e.canonicalName)] = e.canonicalName;
    for(const auto& a : e.aliases) out[NormalizeToken(a)] = e.canonicalName;
  }
  return out;
}

const MaterialCatalogEntry* FindCatalogEntryByCanonical(const std::string& canonical){
  for(const auto& e : ShieldMaterialCatalog())
    if(e.canonicalName==canonical) return &e;
  return nullptr;
}

} // namespace

const std::vector<MaterialCatalogEntry>& ShieldMaterialCatalog(){
  // The order below controls the order printed in --help and --list-materials.
  // Keep the list close to the user-interface grouping shown in the material
  // selection mock-up: metals, hydrogenous polymers, composites, water/regolith.
  static const std::vector<MaterialCatalogEntry> entries = {
    {"Structural Metals", "Al", "Aluminum (Al)", "G4_Al",
     "Geant4/NIST elemental aluminum; density from installed G4 material database.",
     "Refs: [G4-MAT], [NIST-XCOM]", {"Aluminum"}, false},
    {"Structural Metals", "Cu", "Copper (Cu)", "G4_Cu",
     "Geant4/NIST elemental copper; density from installed G4 material database.",
     "Refs: [G4-MAT], [NIST-XCOM]", {"Copper"}, false},
    {"Structural Metals", "W", "Tungsten (W)", "G4_W",
     "Geant4/NIST elemental tungsten; density from installed G4 material database.",
     "Refs: [G4-MAT], [NIST-XCOM]", {"Tungsten"}, false},
    {"Structural Metals", "Ta", "Tantalum (Ta)", "G4_Ta",
     "Geant4/NIST elemental tantalum; density from installed G4 material database.",
     "Refs: [G4-MAT], [NIST-XCOM]", {"Tantalum"}, false},

    {"Hydrogenous Polymers", "HDPE", "HDPE - High-Density Polyethylene", "ShieldSim_HDPE",
     "Custom polyethylene approximation: CH2 repeat unit, density 0.95 g/cm3.",
     "Refs: [G4-PE], [FSRI-HDPE]", {"Polyethylene","HighDensityPolyethylene","High-Density Polyethylene","PE"}, true},
    {"Hydrogenous Polymers", "BPE", "BPE - Borated Polyethylene (5% B)", "ShieldSim_BPE_5pctB",
     "Custom mixture: 95 wt% ShieldSim_HDPE + 5 wt% natural boron, density 0.98 g/cm3.",
     "Refs: [G4-PE], [FSRI-HDPE]; 5 wt% B is the selected catalog definition", {"BoratedPolyethylene","BPE5","BPE5B","Borated Polyethylene 5 B"}, true},
    {"Hydrogenous Polymers", "Kevlar", "Kevlar / Aramid", "ShieldSim_Kevlar",
     "Custom para-aramid approximation: C14H10N2O2 repeat unit, density 1.44 g/cm3.",
     "Ref: [DUPONT-KEVLAR]", {"Aramid","KevlarAramid","ParaAramid"}, true},
    {"Hydrogenous Polymers", "Kapton", "Polyimide / Kapton", "G4_KAPTON",
     "Geant4/NIST Kapton/polyimide material.",
     "Refs: [G4-MAT], [NIST-KAPTON]", {"Polyimide","G4_KAPTON"}, false},

    {"Composites", "CFRP", "CFRP - Carbon Fiber / Epoxy", "ShieldSim_CFRP",
     "Custom trade-study composite: 60 wt% carbon fiber + 40 wt% epoxy, density 1.60 g/cm3.",
     "Engineering placeholder; replace with measured layup/resin data", {"CarbonFiberEpoxy","CarbonFiber","Carbon Fiber Epoxy"}, true},
    {"Composites", "SiCComposite", "SiC Composite Plastic", "ShieldSim_SiCCompositePlastic",
     "Custom trade-study composite: 70 wt% SiC + 30 wt% epoxy, density 2.20 g/cm3.",
     "Refs: [SIC] for SiC phase; matrix/loading are catalog assumptions", {"SiCPlastic","SiCCompositePlastic","SiC Composite"}, true},

    {"Water / Regolith", "Water", "Water (H2O)", "G4_WATER",
     "Geant4/NIST liquid water.",
     "Refs: [G4-MAT], [NIST-XCOM]", {"H2O","G4_WATER"}, false},
    {"Water / Regolith", "LunarRegolith", "Lunar Regolith", "ShieldSim_LunarRegolith",
     "Custom rough bulk lunar-regolith mixture, density 1.60 g/cm3; mass fractions in source.",
     "Ref: [LUNAR-SOURCEBOOK]; site/sample dependent", {"MoonRegolith","Regolith","Lunar/MarsRegolith","Lunar Mars Regolith"}, true},
    {"Water / Regolith", "MarsRegolith", "Mars Regolith", "ShieldSim_MarsRegolith",
     "Custom rough bulk Martian-regolith mixture, density 1.50 g/cm3; mass fractions in source.",
     "Ref: [MARS-APXS]; site/sample dependent", {"MartianRegolith","MarsSoil","MartianSoil"}, true}
  };
  return entries;
}

std::string ResolveShieldMaterialName(const std::string& userName){
  if(userName.empty()) return userName;
  if(StartsWithG4(userName)) return userName;

  static const auto aliasMap = MakeAliasMap();
  const auto key = NormalizeToken(userName);
  const auto it = aliasMap.find(key);
  if(it!=aliasMap.end()) return it->second;

  // Unknown names are passed through.  This preserves a useful extension path:
  // future code can create a G4Material elsewhere before DetectorConstruction is
  // built, then pass that material name through the CLI.
  return userName;
}

G4Material* FindOrBuildShieldMaterial(const std::string& userName){
  const std::string canonical = ResolveShieldMaterialName(userName);

  // First check whether Geant4 already has the material in its global material
  // table.  This catches custom materials created earlier in this run and any
  // NIST materials already requested by another part of the application.
  if(auto* m=G4Material::GetMaterial(canonical,false)) return m;

  // Then try shieldSim custom material builders.
  if(auto* m=BuildCustomMaterial(canonical)) return m;

  // Finally, ask the Geant4/NIST material manager.  This is where ordinary G4_*
  // material names are resolved.
  auto* nist=G4NistManager::Instance();
  auto* m=nist->FindOrBuildMaterial(canonical,false);
  if(!m){
    const std::string msg = "Material '" + userName + "' could not be resolved. "
                          + "Use --list-materials to see built-in aliases, or "
                          + "provide a valid Geant4 NIST material name such as G4_Al.";
    G4Exception("FindOrBuildShieldMaterial","BadMaterial",FatalException,msg.c_str());
  }
  return m;
}

std::string ShieldMaterialCatalogText(){
  std::ostringstream out;
  std::string currentCategory;

  out<<"\n  Allowed shielding-material catalog keys\n";
  out<<"  ---------------------------------------\n";
  for(const auto& e : ShieldMaterialCatalog()){
    if(e.category!=currentCategory){
      currentCategory=e.category;
      out<<"\n  "<<currentCategory<<"\n";
    }
    out<<"    "<<std::left<<std::setw(16)<<e.key
       <<" -> "<<std::setw(32)<<e.displayName
       <<" canonical="<<std::setw(32)<<e.canonicalName
       <<(e.isCustom?" custom":" G4/NIST")<<"\n";
    out<<"       "<<e.description<<"\n";
    out<<"       "<<e.sourceNote<<"\n";
    if(!e.aliases.empty()){
      out<<"       aliases: ";
      for(std::size_t i=0;i<e.aliases.size();++i){
        if(i) out<<", ";
        out<<e.aliases[i];
      }
      out<<"\n";
    }
  }

  out<<"\n  Any Geant4/NIST material name beginning with G4_ is also accepted,\n";
  out<<"  for example G4_Al, G4_Cu, G4_WATER, G4_Si, G4_Pb, or G4_POLYETHYLENE.\n";

  out<<"\n  Material-property reference notes\n";
  out<<"  ---------------------------------\n";
  out<<"    [G4-MAT]          Geant4 material database / installed Geant4 NIST definitions.\n";
  out<<"    [NIST-XCOM]       NIST material constants and elemental-density tables.\n";
  out<<"    [NIST-KAPTON]     NIST STAR Kapton composition/density.\n";
  out<<"    [G4-PE]           Geant4 G4_POLYETHYLENE reference material.\n";
  out<<"    [FSRI-HDPE]       Measured HDPE density scale from FSRI material database.\n";
  out<<"    [DUPONT-KEVLAR]   Kevlar technical-guide density scale.\n";
  out<<"    [SIC]             SiC stoichiometry/density from materials references.\n";
  out<<"    [LUNAR-SOURCEBOOK] Lunar Sourcebook regolith composition/density context.\n";
  out<<"    [MARS-APXS]       Mars APXS literature for Martian soil/regolith chemistry.\n";
  out<<"\n  Full URLs and implementation details are documented in src/MaterialCatalog.cc.\n";
  out<<"  Custom composite/regolith definitions are approximate trade-study materials;\n";
  out<<"  replace them with measured project-specific material data for final analyses.\n";
  out<<"\n";
  return out.str();
}

std::string DescribeShieldMaterial(const std::string& userName){
  const std::string canonical = ResolveShieldMaterialName(userName);
  if(const auto* e = FindCatalogEntryByCanonical(canonical)){
    return e->displayName + " [" + e->key + "; " + e->canonicalName + "]";
  }
  return userName;
}
