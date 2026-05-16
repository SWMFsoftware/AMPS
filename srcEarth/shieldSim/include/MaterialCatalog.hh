#ifndef SHIELDSIM_MATERIAL_CATALOG_HH
#define SHIELDSIM_MATERIAL_CATALOG_HH

/* ============================================================================
 * MaterialCatalog.hh
 *
 * Centralized shielding-material catalog for shieldSim
 * ---------------------------------------------------
 * This header defines the small, user-facing material database used by the
 * shieldSim command line.  The purpose of the catalog is not to replace the
 * Geant4/NIST material database.  Instead, it gives the user short and stable
 * shielding-material names such as
 *
 *     Al, HDPE, BPE, Kevlar, CFRP, LunarRegolith
 *
 * while still allowing any ordinary Geant4/NIST material name beginning with
 * "G4_", for example G4_Al, G4_WATER, G4_Si, or G4_Pb.
 *
 * Why this catalog exists
 * -----------------------
 * Geant4 already contains many standard NIST materials.  However, for a
 * shielding trade study it is useful to expose a smaller list of engineering
 * materials in the CLI and to define a few approximate materials that are not
 * guaranteed to exist as a predefined Geant4 entry, for example borated
 * polyethylene, carbon-fiber/epoxy, and regolith-like mixtures.
 *
 * Where material properties are defined
 * -------------------------------------
 * The actual material definitions are implemented in src/MaterialCatalog.cc.
 * Some catalog entries are simple aliases to Geant4/NIST materials; those are
 * marked with isCustom=false and their canonicalName is a G4_* material name.
 * Other entries are built explicitly by shieldSim; those are marked with
 * isCustom=true and have canonical names beginning with ShieldSim_*.
 *
 * IMPORTANT: custom entries are approximate trade-study materials.  For a final
 * mission or instrument analysis, replace density, composition, filler fraction,
 * composite layup, and regolith chemistry with measured project-specific values.
 *
 * How to add a new shielding material
 * -----------------------------------
 * Case A: the material already exists in Geant4/NIST
 *   1. Add a new MaterialCatalogEntry in ShieldMaterialCatalog() in
 *      src/MaterialCatalog.cc.
 *   2. Set canonicalName to the Geant4 material name, e.g. "G4_Pb".
 *   3. Set isCustom=false.
 *   4. Add aliases that users might type on the command line.
 *   5. Add a clear description and a source/reference note.
 *
 * Case B: the material must be built by shieldSim
 *   1. Add a new MaterialCatalogEntry in ShieldMaterialCatalog() with a unique
 *      canonicalName beginning with "ShieldSim_" and isCustom=true.
 *   2. Implement a builder in src/MaterialCatalog.cc, following the pattern of
 *      BuildHDPE(), BuildBPE(), BuildKevlar(), or BuildLunarRegolith().
 *      Use Geant4 units explicitly, e.g. 1.23*g/cm3.
 *   3. At the beginning of the builder, return an existing material if it has
 *      already been created:
 *
 *          const char* name = "ShieldSim_NewMaterial";
 *          if (auto* m = Existing(name)) return m;
 *
 *      This avoids duplicate material definitions when geometry is rebuilt in
 *      sweep mode.
 *   4. Add the builder to BuildCustomMaterial().
 *   5. For compounds or polymers, use AddElement(element, nAtoms).  For mixtures
 *      and composites, use AddElement/AddMaterial with mass fractions that sum
 *      to 1.0.
 *   6. Document the density, composition convention, and references in the
 *      comment above the builder and in the catalog entry.  The --help and
 *      --list-materials text is generated from the catalog, so no separate help
 *      table needs to be edited.
 *
 * Common pitfalls when adding materials
 * -------------------------------------
 *   - Do not omit Geant4 units for density or thickness.
 *   - Do not create the same G4Material twice under the same name.
 *   - Be explicit about whether fractions are mass fractions or atom counts.
 *   - Avoid pretending that engineering composites/regolith are unique materials;
 *     their properties are sample- and process-dependent.
 * ========================================================================== */

#include <string>
#include <vector>

class G4Material;

struct MaterialCatalogEntry {
  std::string category;       // UI grouping, e.g. "Structural Metals".
  std::string key;            // Recommended CLI key, e.g. "Al" or "HDPE".
  std::string displayName;    // Human-readable label shown in help text.
  std::string canonicalName;  // Geant4/NIST name or ShieldSim_* custom name.
  std::string description;    // Density/composition/meaning note.
  std::string sourceNote;     // Short source/reference note shown in help text.
  std::vector<std::string> aliases; // Alternative CLI spellings.
  bool isCustom = false;      // True if shieldSim builds the material itself.
};

// Full material catalog in the order used by --help, --list-materials, and the
// README.  Keeping one authoritative list avoids inconsistencies between the
// CLI parser, geometry builder, and documentation.
const std::vector<MaterialCatalogEntry>& ShieldMaterialCatalog();

// Resolve a user-provided name or alias to the canonical material name.  Names
// beginning with "G4_" are passed through unchanged so arbitrary Geant4/NIST
// materials remain supported.
std::string ResolveShieldMaterialName(const std::string& userName);

// Build or return the Geant4 material corresponding to a catalog key, alias, or
// ordinary NIST material name.  This function is used for both shield and
// scoring-slab materials so custom catalog materials work anywhere a material
// name is accepted.
G4Material* FindOrBuildShieldMaterial(const std::string& userName);

// User-facing table for --help and --list-materials.  The returned string lists
// allowed catalog keys, aliases, material descriptions, and reference notes.
std::string ShieldMaterialCatalogText();

// One-line material description for configuration printouts and metadata.
std::string DescribeShieldMaterial(const std::string& userName);

#endif // SHIELDSIM_MATERIAL_CATALOG_HH
