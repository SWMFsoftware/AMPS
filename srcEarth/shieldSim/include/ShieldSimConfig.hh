#ifndef SHIELDSIM_CONFIG_HH
#define SHIELDSIM_CONFIG_HH

/* ============================================================================
 * ShieldSimConfig.hh
 *
 * Shared run configuration and simple result containers for shieldSim.
 *
 * All command-line options are parsed once at application startup and stored in
 * Options.  The same Options object, or a local copy of it, is then passed to
 * DetectorConstruction, PrimaryGeneratorAction, and RunAction so that the whole
 * application uses one consistent configuration.
 *
 * User-facing length values in the command line are specified in mm and are
 * converted immediately to Geant4 internal length units during parsing.
 * User-facing energy values are specified in MeV total kinetic energy per
 * particle.  For alpha particles this means total alpha energy, not MeV/n.
 * ========================================================================== */

#include <G4Types.hh>
#include <G4SystemOfUnits.hh>

#include <string>
#include <utility>
#include <vector>

struct Options {
  // ---- physics configuration ----------------------------------------------
  // Geant4 reference physics list created through G4PhysListFactory.
  // Supported CLI values are currently:
  //   FTFP_BERT, FTFP_BERT_HP, Shielding, QGSP_BIC_HP
  std::string physicsList = "FTFP_BERT";

  // ---- source and spectrum configuration ----------------------------------
  // Empty spectrumFile selects the built-in approximate GCR spectrum.
  std::string spectrumFile = "";

  // Supported source modes:
  //   beam      : normal-incidence pencil beam along +z
  //   isotropic : finite upstream plane source with cosine-law directions
  std::string sourceMode   = "beam";

  // ---- energy limits for sampling [MeV total kinetic energy] ---------------
  G4double eMinProton =  10.0;
  G4double eMaxProton = 1.0e5;
  G4double eMinAlpha  =  10.0;
  G4double eMaxAlpha  = 1.0e5;

  // ---- single-run shield geometry -----------------------------------------
  // Material can be either a Geant4/NIST material name (for example G4_Al)
  // or a user-friendly catalog key/alias (for example Al, HDPE, BPE, CFRP,
  // Water, LunarRegolith).  See MaterialCatalog.hh for the built-in list.
  std::string shieldMaterial  = "Al";
  G4double    shieldThickness = 2.0*mm;

  // ---- downstream scoring slabs -------------------------------------------
  // Each entry is (NIST material name, slab thickness in Geant4 length units).
  std::vector<std::pair<std::string,G4double>> scoringMaterials;


  // ---- post-processed computed quantities ---------------------------------
  // The Geant4 transport always scores energy deposition and transmitted
  // spectra.  These switches control which derived radiation-effect quantities
  // are written at the end of each run.  Defaults compute all quantities shown
  // in the UI: TID, DDD, n_eq, LET spectrum, and H100/10.
  bool calcTID      = true;
  bool calcDDD      = true;
  bool calcNEq      = true;
  bool calcLET      = true;
  bool calcHardness = true;

  // ---- statistics ----------------------------------------------------------
  G4int nEvents = 10000;

  // ---- dose-vs-thickness sweep --------------------------------------------
  bool        doSweep       = false;
  std::string sweepMaterial = "";    // default: same as shieldMaterial
  G4double    sweepTmin     = 0.5;   // mm, still user-facing here
  G4double    sweepTmax     = 50.0;  // mm, still user-facing here
  G4int       sweepN        = 10;
  bool        sweepLog      = false;

  bool showHelp = false;
  bool listMaterials = false;       // print shielding-material catalog and exit
  bool listTargetMaterials = false; // print detector/absorber/target-material catalog and exit
  bool listQuantities = false;      // print computed-quantity definitions and exit
};

// One row in the dose-vs-thickness sweep output table.
struct SweepPoint {
  G4double thickMM;                       // shield thickness [mm]
  G4double arealDensity;                  // shield areal density [g/cm^2]
  std::vector<G4double> dose_Gy;          // dose per primary, in G4 dose units
  std::vector<G4double> doseRate_Gy_s;    // source-normalized ionizing dose rate

  // Spectrum-folded radiation-effect outputs.  DDD is accumulated in MeV/g
  // because NIEL response functions are conventionally in MeV cm^2/g and
  // fluence is in cm^-2.  The rad-equivalent conversion is written at output
  // time.  n_eq is the equivalent 1-MeV-neutron fluence in cm^-2.
  std::vector<G4double> ddd_MeV_g_perPrimary;
  std::vector<G4double> dddRate_MeV_g_s;
  std::vector<G4double> neq_cm2_perPrimary;
  std::vector<G4double> neqRate_cm2_s;

  // H100/10 is independent of the target material, but it is repeated in the
  // per-target output rows so each shielding x absorber row is self-contained.
  G4double hardness_H100_10 = 0.0;
};

#endif // SHIELDSIM_CONFIG_HH
