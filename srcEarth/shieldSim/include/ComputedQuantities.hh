#ifndef SHIELDSIM_COMPUTED_QUANTITIES_HH
#define SHIELDSIM_COMPUTED_QUANTITIES_HH

/* ============================================================================
 * ComputedQuantities.hh
 *
 * Post-processing quantities derived from Geant4 scoring and transmitted spectra
 * ---------------------------------------------------------------------------
 * This module contains the small amount of radiation-effects post-processing
 * that is independent of Geant4 tracking itself.  Geant4 transports particles,
 * creates secondaries, and scores ionizing energy deposition in the selected
 * target slabs.  At the end of each run, RunAction uses the helpers declared
 * here to convert those raw Monte Carlo accumulators into the quantities shown
 * in the shieldSim UI:
 *
 *   TID       total ionizing dose in each target material;
 *   DDD       legacy compatibility proxy from analytic NIEL x shield-exit fluence;
 *   n_eq      legacy material proxy derived from aggregate surrogate DDD;
 *   LET       transmitted charged-particle fluence spectrum binned by LET;
 *   H100/10   primary-proton shield-exit fluence(E>100 MeV) divided by
 *             fluence(E>10 MeV), evaluated from exact crossing energies.
 *
 * Important physics boundary between transport and post-processing
 * ----------------------------------------------------------------
 * This module is the compatibility path. TID is scored directly from Geant4
 * deposition; its DDD, n_eq, and LET helpers fold the shield-exit spectrum with
 * analytic or semi-analytic surrogates. Versioned Stage-5 cases do NOT use the
 * DDD/n_eq functions below: RunAction folds exact absorber-entry crossings with
 * the authenticated DamageResponseLibrary and implements neutron-only ASTM E722
 * n_eq. Keeping that boundary explicit prevents legacy output from being
 * mistaken for the enabled Stage-5 observable.
 *
 * References / provenance for equations and conventions
 * -----------------------------------------------------
 * [ICRU-85]      ICRU Report 85a, Fundamental Quantities and Units for Ionizing
 *                Radiation.  Used for dose and LET terminology.
 * [NIST-PSTAR]   NIST PSTAR/ASTAR stopping-power tables.  Used as the reference
 *                standard for charged-particle mass stopping power.  The code
 *                implements a Bethe-Bloch mass-stopping-power approximation for
 *                LET binning so it remains self-contained; replace it with NIST
 *                or Geant4 tabulated stopping powers for precision work.
 * [ASTM-E722]    ASTM E722, standard practice for characterizing neutron fluence
 *                spectra in terms of 1-MeV neutron-equivalent fluence in silicon.
 *                Used for the n_eq convention.
 * [Messenger]    Messenger, Burke, Summers, and related NIEL literature for the
 *                displacement-damage-dose convention D_d = integral Phi(E)
 *                NIEL(E) dE.
 * [SR-NIEL]      Screened-relativistic NIEL literature/databases.  Recommended
 *                source of replacement NIEL(E) tables for quantitative studies.
 *
 * Units used by this module
 * -------------------------
 * Energies are MeV total kinetic energy per particle.  TID is reported in Gy
 * and rad.  DDD is reported internally as MeV/g and additionally converted to
 * rad-equivalent units using 1 rad = 100 erg/g = 6.241509074e7 MeV/g.  The n_eq
 * proxy has units cm^-2.  LET is binned in MeV cm^2/mg, the common space-radiation
 * LET unit; internally the Bethe-Bloch helper returns MeV cm^2/g and divides by
 * 1000 to convert g to mg.
 * ========================================================================== */

#include <G4String.hh>
#include <G4Types.hh>

#include <string>

class G4Material;

namespace ComputedQuantities {

struct Selection {
  bool tid      = true;  // total ionizing dose from Geant4 energy deposition
  bool ddd      = true;  // NIEL-folded displacement-damage-dose proxy
  bool neq      = true;  // 1-MeV neutron-equivalent fluence proxy
  bool let      = true;  // LET spectrum of transmitted charged particles
  bool hardness = true;  // J(>100 MeV)/J(>10 MeV) spectral hardness index
};

// Parse a comma-separated user selection.  Accepted tokens are:
//   all, none, tid, ddd, neq/n_eq, let, hardness/h100/10/h10010
// The parser is case- and punctuation-tolerant so UI-style names are accepted.
Selection ParseSelection(const std::string& text);

// Return a compact human-readable list of enabled quantities.
std::string SelectionSummary(const Selection& s);

// User-facing help text for --help and --list-quantities.
std::string QuantityCatalogText();

// Convert displacement damage energy per unit mass from MeV/g to rad.  This is
// only an energy-density conversion, not a statement that DDD has the same
// biological or device meaning as ionizing dose.
G4double MeVPerGramToRad(G4double mevPerGram);

// Approximate charged-particle LET / electronic mass stopping power in the
// specified target material.  The return value is MeV cm^2/mg.  Neutrons return
// zero because LET is a charged-particle quantity in this simplified output.
G4double ChargedLET_MeV_cm2_mg(const G4String& particleName,
                               G4double kineticEnergyMeV,
                               const G4Material* material);

// Approximate legacy NIEL mass stopping power in MeV cm^2/g. This function is
// used only by compatibility CLI calculations. Production case-v2 Stage-5 DDD
// uses DamageResponseLibrary instead and must never call this surrogate.
G4double NIEL_MeV_cm2_g(const G4String& particleName,
                        G4double kineticEnergyMeV,
                        const G4Material* material);

// Legacy material-dependent reference in MeV cm^2/g. This converts aggregate
// compatibility DDD only; it is not the Stage-5 neutron-only E722 definition:
//   n_eq = DDD / NIEL_1MeV_neutron.
G4double ReferenceNIEL1MeVNeutron_MeV_cm2_g(const G4Material* material);

// Logarithmic LET bins used for shieldSim_let_spectrum.dat.  The range covers
// lightly ionizing high-energy protons through high-LET lower-energy ions for a
// first-order slab-shielding study.
namespace LETBins {
  static const G4int    N    = 100;
  static const G4double Lmin = 1.0e-4; // MeV cm^2/mg
  static const G4double Lmax = 1.0e2;  // MeV cm^2/mg

  G4double Center(G4int i);
  G4double Edge  (G4int i);
  G4double Width (G4int i);
  G4int    Bin   (G4double let_MeV_cm2_mg);
}

} // namespace ComputedQuantities

#endif // SHIELDSIM_COMPUTED_QUANTITIES_HH
