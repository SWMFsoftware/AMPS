#ifndef SHIELDSIM_RUN_ACTION_HH
#define SHIELDSIM_RUN_ACTION_HH

/* ============================================================================
 * RunAction.hh
 *
 * Run-level scoring, normalization, and output for shieldSim.
 *
 * SteppingAction adds energy deposition and transmitted-particle counts during
 * tracking.  RunAction resets these accumulators at the beginning of each run,
 * computes dose and source-normalized rates at the end of the run, and writes
 * the diagnostic CSV and Tecplot spectra files.
 * ========================================================================== */

#include "ShieldSimConfig.hh"

#include <G4UserRunAction.hh>
#include <G4Types.hh>

#include <string>
#include <vector>

class DetectorConstruction;
class G4LogicalVolume;
class G4Run;

class RunAction : public G4UserRunAction {
public:
  RunAction(const Options& opts,
            DetectorConstruction* detector,
            const std::vector<std::string>& scoringNames,
            const std::vector<G4double>&    scoringThick);
  ~RunAction() override;

  // Spectral counters.  Energies are MeV total kinetic energy.
  void AddInP  (G4double E);
  void AddInA  (G4double E);
  void AddOutP (G4double E);
  void AddOutA (G4double E);
  void AddOutN (G4double E);

  // Energy deposition is accumulated in Geant4 internal energy units.  Dose is
  // computed later as Edep/mass and converted to gray only at output time.
  void AddEdep (std::size_t i,G4double edep);

  void SetSourceNormalization(G4double normNoAngular,
                              G4double angularFactor,
                              const std::string& sourceMode);

  G4double GetSourcePlaneAreaCM2() const;
  G4double GetIncidentParticleRate() const;

  void SetSweepMode(bool on,const std::string& mat,G4double tMM);
  void AppendSweepPoint(const SweepPoint& sp);
  const std::vector<SweepPoint>& GetSweepData() const;
  const std::vector<G4double>&   GetLastDose() const;
  const std::vector<G4double>&   GetLastDoseRate() const;
  const std::vector<G4double>&   GetLastDDD() const;
  const std::vector<G4double>&   GetLastDDDRate() const;
  const std::vector<G4double>&   GetLastNEq() const;
  const std::vector<G4double>&   GetLastNEqRate() const;
  G4double                      GetLastHardness() const;

  void RefreshLVPointers();
  G4LogicalVolume*                     GetShieldLV() const;
  const std::vector<G4LogicalVolume*>& GetScoringLVs() const;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run* run) override;

private:
  static std::string SanitiseName(const std::string& s);
  static std::string FormatMM(G4double t);
  void WriteSpectraTecplot(G4int nEv);
  void ComputeSpectrumFoldedQuantities(G4int nEv);
  void WriteComputedQuantitiesTecplot(G4int nEv);
  void WriteLETSpectrumTecplot(G4int nEv);
  G4double ComputeHardnessIndex(G4int nEv) const;

  Options                     fOpts;
  DetectorConstruction*       fDetector=nullptr;
  std::vector<std::string>    fScoringNames;
  std::vector<G4double>       fScoringThick;

  std::vector<G4LogicalVolume*> fScoringLVs;
  G4LogicalVolume*            fShieldLV=nullptr;

  std::vector<G4double>       fEdep;
  std::vector<G4double>       fLastDose;
  std::vector<G4double>       fLastDoseRate;
  std::vector<G4double>       fLastDDD;
  std::vector<G4double>       fLastDDDRate;
  std::vector<G4double>       fLastNEq;
  std::vector<G4double>       fLastNEqRate;
  G4double                    fLastHardness=0.0;
  std::vector<G4double>       fInP,fInA,fOutP,fOutA,fOutN;

  G4double    fSourceNormNoAngular=0.0;
  G4double    fSourceAngularFactor=1.0;
  G4double    fSourceNorm=0.0;
  std::string fSourceMode="beam";

  bool        fSweepMode=false;
  std::string fCurrentMat;
  G4double    fCurrentTmm=0;
  bool        fFirstSpectraWrite=true;
  bool        fFirstQuantitiesWrite=true;
  bool        fFirstLETWrite=true;
  std::vector<SweepPoint> fSweepData;
};

#endif // SHIELDSIM_RUN_ACTION_HH
