#ifndef SHIELDSIM_RUN_ACTION_HH
#define SHIELDSIM_RUN_ACTION_HH

/* ============================================================================
 * RunAction.hh
 *
 * Run-level scoring, normalization, and output for shieldSim.
 *
 * SteppingAction adds separate total/non-ionizing/ionizing deposition ledgers
 * and transmitted-particle counts. Versioned runs retain accumulators across
 * adaptive BeamOn batches, evaluate event-level precision, and emit a canonical
 * result JSON as well as compatibility diagnostic/Tecplot files.
 * ========================================================================== */

#include "ShieldSimConfig.hh"

#include <G4UserRunAction.hh>
#include <G4Types.hh>

#include <fstream>
#include <string>
#include <vector>
#include <memory>

class DetectorConstruction;
class G4LogicalVolume;
class G4Run;
namespace shield_sim { class AccuracyController; }

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

  // Diagnostic sample recorders used by automated tests.  These methods are
  // intentionally run-action methods rather than standalone global streams so
  // that all diagnostic files are opened/closed consistently with each run.
  // Positions are written in mm, directions are unit vectors, and energies are
  // MeV total kinetic energy per particle.
  void RecordSourceSample(const std::string& species,
                          G4double energyMeV,
                          G4double xMM, G4double yMM, G4double zMM,
                          G4double ux, G4double uy, G4double uz);
  void RecordExitParticle(const std::string& species,
                          G4double energyMeV,
                          G4double xGlobalMM, G4double yGlobalMM, G4double zGlobalMM,
                          G4double xLocalMM,  G4double yLocalMM,  G4double zLocalMM,
                          G4double uxLocal,   G4double uyLocal,   G4double uzLocal);

  // Separate ledgers are mandatory: TID is ionizing deposited energy/mass,
  // while total and non-ionizing deposit remain available for closure audits.
  // Geant4 reports total and non-ionizing deposit per step; ionizing deposit is
  // their non-negative difference.
  void AddEdep(std::size_t i,G4double totalEdep,G4double nonIonizingEdep);
  void AddPrimaryProtonExit(G4double kineticEnergyMeV);
  void BeginEvent();
  void EndEvent();
  bool AccuracyStatisticallySatisfied() const;
  G4long GetCumulativeEvents() const;

  void SetSourceNormalization(G4double normNoAngular,
                              G4double angularFactor,
                              const std::string& sourceMode);

  G4double GetSourcePlaneAreaCM2() const;
  G4double GetIncidentParticleRate() const;

  void SetSweepMode(bool on,const std::string& mat,G4double tMM);
  void ConfigureIndependentCase(std::size_t absorberIndex,
                                const std::string& outputPrefix,
                                const std::string& scoringName,
                                G4double scoringThickness);
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
  std::string OutputName(const std::string& suffix) const;
  void OpenDiagnosticFiles();
  void CloseDiagnosticFiles();
  bool DiagnosticLimitReached(G4long rowsWritten) const;
  void WriteSpectraTecplot(G4int nEv);
  void ComputeSpectrumFoldedQuantities(G4int nEv);
  void WriteComputedQuantitiesTecplot(G4int nEv);
  void WriteLETSpectrumTecplot(G4int nEv);
  void WriteRunSummary(G4int nEv);
  void WriteResultJson(G4int nEv);
  G4double ComputeHardnessIndex(G4int nEv) const;

  Options                     fOpts;
  DetectorConstruction*       fDetector=nullptr;
  std::vector<std::string>    fScoringNames;
  std::vector<G4double>       fScoringThick;

  std::vector<G4LogicalVolume*> fScoringLVs;
  std::vector<G4LogicalVolume*> fShieldLVs;
  G4LogicalVolume*            fShieldLV=nullptr;

  std::vector<G4double>       fEdepTotal;
  std::vector<G4double>       fEdepIonizing;
  std::vector<G4double>       fEdepNonIonizing;
  std::vector<G4double>       fLastDoseTotal;
  std::vector<G4double>       fLastDoseNonIonizing;
  std::vector<G4double>       fLastDose;
  std::vector<G4double>       fLastDoseRate;
  std::vector<G4double>       fLastDDD;
  std::vector<G4double>       fLastDDDRate;
  std::vector<G4double>       fLastNEq;
  std::vector<G4double>       fLastNEqRate;
  G4double                    fLastHardness=0.0;
  G4double                    fPrimaryProtonAbove10=0.0;
  G4double                    fPrimaryProtonAbove100=0.0;
  std::vector<G4double>       fEventEdepIonizing;
  bool                        fEventPrimaryProtonAbove10=false;
  bool                        fEventPrimaryProtonAbove100=false;
  std::unique_ptr<shield_sim::AccuracyController> fAccuracyController;
  G4long                      fCumulativeEvents=0;
  bool                        fCumulativeMode=false;
  bool                        fAccumulationStarted=false;
  std::vector<G4double>       fInP,fInA,fOutP,fOutA,fOutN;

  G4double    fSourceNormNoAngular=0.0;
  G4double    fSourceAngularFactor=1.0;
  G4double    fSourceNorm=0.0;
  std::string fSourceMode="beam";

  bool        fSweepMode=false;
  std::string fCurrentMat;
  G4double    fCurrentTmm=0;

  // Optional Layer-2 diagnostic streams.  The source file records what the
  // generator requested from Geant4.  The exit file records only particles that
  // the stepping action accepted as crossing the downstream shield face in the
  // shield-local coordinate system.
  std::ofstream fSourceDump;
  std::ofstream fExitDump;
  G4long        fSourceDumpRows=0;
  G4long        fExitDumpRows=0;

  bool        fFirstSpectraWrite=true;
  bool        fFirstQuantitiesWrite=true;
  bool        fFirstLETWrite=true;
  bool        fFirstSummaryWrite=true;
  std::vector<SweepPoint> fSweepData;
};

#endif // SHIELDSIM_RUN_ACTION_HH
