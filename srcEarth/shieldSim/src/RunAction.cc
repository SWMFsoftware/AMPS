#include "RunAction.hh"

#include "AccuracyController.hh"
#include "BuildProvenance.hh"
#include "DetectorConstruction.hh"
#include "CaseDefinition.hh"
#include "ComputedQuantities.hh"
#include "ContentHash.hh"
#include "MaterialCatalog.hh"
#include "SpecBins.hh"
#include "StrictJson.hh"

#include <G4AnalysisManager.hh>
#include <G4Box.hh>
#include <G4Exception.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4Run.hh>
#include <G4SystemOfUnits.hh>
#include <G4ios.hh>
#include <G4Version.hh>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

RunAction::RunAction(const Options& opts,
                     DetectorConstruction* detector,
                     const std::vector<std::string>& scoringNames,
                     const std::vector<G4double>&    scoringThick)
  : fOpts(opts), fDetector(detector), fScoringNames(scoringNames), fScoringThick(scoringThick)
{
  G4int N=SpecBins::N;
  fEdepTotal.assign(scoringNames.size(),0);
  fEdepIonizing.assign(scoringNames.size(),0);
  fEdepNonIonizing.assign(scoringNames.size(),0);
  fEventEdepIonizing.assign(scoringNames.size(),0);
  fEventEdepTotal.assign(scoringNames.size(),0);
  fEventEdepNonIonizing.assign(scoringNames.size(),0);
  fLastDoseTotal.assign(scoringNames.size(),0);
  fLastDoseNonIonizing.assign(scoringNames.size(),0);
  fLastDose.assign(scoringNames.size(),0);
  fLastDoseRate.assign(scoringNames.size(),0);
  fLastDDD.assign(scoringNames.size(),0);
  fLastDDDRate.assign(scoringNames.size(),0);
  fLastNEq.assign(scoringNames.size(),0);
  fLastNEqRate.assign(scoringNames.size(),0);
  fInP.assign(N,0); fInA.assign(N,0);
  fOutP.assign(N,0); fOutA.assign(N,0); fOutN.assign(N,0);

  fCumulativeMode=static_cast<bool>(fOpts.caseDefinition);
  if(fOpts.caseDefinition){
    fAccuracyController=std::make_unique<shield_sim::AccuracyController>(
        fOpts.caseDefinition->accuracy);
    for(const auto& quantity:fOpts.caseDefinition->quantities){
      const auto& requirement=fOpts.caseDefinition->AccuracyFor(quantity);
      if(requirement.systematicBudgetComplete)
        fAccuracyController->SetSystematicRelativeUncertainty(
            quantity,requirement.systematicRelativeUncertainty);
      else
        fAccuracyController->MarkSystematicBudgetIncomplete(quantity);
      if(requirement.hasKnownSystematicFloor)
        fAccuracyController->SetKnownSystematicFloor(quantity,requirement.knownSystematicFloor);
    }
  }

  // CSV histograms are diagnostic outputs.  The Tecplot spectra file contains
  // the more detailed logarithmic spectra and normalization metadata.
  auto* am=G4AnalysisManager::Instance();
  am->SetVerboseLevel(0);
  am->SetFirstHistoId(0);
  am->SetDefaultFileType("csv");
  G4int nb=500; G4double eH=1e5;
  am->CreateH1("ProtonExitE",  "Proton KE at exit",  nb,0,eH,"MeV");
  am->CreateH1("AlphaExitE",   "Alpha KE at exit",   nb,0,eH,"MeV");
  am->CreateH1("NeutronExitE", "Neutron KE at exit", nb,0,eH,"MeV");
}

RunAction::~RunAction() {}

void RunAction::AddInP  (G4double E){ fInP [SpecBins::Bin(E)]+=1; }
void RunAction::AddInA  (G4double E){ fInA [SpecBins::Bin(E)]+=1; }
void RunAction::AddOutP (G4double E){ fOutP[SpecBins::Bin(E)]+=1; }
void RunAction::AddOutA (G4double E){ fOutA[SpecBins::Bin(E)]+=1; }
void RunAction::AddOutN (G4double E){ fOutN[SpecBins::Bin(E)]+=1; }

bool RunAction::DiagnosticLimitReached(G4long rowsWritten) const {
  return (fOpts.diagnosticMaxRows>0 && rowsWritten>=fOpts.diagnosticMaxRows);
}

void RunAction::RecordSourceSample(const std::string& species,
                                   G4double energyMeV,
                                   G4double xMM, G4double yMM, G4double zMM,
                                   G4double ux, G4double uy, G4double uz){
  // The source diagnostic is designed for automated source-sampling tests.
  // It records the final particle state handed to G4ParticleGun, not the
  // spectrum table itself.  Therefore beam tests can verify that x=y=0 and
  // u=(0,0,1), while isotropic tests can verify p(mu)=2mu by taking mu=uz.
  if(!fSourceDump.is_open() || DiagnosticLimitReached(fSourceDumpRows)) return;
  fSourceDump<<std::scientific<<std::setprecision(12)
             <<fSourceDumpRows<<' '<<species<<' '<<energyMeV<<' '
             <<xMM<<' '<<yMM<<' '<<zMM<<' '
             <<ux<<' '<<uy<<' '<<uz<<'\n';
  ++fSourceDumpRows;
}

void RunAction::RecordExitParticle(const std::string& species,
                                   G4double energyMeV,
                                   G4double xGlobalMM, G4double yGlobalMM, G4double zGlobalMM,
                                   G4double xLocalMM,  G4double yLocalMM,  G4double zLocalMM,
                                   G4double uxLocal,   G4double uyLocal,   G4double uzLocal){
  // The exit diagnostic is designed to validate shield-rear-face scoring.
  // It is called only after SteppingAction has transformed the post-step point
  // into shield-local coordinates and confirmed that the particle crossed the
  // downstream face.  The local z coordinate should therefore be close to the
  // local shield half-thickness, and the local direction cosine uzLocal should
  // be positive.
  if(!fExitDump.is_open() || DiagnosticLimitReached(fExitDumpRows)) return;
  fExitDump<<std::scientific<<std::setprecision(12)
           <<fExitDumpRows<<' '<<species<<' '<<energyMeV<<' '
           <<xGlobalMM<<' '<<yGlobalMM<<' '<<zGlobalMM<<' '
           <<xLocalMM <<' '<<yLocalMM <<' '<<zLocalMM <<' '
           <<uxLocal  <<' '<<uyLocal  <<' '<<uzLocal  <<'\n';
  ++fExitDumpRows;
}

void RunAction::AddEdep(std::size_t i,G4double totalEdep,G4double nonIonizingEdep){
  if(i>=fEdepTotal.size()) return;
  // Protect only against sub-roundoff overshoot in the Geant4 ledgers.  A
  // material negative discrepancy is retained as a fatal scoring error rather
  // than silently changing the energy accounting.
  if(nonIonizingEdep<0.0 || totalEdep<0.0 || nonIonizingEdep>totalEdep+1.0e-12*MeV)
    G4Exception("RunAction::AddEdep","EnergyLedger",FatalException,
                "Invalid total/non-ionizing energy-deposit ledger");
  const G4double nonIon=std::min(totalEdep,nonIonizingEdep);
  fEdepTotal[i]+=totalEdep;
  fEdepNonIonizing[i]+=nonIon;
  fEdepIonizing[i]+=totalEdep-nonIon;
  fEventEdepIonizing[i]+=totalEdep-nonIon;
  fEventEdepTotal[i]+=totalEdep;
  fEventEdepNonIonizing[i]+=nonIon;
}

namespace {
std::string ScoreToken(std::string value){
  for(char& character:value)
    if(!(std::isalnum(static_cast<unsigned char>(character)) || character=='_' || character=='-'))
      character='_';
  return value;
}
}

void RunAction::RecordShieldExit(const std::string& particle,int,int,int,
                                 bool primary,const std::string& rootParticle,
                                 int generation,G4double kineticEnergyMeV){
  const std::string component=primary ? "primary" : "secondary";
  const std::string base="shield_exit/"+component+"/"+ScoreToken(particle);
  fEventRawScores[base+"/crossing_count"]+=1.0;
  fEventRawScores[base+"/kinetic_energy_sum_MeV"]+=kineticEnergyMeV;
  fRawScoreUnits[base+"/crossing_count"]="crossings/primary";
  fRawScoreUnits[base+"/kinetic_energy_sum_MeV"]="MeV/primary";
  // The ancestry-qualified score is deliberately additional to the ordinary
  // species/component score: consumers can aggregate by transported species
  // without losing the root-primary provenance needed for mixed environments.
  const std::string ancestry=base+"/root_"+ScoreToken(rootParticle)+
      "/generation_"+std::to_string(generation)+"/crossing_count";
  fEventRawScores[ancestry]+=1.0;
  fRawScoreUnits[ancestry]="crossings/primary";
}

void RunAction::RecordLayerStep(std::size_t layer,G4double totalEdep,
                                G4double nonIonizingEdep,G4double trackLength){
  if(!fOpts.caseDefinition || layer>=fOpts.caseDefinition->shieldLayers.size()) return;
  const std::string base="shield_layer/"+
      ScoreToken(fOpts.caseDefinition->shieldLayers[layer].id);
  const double total=totalEdep/MeV;
  const double nonIon=std::min(totalEdep,nonIonizingEdep)/MeV;
  fEventRawScores[base+"/total_deposited_energy_MeV"]+=total;
  fEventRawScores[base+"/nonionizing_deposited_energy_MeV"]+=nonIon;
  fEventRawScores[base+"/ionizing_deposited_energy_MeV"]+=total-nonIon;
  fEventRawScores[base+"/track_length_mm"]+=trackLength/mm;
  fEventRawScores[base+"/step_count"]+=1.0;
  fRawScoreUnits[base+"/total_deposited_energy_MeV"]="MeV/primary";
  fRawScoreUnits[base+"/nonionizing_deposited_energy_MeV"]="MeV/primary";
  fRawScoreUnits[base+"/ionizing_deposited_energy_MeV"]="MeV/primary";
  fRawScoreUnits[base+"/track_length_mm"]="mm/primary";
  fRawScoreUnits[base+"/step_count"]="steps/primary";
}

void RunAction::RecordSecondaryBirth(std::size_t layer,const std::string& particle,
                                     int pdg,int atomicNumber,int massNumber,
                                     const std::string& rootParticle,int generation,
                                     const std::string& creatorProcess,
                                     G4double kineticEnergyMeV){
  if(!fOpts.caseDefinition || layer>=fOpts.caseDefinition->shieldLayers.size()) return;
  const std::string key=std::to_string(layer)+"|"+particle+"|"+rootParticle+"|"+
                        std::to_string(generation)+"|"+creatorProcess;
  auto& ledger=fSecondaryLedgers[key];
  ledger.layer=layer;
  ledger.particle=particle;
  ledger.pdg=pdg;
  ledger.atomicNumber=atomicNumber;
  ledger.massNumber=massNumber;
  ledger.rootParticle=rootParticle;
  ledger.generation=generation;
  ledger.creatorProcess=creatorProcess;
  ++ledger.count;
  ledger.initialKineticEnergyMeV+=kineticEnergyMeV;
  const std::string score="shield_layer/"+
      ScoreToken(fOpts.caseDefinition->shieldLayers[layer].id)+
      "/secondary_birth/"+ScoreToken(particle)+"/count";
  fEventRawScores[score]+=1.0;
  fRawScoreUnits[score]="particles/primary";
}

void RunAction::AddPrimaryProtonExit(G4double kineticEnergyMeV){
  // Thresholds are evaluated at the exact crossing energy, not at an output
  // histogram bin center. The approved definition uses strict E > threshold.
  // Count a primary history at most once at each nested threshold: a rare
  // backscatter/re-entry of the same primary must not masquerade as another
  // independent incident particle in the ratio or its binomial covariance.
  if(kineticEnergyMeV>10.0 && !fEventPrimaryProtonAbove10){
    ++fPrimaryProtonAbove10;
    fEventPrimaryProtonAbove10=true;
  }
  if(kineticEnergyMeV>100.0 && !fEventPrimaryProtonAbove100){
    ++fPrimaryProtonAbove100;
    fEventPrimaryProtonAbove100=true;
  }
}

void RunAction::BeginEvent(){
  std::fill(fEventEdepIonizing.begin(),fEventEdepIonizing.end(),0.0);
  std::fill(fEventEdepTotal.begin(),fEventEdepTotal.end(),0.0);
  std::fill(fEventEdepNonIonizing.begin(),fEventEdepNonIonizing.end(),0.0);
  fEventRawScores.clear();
  fEventPrimaryProtonAbove10=false;
  fEventPrimaryProtonAbove100=false;
}

void RunAction::EndEvent(){
  // Absorber deposition belongs in the same event vector as boundary/layer
  // scores. This preserves cross-score covariance needed by future nonlinear
  // quantities (which must fold primitive numerator/denominator scores).
  for(std::size_t i=0;i<fEventEdepTotal.size();++i){
    const std::string absorberId=(fOpts.caseDefinition &&
        fOpts.activeAbsorberIndex<fOpts.caseDefinition->absorbers.size())
        ? fOpts.caseDefinition->absorbers[fOpts.activeAbsorberIndex].id
        : std::to_string(i);
    const std::string base="absorber/"+ScoreToken(absorberId);
    const double values[3]={fEventEdepTotal[i]/MeV,fEventEdepIonizing[i]/MeV,
                            fEventEdepNonIonizing[i]/MeV};
    const char* names[3]={"total_deposited_energy_MeV","ionizing_deposited_energy_MeV",
                          "nonionizing_deposited_energy_MeV"};
    for(int j=0;j<3;++j){
      fEventRawScores[base+"/"+names[j]]+=values[j];
      fRawScoreUnits[base+"/"+names[j]]="MeV/primary";
    }
  }
  // Every response accumulator uses the primary history as its independent
  // sample. Scores absent from this map for a history contribute zero; storing
  // only nonzero values is exact because both their sum and sum of squares are
  // unchanged by explicit zeros.
  for(const auto& score:fEventRawScores){
    auto& moments=fRawScores[score.first];
    moments.sum+=score.second;
    moments.sumSquares+=score.second*score.second;
  }
  for(auto first=fEventRawScores.begin();first!=fEventRawScores.end();++first){
    for(auto second=first;second!=fEventRawScores.end();++second){
      fRawCrossProducts[{first->first,second->first}]+=first->second*second->second;
    }
  }
  if(!fAccuracyController) return;
  fAccuracyController->RecordHistory();
  if(fOpts.calcTID){
    G4double sum=0.0;
    for(const G4double value:fEventEdepIonizing) sum+=value;
    // Scaling all samples by target mass leaves relative uncertainty and the
    // adaptive stopping decision unchanged; final dose is computed with the
    // exact Geant4 mass at EndOfRunAction.
    fAccuracyController->AddSample("TID",sum/MeV);
  }
  if(fOpts.calcHardness && fEventPrimaryProtonAbove10)
    fAccuracyController->AddSample("H100/10",fEventPrimaryProtonAbove100?1.0:0.0);
  if(fOpts.caseDefinition &&
     std::find(fOpts.caseDefinition->quantities.begin(),fOpts.caseDefinition->quantities.end(),
               "transport_raw")!=fOpts.caseDefinition->quantities.end()){
    double total=0.0;
    for(double value:fEventEdepTotal) total+=value/MeV;
    fAccuracyController->AddSample("transport_raw",total);
  }
}

bool RunAction::AccuracyStatisticallySatisfied() const{
  if(!fOpts.caseDefinition || !fAccuracyController) return true;
  if(fCumulativeEvents<fOpts.caseDefinition->accuracy.minHistories) return false;
  for(const auto& quantity:fOpts.caseDefinition->quantities){
    const auto& requirement=fOpts.caseDefinition->AccuracyFor(quantity);
    // A documented systematic floor above a combined target is already
    // unachievable. Stop after the required minimum histories and report
    // not_met; spending more histories cannot change that conclusion.
    if(requirement.scope=="combined" && requirement.hasKnownSystematicFloor &&
       requirement.knownSystematicFloor>requirement.target) continue;
    const auto* estimate=fAccuracyController->Estimate(quantity);
    if(!estimate || estimate->count<2 || estimate->mean==0.0 ||
       !std::isfinite(estimate->StandardError())) return false;
    if(estimate->StandardError()/std::abs(estimate->mean)>requirement.target) return false;
  }
  return true;
}

G4long RunAction::GetCumulativeEvents() const { return fCumulativeEvents; }

void RunAction::SetSourceNormalization(G4double normNoAngular,
                                       G4double angularFactor,
                                       const std::string& sourceMode){
  fSourceNormNoAngular = normNoAngular;
  fSourceAngularFactor = angularFactor;
  fSourceNorm          = normNoAngular*angularFactor;
  fSourceMode          = sourceMode;
}

G4double RunAction::GetSourcePlaneAreaCM2() const {
  if(!fDetector) return 0.0;
  return (2.0*fDetector->GetSourcePlaneHalfX()/cm)*
         (2.0*fDetector->GetSourcePlaneHalfY()/cm);
}

G4double RunAction::GetIncidentParticleRate() const {
  if(fSourceMode=="isotropic") return fSourceNorm*GetSourcePlaneAreaCM2();
  return fSourceNorm;
}

void RunAction::SetSweepMode(bool on,const std::string& mat,G4double tMM){
  fSweepMode=on; fCurrentMat=mat; fCurrentTmm=tMM;
}

void RunAction::ConfigureIndependentCase(std::size_t absorberIndex,
                                         const std::string& outputPrefix,
                                         const std::string& scoringName,
                                         G4double scoringThickness){
  if(!fOpts.caseDefinition)
    G4Exception("RunAction::ConfigureIndependentCase","Mode",FatalException,
                "Independent-case reconfiguration requires a versioned case");
  fOpts.activeAbsorberIndex=absorberIndex;
  fOpts.outputPrefix=outputPrefix;
  fOpts.dumpRunSummaryFile=outputPrefix+"_run_summary.dat";
  fScoringNames.assign(1,scoringName);
  fScoringThick.assign(1,scoringThickness);
  fEdepTotal.assign(1,0.0);
  fEdepIonizing.assign(1,0.0);
  fEdepNonIonizing.assign(1,0.0);
  fEventEdepIonizing.assign(1,0.0);
  fEventEdepTotal.assign(1,0.0);
  fEventEdepNonIonizing.assign(1,0.0);
  fLastDose.assign(1,0.0);
  fLastDoseTotal.assign(1,0.0);
  fLastDoseNonIonizing.assign(1,0.0);
  fLastDoseRate.assign(1,0.0);
  fLastDDD.assign(1,0.0);
  fLastDDDRate.assign(1,0.0);
  fLastNEq.assign(1,0.0);
  fLastNEqRate.assign(1,0.0);
  fAccuracyController=std::make_unique<shield_sim::AccuracyController>(
      fOpts.caseDefinition->accuracy);
  for(const auto& quantity:fOpts.caseDefinition->quantities){
    const auto& requirement=fOpts.caseDefinition->AccuracyFor(quantity);
    if(requirement.systematicBudgetComplete)
      fAccuracyController->SetSystematicRelativeUncertainty(
          quantity,requirement.systematicRelativeUncertainty);
    else
      fAccuracyController->MarkSystematicBudgetIncomplete(quantity);
    if(requirement.hasKnownSystematicFloor)
      fAccuracyController->SetKnownSystematicFloor(quantity,requirement.knownSystematicFloor);
  }
  fCumulativeEvents=0;
  fAccumulationStarted=false;
  fPrimaryProtonAbove10=0.0;
  fPrimaryProtonAbove100=0.0;
  fLastHardness=0.0;
  fEventRawScores.clear();
  fRawScores.clear();
  fRawCrossProducts.clear();
  fRawScoreUnits.clear();
  fSecondaryLedgers.clear();
}

void RunAction::AppendSweepPoint(const SweepPoint& sp){ fSweepData.push_back(sp); }
const std::vector<SweepPoint>& RunAction::GetSweepData() const { return fSweepData; }
const std::vector<G4double>&   RunAction::GetLastDose() const { return fLastDose; }
const std::vector<G4double>&   RunAction::GetLastDoseRate() const { return fLastDoseRate; }
const std::vector<G4double>&   RunAction::GetLastDDD() const { return fLastDDD; }
const std::vector<G4double>&   RunAction::GetLastDDDRate() const { return fLastDDDRate; }
const std::vector<G4double>&   RunAction::GetLastNEq() const { return fLastNEq; }
const std::vector<G4double>&   RunAction::GetLastNEqRate() const { return fLastNEqRate; }
G4double RunAction::GetLastHardness() const { return fLastHardness; }

void RunAction::RefreshLVPointers(){
  fShieldLV = fDetector ? fDetector->GetShieldLV() : nullptr;
  fShieldLVs.clear();
  fScoringLVs.clear();
  if(fDetector){
    fShieldLVs=fDetector->GetShieldLVs();
    fScoringLVs = fDetector->GetScoringLVs();
  }
}

G4LogicalVolume* RunAction::GetShieldLV() const{ return fShieldLV; }
const std::vector<G4LogicalVolume*>& RunAction::GetShieldLVs() const{ return fShieldLVs; }
const std::vector<G4LogicalVolume*>& RunAction::GetScoringLVs() const{ return fScoringLVs; }

void RunAction::BeginOfRunAction(const G4Run*) {
  RefreshLVPointers();
  OpenDiagnosticFiles();

  auto* am=G4AnalysisManager::Instance();
  am->Reset();
  std::string csvName=fOpts.outputPrefix+"Output";
  if(fSweepMode)
    csvName+="_"+SanitiseName(fCurrentMat)+"_"+FormatMM(fCurrentTmm)+"mm";
  csvName+=".csv";
  am->OpenFile(csvName);

  const bool reset=!fCumulativeMode || !fAccumulationStarted;
  if(!reset) return;
  fAccumulationStarted=true;
  fCumulativeEvents=0;
  std::fill(fEdepTotal.begin(),fEdepTotal.end(),0);
  std::fill(fEdepIonizing.begin(),fEdepIonizing.end(),0);
  std::fill(fEdepNonIonizing.begin(),fEdepNonIonizing.end(),0);
  fEventRawScores.clear();
  fRawScores.clear();
  fRawCrossProducts.clear();
  fRawScoreUnits.clear();
  fSecondaryLedgers.clear();
  std::fill(fLastDoseTotal.begin(),fLastDoseTotal.end(),0);
  std::fill(fLastDoseNonIonizing.begin(),fLastDoseNonIonizing.end(),0);
  std::fill(fLastDoseRate.begin(),fLastDoseRate.end(),0);
  std::fill(fLastDDD.begin(),fLastDDD.end(),0);
  std::fill(fLastDDDRate.begin(),fLastDDDRate.end(),0);
  std::fill(fLastNEq.begin(),fLastNEq.end(),0);
  std::fill(fLastNEqRate.begin(),fLastNEqRate.end(),0);
  fLastHardness=0.0;
  fPrimaryProtonAbove10=0.0;
  fPrimaryProtonAbove100=0.0;
  std::fill(fInP .begin(),fInP .end(),0);
  std::fill(fInA .begin(),fInA .end(),0);
  std::fill(fOutP.begin(),fOutP.end(),0);
  std::fill(fOutA.begin(),fOutA.end(),0);
  std::fill(fOutN.begin(),fOutN.end(),0);
}

void RunAction::EndOfRunAction(const G4Run* run) {
  fCumulativeEvents+=run->GetNumberOfEvent();
  G4int nEv=fCumulativeMode ? static_cast<G4int>(fCumulativeEvents)
                            : run->GetNumberOfEvent();
  if(nEv==0){ CloseDiagnosticFiles(); return; }

  // Dose per primary. TID is the ionizing ledger. Total absorbed dose and
  // non-ionizing deposited-energy dose are reported independently so the
  // identity total = ionizing + non-ionizing can be audited.
  for(std::size_t i=0;i<fScoringLVs.size();++i){
    auto* lv=fScoringLVs[i];
    G4double dose=0.0;
    if(lv){
      G4double mass = lv->GetMass(true,false);
      if(mass>0.0){
        dose = fEdepIonizing[i]/mass/nEv;
        fLastDoseTotal[i]=fEdepTotal[i]/mass/nEv;
        fLastDoseNonIonizing[i]=fEdepNonIonizing[i]/mass/nEv;
      }
    }
    fLastDose[i]=dose;
  }

  const G4double incidentRate = GetIncidentParticleRate();
  for(std::size_t i=0;i<fLastDose.size();++i)
    fLastDoseRate[i] = fLastDose[i]*incidentRate;

  G4cout<<"\nDose per primary ["
        <<(fSweepMode ? DescribeShieldMaterial(fCurrentMat)+" "+FormatMM(fCurrentTmm)+"mm" : "single run")
        <<"]:"<<G4endl;
  for(std::size_t i=0;i<fScoringNames.size();++i)
    G4cout<<"  "<<fScoringNames[i]<<": "
          <<std::setprecision(5)<<std::scientific
          <<fLastDose[i]/gray<<" Gy/primary"
          <<"   (source-normalized: "
          <<fLastDoseRate[i]/gray<<" Gy/s, if input spectrum units are physical)"
          <<G4endl;

  G4cout<<"Source normalization: mode="<<fSourceMode
        <<", integral(no angular factor)="<<fSourceNormNoAngular
        <<", angular factor="<<fSourceAngularFactor
        <<", normalized source="<<fSourceNorm;
  if(fSourceMode=="isotropic")
    G4cout<<" per cm2/s, source area="<<GetSourcePlaneAreaCM2()
          <<" cm2, incident particle rate="<<incidentRate<<" 1/s";
  else
    G4cout<<" 1/s for the pencil beam";
  G4cout<<G4endl;

  // Compute spectrum-folded quantities after the Geant4 dose has been
  // normalized.  These quantities use the transmitted proton/alpha/neutron
  // spectra scored at the downstream shield face, not additional Geant4 energy
  // deposition inside the target slabs.
  ComputeSpectrumFoldedQuantities(nEv);

  auto* am=G4AnalysisManager::Instance();
  am->Write();
  am->CloseFile();

  WriteSpectraTecplot(nEv);
  WriteComputedQuantitiesTecplot(nEv);
  if(fOpts.calcLET) WriteLETSpectrumTecplot(nEv);

  // Optional machine-readable scalar summary. Automated compatibility tests
  // need stable key/value rows independent of Tecplot variable ordering. For a
  // versioned case, shieldSim.result.v1 is the authoritative result contract.
  WriteRunSummary(nEv);
  WriteResultJson(nEv);
  if(fOpts.caseDefinition && fOpts.caseDefinition->responseIncident.enabled)
    WriteResponseCellJson(nEv);

  CloseDiagnosticFiles();
}

std::string RunAction::SanitiseName(const std::string& s){
  std::string r=s;
  for(auto& c:r) if(c=='/'||c=='\\'||c==' ') c='_';
  return r;
}

std::string RunAction::FormatMM(G4double t){
  std::ostringstream ss; ss<<std::fixed<<std::setprecision(2)<<t;
  return ss.str();
}

std::string RunAction::OutputName(const std::string& suffix) const {
  // Preserve the historical names when outputPrefix is the default
  // "shieldSim".  For example, OutputName("_spectra.dat") becomes
  // "shieldSim_spectra.dat".  Test scripts can set --output-prefix=caseA to
  // obtain caseA_spectra.dat, caseA_quantities.dat, and so on.
  return fOpts.outputPrefix + suffix;
}

void RunAction::OpenDiagnosticFiles(){
  // Diagnostic files are opened at BeginOfRunAction so a geometry sweep can
  // reset them for each run if desired.  The Layer-2 test script uses one run
  // per directory, so the non-append behavior is simplest and avoids stale
  // rows from earlier tests.
  fSourceDumpRows=0;
  fExitDumpRows=0;

  if(!fOpts.dumpSourceSamplesFile.empty()){
    fSourceDump.open(fOpts.dumpSourceSamplesFile);
    if(fSourceDump){
      fSourceDump<<"# shieldSim source-sample diagnostic\n";
      fSourceDump<<"# Columns:\n";
      fSourceDump<<"#   row species E_MeV x_mm y_mm z_mm ux uy uz\n";
      fSourceDump<<"# Positions are global coordinates in mm; u is a unit direction vector.\n";
    } else {
      G4cerr<<"Cannot open source diagnostic file "<<fOpts.dumpSourceSamplesFile<<G4endl;
    }
  }

  if(!fOpts.dumpExitParticlesFile.empty()){
    fExitDump.open(fOpts.dumpExitParticlesFile);
    if(fExitDump){
      fExitDump<<"# shieldSim shield-exit diagnostic\n";
      fExitDump<<"# Columns:\n";
      fExitDump<<"#   row species E_MeV xg_mm yg_mm zg_mm xl_mm yl_mm zl_mm uxl uyl uzl\n";
      fExitDump<<"# Global positions are Geant4 world coordinates.  Local coordinates and directions are in the shield frame.\n";
    } else {
      G4cerr<<"Cannot open exit diagnostic file "<<fOpts.dumpExitParticlesFile<<G4endl;
    }
  }
}

void RunAction::CloseDiagnosticFiles(){
  if(fSourceDump.is_open()) fSourceDump.close();
  if(fExitDump.is_open())   fExitDump.close();
}

void RunAction::WriteSpectraTecplot(G4int nEv){
  const std::string fname=OutputName("_spectra.dat");
  std::ios::openmode mode= (fSweepMode && !fFirstSpectraWrite)
                           ? (std::ios::out|std::ios::app)
                           :  std::ios::out;
  std::ofstream out(fname,mode);
  if(!out){ G4cerr<<"Cannot write "<<fname<<G4endl; return; }

  if(!fSweepMode || fFirstSpectraWrite){
    out<<"TITLE = \"Energetic Particle Shielding Simulation Spectra - Geant4/"<<fOpts.physicsList<<"\"\n";
    out<<"# Physics list: "<<fOpts.physicsList<<"\n";
    out<<"# Energy unit: MeV total kinetic energy per particle.\n";
    out<<"# MC spectra: raw Monte Carlo counts/(MeV primary).\n";
    out<<"# Source-normalized spectra: MC spectra multiplied by the integrated source normalization.\n";
    out<<"# Beam mode normalized units: particles/s/MeV if input spectrum columns are particles/s/MeV.\n";
    out<<"# Isotropic mode normalized units: particles/(cm2 s MeV) if input spectrum columns are particles/(cm2 s sr MeV); pi angular factor applied.\n";
    out<<"# Source mode: "<<fSourceMode<<"\n";
    out<<"# Source normalization integral(no angular factor): "<<fSourceNormNoAngular<<"\n";
    out<<"# Source angular factor: "<<fSourceAngularFactor<<"\n";
    out<<"# Source normalization used for spectra: "<<fSourceNorm<<"\n";
    out<<"# Shield: "<<DescribeShieldMaterial(fSweepMode?fCurrentMat:fOpts.shieldMaterial)<<"\n";
    if(fOpts.caseDefinition && (fOpts.caseDefinition->responseIncident.enabled ||
                                !fOpts.caseDefinition->directIncidents.empty()))
      out<<"# Spectrum: "<<fOpts.caseDefinition->sourceFormat
         <<" Stage-3 verification source; one generated primary per history.\n";
    else if(fOpts.spectrumFile.empty())
      out<<"# Spectrum: built-in approximate Badhwar-O'Neill, phi = 550 MV.\n";
    else
      out<<"# Spectrum file: "<<fOpts.spectrumFile<<"\n";
    out<<"VARIABLES = \"Energy [MeV]\""
       <<" \"Input_Proton_MC [1/(MeV primary)]\""
       <<" \"Input_Alpha_MC [1/(MeV primary)]\""
       <<" \"Output_Proton_MC [1/(MeV primary)]\""
       <<" \"Output_Alpha_MC [1/(MeV primary)]\""
       <<" \"Output_Neutron_MC [1/(MeV primary)]\""
       <<" \"Input_Proton_Norm [source/(MeV)]\""
       <<" \"Input_Alpha_Norm [source/(MeV)]\""
       <<" \"Output_Proton_Norm [source/(MeV)]\""
       <<" \"Output_Alpha_Norm [source/(MeV)]\""
       <<" \"Output_Neutron_Norm [source/(MeV)]\"\n";
    if(fSweepMode) fFirstSpectraWrite=false;
  }

  std::ostringstream zt;
  if(fSweepMode)
    zt<<SanitiseName(DescribeShieldMaterial(fCurrentMat))<<"_"<<FormatMM(fCurrentTmm)<<"mm";
  else
    zt<<SanitiseName(DescribeShieldMaterial(fOpts.shieldMaterial))<<"_"<<FormatMM(fOpts.shieldThickness/mm)<<"mm";
  zt<<" | source-mode="<<fSourceMode<<" | "<<nEv<<" events";

  out<<"ZONE T=\""<<zt.str()<<"\", I="<<SpecBins::N<<", DATAPACKING=POINT\n";
  out<<std::scientific<<std::setprecision(5);
  for(G4int i=0;i<SpecBins::N;++i){
    G4double dE=SpecBins::Width(i);
    G4double nr=dE*nEv;

    const G4double inPmc  = fInP [i]/nr;
    const G4double inAmc  = fInA [i]/nr;
    const G4double outPmc = fOutP[i]/nr;
    const G4double outAmc = fOutA[i]/nr;
    const G4double outNmc = fOutN[i]/nr;

    out<<std::setw(13)<<SpecBins::Center(i)
       <<std::setw(13)<<inPmc
       <<std::setw(13)<<inAmc
       <<std::setw(13)<<outPmc
       <<std::setw(13)<<outAmc
       <<std::setw(13)<<outNmc
       <<std::setw(13)<<inPmc *fSourceNorm
       <<std::setw(13)<<inAmc *fSourceNorm
       <<std::setw(13)<<outPmc*fSourceNorm
       <<std::setw(13)<<outAmc*fSourceNorm
       <<std::setw(13)<<outNmc*fSourceNorm<<"\n";
  }
  out.close();
  G4cout<<"Spectra written: "<<fname
        <<(fSweepMode?" (appended zone)":"")<<G4endl;
}

G4double RunAction::ComputeHardnessIndex(G4int nEv) const{
  // Approved Stage-0 definition: primary protons only, at the shield rear
  // face, using exact crossing kinetic energies and strict > thresholds.
  // A zero denominator is undefined and remains NaN in machine output.
  if(nEv<=0 || fPrimaryProtonAbove10<=0.0)
    return std::numeric_limits<G4double>::quiet_NaN();
  return fPrimaryProtonAbove100/fPrimaryProtonAbove10;
}

void RunAction::ComputeSpectrumFoldedQuantities(G4int nEv){
  // Convert transmitted spectra into DDD and n_eq for each selected target.
  //
  // Step 1: convert raw bin counts to area-averaged fluence per primary.
  //   The transmitted counters fOutP/fOutA/fOutN store particles per energy bin.
  //   For a bin i, count_i/nEv is the number of transmitted particles per
  //   source primary in that energy bin.  Dividing by the finite source/scoring
  //   plane area gives an area-averaged fluence per primary [cm^-2 primary^-1].
  //
  // Step 2: fold the fluence with a material-dependent NIEL response.
  //   D_d = sum_species sum_bins Phi_i * NIEL_species(E_i,material).
  //   With Phi_i in cm^-2 and NIEL in MeV cm^2/g, D_d is MeV/g.
  //
  // Step 3: multiply by the incident source-primary rate to obtain a rate.
  //   In isotropic mode the source rate is the plane-crossing flux times the
  //   finite source-plane area.  In beam mode it is the pencil-beam particle
  //   rate.  This keeps beam-mode DDD an area-averaged estimate over the finite
  //   5 cm x 5 cm scoring plane, matching the dose-geometry convention.
  //
  // Step 4: convert DDD to 1-MeV neutron-equivalent fluence.
  //   n_eq = D_d / NIEL_1MeV_neutron(material).
  //   This is the ASTM-E722-style convention used most commonly for silicon;
  //   for non-silicon materials it is a proxy unless a proper material-specific
  //   damage function is supplied.
  const G4double areaCM2 = std::max(1.0e-30,GetSourcePlaneAreaCM2());
  const G4double incidentRate = GetIncidentParticleRate();

  fLastHardness = fOpts.calcHardness ? ComputeHardnessIndex(nEv) : 0.0;

  for(std::size_t t=0;t<fScoringLVs.size();++t){
    const G4Material* mat = fScoringLVs[t] ? fScoringLVs[t]->GetMaterial() : nullptr;
    G4double dddPerPrimary = 0.0;

    if(fOpts.calcDDD || fOpts.calcNEq){
      for(G4int i=0;i<SpecBins::N;++i){
        const G4double E = SpecBins::Center(i);

        // Bin-integrated fluence per primary, area averaged over the finite
        // downstream plane.  We do not multiply by dE here because fOut* already
        // stores bin-integrated counts; fOut/(nEv*dE) would be the differential
        // spectrum and the dE factor in the integral would cancel.
        const G4double pFluence = fOutP[i]/nEv/areaCM2;
        const G4double aFluence = fOutA[i]/nEv/areaCM2;
        const G4double nFluence = fOutN[i]/nEv/areaCM2;

        dddPerPrimary += pFluence * ComputedQuantities::NIEL_MeV_cm2_g("proton", E, mat);
        dddPerPrimary += aFluence * ComputedQuantities::NIEL_MeV_cm2_g("alpha",  E, mat);
        dddPerPrimary += nFluence * ComputedQuantities::NIEL_MeV_cm2_g("neutron",E, mat);
      }
    }

    const G4double refNIEL = ComputedQuantities::ReferenceNIEL1MeVNeutron_MeV_cm2_g(mat);

    fLastDDD[t]     = fOpts.calcDDD ? dddPerPrimary : 0.0;
    fLastDDDRate[t] = fOpts.calcDDD ? dddPerPrimary*incidentRate : 0.0;
    fLastNEq[t]     = (fOpts.calcNEq && refNIEL>0.0) ? dddPerPrimary/refNIEL : 0.0;
    fLastNEqRate[t] = (fOpts.calcNEq && refNIEL>0.0) ? dddPerPrimary*incidentRate/refNIEL : 0.0;
  }

  if(fOpts.calcDDD || fOpts.calcNEq || fOpts.calcHardness){
    G4cout<<"Computed quantities:"<<G4endl;
    for(std::size_t i=0;i<fScoringNames.size();++i){
      G4cout<<"  "<<fScoringNames[i];
      if(fOpts.calcDDD)
        G4cout<<"  DDD="<<std::setprecision(5)<<std::scientific
              <<fLastDDD[i]<<" MeV/g/primary, "
              <<fLastDDDRate[i]<<" MeV/g/s";
      if(fOpts.calcNEq)
        G4cout<<"  n_eq="<<std::setprecision(5)<<std::scientific
              <<fLastNEq[i]<<" cm^-2/primary, "
              <<fLastNEqRate[i]<<" cm^-2/s";
      G4cout<<G4endl;
    }
    if(fOpts.calcHardness)
      G4cout<<"  H100/10="<<std::setprecision(5)<<std::scientific
            <<fLastHardness<<G4endl;
  }
}

void RunAction::WriteComputedQuantitiesTecplot(G4int nEv){
  // Scalar computed quantities are written as one row per target material for
  // the current shielding configuration.  In sweep mode a new Tecplot zone is
  // appended for each thickness, so the file represents all
  // (shielding x absorber) combinations requested by the run.
  if(!(fOpts.calcTID || fOpts.calcDDD || fOpts.calcNEq || fOpts.calcHardness)) return;

  const std::string fname=OutputName("_quantities.dat");
  const bool append = fSweepMode && !fFirstQuantitiesWrite;
  std::ofstream out(fname, append ? (std::ios::out|std::ios::app) : std::ios::out);
  if(!out){ G4cerr<<"Cannot write "<<fname<<G4endl; return; }

  const G4double shieldThicknessMM = fSweepMode ? fCurrentTmm : fOpts.shieldThickness/mm;
  G4double areal = 0.0;
  if(fOpts.caseDefinition){
    areal=fOpts.caseDefinition->ShieldTotalArealDensityGPerCm2();
  } else if(fShieldLV && fShieldLV->GetMaterial()){
    const G4double rho_gcc = fShieldLV->GetMaterial()->GetDensity()/(g/cm3);
    areal = rho_gcc*shieldThicknessMM*0.1; // 1 mm = 0.1 cm
  }

  if(!append){
    out<<"TITLE = \"ShieldSim Computed Quantities - Geant4/"<<fOpts.physicsList<<"\"\n";
    out<<"# Physics list: "<<fOpts.physicsList<<"\n";
    out<<"# Source mode: "<<fSourceMode<<"\n";
    out<<"# Quantities enabled: "
       <<(fOpts.calcTID?"TID ":"")
       <<(fOpts.calcDDD?"DDD ":"")
       <<(fOpts.calcNEq?"n_eq ":"")
       <<(fOpts.calcLET?"LET ":"")
       <<(fOpts.calcHardness?"H100/10 ":"")<<"\n";
    out<<"# TID equation: D_ion = (E_dep,total - E_dep,nonionizing)/m in the finite target.\n";
    out<<"# Total and non-ionizing dose ledgers are preserved in the run summary for closure.\n";
    out<<"# DDD equation: D_d = integral Phi(E) NIEL(E) dE.  Phi is area-averaged transmitted fluence.\n";
    out<<"# n_eq equation: n_eq = D_d / NIEL_1MeV_neutron(material).  Standard convention is most meaningful for silicon.\n";
    out<<"# H100/10 equation: primary-proton shield-exit fluence(E>100 MeV) / fluence(E>10 MeV); exact crossing energy.\n";
    out<<"# DDD NIEL and LET response functions are documented in src/ComputedQuantities.cc.\n";
    out<<"# Target index mapping for each zone follows the row order below; detailed material definitions are in MaterialCatalog.cc.\n";
    for(std::size_t i=0;i<fScoringNames.size();++i)
      out<<"#   TargetIndex "<<i<<": "<<DescribeDetectorMaterial(fScoringNames[i])
         <<", thickness="<<fScoringThick[i]/mm<<" mm\n";

    out<<"VARIABLES = \"TargetIndex\" \"ShieldThickness [mm]\" \"ArealDensity [g/cm2]\" \"TargetThickness [mm]\"";
    if(fOpts.calcTID)
      out<<" \"TID [Gy/primary]\" \"TID [rad/primary]\" \"TIDRate [Gy/s]\" \"TIDRate [rad/s]\"";
    if(fOpts.calcDDD)
      out<<" \"DDD [MeV/g/primary]\" \"DDDRate [MeV/g/s]\" \"DDD [radEq/primary]\" \"DDDRate [radEq/s]\"";
    if(fOpts.calcNEq)
      out<<" \"n_eq [cm-2/primary]\" \"n_eq_rate [cm-2/s]\"";
    if(fOpts.calcHardness)
      out<<" \"H100_10\"";
    out<<"\n";
  }

  std::ostringstream zt;
  if(fSweepMode)
    zt<<SanitiseName(DescribeShieldMaterial(fCurrentMat))<<"_"<<FormatMM(fCurrentTmm)<<"mm";
  else
    zt<<SanitiseName(DescribeShieldMaterial(fOpts.shieldMaterial))<<"_"<<FormatMM(fOpts.shieldThickness/mm)<<"mm";
  zt<<" | source-mode="<<fSourceMode<<" | "<<nEv<<" events";

  out<<"ZONE T=\""<<zt.str()<<"\", I="<<fScoringNames.size()<<", DATAPACKING=POINT\n";
  out<<std::scientific<<std::setprecision(6);
  for(std::size_t i=0;i<fScoringNames.size();++i){
    out<<std::setw(14)<<i
       <<std::setw(14)<<shieldThicknessMM
       <<std::setw(14)<<areal
       <<std::setw(14)<<fScoringThick[i]/mm;
    if(fOpts.calcTID){
      out<<std::setw(14)<<fLastDose[i]/gray
         <<std::setw(14)<<100.0*(fLastDose[i]/gray)
         <<std::setw(14)<<fLastDoseRate[i]/gray
         <<std::setw(14)<<100.0*(fLastDoseRate[i]/gray);
    }
    if(fOpts.calcDDD){
      out<<std::setw(14)<<fLastDDD[i]
         <<std::setw(14)<<fLastDDDRate[i]
         <<std::setw(14)<<ComputedQuantities::MeVPerGramToRad(fLastDDD[i])
         <<std::setw(14)<<ComputedQuantities::MeVPerGramToRad(fLastDDDRate[i]);
    }
    if(fOpts.calcNEq){
      out<<std::setw(14)<<fLastNEq[i]
         <<std::setw(14)<<fLastNEqRate[i];
    }
    if(fOpts.calcHardness){
      out<<std::setw(14)<<fLastHardness;
    }
    out<<"\n";
  }
  out.close();
  if(fSweepMode) fFirstQuantitiesWrite=false;
  G4cout<<"Computed quantities written: "<<fname
        <<(fSweepMode?" (appended zone)":"")<<G4endl;
}

void RunAction::WriteLETSpectrumTecplot(G4int nEv){
  // LET spectra are computed from the transmitted charged-particle spectra.
  // For each target material, every transmitted proton/alpha energy bin is
  // mapped to an LET bin using ChargedLET_MeV_cm2_mg().  The LET spectrum is
  // therefore not an additional Geant4 detector tally; it is the standard
  // response-function folding operation
  //
  //   dPhi/dLET = integral dPhi/dE * delta(LET - LET(E,material)) dE .
  //
  // The output contains an area-averaged per-primary fluence spectrum and a
  // source-normalized fluence-rate spectrum.  This is useful for comparing the
  // radiation quality of different shielding/target configurations.
  const std::string fname=OutputName("_let_spectrum.dat");
  const bool append = fSweepMode && !fFirstLETWrite;
  std::ofstream out(fname, append ? (std::ios::out|std::ios::app) : std::ios::out);
  if(!out){ G4cerr<<"Cannot write "<<fname<<G4endl; return; }

  const G4double areaCM2 = std::max(1.0e-30,GetSourcePlaneAreaCM2());
  const G4double incidentRate = GetIncidentParticleRate();

  if(!append){
    out<<"TITLE = \"ShieldSim LET Spectra - Geant4/"<<fOpts.physicsList<<"\"\n";
    out<<"# LET unit: MeV cm2/mg.\n";
    out<<"# Equation: dPhi/dLET = integral dPhi/dE delta(LET-LET(E,material)) dE.\n";
    out<<"# LET(E,material) is approximated by a Bethe-Bloch electronic mass stopping power in src/ComputedQuantities.cc.\n";
    out<<"# PerPrimary columns are area-averaged fluence per source primary per LET unit [cm-2 primary-1 /(MeV cm2 mg-1)].\n";
    out<<"# Rate columns multiply the per-primary fluence by the source-primary rate [cm-2 s-1 /(MeV cm2 mg-1)].\n";
    out<<"VARIABLES = \"LET [MeV cm2/mg]\""
       <<" \"Proton_PerPrimary\" \"Alpha_PerPrimary\" \"Charged_PerPrimary\""
       <<" \"Proton_Rate\" \"Alpha_Rate\" \"Charged_Rate\"\n";
  }

  for(std::size_t t=0;t<fScoringLVs.size();++t){
    const G4Material* mat = fScoringLVs[t] ? fScoringLVs[t]->GetMaterial() : nullptr;
    std::vector<G4double> p(ComputedQuantities::LETBins::N,0.0);
    std::vector<G4double> a(ComputedQuantities::LETBins::N,0.0);

    for(G4int i=0;i<SpecBins::N;++i){
      const G4double E = SpecBins::Center(i);
      const G4double pFluence = fOutP[i]/nEv/areaCM2;
      const G4double aFluence = fOutA[i]/nEv/areaCM2;

      const G4double pLET = ComputedQuantities::ChargedLET_MeV_cm2_mg("proton",E,mat);
      const G4double aLET = ComputedQuantities::ChargedLET_MeV_cm2_mg("alpha", E,mat);
      if(pLET>0.0) p[ComputedQuantities::LETBins::Bin(pLET)] += pFluence;
      if(aLET>0.0) a[ComputedQuantities::LETBins::Bin(aLET)] += aFluence;
    }

    std::ostringstream zt;
    if(fSweepMode)
      zt<<SanitiseName(DescribeShieldMaterial(fCurrentMat))<<"_"<<FormatMM(fCurrentTmm)<<"mm";
    else
      zt<<SanitiseName(DescribeShieldMaterial(fOpts.shieldMaterial))<<"_"<<FormatMM(fOpts.shieldThickness/mm)<<"mm";
    zt<<" | target="<<SanitiseName(fScoringNames[t])<<" | source-mode="<<fSourceMode;

    out<<"ZONE T=\""<<zt.str()<<"\", I="<<ComputedQuantities::LETBins::N<<", DATAPACKING=POINT\n";
    out<<std::scientific<<std::setprecision(6);
    for(G4int j=0;j<ComputedQuantities::LETBins::N;++j){
      const G4double dL = ComputedQuantities::LETBins::Width(j);
      const G4double pDiff = p[j]/dL;
      const G4double aDiff = a[j]/dL;
      out<<std::setw(14)<<ComputedQuantities::LETBins::Center(j)
         <<std::setw(14)<<pDiff
         <<std::setw(14)<<aDiff
         <<std::setw(14)<<pDiff+aDiff
         <<std::setw(14)<<pDiff*incidentRate
         <<std::setw(14)<<aDiff*incidentRate
         <<std::setw(14)<<(pDiff+aDiff)*incidentRate<<"\n";
    }
  }

  out.close();
  if(fSweepMode) fFirstLETWrite=false;
  G4cout<<"LET spectra written: "<<fname
        <<(fSweepMode?" (appended zones)":"")<<G4endl;
}


void RunAction::WriteRunSummary(G4int nEv){
  // The run-summary file is a deliberately simple, machine-readable diagnostic
  // for Layer-3 physics/numerics tests and Layer-4 regression tests.  It
  // duplicates integrated scalar quantities that are also present in the normal
  // Tecplot outputs, but writes them as fixed keyword rows so a bash/Python
  // test script can parse them robustly.  This file is not intended to replace
  // the science output files.
  //
  // Layer-4 regression tests use this summary as the primary comparison source
  // because it avoids brittle parsing of column-oriented Tecplot files.  For
  // that reason, the first few metadata rows below form a tiny schema.  When new
  // rows are added in the future, keep existing row names and column ordering
  // stable whenever possible; this lets old regression baselines remain useful.
  //
  // Format notes:
  //   meta   <name> <string-value>
  //   scalar <name> <floating-point-value>
  //   count  <name> <floating-point-count>
  //   target <index> <name> <thick_mm> <TID_Gy/primary> <TIDRate_Gy/s>
  //          <DDD_MeV/g/primary> <DDDRate_MeV/g/s> <n_eq/primary> <n_eq/s>
  //
  // The same file can contain several blocks in sweep mode.  The first sweep
  // point overwrites any old file, while subsequent points append.
  if(fOpts.dumpRunSummaryFile.empty()) return;

  const bool append = fSweepMode && !fFirstSummaryWrite;
  std::ofstream out(fOpts.dumpRunSummaryFile, append ? (std::ios::out|std::ios::app) : std::ios::out);
  if(!out){
    G4cerr<<"Cannot write run-summary diagnostic "<<fOpts.dumpRunSummaryFile<<G4endl;
    return;
  }

  auto sumCounts=[](const std::vector<G4double>& v){
    G4double s=0.0;
    for(G4double x:v) s+=x;
    return s;
  };

  const G4double shieldThicknessMM = fSweepMode ? fCurrentTmm : fOpts.shieldThickness/mm;
  G4double areal = 0.0;
  if(fOpts.caseDefinition){
    areal=fOpts.caseDefinition->ShieldTotalArealDensityGPerCm2();
  } else if(fShieldLV && fShieldLV->GetMaterial()){
    const G4double rho_gcc = fShieldLV->GetMaterial()->GetDensity()/(g/cm3);
    areal = rho_gcc*shieldThicknessMM*0.1; // 1 mm = 0.1 cm
  }

  out<<std::scientific<<std::setprecision(12);
  out<<"# shieldSim machine-readable run summary\n";
  out<<"# The target row format is:\n";
  out<<"# target index name thickness_mm TID_Gy_perPrimary TIDRate_Gy_s DDD_MeV_g_perPrimary DDDRate_MeV_g_s n_eq_cm2_perPrimary n_eq_rate_cm2_s\n";
  out<<"begin_run\n";
  out<<"meta summary_schema layer4_v1\n";
  out<<"meta output_prefix "<<fOpts.outputPrefix<<"\n";
  out<<"meta physics_list "<<fOpts.physicsList<<"\n";
  out<<"meta source_mode "<<fSourceMode<<"\n";
  out<<"meta spectrum_file "
     <<(fOpts.caseDefinition && (fOpts.caseDefinition->responseIncident.enabled ||
                                 !fOpts.caseDefinition->directIncidents.empty())
          ? fOpts.caseDefinition->sourceFormat
          : fOpts.spectrumFile.empty()?"builtin":fOpts.spectrumFile)<<"\n";
  out<<"meta shield_material "<<(fSweepMode?fCurrentMat:fOpts.shieldMaterial)<<"\n";
  out<<"meta result_capability_status "
     <<(fOpts.caseDefinition
          ? (fOpts.caseDefinition->schemaVersion=="shieldSim.case.v2"
               ? "stage3_transport_engine_implemented_pending_independent_review_no_science_quantity_claim"
               : fOpts.caseDefinition->sourceFormat=="ams02_hapi_component"
               ? "stage2_published_spectrum_transport_test_not_physical_validation"
               : "stage1_core_not_yet_quantity_validated")
          : "legacy_prototype")<<"\n";
  out<<"scalar events "<<nEv<<"\n";
  out<<"scalar target_count "<<fScoringNames.size()<<"\n";
  out<<"scalar shield_thickness_mm "<<shieldThicknessMM<<"\n";
  out<<"scalar shield_areal_density_g_cm2 "<<areal<<"\n";
  out<<"scalar source_norm_no_angular "<<fSourceNormNoAngular<<"\n";
  out<<"scalar source_angular_factor "<<fSourceAngularFactor<<"\n";
  out<<"scalar source_norm "<<fSourceNorm<<"\n";
  out<<"scalar source_plane_area_cm2 "<<GetSourcePlaneAreaCM2()<<"\n";
  out<<"scalar incident_particle_rate_s "<<GetIncidentParticleRate()<<"\n";
  out<<"scalar H100_10 "<<fLastHardness<<"\n";
  out<<"count H100_10_primary_proton_numerator "<<fPrimaryProtonAbove100<<"\n";
  out<<"count H100_10_primary_proton_denominator "<<fPrimaryProtonAbove10<<"\n";
  if(fOpts.caseDefinition && fAccuracyController){
    for(const auto& quantity:fOpts.caseDefinition->quantities){
      const auto assessment=fAccuracyController->Assess(quantity);
      const auto& requirement=fOpts.caseDefinition->AccuracyFor(quantity);
      out<<"accuracy "<<SanitiseName(quantity)<<' '
         <<requirement.target<<' '<<requirement.scope<<' '
         <<shield_sim::AccuracyStatusName(assessment.status)<<' '
         <<assessment.statisticalRelativeUncertainty<<' '
         <<assessment.combinedRelativeUncertainty<<' '
         <<assessment.interval95Low<<' '<<assessment.interval95High<<' '
         <<'"'<<assessment.reason<<'"'<<"\n";
    }
  }
  out<<"count input_proton "<<sumCounts(fInP)<<"\n";
  out<<"count input_alpha "<<sumCounts(fInA)<<"\n";
  out<<"count output_proton "<<sumCounts(fOutP)<<"\n";
  out<<"count output_alpha "<<sumCounts(fOutA)<<"\n";
  out<<"count output_neutron "<<sumCounts(fOutN)<<"\n";
  out<<"scalar production_cut_mm "<<(fOpts.productionCut>0.0 ? fOpts.productionCut/mm : -1.0)<<"\n";
  out<<"scalar max_step_mm "<<(fOpts.maxStepLength>0.0 ? fOpts.maxStepLength/mm : -1.0)<<"\n";

  for(std::size_t i=0;i<fShieldLVs.size();++i){
    const auto* lv=fShieldLVs[i];
    out<<"shield_layer_material "<<i<<' '
       <<(lv&&lv->GetMaterial()?SanitiseName(lv->GetMaterial()->GetName()):"null")<<' '
       <<(lv&&lv->GetMaterial()?lv->GetMaterial()->GetDensity()/(g/cm3):0.0)<<"\n";
  }

  for(std::size_t i=0;i<fScoringNames.size();++i){
    out<<"target "<<i<<' '<<SanitiseName(fScoringNames[i])<<' '
       <<fScoringThick[i]/mm<<' '
       <<(i<fLastDose.size()?fLastDose[i]/gray:0.0)<<' '
       <<(i<fLastDoseRate.size()?fLastDoseRate[i]/gray:0.0)<<' '
       <<(i<fLastDDD.size()?fLastDDD[i]:0.0)<<' '
       <<(i<fLastDDDRate.size()?fLastDDDRate[i]:0.0)<<' '
       <<(i<fLastNEq.size()?fLastNEq[i]:0.0)<<' '
       <<(i<fLastNEqRate.size()?fLastNEqRate[i]:0.0)<<"\n";
    auto* targetLV=i<fScoringLVs.size()?fScoringLVs[i]:nullptr;
    out<<"target_material "<<i<<' '
       <<(targetLV&&targetLV->GetMaterial()?SanitiseName(targetLV->GetMaterial()->GetName()):"null")<<' '
       <<(targetLV&&targetLV->GetMaterial()?targetLV->GetMaterial()->GetDensity()/(g/cm3):0.0)<<' '
       <<(targetLV?targetLV->GetMass(true,false)/kg:0.0)<<"\n";
    out<<"deposit_ledger "<<i<<' '
       <<(i<fLastDoseTotal.size()?fLastDoseTotal[i]/gray:0.0)<<' '
       <<(i<fLastDose.size()?fLastDose[i]/gray:0.0)<<' '
       <<(i<fLastDoseNonIonizing.size()?fLastDoseNonIonizing[i]/gray:0.0)<<"\n";
  }
  out<<"end_run\n";

  if(fSweepMode) fFirstSummaryWrite=false;
  G4cout<<"Run summary written: "<<fOpts.dumpRunSummaryFile
        <<(append?" (appended block)":"")<<G4endl;
}

shield_sim::Json RunAction::BuildRawScoreCovarianceJson(G4int nEv) const {
  // Covariance is emitted for every upper-triangular score pair, including
  // exact zeros. A sparse event-product accumulator keeps runtime proportional
  // to scores present in an event; missing cross-products are mathematically
  // zero because absent scores contribute zero for that primary history.
  shield_sim::Json::Array covariance;
  if(nEv<=0) return shield_sim::Json(std::move(covariance));
  for(auto first=fRawScores.begin();first!=fRawScores.end();++first){
    for(auto second=first;second!=fRawScores.end();++second){
      const auto product=fRawCrossProducts.find({first->first,second->first});
      const double sumProducts=product==fRawCrossProducts.end()?0.0:product->second;
      double sampleCovariance=nEv>1
          ? (sumProducts-first->second.sum*second->second.sum/nEv)/(nEv-1.0) : 0.0;
      if(first==second) sampleCovariance=std::max(0.0,sampleCovariance);
      const auto firstUnit=fRawScoreUnits.find(first->first);
      const auto secondUnit=fRawScoreUnits.find(second->first);
      covariance.emplace_back(shield_sim::Json::Object{
        {"covariance_of_mean",sampleCovariance/nEv},
        {"histories",static_cast<double>(nEv)},
        {"sample_covariance",sampleCovariance},
        {"score_id_i",first->first},{"score_id_j",second->first},
        {"sum_products",sumProducts},
        {"unit",(firstUnit==fRawScoreUnits.end()?"unknown":firstUnit->second)+
                    std::string(" * ")+
                    (secondUnit==fRawScoreUnits.end()?"unknown":secondUnit->second)}});
    }
  }
  return shield_sim::Json(std::move(covariance));
}

void RunAction::WriteResultJson(G4int nEv){
  if(!fOpts.caseDefinition || !fAccuracyController) return;
  const auto& c=*fOpts.caseDefinition;
  const std::string validationStatus = c.schemaVersion=="shieldSim.case.v2"
      ? "implemented_stage3_transport_raw_not_physically_validated"
      : c.sourceFormat=="ams02_hapi_component"
      ? "stage2_published_spectrum_transport_test_not_physical_validation"
      : "implemented_not_physically_validated_stage1";
  const auto& absorber=c.absorbers.at(fOpts.activeAbsorberIndex);
  auto finiteOrNull=[](double value)->shield_sim::Json{
    return std::isfinite(value) ? shield_sim::Json(value) : shield_sim::Json(nullptr);
  };
  auto accuracyJson=[&](const std::string& quantity,
                        shield_sim::AccuracyAssessment assessment){
    const auto& requirement=c.AccuracyFor(quantity);
    shield_sim::Json::Object object;
    object["requested_target"]=requirement.target;
    object["scope"]=requirement.scope;
    object["standard_k"]=requirement.standardK;
    object["standard_probability"]=requirement.standardProbability;
    object["interval_probability"]=requirement.intervalProbability;
    object["interval_method"]=requirement.intervalMethod;
    object["status"]=shield_sim::AccuracyStatusName(assessment.status);
    object["reason"]=assessment.reason;
    object["statistical_relative_standard_uncertainty"]=finiteOrNull(
        assessment.statisticalRelativeUncertainty);
    object["combined_relative_standard_uncertainty"]=finiteOrNull(
        assessment.combinedRelativeUncertainty);
    object["systematic_budget_status"]=requirement.systematicBudgetComplete
        ? "complete" : "incomplete";
    if(requirement.systematicBudgetComplete)
      object["systematic_relative_standard_uncertainty"]=
          requirement.systematicRelativeUncertainty;
    if(requirement.hasKnownSystematicFloor)
      object["known_systematic_relative_floor"]=requirement.knownSystematicFloor;
    return shield_sim::Json(std::move(object));
  };

  shield_sim::Json::Array outputs;
  if(fOpts.calcTID){
    for(std::size_t i=0;i<fScoringLVs.size();++i){
      auto* lv=fScoringLVs[i];
      const G4double mass=lv?lv->GetMass(true,false):0.0;
      auto assessment=fAccuracyController->Assess("TID");
      shield_sim::Json::Object interval;
      if(mass>0.0 && std::isfinite(assessment.interval95Low) &&
         std::isfinite(assessment.interval95High)){
        interval["low"]=assessment.interval95Low*MeV/mass/gray;
        interval["high"]=assessment.interval95High*MeV/mass/gray;
        interval["unit"]="Gy/primary";
      } else {
        interval["low"]=shield_sim::Json(nullptr);
        interval["high"]=shield_sim::Json(nullptr);
        interval["unit"]="Gy/primary";
      }
      outputs.emplace_back(shield_sim::Json::Object{
        {"accuracy",accuracyJson("TID",assessment)},
        {"absorber_id",absorber.id},
        {"absorber_mass",shield_sim::Json(shield_sim::Json::Object{{"unit","kg"},{"value",mass/kg}})},
        {"definition","ionizing deposited energy divided by finite absorber mass"},
        {"deposited_energy_ledger",shield_sim::Json(shield_sim::Json::Object{
          {"identity","total = ionizing + nonionizing"},
          {"ionizing_dose_Gy_per_primary",fLastDose[i]/gray},
          {"nonionizing_dose_Gy_per_primary",fLastDoseNonIonizing[i]/gray},
          {"total_dose_Gy_per_primary",fLastDoseTotal[i]/gray}})},
        {"interval_95",shield_sim::Json(std::move(interval))},
        {"name","TID"}, {"normalization","per incident primary"},
        {"source_normalized_unit","Gy/s"},
        {"source_normalized_value",fLastDoseRate[i]/gray},
        {"unit","Gy/primary"}, {"validation_status",validationStatus},
        {"value",fLastDose[i]/gray}
      });
    }
  }
  if(fOpts.calcHardness){
    auto assessment=fAccuracyController->Assess("H100/10");
    // H100/10 is a conditional binomial proportion because every >100 MeV
    // crossing is also a >10 MeV crossing. Use a bounded Wilson 95% interval,
    // not an unbounded normal interval.
    shield_sim::Json low(nullptr),high(nullptr);
    if(fPrimaryProtonAbove10>0.0){
      const double n=fPrimaryProtonAbove10;
      const double p=fPrimaryProtonAbove100/n;
      const double z=1.959963984540054;
      const double denominator=1.0+z*z/n;
      const double center=(p+z*z/(2.0*n))/denominator;
      const double half=z*std::sqrt(p*(1.0-p)/n+z*z/(4.0*n*n))/denominator;
      low=std::max(0.0,center-half);
      high=std::min(1.0,center+half);
    }
    shield_sim::Json accuracy=accuracyJson("H100/10",assessment);
    accuracy.AsObject()["interval_95_method"]="Wilson score interval for bounded conditional proportion";
    outputs.emplace_back(shield_sim::Json::Object{
      {"accuracy",std::move(accuracy)}, {"absorber_id",absorber.id},
      {"covariance_note","event indicators are nested: I(E>100) is a subset of I(E>10)"},
      {"indicator_sample_covariance",nEv>1
          ? (fPrimaryProtonAbove100-
             fPrimaryProtonAbove100*fPrimaryProtonAbove10/nEv)/(nEv-1.0)
          : 0.0},
      {"estimated_fluence_covariance_per_primary_squared",nEv>1
          ? (fPrimaryProtonAbove100-
             fPrimaryProtonAbove100*fPrimaryProtonAbove10/nEv)/((nEv-1.0)*nEv)
          : 0.0},
      {"definition","primary-proton shield-exit fluence(E>100 MeV) / fluence(E>10 MeV)"},
      {"denominator_count",fPrimaryProtonAbove10},
      {"interval_95",shield_sim::Json(shield_sim::Json::Object{{"high",high},{"low",low}})},
      {"name","H100/10"}, {"numerator_count",fPrimaryProtonAbove100},
      {"unit","1"}, {"validation_status",validationStatus},
      {"value",finiteOrNull(fLastHardness)}
    });
  }

  shield_sim::Json::Array rawScores;
  for(const auto& item:fRawScores){
    const double n=nEv;
    const double mean=n>0.0?item.second.sum/n:0.0;
    const double variance=n>1.0?std::max(0.0,(item.second.sumSquares-
        item.second.sum*item.second.sum/n)/(n-1.0)):0.0;
    const auto unit=fRawScoreUnits.find(item.first);
    rawScores.emplace_back(shield_sim::Json::Object{
      {"histories",n},{"mean_per_primary",mean},{"sample_variance",variance},
      {"score_id",item.first},{"standard_error",n>0.0?std::sqrt(variance/n):0.0},
      {"sum",item.second.sum},{"sum_squares",item.second.sumSquares},
      {"unit",unit==fRawScoreUnits.end()?"declared_by_score":unit->second}});
  }

  shield_sim::Json::Object result{
    {"absorber_geometry_hash",absorber.geometryHash}, {"absorber_id",absorber.id},
    {"case_hash",c.caseHash}, {"code_revision",shield_sim_build::GitRevision},
    {"events",static_cast<double>(nEv)}, {"geant4_version",G4Version},
    {"geometry_hash",c.geometryHash}, {"material_library_hash",c.materials->ContentHash()},
    {"outputs",shield_sim::Json(std::move(outputs))}, {"physics_list",fOpts.physicsList},
    {"ion_physics",c.ionPhysics},
    {"energy_ledger_scope",shield_sim::Json(shield_sim::Json::Object{
      {"absorber_deposition","Geant4 total, ionizing, and nonionizing deposited energy per primary"},
      {"shield_crossings","kinetic energy only at the downstream stack boundary"},
      {"secondary_births","count and initial kinetic energy grouped by production layer, particle, ancestry generation, and creator process"},
      {"nuclear_rest_mass_and_q_values","not included in the kinetic/deposition closure; no global conservation claim is made"}})},
    {"raw_transport_scores",shield_sim::Json(std::move(rawScores))},
    {"raw_transport_score_covariance",BuildRawScoreCovarianceJson(nEv)},
    {"random_seed",static_cast<double>(fOpts.randomSeed)},
    {"schema_version",c.schemaVersion=="shieldSim.case.v2"
        ? "shieldSim.result.v2" : "shieldSim.result.v1"},
    {"source_component_id",c.sourceComponentId}, {"source_format",c.sourceFormat},
    {"source_hash",c.spectrum.contentHash},
    {"source_id",c.spectrum.id},
    {"source_normalization",shield_sim::Json(shield_sim::Json::Object{
      {"angular_factor",fSourceAngularFactor},
      {"incident_particle_rate",GetIncidentParticleRate()},
      {"incident_particle_rate_unit","1/s"},
      {"integral_no_angular",fSourceNormNoAngular},
      {"integral_no_angular_unit",c.spectrum.valueType=="differential_intensity"
          ? "1/(cm2 sr s)" : "source_unit_integrated_over_energy"},
      {"source_plane_area_cm2",GetSourcePlaneAreaCM2()}})},
    {"source_science_status",c.spectrum.scienceStatus}
  };
  const std::string path=fOpts.outputPrefix+"_result.json";
  std::ofstream out(path);
  if(!out)
    G4Exception("RunAction::WriteResultJson","Output",FatalException,
                ("Cannot write result JSON '"+path+"'").c_str());
  out<<shield_sim::Json(std::move(result)).DumpCanonical()<<'\n';
  G4cout<<"Result JSON written: "<<path<<G4endl;
}

void RunAction::WriteResponseCellJson(G4int nEv){
  if(!fOpts.caseDefinition || !fOpts.caseDefinition->responseIncident.enabled) return;
  const auto& c=*fOpts.caseDefinition;
  const auto& incident=c.responseIncident;
  const auto& absorber=c.absorbers.at(fOpts.activeAbsorberIndex);

  // The response key intentionally excludes incident energy/direction, random
  // seed, history count, and output path: those identify a cell or realization,
  // not the reusable transport configuration.  Conversely it includes every
  // geometry/material/scoring/physics/build field whose change can alter a
  // response.  This errs on the safe side: a different executable revision or
  // data-library path invalidates reuse even if a later review finds it benign.
  shield_sim::Json::Object dataLibraries;
  for(const char* name:{"G4LEDATA","G4NEUTRONHPDATA","G4PARTICLEXSDATA","G4SAIDXSDATA",
                        "G4ABLADATA","G4ENSDFSTATEDATA","G4INCLDATA","G4PIIDATA"}){
    const char* value=std::getenv(name);
    dataLibraries[name]=value?value:"not_set";
  }
  shield_sim::Json::Object keyDocument{
    {"absorber_geometry_hash",absorber.geometryHash},
    {"absorber_geometry_provider","slab.v1"},
    {"build_type",shield_sim_build::BuildType},
    {"code_revision",shield_sim_build::GitRevision},
    {"compiler",shield_sim_build::Compiler},
    {"compiler_flags",shield_sim_build::CompilerFlags},
    {"geant4_data_libraries",shield_sim::Json(dataLibraries)},
    {"geant4_version",G4Version},{"geometry_hash",c.geometryHash},
    {"ion_physics",c.ionPhysics},{"material_library_hash",c.materials->ContentHash()},
    {"max_step_mm",c.maxStepMm},{"physics_list",c.physicsList},
    {"production_cut_mm",c.productionCutMm},
    {"scoring_contract","shieldSim.transport_raw.v1"},
    {"source_tree_state",shield_sim_build::GitState}
  };
  const std::string responseKey=shield_sim::Sha256(
      shield_sim::Json(keyDocument).DumpCanonical());

  shield_sim::Json::Array scores;
  for(const auto& item:fRawScores){
    const double n=nEv;
    const double mean=n>0.0?item.second.sum/n:0.0;
    const double variance=n>1.0?std::max(0.0,(item.second.sumSquares-
        item.second.sum*item.second.sum/n)/(n-1.0)):0.0;
    const auto unit=fRawScoreUnits.find(item.first);
    scores.emplace_back(shield_sim::Json::Object{
      {"histories",n},{"mean_per_primary",mean},{"sample_variance",variance},
      {"score_id",item.first},{"standard_error",n>0.0?std::sqrt(variance/n):0.0},
      {"sum",item.second.sum},{"sum_squares",item.second.sumSquares},
      {"unit",unit==fRawScoreUnits.end()?"declared_by_score":unit->second}});
  }

  shield_sim::Json::Array secondaryLedgers;
  for(const auto& item:fSecondaryLedgers){
    const auto& ledger=item.second;
    secondaryLedgers.emplace_back(shield_sim::Json::Object{
      {"A",ledger.massNumber},{"Z",ledger.atomicNumber},
      {"count",static_cast<double>(ledger.count)},
      {"creator_process",ledger.creatorProcess},{"generation",ledger.generation},
      {"initial_kinetic_energy_sum_MeV",ledger.initialKineticEnergyMeV},
      {"layer_id",c.shieldLayers.at(ledger.layer).id},{"particle",ledger.particle},
      {"pdg",ledger.pdg},{"root_particle",ledger.rootParticle}});
  }

  shield_sim::Json cell(shield_sim::Json::Object{
    {"absorber_id",absorber.id},{"case_hash",c.caseHash},
    {"histories",static_cast<double>(nEv)},
    {"incident",shield_sim::Json(shield_sim::Json::Object{
      {"A",incident.massNumber},{"Z",incident.atomicNumber},
      {"energy_coordinate",incident.energyCoordinate},
      {"energy_coordinate_unit",incident.energyCoordinateUnit},
      {"energy_coordinate_value",incident.energyCoordinateValue},
      {"incident_id",incident.incidentId},{"mu",incident.mu},
      {"particle",incident.particleName},{"phi_rad",incident.phiRad},
      {"total_kinetic_energy_MeV",incident.totalKineticEnergyMeV}})},
    {"normalization","mean response per generated primary crossing the declared finite source plane"},
    {"random_seed",static_cast<double>(c.randomSeed)},
    {"response_key",responseKey},{"response_key_document",shield_sim::Json(std::move(keyDocument))},
    {"schema_version","shieldSim.response_cell.v1"},
    {"score_covariance",BuildRawScoreCovarianceJson(nEv)},
    {"scores",shield_sim::Json(std::move(scores))},
    {"secondary_production_by_layer",shield_sim::Json(std::move(secondaryLedgers))},
    {"validation_status","implemented_stage3_transport_raw_not_physically_validated"}
  });
  const std::string path=fOpts.outputPrefix+"_response_cell.json";
  std::ofstream out(path);
  if(!out)
    G4Exception("RunAction::WriteResponseCellJson","Output",FatalException,
                ("Cannot write response cell JSON '"+path+"'").c_str());
  out<<cell.DumpCanonical()<<'\n';
  G4cout<<"Response cell written: "<<path<<G4endl;
}
