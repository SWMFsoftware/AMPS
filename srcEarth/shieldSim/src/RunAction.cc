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
#include "Stage6Physics.hh"

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
  if(fOpts.caseDefinition && fOpts.caseDefinition->letResponse.enabled){
    const std::size_t bins=fOpts.caseDefinition->letResponse.edgesKeVPerUm.size()-1;
    fEventLETBins.assign(bins,0.0);
    fLETBins.assign(bins,RawMoments{});
  }

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

void RunAction::RecordAbsorberStep(std::size_t absorber,
                                   const std::string& particle,bool primary,
                                   const std::string& rootParticle,int generation,
                                   const std::string& creatorProcess,
                                   const std::string& stepDefiningProcess,
                                   G4double totalEdep,G4double nonIonizingEdep){
  if(absorber>=fEdepTotal.size() || (totalEdep==0.0 && nonIonizingEdep==0.0)) return;
  AddEdep(absorber,totalEdep,nonIonizingEdep);

  const std::string absorberId=(fOpts.caseDefinition &&
      fOpts.activeAbsorberIndex<fOpts.caseDefinition->absorbers.size())
      ? fOpts.caseDefinition->absorbers[fOpts.activeAbsorberIndex].id
      : std::to_string(absorber);
  const double total=totalEdep/MeV;
  const double nonIon=std::min(totalEdep,nonIonizingEdep)/MeV;
  const double values[3]={total,total-nonIon,nonIon};
  const char* names[3]={"total_deposited_energy_MeV",
                        "ionizing_deposited_energy_MeV",
                        "nonionizing_deposited_energy_MeV"};
  auto addTriplet=[&](const std::string& base){
    for(int i=0;i<3;++i){
      fEventRawScores[base+"/"+names[i]]+=values[i];
      fRawScoreUnits[base+"/"+names[i]]="MeV/primary";
    }
  };

  // These two views are authoritative classifications of the transported
  // particle depositing the energy. They are deliberately based on track
  // identity and propagated ancestry, not on a process-name interpretation.
  const std::string component=primary ? "primary" : "secondary";
  const std::string authoritative="absorber/"+ScoreToken(absorberId)+
      "/deposition_authoritative/by_particle/"+ScoreToken(particle)+
      "/"+component;
  addTriplet(authoritative);
  addTriplet("absorber/"+ScoreToken(absorberId)+
      "/deposition_authoritative/by_ancestry/root_"+ScoreToken(rootParticle)+
      "/generation_"+std::to_string(generation)+"/particle_"+ScoreToken(particle));

  // Process names answer diagnostic questions about which process created a
  // track or limited a step. They do not provide a unique causal allocation of
  // deposited energy, so these scores are explicitly segregated from the
  // authoritative species/ancestry breakdown above.
  addTriplet("absorber/"+ScoreToken(absorberId)+
      "/deposition_diagnostic_by_creator_process/"+ScoreToken(creatorProcess));
  addTriplet("absorber/"+ScoreToken(absorberId)+
      "/deposition_diagnostic_by_step_defining_process/"+
      ScoreToken(stepDefiningProcess));
}

void RunAction::RecordAbsorberEntry(std::size_t absorber,
                                    const std::string& particle,
                                    G4double kineticEnergyMeV,
                                    G4double fluenceWeightPerCm2){
  if(!fOpts.caseDefinition || !fOpts.caseDefinition->damageResponse ||
     absorber!=0 || absorber>=fScoringLVs.size()) return;
  const auto& c=*fOpts.caseDefinition;
  const auto& definition=c.absorbers.at(fOpts.activeAbsorberIndex);
  const std::string base="absorber/"+ScoreToken(definition.id)+
      "/entry/by_particle/"+ScoreToken(particle);
  fEventRawScores[base+"/crossing_count"]+=1.0;
  fEventRawScores[base+"/kinetic_energy_sum_MeV"]+=kineticEnergyMeV;
  fEventRawScores[base+"/scalar_fluence_cm2"]+=fluenceWeightPerCm2;
  fRawScoreUnits[base+"/crossing_count"]="crossings/primary";
  fRawScoreUnits[base+"/kinetic_energy_sum_MeV"]="MeV/primary";
  fRawScoreUnits[base+"/scalar_fluence_cm2"]="cm^-2/primary";

  const bool dddRequested=fOpts.calcDDD;
  const bool neqRequested=fOpts.calcNEq;
  const auto response=c.damageResponse->Evaluate(particle,kineticEnergyMeV);
  if(response.status==shield_sim::DamageResponseStatus::UnsupportedSpecies){
    // n_eq is explicitly a neutron-only quantity, so a proton, ion, electron,
    // or photon crossing is not an omitted n_eq term. For total supported-
    // particle DDD, however, an unrepresented species makes the requested sum
    // incomplete and therefore not assessable.
    if(dddRequested){
      ++fDddUnsupportedCrossings;
      fEventRawScores[base+"/unsupported_DDD_crossing_count"]+=1.0;
      fRawScoreUnits[base+"/unsupported_DDD_crossing_count"]="crossings/primary";
    }
    return;
  }
  if(response.status==shield_sim::DamageResponseStatus::OutOfDomain){
    if(dddRequested){
      ++fDddOutOfDomainCrossings;
      fEventRawScores[base+"/out_of_domain_DDD_crossing_count"]+=1.0;
      fRawScoreUnits[base+"/out_of_domain_DDD_crossing_count"]="crossings/primary";
    }
    if(neqRequested && particle=="neutron"){
      ++fNeqOutOfDomainCrossings;
      fEventRawScores[base+"/out_of_domain_n_eq_crossing_count"]+=1.0;
      fRawScoreUnits[base+"/out_of_domain_n_eq_crossing_count"]="crossings/primary";
    }
    return;
  }

  // The geometry provider supplies scalar-fluence weight. For slab.v1 it is
  // 1/(entryArea*localDirection.z), including the required projected-cosine
  // correction. A future curved/segmented provider owns its local surface
  // normal and patch weighting; the quantity layer remains shape-neutral.
  if(!(fluenceWeightPerCm2>0.0) || !std::isfinite(fluenceWeightPerCm2))
    G4Exception("RunAction::RecordAbsorberEntry","Geometry",FatalException,
                "Stage-5 absorber-entry fluence weight must be finite and positive");
  const double ddd=response.nielMeVCm2PerG*fluenceWeightPerCm2;
  if(dddRequested){
    fEventDDDMeVPerG+=ddd;
    fEventRawScores[base+"/DDD_MeV_g"]+=ddd;
    fRawScoreUnits[base+"/DDD_MeV_g"]="MeV/g/primary";
  }
  if(neqRequested && particle=="neutron"){
    const double equivalent=response.neutronHardnessFactor*fluenceWeightPerCm2;
    fEventNeutronEquivalentPerCm2+=equivalent;
    fEventRawScores[base+"/n_eq_cm2"]+=equivalent;
    fRawScoreUnits[base+"/n_eq_cm2"]="cm^-2/primary";
  }
}

void RunAction::RecordAbsorberLETStep(std::size_t absorber,
                                      const std::string& particle,int pdg,
                                      int atomicNumber,int massNumber,
                                      G4double charge,bool primary,
                                      const std::string& rootParticle,
                                      int generation,
                                      G4double linearLETKeVPerUm,
                                      G4double stepLength){
  if(!LETScoringEnabled() || absorber>=fScoringLVs.size()) return;
  std::ostringstream keyStream;
  keyStream<<ScoreToken(particle)<<'|'<<pdg<<'|'<<atomicNumber<<'|'<<massNumber
           <<'|'<<std::setprecision(12)<<charge<<'|'<<(primary?1:0)<<'|'
           <<ScoreToken(rootParticle)<<'|'<<generation;
  const std::string key=keyStream.str();
  auto found=fEventLETComponents.find(key);
  if(found==fEventLETComponents.end()){
    EventLETComponent component;
    component.metadata.particle=particle;
    component.metadata.pdg=pdg;
    component.metadata.atomicNumber=atomicNumber;
    component.metadata.massNumber=massNumber;
    component.metadata.charge=charge;
    component.metadata.primary=primary;
    component.metadata.rootParticle=rootParticle;
    component.metadata.generation=generation;
    component.bins.assign(fEventLETBins.size(),0.0);
    found=fEventLETComponents.emplace(key,std::move(component)).first;
  }
  if(!(stepLength>0.0)){
    ++fLETZeroLengthSteps;
    ++found->second.metadata.zeroLengthSteps;
    return;
  }
  if(!std::isfinite(linearLETKeVPerUm) || linearLETKeVPerUm<0.0)
    G4Exception("RunAction::RecordAbsorberLETStep","LET",FatalException,
                "Electronic LET must be finite and non-negative");
  auto* logical=fScoringLVs[absorber];
  const double volumeCm3=(logical && logical->GetSolid())
      ? logical->GetSolid()->GetCubicVolume()/cm3 : 0.0;
  if(!(volumeCm3>0.0) || !std::isfinite(volumeCm3))
    G4Exception("RunAction::RecordAbsorberLETStep","Geometry",FatalException,
                "LET scoring requires a finite positive absorber volume");
  const double trackFluence=(stepLength/cm)/volumeCm3;
  ++fLETPositiveLengthSteps;
  ++found->second.metadata.positiveLengthSteps;
  fEventLETT0+=trackFluence;
  fEventLETT1+=linearLETKeVPerUm*trackFluence;
  fEventLETT2+=linearLETKeVPerUm*linearLETKeVPerUm*trackFluence;
  found->second.t0+=trackFluence;
  found->second.t1+=linearLETKeVPerUm*trackFluence;
  found->second.t2+=linearLETKeVPerUm*linearLETKeVPerUm*trackFluence;

  const shield_sim::LogarithmicLETBinning binning(
      fOpts.caseDefinition->letResponse.minimumKeVPerUm,
      fOpts.caseDefinition->letResponse.maximumKeVPerUm,
      fOpts.caseDefinition->letResponse.binsPerDecade);
  const long bin=binning.Locate(linearLETKeVPerUm);
  if(bin<0){
    fEventLETUnderflow+=trackFluence;
    found->second.underflow+=trackFluence;
  } else if(static_cast<std::size_t>(bin)>=fEventLETBins.size()){
    fEventLETOverflow+=trackFluence;
    found->second.overflow+=trackFluence;
  } else {
    fEventLETBins[static_cast<std::size_t>(bin)]+=trackFluence;
    found->second.bins[static_cast<std::size_t>(bin)]+=trackFluence;
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

void RunAction::AddPrimaryProtonExit(G4double kineticEnergyMeV,
                                     G4double fluenceWeightPerCm2){
  // Thresholds are evaluated at the exact crossing energy, not at an output
  // histogram bin center. The approved definition uses strict E > threshold.
  // Every crossing is a surface-fluence contribution. Multiple crossings in
  // one primary history remain correlated and are therefore accumulated into
  // the event numerator/denominator before their covariance is committed.
  if(!(fluenceWeightPerCm2>0.0) || !std::isfinite(fluenceWeightPerCm2))
    G4Exception("RunAction::AddPrimaryProtonExit","Geometry",FatalException,
                "Shield-exit fluence weight must be finite and positive");
  const bool stage6=fOpts.caseDefinition &&
      fOpts.caseDefinition->schemaVersion=="shieldSim.case.v2";
  if(kineticEnergyMeV>10.0 && (stage6 || fEventPrimaryProtonAbove10==0.0)){
    ++fPrimaryProtonAbove10;
    fEventPrimaryProtonAbove10+=stage6?fluenceWeightPerCm2:1.0;
  }
  if(kineticEnergyMeV>100.0 && (stage6 || fEventPrimaryProtonAbove100==0.0)){
    ++fPrimaryProtonAbove100;
    fEventPrimaryProtonAbove100+=stage6?fluenceWeightPerCm2:1.0;
  }
}

void RunAction::BeginEvent(){
  std::fill(fEventEdepIonizing.begin(),fEventEdepIonizing.end(),0.0);
  std::fill(fEventEdepTotal.begin(),fEventEdepTotal.end(),0.0);
  std::fill(fEventEdepNonIonizing.begin(),fEventEdepNonIonizing.end(),0.0);
  fEventRawScores.clear();
  fEventDDDMeVPerG=0.0;
  fEventNeutronEquivalentPerCm2=0.0;
  fEventPrimaryProtonAbove10=0.0;
  fEventPrimaryProtonAbove100=0.0;
  std::fill(fEventLETBins.begin(),fEventLETBins.end(),0.0);
  fEventLETComponents.clear();
  fEventLETUnderflow=0.0;
  fEventLETOverflow=0.0;
  fEventLETT0=fEventLETT1=fEventLETT2=0.0;
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
  auto addMoment=[](RawMoments& moments,double value){
    moments.sum+=value;
    moments.sumSquares+=value*value;
  };
  for(std::size_t bin=0;bin<fEventLETBins.size();++bin)
    addMoment(fLETBins[bin],fEventLETBins[bin]);
  addMoment(fLETUnderflow,fEventLETUnderflow);
  addMoment(fLETOverflow,fEventLETOverflow);
  auto addRatioEvent=[](shield_sim::RatioMoments& moments,double numerator,
                        double denominator){
    ++moments.histories;
    moments.numeratorSum+=numerator;
    moments.denominatorSum+=denominator;
    moments.numeratorSumSquares+=numerator*numerator;
    moments.denominatorSumSquares+=denominator*denominator;
    moments.numeratorDenominatorSumProducts+=numerator*denominator;
  };
  for(const auto& item:fEventLETComponents){
    auto inserted=fLETComponents.emplace(item.first,item.second.metadata);
    auto& ledger=inserted.first->second;
    if(inserted.second){
      ledger.bins.assign(fEventLETBins.size(),RawMoments{});
      ledger.positiveLengthSteps=0;
      ledger.zeroLengthSteps=0;
      // Components first created after earlier primary histories receive the
      // mathematically required implicit-zero histories.
      ledger.trackMeanMoments.histories=fLETTrackMeanMoments.histories;
      ledger.doseMeanMoments.histories=fLETDoseMeanMoments.histories;
    }
    for(std::size_t bin=0;bin<item.second.bins.size();++bin)
      addMoment(ledger.bins[bin],item.second.bins[bin]);
    addMoment(ledger.underflow,item.second.underflow);
    addMoment(ledger.overflow,item.second.overflow);
    // Counts are cumulative diagnostics, not stochastic scores.
    ledger.positiveLengthSteps+=item.second.metadata.positiveLengthSteps;
    ledger.zeroLengthSteps+=item.second.metadata.zeroLengthSteps;
  }
  for(auto& item:fLETComponents){
    const auto event=fEventLETComponents.find(item.first);
    const double t0=event==fEventLETComponents.end()?0.0:event->second.t0;
    const double t1=event==fEventLETComponents.end()?0.0:event->second.t1;
    const double t2=event==fEventLETComponents.end()?0.0:event->second.t2;
    addRatioEvent(item.second.trackMeanMoments,t1,t0);
    addRatioEvent(item.second.doseMeanMoments,t2,t1);
  }
  if(LETScoringEnabled()){
    addRatioEvent(fLETTrackMeanMoments,fEventLETT1,fEventLETT0);
    addRatioEvent(fLETDoseMeanMoments,fEventLETT2,fEventLETT1);
  }
  if(fOpts.calcHardness && fOpts.caseDefinition &&
     fOpts.caseDefinition->schemaVersion=="shieldSim.case.v2")
    addRatioEvent(fHardnessMoments,fEventPrimaryProtonAbove100,
                  fEventPrimaryProtonAbove10);
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
  if(fOpts.calcDDD) fAccuracyController->AddSample("DDD",fEventDDDMeVPerG);
  if(fOpts.calcNEq) fAccuracyController->AddSample("n_eq",fEventNeutronEquivalentPerCm2);
  // LET means and H100/10 are correlated ratios. Their precision is evaluated
  // directly from event-level sufficient moments in AccuracyStatisticallySatisfied;
  // feeding either into the scalar mean controller would discard covariance.
  if(fOpts.calcHardness && fOpts.caseDefinition &&
     fOpts.caseDefinition->schemaVersion!="shieldSim.case.v2" &&
     fEventPrimaryProtonAbove10>0.0)
    fAccuracyController->AddSample("H100/10",
        fEventPrimaryProtonAbove100>0.0?1.0:0.0);
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
    if(quantity=="H100/10" && fOpts.caseDefinition->schemaVersion=="shieldSim.case.v2"){
      const auto estimate=shield_sim::EstimateRatioOfMeans(fHardnessMoments);
      if(estimate.status!=shield_sim::RatioStatus::Available || estimate.value==0.0 ||
         estimate.standardError/std::abs(estimate.value)>requirement.target) return false;
      continue;
    }
    if(quantity=="LET"){
      const auto track=shield_sim::EstimateRatioOfMeans(fLETTrackMeanMoments);
      const auto dose=shield_sim::EstimateRatioOfMeans(fLETDoseMeanMoments);
      if(track.status!=shield_sim::RatioStatus::Available ||
         dose.status!=shield_sim::RatioStatus::Available || track.value==0.0 ||
         dose.value==0.0 || track.standardError/std::abs(track.value)>requirement.target ||
         dose.standardError/std::abs(dose.value)>requirement.target) return false;
      continue;
    }
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
  fDddOutOfDomainCrossings=0;
  fDddUnsupportedCrossings=0;
  fNeqOutOfDomainCrossings=0;
  std::fill(fEventLETBins.begin(),fEventLETBins.end(),0.0);
  std::fill(fLETBins.begin(),fLETBins.end(),RawMoments{});
  fEventLETComponents.clear();
  fLETComponents.clear();
  fLETUnderflow=RawMoments{};
  fLETOverflow=RawMoments{};
  fLETPositiveLengthSteps=0;
  fLETZeroLengthSteps=0;
  fLETTrackMeanMoments=shield_sim::RatioMoments{};
  fLETDoseMeanMoments=shield_sim::RatioMoments{};
  fHardnessMoments=shield_sim::RatioMoments{};
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
bool RunAction::AcceptsAbsorberEntry(
    std::size_t absorber,const G4StepPoint& postStepPoint) const {
  return fDetector && fDetector->IsAbsorberEntry(absorber,postStepPoint);
}
G4double RunAction::AbsorberEntryFluenceWeightPerCm2(
    std::size_t absorber,const G4StepPoint& postStepPoint) const {
  return fDetector
      ? fDetector->GetAbsorberEntryFluenceWeightPerCM2(absorber,postStepPoint) : 0.0;
}

G4double RunAction::ShieldExitFluenceWeightPerCm2(
    G4double localDirectionCosine) const {
  if(!(localDirectionCosine>0.0) || !fShieldLV) return 0.0;
  const auto* solid=dynamic_cast<const G4Box*>(fShieldLV->GetSolid());
  if(!solid) return 0.0;
  const double areaCm2=(2.0*solid->GetXHalfLength()/cm)*
                       (2.0*solid->GetYHalfLength()/cm);
  return areaCm2>0.0 ? 1.0/(areaCm2*std::abs(localDirectionCosine)) : 0.0;
}

bool RunAction::LETScoringEnabled() const {
  return fOpts.calcLET && fOpts.caseDefinition &&
         fOpts.caseDefinition->letResponse.enabled;
}

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
  fDddOutOfDomainCrossings=0;
  fDddUnsupportedCrossings=0;
  fNeqOutOfDomainCrossings=0;
  std::fill(fEventLETBins.begin(),fEventLETBins.end(),0.0);
  std::fill(fLETBins.begin(),fLETBins.end(),RawMoments{});
  fEventLETComponents.clear();
  fLETComponents.clear();
  fLETUnderflow=RawMoments{};
  fLETOverflow=RawMoments{};
  fLETPositiveLengthSteps=0;
  fLETZeroLengthSteps=0;
  fLETTrackMeanMoments=shield_sim::RatioMoments{};
  fLETDoseMeanMoments=shield_sim::RatioMoments{};
  fHardnessMoments=shield_sim::RatioMoments{};
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
         <<" Stage-3/4/5/6 verification source; one generated primary per history.\n";
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
  // Versioned Stage-6 cases use scalar-fluence crossing weights and a ratio of
  // event means. The point estimate is emitted only when Fieller's denominator
  // confidence set excludes zero; unresolved/zero denominators remain NaN.
  if(fOpts.caseDefinition && fOpts.caseDefinition->schemaVersion=="shieldSim.case.v2"){
    const auto estimate=shield_sim::EstimateRatioOfMeans(fHardnessMoments);
    return estimate.status==shield_sim::RatioStatus::Available
        ? estimate.value : std::numeric_limits<G4double>::quiet_NaN();
  }
  // Legacy CLI compatibility retains its historical unweighted count ratio.
  if(nEv<=0 || fPrimaryProtonAbove10<=0.0)
    return std::numeric_limits<G4double>::quiet_NaN();
  return fPrimaryProtonAbove100/fPrimaryProtonAbove10;
}

void RunAction::ComputeSpectrumFoldedQuantities(G4int nEv){
  // Versioned Stage-5 cases never enter the legacy shield-exit histogram
  // approximation below. Their authoritative response is accumulated at exact
  // absorber-entry crossing energies and uses only the authenticated tables in
  // DamageResponseLibrary. This branch is intentionally early so the analytic
  // compatibility surrogate cannot overwrite a Stage-5 result.
  if(fOpts.caseDefinition && fOpts.caseDefinition->damageResponse){
    const G4double incidentRate=GetIncidentParticleRate();
    fLastHardness=fOpts.calcHardness?ComputeHardnessIndex(nEv):0.0;
    const auto* dddEstimate=fAccuracyController?fAccuracyController->Estimate("DDD"):nullptr;
    const auto* neqEstimate=fAccuracyController?fAccuracyController->Estimate("n_eq"):nullptr;
    const bool dddAssessable=fDddOutOfDomainCrossings==0 && fDddUnsupportedCrossings==0;
    const bool neqAssessable=fNeqOutOfDomainCrossings==0;
    for(std::size_t target=0;target<fScoringLVs.size();++target){
      fLastDDD[target]=(fOpts.calcDDD && dddEstimate && dddAssessable)
          ? dddEstimate->mean : (fOpts.calcDDD
              ? std::numeric_limits<G4double>::quiet_NaN() : 0.0);
      fLastDDDRate[target]=fOpts.calcDDD?fLastDDD[target]*incidentRate:0.0;
      fLastNEq[target]=(fOpts.calcNEq && neqEstimate && neqAssessable)
          ? neqEstimate->mean : (fOpts.calcNEq
              ? std::numeric_limits<G4double>::quiet_NaN() : 0.0);
      fLastNEqRate[target]=fOpts.calcNEq?fLastNEq[target]*incidentRate:0.0;
    }
    G4cout<<"Stage-5 absorber-entry quantities:"<<G4endl;
    if(fOpts.calcDDD)
      G4cout<<"  DDD="<<fLastDDD.front()<<" MeV/g/primary; out_of_domain_crossings="
            <<fDddOutOfDomainCrossings<<"; unsupported_species_crossings="
            <<fDddUnsupportedCrossings<<G4endl;
    if(fOpts.calcNEq)
      G4cout<<"  neutron-only ASTM n_eq="<<fLastNEq.front()
            <<" cm^-2/primary; out_of_domain_neutrons="
            <<fNeqOutOfDomainCrossings<<G4endl;
    return;
  }

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
    if(fOpts.caseDefinition && fOpts.caseDefinition->damageResponse){
      out<<"# Stage-5 DDD equation: absorber-entry fluence at exact crossing energy times pinned natural-silicon NIEL.\n";
      out<<"# Stage-5 n_eq: ASTM E722-19 neutron contribution only, divided by the 95 MeV-mb silicon reference; charged-particle DDD is excluded.\n";
      out<<"# Any nonzero unsupported/out-of-domain entry makes the affected value not assessable; no extrapolation or clipping is used.\n";
    } else {
      out<<"# Legacy DDD equation: D_d = integral Phi(E) NIEL(E) dE using the compatibility-only shield-exit analytic surrogate.\n";
      out<<"# Legacy n_eq divides aggregate surrogate DDD by a material proxy; it is not the Stage-5 ASTM E722 result.\n";
    }
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
  if(LETScoringEnabled()){
    const std::string fname=OutputName("_let_spectrum.dat");
    std::ofstream out(fname);
    if(!out)
      G4Exception("RunAction::WriteLETSpectrumTecplot","Output",FatalException,
                  ("Cannot write LET spectrum '"+fname+"'").c_str());
    const auto& response=fOpts.caseDefinition->letResponse;
    const auto& absorber=fOpts.caseDefinition->absorbers.at(fOpts.activeAbsorberIndex);
    const shield_sim::LogarithmicLETBinning bins(response.minimumKeVPerUm,
        response.maximumKeVPerUm,response.binsPerDecade);
    out<<"TITLE = \"ShieldSim Stage-6 Step-Fluence LET Spectrum\"\n";
    out<<"# Definition: unrestricted Geant4 electronic stopping at transported-step midpoint.\n";
    out<<"# Weight: step track length / absorber volume; delta electrons are separate components.\n";
    out<<"# Out-of-range policy: not_assessable_no_clamping.\n";
    out<<"# Underflow fluence sum [cm-2 primary-1]: "<<fLETUnderflow.sum/nEv<<"\n";
    out<<"# Overflow fluence sum [cm-2 primary-1]: "<<fLETOverflow.sum/nEv<<"\n";
    out<<"VARIABLES = \"LET_low [keV/um]\" \"LET_high [keV/um]\" "
          "\"LET_center [keV/um]\" \"mass_LET_center [MeV cm2/g]\" "
          "\"fluence [cm-2 primary-1]\" \"fluence_SE [cm-2 primary-1]\" "
          "\"dPhi_dLET [cm-2 primary-1 (keV/um)-1]\" "
          "\"fluence_rate [cm-2 s-1]\" \"dPhi_dLET_rate [cm-2 s-1 (keV/um)-1]\"\n";
    auto writeZone=[&](const std::string& title,const std::vector<RawMoments>& moments){
      out<<"ZONE T=\""<<title<<"\", I="<<bins.BinCount()<<", DATAPACKING=POINT\n";
      for(std::size_t i=0;i<bins.BinCount();++i){
        const double mean=moments[i].sum/nEv;
        const double variance=nEv>1?std::max(0.0,(moments[i].sumSquares-
            moments[i].sum*moments[i].sum/nEv)/(nEv-1.0)):0.0;
        const double se=std::sqrt(variance/nEv);
        const double center=bins.GeometricCenter(i);
        out<<std::scientific<<std::setprecision(12)
           <<bins.Lower(i)<<' '<<bins.Upper(i)<<' '<<center<<' '
           <<10.0*center/absorber.densityGPerCm3<<' '
           <<mean<<' '<<se<<' '<<mean/bins.Width(i)<<' '
           <<mean*GetIncidentParticleRate()<<' '
           <<mean*GetIncidentParticleRate()/bins.Width(i)<<'\n';
      }
    };
    writeZone("total",fLETBins);
    for(const auto& item:fLETComponents){
      const auto& component=item.second;
      std::ostringstream title;
      title<<component.particle<<";q="<<component.charge
           <<";"<<(component.primary?"primary":"secondary")
           <<";root="<<component.rootParticle<<";generation="<<component.generation;
      writeZone(title.str(),component.bins);
    }
    G4cout<<"Stage-6 LET spectrum written: "<<fname<<G4endl;
    return;
  }
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
  if(fOpts.caseDefinition){
    out<<"meta em_physics "<<fOpts.caseDefinition->emPhysics<<"\n";
    out<<"meta ion_physics "<<fOpts.caseDefinition->ionPhysics<<"\n";
  }
  out<<"meta source_mode "<<fSourceMode<<"\n";
  out<<"meta spectrum_file "
     <<(fOpts.caseDefinition && (fOpts.caseDefinition->responseIncident.enabled ||
                                 !fOpts.caseDefinition->directIncidents.empty())
          ? fOpts.caseDefinition->sourceFormat
          : fOpts.spectrumFile.empty()?"builtin":fOpts.spectrumFile)<<"\n";
  out<<"meta shield_material "<<(fSweepMode?fCurrentMat:fOpts.shieldMaterial)<<"\n";
  out<<"meta result_capability_status "
     <<(fOpts.caseDefinition
          ? (fOpts.caseDefinition->absorbers.at(fOpts.activeAbsorberIndex).provider!=
                 "slab.v1"
               ? "stage10_geometry_providers_verification_only_incomplete_no_independent_nontrivial_shape_validation"
             : fOpts.caseDefinition->schemaVersion=="shieldSim.case.v2" &&
             std::find(fOpts.caseDefinition->quantities.begin(),
                       fOpts.caseDefinition->quantities.end(),"TID")!=
                 fOpts.caseDefinition->quantities.end()
               ? "stage4_finite_slab_tid_implemented_and_acceptance_tested_pending_independent_review"
             : fOpts.caseDefinition->schemaVersion=="shieldSim.case.v2"
               ? "stage3_transport_engine_implemented_and_acceptance_tested_pending_independent_review_no_science_quantity_claim"
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
         <<assessment.combinedInterval95Low<<' '<<assessment.combinedInterval95High<<' '
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
  const bool stage4Tid=c.schemaVersion=="shieldSim.case.v2" &&
      std::find(c.quantities.begin(),c.quantities.end(),"TID")!=c.quantities.end();
  const bool stage5Damage=c.schemaVersion=="shieldSim.case.v2" && c.damageResponse &&
      (std::find(c.quantities.begin(),c.quantities.end(),"DDD")!=c.quantities.end() ||
       std::find(c.quantities.begin(),c.quantities.end(),"n_eq")!=c.quantities.end());
  const bool stage6=c.schemaVersion=="shieldSim.case.v2" &&
      (std::find(c.quantities.begin(),c.quantities.end(),"LET")!=c.quantities.end() ||
       std::find(c.quantities.begin(),c.quantities.end(),"H100/10")!=c.quantities.end());
  const std::string stage4ValidationStatus=
      "stage4_finite_slab_tid_implemented_and_acceptance_tested_pending_independent_review";
  const std::string stage5ValidationStatus=
      "stage5_silicon_DDD_and_neutron_only_n_eq_implemented_and_reference_tested_pending_accelerator_device_validation_and_independent_review";
  const std::string stage6ValidationStatus=
      "stage6_LET_and_H100_10_implemented_and_reference_tested_pending_independent_review";
  const auto& absorber=c.absorbers.at(fOpts.activeAbsorberIndex);
  const bool stage10Geometry=absorber.provider!="slab.v1";
  const std::string stage10ValidationStatus=
      "stage10_geometry_providers_verification_only_incomplete_no_independent_nontrivial_shape_validation";
  const std::string validationStatus = stage10Geometry
      ? stage10ValidationStatus
      : stage6
      ? stage6ValidationStatus
      : stage5Damage
      ? stage5ValidationStatus
      : stage4Tid
      ? stage4ValidationStatus
      : c.schemaVersion=="shieldSim.case.v2"
      ? "implemented_stage3_transport_raw_acceptance_tested_pending_independent_review"
      : c.sourceFormat=="ams02_hapi_component"
      ? "stage2_published_spectrum_transport_test_not_physical_validation"
      : "implemented_not_physically_validated_stage1";
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
    if(!requirement.systematicComponents.empty()){
      shield_sim::Json::Array components;
      for(const auto& component:requirement.systematicComponents)
        components.emplace_back(shield_sim::Json::Object{
          {"category",component.category},{"id",component.id},
          {"provenance",component.provenance},
          {"relative_standard_uncertainty",component.relativeStandardUncertainty}});
      shield_sim::Json::Array correlations;
      for(const auto& correlation:requirement.systematicCorrelations)
        correlations.emplace_back(shield_sim::Json::Object{
          {"component_i",correlation.firstComponent},
          {"component_j",correlation.secondComponent},
          {"correlation_coefficient",correlation.coefficient},
          {"rationale",correlation.rationale}});
      object["systematic_components"]=shield_sim::Json(std::move(components));
      object["systematic_correlations"]=shield_sim::Json(std::move(correlations));
      object["systematic_pair_policy"]=requirement.systematicPairPolicy;
    }
    if(requirement.hasKnownSystematicFloor)
      object["known_systematic_relative_floor"]=requirement.knownSystematicFloor;
    return shield_sim::Json(std::move(object));
  };
  auto ratioAccuracyJson=[&](const std::string& quantity,
      const std::vector<shield_sim::RatioEstimate>& estimates){
    shield_sim::Json result=accuracyJson(quantity,fAccuracyController->Assess(quantity));
    auto& object=result.AsObject();
    const auto& requirement=c.AccuracyFor(quantity);
    bool available=!estimates.empty();
    double maximumRelative=0.0;
    std::string reason;
    for(const auto& estimate:estimates){
      if(estimate.status!=shield_sim::RatioStatus::Available || estimate.value==0.0){
        available=false;
        if(!reason.empty()) reason+="; ";
        reason+=estimate.reason;
      } else {
        maximumRelative=std::max(maximumRelative,
            estimate.standardError/std::abs(estimate.value));
      }
    }
    object["statistical_relative_standard_uncertainty"]=available
        ? shield_sim::Json(maximumRelative) : shield_sim::Json(nullptr);
    object["ratio_interval_method"]="two-sided 95% Fieller interval using event-level numerator/denominator covariance";
    if(!available){
      object["status"]="not_assessable";
      object["reason"]=reason.empty()?"ratio estimate is unavailable":reason;
      object["combined_relative_standard_uncertainty"]=shield_sim::Json(nullptr);
      return result;
    }
    if(requirement.scope=="combined" && !requirement.systematicBudgetComplete){
      object["status"]="not_assessable";
      object["reason"]="combined accuracy requested but the systematic uncertainty budget is incomplete";
      object["combined_relative_standard_uncertainty"]=shield_sim::Json(nullptr);
      return result;
    }
    double assessed=maximumRelative;
    if(requirement.scope=="combined"){
      const double systematic=std::max(requirement.systematicRelativeUncertainty,
          requirement.hasKnownSystematicFloor?requirement.knownSystematicFloor:0.0);
      assessed=std::hypot(maximumRelative,systematic);
      object["combined_relative_standard_uncertainty"]=assessed;
    } else {
      object["combined_relative_standard_uncertainty"]=shield_sim::Json(nullptr);
    }
    object["status"]=assessed<=requirement.target?"met":"not_met";
    object["reason"]=assessed<=requirement.target
        ? "all required ratio means meet the requested relative standard uncertainty"
        : "at least one required ratio mean exceeds the requested relative standard uncertainty";
    return result;
  };
  auto ratioMomentsJson=[&](const shield_sim::RatioMoments& moments,
                            const shield_sim::RatioEstimate& estimate){
    return shield_sim::Json(shield_sim::Json::Object{
      {"covariance_of_means",estimate.covarianceOfMeans},
      {"denominator_mean",estimate.denominatorMean},
      {"denominator_sum",moments.denominatorSum},
      {"denominator_sum_squares",moments.denominatorSumSquares},
      {"denominator_variance_of_mean",estimate.denominatorVarianceOfMean},
      {"histories",static_cast<double>(moments.histories)},
      {"numerator_denominator_sum_products",moments.numeratorDenominatorSumProducts},
      {"numerator_mean",estimate.numeratorMean},
      {"numerator_sum",moments.numeratorSum},
      {"numerator_sum_squares",moments.numeratorSumSquares},
      {"numerator_variance_of_mean",estimate.numeratorVarianceOfMean}});
  };

  shield_sim::Json::Array outputs;
  if(fOpts.calcTID){
    for(std::size_t i=0;i<fScoringLVs.size();++i){
      auto* lv=fScoringLVs[i];
      const G4double mass=lv?lv->GetMass(true,false):0.0;
      auto assessment=fAccuracyController->Assess("TID");
      const auto& requirement=c.AccuracyFor("TID");
      shield_sim::Json::Object interval;
      const double intervalLow=requirement.scope=="combined"
          ? assessment.combinedInterval95Low : assessment.interval95Low;
      const double intervalHigh=requirement.scope=="combined"
          ? assessment.combinedInterval95High : assessment.interval95High;
      if(mass>0.0 && std::isfinite(intervalLow) && std::isfinite(intervalHigh)){
        interval["low"]=intervalLow*MeV/mass/gray;
        interval["high"]=intervalHigh*MeV/mass/gray;
        interval["unit"]="Gy/primary";
      } else {
        interval["low"]=shield_sim::Json(nullptr);
        interval["high"]=shield_sim::Json(nullptr);
        interval["unit"]="Gy/primary";
      }
      shield_sim::Json::Object statisticalInterval;
      if(mass>0.0 && std::isfinite(assessment.interval95Low) &&
         std::isfinite(assessment.interval95High)){
        statisticalInterval["low"]=assessment.interval95Low*MeV/mass/gray;
        statisticalInterval["high"]=assessment.interval95High*MeV/mass/gray;
      } else {
        statisticalInterval["low"]=shield_sim::Json(nullptr);
        statisticalInterval["high"]=shield_sim::Json(nullptr);
      }
      statisticalInterval["unit"]="Gy/primary";
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
        {"statistical_interval_95",shield_sim::Json(std::move(statisticalInterval))},
        {"name","TID"}, {"normalization","per incident primary"},
        {"rad_per_primary",100.0*fLastDose[i]/gray},
        {"source_normalized_unit","Gy/s"},
        {"source_normalized_value",fLastDoseRate[i]/gray},
        {"unit","Gy/primary"}, {"validation_status",stage10Geometry
            ? stage10ValidationStatus : stage4Tid
            ? stage4ValidationStatus : validationStatus},
        {"value",fLastDose[i]/gray}
      });
    }
  }
  if(fOpts.calcDDD){
    const auto* estimate=fAccuracyController->Estimate("DDD");
    auto assessment=fAccuracyController->Assess("DDD");
    const bool domainComplete=fDddOutOfDomainCrossings==0 &&
                              fDddUnsupportedCrossings==0;
    shield_sim::Json accuracy=accuracyJson("DDD",assessment);
    if(!domainComplete){
      accuracy.AsObject()["status"]="not_assessable";
      accuracy.AsObject()["reason"]="nonzero absorber-entry fluence was outside a pinned table domain or belonged to an unsupported species; no contribution was clamped, extrapolated, or silently omitted";
    }
    shield_sim::Json::Object interval{{"unit","MeV/g/primary"}};
    const auto& requirement=c.AccuracyFor("DDD");
    const double low=requirement.scope=="combined"
        ? assessment.combinedInterval95Low : assessment.interval95Low;
    const double high=requirement.scope=="combined"
        ? assessment.combinedInterval95High : assessment.interval95High;
    interval["low"]=(domainComplete && std::isfinite(low))
        ? shield_sim::Json(low) : shield_sim::Json(nullptr);
    interval["high"]=(domainComplete && std::isfinite(high))
        ? shield_sim::Json(high) : shield_sim::Json(nullptr);
    const double partial=estimate?estimate->mean:0.0;
    auto* lv=fScoringLVs.empty()?nullptr:fScoringLVs.front();
    const double massG=lv?lv->GetMass(true,false)/g:0.0;
    const double diagnostic=(massG>0.0 && nEv>0)
        ? fEdepNonIonizing.front()/MeV/massG/nEv : 0.0;
    outputs.emplace_back(shield_sim::Json::Object{
      {"accuracy",std::move(accuracy)},{"absorber_id",absorber.id},
      {"authoritative_method","absorber-entry scalar fluence folded at exact crossing energy with pinned species/target NIEL"},
      {"definition","sum over proton, He-4, and neutron crossings of natural-silicon NIEL divided by entry area and absolute direction/surface-normal cosine"},
      {"diagnostic_Geant4_nonionizing_deposition",shield_sim::Json(shield_sim::Json::Object{
        {"definition","Geant4 nonionizing energy deposited in the finite sensitive region divided by its mass; separate diagnostic, never added to table-folded DDD"},
        {"unit","MeV/g/primary"},{"value",diagnostic}})},
      {"interval_95",shield_sim::Json(std::move(interval))},
      {"name","DDD"},{"normalization","per incident primary"},
      {"out_of_domain_crossings",static_cast<double>(fDddOutOfDomainCrossings)},
      {"partial_supported_species_value",partial},
      {"partial_supported_species_value_unit","MeV/g/primary"},
      {"source_normalized_unit","MeV/g/s"},
      {"source_normalized_value",domainComplete
          ? shield_sim::Json(partial*GetIncidentParticleRate()) : shield_sim::Json(nullptr)},
      {"table_domain_status",domainComplete?"complete":"not_assessable"},
      {"unit","MeV/g/primary"},{"unsupported_species_crossings",
          static_cast<double>(fDddUnsupportedCrossings)},
      {"validation_status",stage10Geometry ? stage10ValidationStatus : stage5ValidationStatus},
      {"value",domainComplete?shield_sim::Json(partial):shield_sim::Json(nullptr)}
    });
  }
  if(fOpts.calcNEq){
    const auto* estimate=fAccuracyController->Estimate("n_eq");
    auto assessment=fAccuracyController->Assess("n_eq");
    const bool domainComplete=fNeqOutOfDomainCrossings==0;
    shield_sim::Json accuracy=accuracyJson("n_eq",assessment);
    if(!domainComplete){
      accuracy.AsObject()["status"]="not_assessable";
      accuracy.AsObject()["reason"]="nonzero neutron absorber-entry fluence was outside the ASTM E722-19 tabulated domain; no contribution was clamped or extrapolated";
    }
    const double partial=estimate?estimate->mean:0.0;
    shield_sim::Json::Object interval{{"unit","cm^-2/primary"}};
    const auto& requirement=c.AccuracyFor("n_eq");
    const double low=requirement.scope=="combined"
        ? assessment.combinedInterval95Low : assessment.interval95Low;
    const double high=requirement.scope=="combined"
        ? assessment.combinedInterval95High : assessment.interval95High;
    interval["low"]=(domainComplete && std::isfinite(low))
        ? shield_sim::Json(low) : shield_sim::Json(nullptr);
    interval["high"]=(domainComplete && std::isfinite(high))
        ? shield_sim::Json(high) : shield_sim::Json(nullptr);
    outputs.emplace_back(shield_sim::Json::Object{
      {"accuracy",std::move(accuracy)},{"absorber_id",absorber.id},
      {"definition","ASTM E722-19 natural-silicon neutron contribution divided groupwise by the prescribed 95 MeV-mb 1-MeV reference"},
      {"excluded_from_definition","proton and He-4 DDD are not converted to E722 neutron-equivalent fluence"},
      {"fluence_normalization","each accepted neutron crossing is weighted by 1/(entry area times absolute direction/surface-normal cosine)"},
      {"interval_95",shield_sim::Json(std::move(interval))},
      {"name","n_eq"},{"normalization","per incident primary"},
      {"out_of_domain_neutron_crossings",static_cast<double>(fNeqOutOfDomainCrossings)},
      {"reference_damage_function_MeV_mb",95.0},
      {"reference_niel_MeV_cm2_g",c.damageResponse->ReferenceNeutronNielMeVCm2PerG()},
      {"source_normalized_unit","cm^-2/s"},
      {"source_normalized_value",domainComplete
          ? shield_sim::Json(partial*GetIncidentParticleRate()) : shield_sim::Json(nullptr)},
      {"table_domain_status",domainComplete?"complete":"not_assessable"},
      {"unit","cm^-2/primary"},{"validation_status",stage10Geometry
          ? stage10ValidationStatus : stage5ValidationStatus},
      {"value",domainComplete?shield_sim::Json(partial):shield_sim::Json(nullptr)}
    });
  }
  if(fOpts.calcLET){
    const auto track=shield_sim::EstimateRatioOfMeans(fLETTrackMeanMoments);
    const auto dose=shield_sim::EstimateRatioOfMeans(fLETDoseMeanMoments);
    const bool completeRange=fLETUnderflow.sum==0.0 && fLETOverflow.sum==0.0;
    shield_sim::Json accuracy=ratioAccuracyJson("LET",{track,dose});
    if(!completeRange){
      accuracy.AsObject()["status"]="not_assessable";
      accuracy.AsObject()["reason"]="nonzero track-length fluence lies outside the declared LET response; no contribution was clamped";
    }
    auto ratioResult=[&](const shield_sim::RatioEstimate& estimate,
                         const shield_sim::RatioMoments& moments){
      const bool available=estimate.status==shield_sim::RatioStatus::Available;
      const double massScale=10.0/absorber.densityGPerCm3;
      return shield_sim::Json(shield_sim::Json::Object{
        {"event_moments",ratioMomentsJson(moments,estimate)},
        {"interval_95_keV_per_um",shield_sim::Json(shield_sim::Json::Object{
          {"high",available?finiteOrNull(estimate.interval95High):shield_sim::Json(nullptr)},
          {"low",available?finiteOrNull(estimate.interval95Low):shield_sim::Json(nullptr)}})},
        {"interval_95_MeV_cm2_per_g",shield_sim::Json(shield_sim::Json::Object{
          {"high",available?finiteOrNull(estimate.interval95High*massScale):shield_sim::Json(nullptr)},
          {"low",available?finiteOrNull(estimate.interval95Low*massScale):shield_sim::Json(nullptr)}})},
        {"reason",estimate.reason},{"standard_error_keV_per_um",available
            ? finiteOrNull(estimate.standardError):shield_sim::Json(nullptr)},
        {"status",shield_sim::RatioStatusName(estimate.status)},
        {"value_keV_per_um",available?finiteOrNull(estimate.value):shield_sim::Json(nullptr)},
        {"value_MeV_cm2_per_g",available
            ? finiteOrNull(estimate.value*massScale):shield_sim::Json(nullptr)}});
    };
    auto binArray=[&](const std::vector<RawMoments>& moments){
      shield_sim::Json::Array result;
      const shield_sim::LogarithmicLETBinning binning(c.letResponse.minimumKeVPerUm,
          c.letResponse.maximumKeVPerUm,c.letResponse.binsPerDecade);
      for(std::size_t i=0;i<moments.size();++i){
        const double mean=moments[i].sum/nEv;
        const double variance=nEv>1?std::max(0.0,(moments[i].sumSquares-
            moments[i].sum*moments[i].sum/nEv)/(nEv-1.0)):0.0;
        result.emplace_back(shield_sim::Json::Object{
          {"fluence_cm2_per_primary",mean},
          {"fluence_rate_cm2_per_s",mean*GetIncidentParticleRate()},
          {"fluence_standard_error_cm2_per_primary",std::sqrt(variance/nEv)},
          {"fluence_rate_standard_error_cm2_per_s",
              std::sqrt(variance/nEv)*GetIncidentParticleRate()},
          {"lower_keV_per_um",binning.Lower(i)},
          {"upper_keV_per_um",binning.Upper(i)}});
      }
      return shield_sim::Json(std::move(result));
    };
    shield_sim::Json::Array components;
    for(const auto& item:fLETComponents){
      const auto& component=item.second;
      const auto componentTrack=shield_sim::EstimateRatioOfMeans(component.trackMeanMoments);
      const auto componentDose=shield_sim::EstimateRatioOfMeans(component.doseMeanMoments);
      components.emplace_back(shield_sim::Json::Object{
        {"A",component.massNumber},{"Z",component.atomicNumber},
        {"bins",binArray(component.bins)},{"charge_e",component.charge},
        {"dose_mean",ratioResult(componentDose,component.doseMeanMoments)},
        {"generation",component.generation},{"particle",component.particle},
        {"pdg",component.pdg},{"positive_length_steps",static_cast<double>(component.positiveLengthSteps)},
        {"primary",component.primary},{"root_particle",component.rootParticle},
        {"track_mean",ratioResult(componentTrack,component.trackMeanMoments)},
        {"underflow_fluence_sum_cm2",component.underflow.sum},
        {"overflow_fluence_sum_cm2",component.overflow.sum},
        {"zero_length_steps",static_cast<double>(component.zeroLengthSteps)}});
    }
    outputs.emplace_back(shield_sim::Json::Object{
      {"accuracy",std::move(accuracy)},{"absorber_id",absorber.id},
      {"bins",binArray(fLETBins)},{"components",shield_sim::Json(std::move(components))},
      {"definition","unrestricted electronic stopping at each positive-length charged-particle absorber-step midpoint, weighted by step track length divided by absorber volume"},
      {"dose_mean",ratioResult(dose,fLETDoseMeanMoments)},
      {"histogram_domain_status",completeRange?"complete":"not_assessable"},
      {"name","LET"},{"out_of_range_policy",c.letResponse.outOfRangePolicy},
      {"overflow_fluence_cm2_per_primary",fLETOverflow.sum/nEv},
      {"overflow_fluence_rate_cm2_per_s",fLETOverflow.sum/nEv*GetIncidentParticleRate()},
      {"positive_length_steps",static_cast<double>(fLETPositiveLengthSteps)},
      {"track_mean",ratioResult(track,fLETTrackMeanMoments)},
      {"underflow_fluence_cm2_per_primary",fLETUnderflow.sum/nEv},
      {"underflow_fluence_rate_cm2_per_s",fLETUnderflow.sum/nEv*GetIncidentParticleRate()},
      {"validation_status",stage10Geometry ? stage10ValidationStatus : stage6ValidationStatus},
      {"zero_length_charged_steps",static_cast<double>(fLETZeroLengthSteps)}});
  }
  if(fOpts.calcHardness && stage6){
    const auto estimate=shield_sim::EstimateRatioOfMeans(fHardnessMoments);
    const bool available=estimate.status==shield_sim::RatioStatus::Available;
    shield_sim::Json accuracy=ratioAccuracyJson("H100/10",{estimate});
    outputs.emplace_back(shield_sim::Json::Object{
      {"accuracy",std::move(accuracy)}, {"absorber_id",absorber.id},
      {"event_moments",ratioMomentsJson(fHardnessMoments,estimate)},
      {"definition","ratio of mean scalar fluence from every forward primary-proton shield-rear-face crossing with strict E>100 MeV and E>10 MeV thresholds"},
      {"denominator_count",fPrimaryProtonAbove10},
      {"fluence_weight","1/(shield rear-face area times absolute local direction/surface-normal cosine)"},
      {"interval_95",shield_sim::Json(shield_sim::Json::Object{
        {"high",available?finiteOrNull(estimate.interval95High):shield_sim::Json(nullptr)},
        {"low",available?finiteOrNull(estimate.interval95Low):shield_sim::Json(nullptr)},
        {"method","Fieller, two-sided 95%, event-level covariance"}})},
      {"name","H100/10"}, {"numerator_count",fPrimaryProtonAbove100},
      {"ratio_status",shield_sim::RatioStatusName(estimate.status)},
      {"reason",estimate.reason},{"unit","1"},
      {"validation_status",stage10Geometry ? stage10ValidationStatus : stage6ValidationStatus},
      {"value",available?finiteOrNull(estimate.value):shield_sim::Json(nullptr)}
    });
  }
  if(fOpts.calcHardness && !stage6){
    auto assessment=fAccuracyController->Assess("H100/10");
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
    accuracy.AsObject()["interval_95_method"]=
        "Wilson score interval for legacy bounded conditional proportion";
    outputs.emplace_back(shield_sim::Json::Object{
      {"accuracy",std::move(accuracy)},{"absorber_id",absorber.id},
      {"covariance_note","legacy event indicators are nested: I(E>100) is a subset of I(E>10)"},
      {"indicator_sample_covariance",nEv>1
          ? (fPrimaryProtonAbove100-
             fPrimaryProtonAbove100*fPrimaryProtonAbove10/nEv)/(nEv-1.0)
          : 0.0},
      {"estimated_fluence_covariance_per_primary_squared",nEv>1
          ? (fPrimaryProtonAbove100-
             fPrimaryProtonAbove100*fPrimaryProtonAbove10/nEv)/((nEv-1.0)*nEv)
          : 0.0},
      {"definition","legacy primary-proton shield-exit count(E>100 MeV) / count(E>10 MeV)"},
      {"denominator_count",fPrimaryProtonAbove10},
      {"interval_95",shield_sim::Json(shield_sim::Json::Object{{"high",high},{"low",low}})},
      {"name","H100/10"},{"numerator_count",fPrimaryProtonAbove100},
      {"unit","1"},{"validation_status",validationStatus},
      {"value",finiteOrNull(fLastHardness)}});
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
    {"em_physics",c.emPhysics},
    {"ion_physics",c.ionPhysics},
    {"energy_ledger_scope",shield_sim::Json(shield_sim::Json::Object{
      {"absorber_deposition","Geant4 total, ionizing, and nonionizing deposited energy per primary"},
      {"authoritative_particle_breakdown","authoritative classification by transported particle species, primary/secondary identity, propagated root particle, and generation; these partitions each close independently to the absorber aggregate"},
      {"diagnostic_process_breakdown","creator-process and step-defining-process views are diagnostics only and are not a unique physical-process attribution"},
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
  if(c.damageResponse){
    result["damage_response"]=shield_sim::Json(shield_sim::Json::Object{
      {"charged_particle_interpolation","log_log_between_positive_nodes"},
      {"library_id",c.damageResponse->LibraryId()},
      {"library_version",c.damageResponse->LibraryVersion()},
      {"manifest_hash",c.damageResponse->ManifestHash()},
      {"manifest_path",c.damageResponse->ManifestPath()},
      {"n_eq_scope","ASTM E722-19 silicon neutron contribution only"},
      {"neutron_integration","ASTM groupwise on explicit edges"},
      {"out_of_domain_policy","not_assessable_no_extrapolation"},
      {"surface_fluence_normalization","provider-defined entry surface; scalar crossing fluence includes the inverse projected-direction cosine"},
      {"target_material_id",c.damageResponse->TargetMaterialId()}});
  }
  if(c.letResponse.enabled){
    result["let_response"]=c.resolvedJson.At("let_response","resolved case");
  }
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
    {"absorber_geometry_provider",absorber.provider},
    {"absorber_geometry_provider_version",absorber.providerVersion},
    {"build_type",shield_sim_build::BuildType},
    {"code_revision",shield_sim_build::GitRevision},
    {"compiler",shield_sim_build::Compiler},
    {"compiler_flags",shield_sim_build::CompilerFlags},
    {"geant4_data_libraries",shield_sim::Json(dataLibraries)},
    {"geant4_version",G4Version},{"geometry_hash",c.geometryHash},
    {"em_physics",c.emPhysics},{"ion_physics",c.ionPhysics},
    {"material_library_hash",c.materials->ContentHash()},
    {"max_step_mm",c.maxStepMm},{"physics_list",c.physicsList},
    {"production_cut_mm",c.productionCutMm},
    {"scoring_contract",c.letResponse.enabled
        ? "shieldSim.transport_raw_plus_step_fluence_LET_and_weighted_hardness.v1"
        : c.damageResponse
        ? "shieldSim.transport_raw_plus_silicon_damage.v1"
        : "shieldSim.transport_raw.v1"},
    {"source_tree_state",shield_sim_build::GitState}
  };
  // Stage-5 entry scores depend on the exact damage-function package even when
  // geometry and Geant4 settings are unchanged. Include its canonical manifest
  // hash in the cache key so a response cell can never be silently reused after
  // a NIEL table, interpolation convention, or normalization change.
  if(c.damageResponse)
    keyDocument["damage_response"]=shield_sim::Json(shield_sim::Json::Object{
      {"library_id",c.damageResponse->LibraryId()},
      {"library_version",c.damageResponse->LibraryVersion()},
      {"manifest_hash",c.damageResponse->ManifestHash()},
      {"target_material_id",c.damageResponse->TargetMaterialId()}});
  if(c.letResponse.enabled)
    keyDocument["let_response"]=c.resolvedJson.At("let_response","resolved case");
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
    {"validation_status",absorber.provider!="slab.v1"
        ? "stage10_geometry_providers_verification_only_incomplete_no_independent_nontrivial_shape_validation"
        : c.damageResponse
        ? "stage5_silicon_DDD_and_neutron_only_n_eq_implemented_and_reference_tested_pending_accelerator_device_validation_and_independent_review"
        : "implemented_stage3_transport_raw_not_physically_validated"}
  });
  const std::string path=fOpts.outputPrefix+"_response_cell.json";
  std::ofstream out(path);
  if(!out)
    G4Exception("RunAction::WriteResponseCellJson","Output",FatalException,
                ("Cannot write response cell JSON '"+path+"'").c_str());
  out<<cell.DumpCanonical()<<'\n';
  G4cout<<"Response cell written: "<<path<<G4endl;
}
