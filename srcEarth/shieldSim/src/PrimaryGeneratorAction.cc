#include "PrimaryGeneratorAction.hh"

#include "CaseDefinition.hh"
#include "DetectorConstruction.hh"
#include "GCRSpectrum.hh"
#include "RunAction.hh"
#include "SpectrumDataset.hh"

#include <G4Box.hh>
#include <G4Event.hh>
#include <G4Exception.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4IonTable.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4ios.hh>
#include <Randomize.hh>
#include <CLHEP/Random/RandGeneral.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <tuple>

PrimaryGeneratorAction::PrimaryGeneratorAction(const Options& opts,
                                               RunAction* ra,
                                               DetectorConstruction* det)
  : fOpts(opts), fRunAction(ra), fDetector(det)
{
  fGun = new G4ParticleGun(1);
  if(fOpts.caseDefinition && (fOpts.caseDefinition->responseIncident.enabled ||
                              !fOpts.caseDefinition->directIncidents.empty())){
    // A response cell is exactly one declared species, total kinetic energy,
    // and direction per history.  Its normalization is therefore one primary,
    // not a physical environmental rate.
    fSourceNormNoAngular=1.0;
    fSourceAngularFactor=1.0;
  } else if(fOpts.caseDefinition)
    ReadDataset(fOpts.caseDefinition->spectrum,
                fOpts.caseDefinition->spectrum.spectrumPath);
  else if(fOpts.spectrumFile.empty()) BuildGCR();
  else                                ReadFile(fOpts.spectrumFile);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fGun;
  delete fPRand;
  delete fARand;
}

G4double PrimaryGeneratorAction::GetSourceNormNoAngular() const { return fSourceNormNoAngular; }
G4double PrimaryGeneratorAction::GetSourceAngularFactor() const { return fSourceAngularFactor; }
G4double PrimaryGeneratorAction::GetSourceNorm() const { return fSourceNormNoAngular*fSourceAngularFactor; }
const std::string& PrimaryGeneratorAction::GetSourceMode() const { return fOpts.sourceMode; }

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* ev) {
  if(fOpts.caseDefinition && fOpts.caseDefinition->responseIncident.enabled){
    const auto& incident=fOpts.caseDefinition->responseIncident;
    fActiveIncident=&incident;
    ConfigureSourcePositionAndDirection();
    GenerateDeclaredIncident(incident,ev);
    return;
  }
  if(fOpts.caseDefinition && !fOpts.caseDefinition->directIncidents.empty()){
    const double draw=G4UniformRand();
    double cumulative=0.0;
    const shield_sim::WeightedResponseIncident* selected=&fOpts.caseDefinition->directIncidents.back();
    for(const auto& candidate:fOpts.caseDefinition->directIncidents){
      cumulative+=candidate.weight;
      if(draw<cumulative){ selected=&candidate; break; }
    }
    fActiveIncident=&selected->incident;
    ConfigureSourcePositionAndDirection();
    GenerateDeclaredIncident(selected->incident,ev);
    return;
  }

  // Query the current world dimensions every event because sweep mode can
  // rebuild the geometry between runs.
  fActiveIncident=nullptr;
  ConfigureSourcePositionAndDirection();

  // Species selection is proportional to the species-integrated source weights.
  // The angular factor is common to both species and cancels in this choice.
  G4double tot=fTotP+fTotA;
  bool useP=(tot>0) ? (G4UniformRand()<fTotP/tot) : true;

  G4double E; G4ParticleDefinition* pd;
  std::string species;
  if(useP){
    E=Sample(fPRand,fEP,fPLower,fPUpper);
    pd=G4ParticleTable::GetParticleTable()->FindParticle("proton");
    species="proton";
    RecordInput(true,E);
  } else {
    E=Sample(fARand,fEA,fALower,fAUpper);
    pd=G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    species="alpha";
    RecordInput(false,E);
  }

  // Record the source diagnostic after species and energy are selected but
  // before Geant4 tracking begins.  This diagnostic is optional and is used by
  // the Layer-2 tests to verify the source algorithm independently from the
  // transport physics.
  RecordSourceDiagnostic(species,E);

  // Geant4 expects total kinetic energy for the particle definition.  For alpha
  // particles this is total alpha energy, not energy per nucleon.
  fGun->SetParticleDefinition(pd);
  fGun->SetParticleEnergy(E*MeV);
  fGun->GeneratePrimaryVertex(ev);
}

void PrimaryGeneratorAction::GenerateDeclaredIncident(
    const shield_sim::ResponseIncidentDefinition& incident,G4Event* event){
  G4ParticleDefinition* particle=nullptr;
  if(incident.particleName=="proton" || incident.particleName=="alpha" ||
     incident.particleName=="neutron" || incident.particleName=="gamma")
    particle=G4ParticleTable::GetParticleTable()->FindParticle(incident.particleName);
  else
    particle=G4IonTable::GetIonTable()->GetIon(incident.atomicNumber,incident.massNumber,0.0);
  if(!particle)
    G4Exception("PrimaryGeneratorAction","ResponseParticle",FatalException,
                "Cannot construct the declared Stage-3/4/5/6 incident particle");
  const G4double energy=incident.totalKineticEnergyMeV;
  RecordSourceDiagnostic(particle->GetParticleName(),energy);
  if(incident.particleName=="proton") RecordInput(true,energy);
  else if(incident.particleName=="alpha") RecordInput(false,energy);
  fGun->SetParticleDefinition(particle);
  fGun->SetParticleEnergy(energy*MeV);
  fGun->GeneratePrimaryVertex(event);
}

G4double PrimaryGeneratorAction::AngularFactorForMode() const {
  // Isotropic intensity crossing a plane has integral J cos(theta)dOmega = pi J
  // over the inward hemisphere.  Beam mode is already a pencil-beam rate.
  return (fOpts.sourceMode=="isotropic") ? pi : 1.0;
}

void PrimaryGeneratorAction::BuildGCR(){
  // The historical implementation clamped requests into the approximation's
  // domain.  That changed a scientific input silently.  Stage 1 requires an
  // explicit failure instead, so the archived request always matches the run.
  if(fOpts.eMinProton<1.0 || fOpts.eMaxProton>1.0e5 ||
     fOpts.eMinAlpha <1.0 || fOpts.eMaxAlpha >1.0e5)
    G4Exception("PrimaryGeneratorAction","GCRDomain",FatalException,
                "Built-in legacy GCR approximation supports only 1--100000 MeV total kinetic energy; request was not clamped");
  G4double epLo=fOpts.eMinProton;
  G4double epHi=fOpts.eMaxProton;
  G4double eaLo=fOpts.eMinAlpha;
  G4double eaHi=fOpts.eMaxAlpha;

  std::vector<G4double> pF,aF;
  fTotP = BuildLogWeightedTable(epLo,epHi,fEP,pF,GCR::ProtonWeight);
  fTotA = BuildLogWeightedTable(eaLo,eaHi,fEA,aF,GCR::AlphaWeight);

  // RandGeneral requires a nonzero distribution even for an unused species.
  std::vector<G4double> pRandW=pF, aRandW=aF;
  if(fTotP<=0.0 && !pRandW.empty()) pRandW[0]=1.0;
  if(fTotA<=0.0 && !aRandW.empty()) aRandW[0]=1.0;
  fPRand=new CLHEP::RandGeneral(pRandW.data(),pRandW.size(),1);
  fARand=new CLHEP::RandGeneral(aRandW.data(),aRandW.size(),1);

  fSourceNormNoAngular = fTotP+fTotA;
  fSourceAngularFactor = AngularFactorForMode();

  G4cout<<"[GCR] Badhwar-O'Neill phi=550MV  "
        <<"proton ["<<epLo<<","<<epHi<<"] MeV  "
        <<"alpha ["<<eaLo<<","<<eaHi<<"] MeV  "
        <<"source-mode="<<fOpts.sourceMode<<G4endl;
}

std::vector<G4double> PrimaryGeneratorAction::EstimateTableBinWidths(
    const std::vector<G4double>& E)
{
  const std::size_t n=E.size();
  std::vector<G4double> dE(n,1.0);
  if(n==0) return dE;
  if(n==1){
    dE[0]=std::max(1.0,1.0e-6*E[0]);
    return dE;
  }

  std::vector<G4double> edge(n+1,0.0);
  for(std::size_t i=1;i<n;++i){
    if(E[i]>0.0 && E[i-1]>0.0 && E[i]!=E[i-1])
      edge[i]=std::sqrt(E[i-1]*E[i]);
    else
      edge[i]=0.5*(E[i-1]+E[i]);
  }
  edge[0] = (edge[1]>0.0) ? E[0]*E[0]/edge[1]
                          : std::max(0.0,E[0]-0.5);
  edge[n] = (edge[n-1]>0.0) ? E[n-1]*E[n-1]/edge[n-1]
                            : E[n-1]+0.5;

  for(std::size_t i=0;i<n;++i){
    dE[i]=edge[i+1]-edge[i];
    if(dE[i]<=0.0) dE[i]=std::max(1.0,1.0e-6*E[i]);
  }
  return dE;
}

void PrimaryGeneratorAction::ReadFile(const std::string& fn){
  // Numeric legacy files are accepted only through the strict adapter and a
  // mandatory metadata sidecar.  The adapter rejects malformed rows, negative
  // or non-finite values, duplicate/out-of-order energies, and undeclared
  // single-bin widths instead of skipping, sorting, or clamping them.
  shield_sim::SpectrumDataset dataset;
  try {
    dataset=shield_sim::SpectrumDataset::LoadLegacyThreeColumn(
        fn,fOpts.legacySpectrumMetadataFile);
  } catch(const std::exception& error) {
    G4Exception("PrimaryGeneratorAction","SpectrumFile",FatalException,error.what());
  }
  ReadDataset(dataset,fn);
}

void PrimaryGeneratorAction::ReadDataset(const shield_sim::SpectrumDataset& dataset,
                                         const std::string& displayPath){
  std::vector<G4double> pF,aF;
  fEP.reserve(dataset.energyMeV.size()); fEA.reserve(dataset.energyMeV.size());
  pF.reserve(dataset.energyMeV.size());  aF.reserve(dataset.energyMeV.size());

  for(std::size_t i=0;i<dataset.energyMeV.size();++i){
    if(dataset.transportBinLowerMeV.size()==dataset.energyMeV.size() &&
       dataset.transportBinUpperMeV.size()==dataset.energyMeV.size()){
      // Published bin values are differential intensities averaged over the
      // stated interval. Clip a boundary bin to the requested energy range,
      // weight it by its retained width, and sample uniformly inside it. This
      // is the unique piecewise-constant interpretation and never invents an
      // effective bin-center energy.
      const G4double pLower=std::max(dataset.transportBinLowerMeV[i],fOpts.eMinProton);
      const G4double pUpper=std::min(dataset.transportBinUpperMeV[i],fOpts.eMaxProton);
      if(pUpper>pLower && dataset.protonValue[i]>0.0){
        fEP.push_back(0.5*(pLower+pUpper));
        fPLower.push_back(pLower); fPUpper.push_back(pUpper);
        pF.push_back(dataset.protonValue[i]*(pUpper-pLower));
      }
      const G4double aLower=std::max(dataset.transportBinLowerMeV[i],fOpts.eMinAlpha);
      const G4double aUpper=std::min(dataset.transportBinUpperMeV[i],fOpts.eMaxAlpha);
      if(aUpper>aLower && dataset.alphaValue[i]>0.0){
        fEA.push_back(0.5*(aLower+aUpper));
        fALower.push_back(aLower); fAUpper.push_back(aUpper);
        aF.push_back(dataset.alphaValue[i]*(aUpper-aLower));
      }
    } else {
      const G4double E=dataset.energyMeV[i];
      G4double p=dataset.protonValue[i];
      G4double a=dataset.alphaValue[i];
      if(E<fOpts.eMinProton||E>fOpts.eMaxProton) p=0;
      if(E<fOpts.eMinAlpha ||E>fOpts.eMaxAlpha)  a=0;
      fEP.push_back(E); fEA.push_back(E);
      pF.push_back(p*dataset.binWidthMeV[i]);
      aF.push_back(a*dataset.binWidthMeV[i]);
    }
  }

  fTotP=0; for(auto v:pF) fTotP+=v;
  fTotA=0; for(auto v:aF) fTotA+=v;

  // RandGeneral requires a nonempty, nonzero table even for a species whose
  // physical integral is zero. A one-bin dummy is safe because species choice
  // uses fTotP/fTotA and therefore never selects a zero-integral species.
  std::vector<G4double> pRandW=pF, aRandW=aF;
  if(pRandW.empty()){ pRandW.push_back(1.0); fEP.push_back(fOpts.eMinProton); }
  else if(fTotP<=0.0) pRandW[0]=1.0;
  if(aRandW.empty()){ aRandW.push_back(1.0); fEA.push_back(fOpts.eMinAlpha); }
  else if(fTotA<=0.0) aRandW[0]=1.0;
  fPRand=new CLHEP::RandGeneral(pRandW.data(),pRandW.size(),1);
  fARand=new CLHEP::RandGeneral(aRandW.data(),aRandW.size(),1);

  if(fTotP<=0 && fTotA<=0)
    G4Exception("PrimaryGeneratorAction","SpectrumFile",FatalException,"Spectrum weights are zero for both species after filtering.");

  fSourceNormNoAngular = fTotP+fTotA;
  fSourceAngularFactor = AngularFactorForMode();

  G4cout<<"[Spectrum file] "<<displayPath
        <<"  dataset="<<dataset.id
        <<"  transport-component="
        <<(dataset.transportComponentId.empty()?"legacy_shared_grid":dataset.transportComponentId)
        <<"  hash="<<dataset.contentHash
        <<"  source-mode="<<fOpts.sourceMode
        <<"  integral(no angular factor)="<<fSourceNormNoAngular
        <<"  angular factor="<<fSourceAngularFactor<<G4endl;
}

G4double PrimaryGeneratorAction::Sample(CLHEP::RandGeneral* rng,
                                        const std::vector<G4double>& E,
                                        const std::vector<G4double>& lower,
                                        const std::vector<G4double>& upper){
  if(!rng||E.empty()) return 0;
  const G4double selected=rng->shoot()*E.size();
  G4int i=static_cast<G4int>(selected);
  i=std::max(0,std::min(i,(G4int)E.size()-1));
  if(lower.size()==E.size() && upper.size()==E.size())
    return lower[i]+G4UniformRand()*(upper[i]-lower[i]);
  return E[i];
}

void PrimaryGeneratorAction::ConfigureSourcePositionAndDirection(){
  auto* worldBox = GetCurrentWorldBox();
  if(!worldBox){
    fGun->SetParticlePosition(G4ThreeVector(0,0,0));
    fGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    return;
  }

  const G4double eps = 1.0*um;
  const G4ThreeVector origin=fDetector ? fDetector->GetSourcePlaneOrigin()
      : G4ThreeVector(0,0,-worldBox->GetZHalfLength()+eps);
  const G4ThreeVector axisX=fDetector ? fDetector->GetSourcePlaneAxisX()
                                      : G4ThreeVector(1,0,0);
  const G4ThreeVector axisY=fDetector ? fDetector->GetSourcePlaneAxisY()
                                      : G4ThreeVector(0,1,0);
  const G4ThreeVector normal=fDetector ? fDetector->GetSourcePlaneNormal()
                                       : G4ThreeVector(0,0,1);

  if(fOpts.sourceMode=="directional_plane"){
    const G4double halfX=fDetector?fDetector->GetSourcePlaneHalfX():worldBox->GetXHalfLength();
    const G4double halfY=fDetector?fDetector->GetSourcePlaneHalfY():worldBox->GetYHalfLength();
    const G4double x=(2.0*G4UniformRand()-1.0)*std::max(0.0,halfX-eps);
    const G4double y=(2.0*G4UniformRand()-1.0)*std::max(0.0,halfY-eps);
    if(!fActiveIncident)
      G4Exception("PrimaryGeneratorAction","DirectionalIncident",FatalException,
                  "directional_plane source has no selected incident node");
    const auto& incident=*fActiveIncident;
    const G4double transverse=std::sqrt(std::max(0.0,1.0-incident.mu*incident.mu));
    fLastSourcePosition=origin+x*axisX+y*axisY;
    fLastSourceDirection=(transverse*std::cos(incident.phiRad))*axisX+
                         (transverse*std::sin(incident.phiRad))*axisY+
                         incident.mu*normal;
    fLastSourceDirection=fLastSourceDirection.unit();
    fGun->SetParticlePosition(fLastSourcePosition);
    fGun->SetParticleMomentumDirection(fLastSourceDirection);
    return;
  }

  if(fOpts.sourceMode=="isotropic"){
    const G4double halfX = fDetector ? fDetector->GetSourcePlaneHalfX()
                                     : worldBox->GetXHalfLength();
    const G4double halfY = fDetector ? fDetector->GetSourcePlaneHalfY()
                                     : worldBox->GetYHalfLength();
    const G4double x = (2.0*G4UniformRand()-1.0)*std::max(0.0,halfX-eps);
    const G4double y = (2.0*G4UniformRand()-1.0)*std::max(0.0,halfY-eps);

    // Isotropic intensity crossing a plane has p(mu)=2*mu for mu in [0,1].
    const G4double mu  = std::sqrt(G4UniformRand());
    const G4double phi = twopi*G4UniformRand();
    const G4double sinTheta = std::sqrt(std::max(0.0,1.0-mu*mu));
    fLastSourcePosition = origin+x*axisX+y*axisY;
    fLastSourceDirection = ((sinTheta*std::cos(phi))*axisX+
                            (sinTheta*std::sin(phi))*axisY+mu*normal).unit();
    fGun->SetParticlePosition(fLastSourcePosition);
    fGun->SetParticleMomentumDirection(fLastSourceDirection);
  } else {
    fLastSourcePosition = origin;
    fLastSourceDirection = normal.unit();
    fGun->SetParticlePosition(fLastSourcePosition);
    fGun->SetParticleMomentumDirection(fLastSourceDirection);
  }
}

void PrimaryGeneratorAction::RecordInput(bool isProton,G4double E){
  if(isProton) fRunAction->AddInP(E);
  else         fRunAction->AddInA(E);
}

void PrimaryGeneratorAction::RecordSourceDiagnostic(const std::string& species,G4double E){
  if(!fRunAction) return;
  fRunAction->RecordSourceSample(species,E,
                                 fLastSourcePosition.x()/mm,
                                 fLastSourcePosition.y()/mm,
                                 fLastSourcePosition.z()/mm,
                                 fLastSourceDirection.x(),
                                 fLastSourceDirection.y(),
                                 fLastSourceDirection.z());
}

G4Box* PrimaryGeneratorAction::GetCurrentWorldBox() const {
  return fDetector ? fDetector->GetWorldBox() : nullptr;
}
