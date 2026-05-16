#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "MaterialCatalog.hh"
#include "SpecBins.hh"

#include <G4AnalysisManager.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4Run.hh>
#include <G4SystemOfUnits.hh>
#include <G4ios.hh>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>

RunAction::RunAction(const Options& opts,
                     DetectorConstruction* detector,
                     const std::vector<std::string>& scoringNames,
                     const std::vector<G4double>&    scoringThick)
  : fOpts(opts), fDetector(detector), fScoringNames(scoringNames), fScoringThick(scoringThick)
{
  G4int N=SpecBins::N;
  fEdep.assign(scoringNames.size(),0);
  fLastDose.assign(scoringNames.size(),0);
  fLastDoseRate.assign(scoringNames.size(),0);
  fInP.assign(N,0); fInA.assign(N,0);
  fOutP.assign(N,0); fOutA.assign(N,0); fOutN.assign(N,0);

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

void RunAction::AddEdep(std::size_t i,G4double edep){
  if(i<fEdep.size()) fEdep[i]+=edep;
}

void RunAction::SetSourceNormalization(G4double normNoAngular,
                                       G4double angularFactor,
                                       const std::string& sourceMode){
  fSourceNormNoAngular = normNoAngular;
  fSourceAngularFactor = angularFactor;
  fSourceNorm          = normNoAngular*angularFactor;
  fSourceMode          = sourceMode;
}

G4double RunAction::GetSourcePlaneAreaCM2() const {
  auto* wb = fDetector ? fDetector->GetWorldBox() : nullptr;
  if(!wb) return 0.0;
  return (2.0*wb->GetXHalfLength()/cm)*(2.0*wb->GetYHalfLength()/cm);
}

G4double RunAction::GetIncidentParticleRate() const {
  if(fSourceMode=="isotropic") return fSourceNorm*GetSourcePlaneAreaCM2();
  return fSourceNorm;
}

void RunAction::SetSweepMode(bool on,const std::string& mat,G4double tMM){
  fSweepMode=on; fCurrentMat=mat; fCurrentTmm=tMM;
}

void RunAction::AppendSweepPoint(const SweepPoint& sp){ fSweepData.push_back(sp); }
const std::vector<SweepPoint>& RunAction::GetSweepData() const { return fSweepData; }
const std::vector<G4double>&   RunAction::GetLastDose() const { return fLastDose; }
const std::vector<G4double>&   RunAction::GetLastDoseRate() const { return fLastDoseRate; }

void RunAction::RefreshLVPointers(){
  fShieldLV = fDetector ? fDetector->GetShieldLV() : nullptr;
  fScoringLVs.clear();
  if(fDetector) fScoringLVs = fDetector->GetScoringLVs();
}

G4LogicalVolume* RunAction::GetShieldLV() const{ return fShieldLV; }
const std::vector<G4LogicalVolume*>& RunAction::GetScoringLVs() const{ return fScoringLVs; }

void RunAction::BeginOfRunAction(const G4Run*) {
  RefreshLVPointers();

  auto* am=G4AnalysisManager::Instance();
  am->Reset();
  std::string csvName="shieldSimOutput";
  if(fSweepMode)
    csvName+="_"+SanitiseName(fCurrentMat)+"_"+FormatMM(fCurrentTmm)+"mm";
  csvName+=".csv";
  am->OpenFile(csvName);

  std::fill(fEdep.begin(),fEdep.end(),0);
  std::fill(fLastDoseRate.begin(),fLastDoseRate.end(),0);
  std::fill(fInP .begin(),fInP .end(),0);
  std::fill(fInA .begin(),fInA .end(),0);
  std::fill(fOutP.begin(),fOutP.end(),0);
  std::fill(fOutA.begin(),fOutA.end(),0);
  std::fill(fOutN.begin(),fOutN.end(),0);
}

void RunAction::EndOfRunAction(const G4Run* run) {
  G4int nEv=run->GetNumberOfEvent();
  if(nEv==0) return;

  // Dose per primary.  fEdep is in Geant4 energy units, mass is in Geant4 mass
  // units, therefore fEdep/mass is already in Geant4 dose units.
  for(std::size_t i=0;i<fScoringLVs.size();++i){
    auto* lv=fScoringLVs[i];
    G4double dose=0.0;
    if(lv){
      G4double mass = lv->GetMass(true,false);
      if(mass>0.0) dose = fEdep[i]/mass/nEv;
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

  auto* am=G4AnalysisManager::Instance();
  am->Write();
  am->CloseFile();

  WriteSpectraTecplot(nEv);
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

void RunAction::WriteSpectraTecplot(G4int nEv){
  const std::string fname="shieldSim_spectra.dat";
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
    if(fOpts.spectrumFile.empty())
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
