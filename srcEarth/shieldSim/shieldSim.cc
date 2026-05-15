// shieldSim.cc
//
// Geant4 GCR shielding simulation with two operating modes:
//
//  Single-run mode (default)
//    One shield thickness, one BeamOn call.  Writes a Tecplot spectrum file
//    (shieldSim_spectra.dat) and a CSV histogram file (shieldSimOutput.csv).
//
//  Sweep mode  (--sweep)
//    Iterates over a range of shield thicknesses using G4RunManager::
//    ReinitializeGeometry() between runs.  For each thickness, dose in every
//    scoring volume is recorded.  Results are written to two Tecplot files:
//      shieldSim_dose_sweep.dat   — dose [Gy/primary] vs thickness [mm]
//      shieldSim_spectra.dat      — transmitted spectra, one ZONE per thickness
//
//  Primary spectrum: Badhwar-O'Neill GCR force-field approximation at 1 AU,
//  phi = 550 MV (solar minimum).  A tabulated file may be supplied instead.
//  Energy sampling range is configurable per species via CLI.

#include <G4RunManager.hh>
#include <G4Run.hh>
#include <G4NistManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4UserRunAction.hh>
#include <G4UserEventAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4VModularPhysicsList.hh>
#include <G4PhysListFactory.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PVPlacement.hh>
#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4Event.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4ThreeVector.hh>
#include <G4AnalysisManager.hh>
#include <Randomize.hh>
#include <CLHEP/Random/RandGeneral.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <map>

// =============================================================================
// GCR spectrum — Badhwar-O'Neill force-field approximation
// =============================================================================
namespace GCR {
// Proton  Z=1, m_p=938.272 MeV, A_LIS=1.9e4, gamma=2.77
inline G4double ProtonWeight(G4double E_MeV) {
  const G4double mp=938.272, phi=550.0;
  G4double EL=E_MeV+phi;
  G4double p2=E_MeV*(E_MeV+2*mp), pL2=EL*(EL+2*mp);
  if(p2<=0||pL2<=0) return 0;
  G4double R=std::sqrt(pL2)/1000., b=std::sqrt(pL2)/(EL+mp);
  return 1.9e4*b*b*std::pow(R,-2.77)*p2/pL2;
}
// Alpha  Z=2, m_a=3727.38 MeV, A_LIS=6.2e3, gamma=2.66
inline G4double AlphaWeight(G4double E_MeV) {
  const G4double ma=3727.38, Z=2, phi=550.0;
  G4double EL=E_MeV+Z*phi;
  G4double p2=E_MeV*(E_MeV+2*ma), pL2=EL*(EL+2*ma);
  if(p2<=0||pL2<=0) return 0;
  G4double R=std::sqrt(pL2)/(Z*1000.), b=std::sqrt(pL2)/(EL+ma);
  return 6.2e3*b*b*std::pow(R,-2.66)*p2/pL2;
}
} // namespace GCR

// =============================================================================
// Shared spectral binning  (120 log-spaced bins, 1 MeV – 100 GeV)
// =============================================================================
namespace SpecBins {
  static const G4int    N=120;
  static const G4double Emin=1.0, Emax=1.0e5;
  static const G4double logR=std::log(Emax/Emin);
  inline G4double Center(G4int i){ return Emin*std::exp((i+.5)/N*logR); }
  inline G4double Edge  (G4int i){ return Emin*std::exp(G4double(i)/N*logR); }
  inline G4double Width (G4int i){ return Edge(i+1)-Edge(i); }
  inline G4int    Bin   (G4double E){
    if(E<=Emin) return 0;
    if(E>=Emax) return N-1;
    return static_cast<G4int>(std::log(E/Emin)/logR*N);
  }
} // namespace SpecBins

// =============================================================================
// Forward declarations
// =============================================================================
class RunAction;
static void PrintHelp();

// =============================================================================
// Options
// =============================================================================
struct Options {
  // ---- spectrum source ----
  std::string spectrumFile = "";      // empty → built-in GCR

  // ---- energy limits for sampling (MeV) ----
  G4double eMinProton =  10.0;
  G4double eMaxProton = 1.0e5;
  G4double eMinAlpha  =  10.0;
  G4double eMaxAlpha  = 1.0e5;

  // ---- single-run geometry ----
  std::string shieldMaterial  = "G4_Al";
  G4double    shieldThickness = 2.0*mm;

  // ---- scoring slabs ----
  std::vector<std::pair<std::string,G4double>> scoringMaterials;

  // ---- statistics ----
  G4int nEvents = 10000;

  // ---- dose-vs-thickness sweep ----
  bool        doSweep       = false;
  std::string sweepMaterial = "";    // default: same as shieldMaterial
  G4double    sweepTmin     = 0.5;   // mm
  G4double    sweepTmax     = 50.0;  // mm
  G4int       sweepN        = 10;
  bool        sweepLog      = false; // log-spaced thicknesses

  bool showHelp = false;
};

// =============================================================================
// DetectorConstruction
// =============================================================================
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  explicit DetectorConstruction(const Options& opts) : fOpts(opts) {}
  ~DetectorConstruction() override {}

  // Called by main() before ReinitializeGeometry() in sweep mode
  void SetShieldThickness(G4double t){ fOpts.shieldThickness=t; }
  const Options& GetOptions() const  { return fOpts; }
  G4LogicalVolume* GetShieldLV() const { return fShieldLV; }
  const std::vector<G4LogicalVolume*>& GetScoringLVs() const { return fScoringLVs; }
  G4Box* GetWorldBox() const { return fWorldBox; }

  G4VPhysicalVolume* Construct() override {
    fWorldBox = nullptr;
    fShieldLV = nullptr;
    fScoringLVs.clear();

    auto* nist=G4NistManager::Instance();
    G4double wXY=5.*cm;
    G4double scSum=0;
    for(const auto& s:fOpts.scoringMaterials) scSum+=s.second;
    G4double wZ=2.*mm+fOpts.shieldThickness+scSum+2.*mm;

    auto* vac  =nist->FindOrBuildMaterial("G4_Galactic");
    auto* wSol =new G4Box("World",wXY/2,wXY/2,wZ/2);
    auto* wLog =new G4LogicalVolume(wSol,vac,"World");
    fWorldBox = wSol;
    auto* wPhys=new G4PVPlacement(nullptr,G4ThreeVector(),wLog,"World",nullptr,false,0);

    auto* shMat=nist->FindOrBuildMaterial(fOpts.shieldMaterial);
    if(!shMat) G4Exception("Detector","ShieldMat",FatalException,
        ("Shield material "+fOpts.shieldMaterial+" not found").c_str());
    G4double shHz=fOpts.shieldThickness/2;
    auto* shS=new G4Box("Shield",wXY/2,wXY/2,shHz);
    auto* shL=new G4LogicalVolume(shS,shMat,"Shield");
    fShieldLV = shL;
    G4double zSh=1.*mm+shHz;
    new G4PVPlacement(nullptr,G4ThreeVector(0,0,zSh-wZ/2),shL,"Shield",wLog,false,0,true);

    G4double curZ=1.*mm+fOpts.shieldThickness;
    for(std::size_t i=0;i<fOpts.scoringMaterials.size();++i){
      const auto& sm=fOpts.scoringMaterials[i];
      auto* mat=nist->FindOrBuildMaterial(sm.first);
      if(!mat) G4Exception("Detector","ScoreMat",FatalException,
          ("Scoring material "+sm.first+" not found").c_str());
      G4double hz=sm.second/2;
      std::string idx=std::to_string(i);
      std::string solidName="ScoringSolid_"+idx+"_"+sm.first;
      std::string lvName   ="ScoringLV_"+idx+"_"+sm.first;
      std::string pvName   ="ScoringPV_"+idx+"_"+sm.first;
      auto* ds=new G4Box(solidName.c_str(),wXY/2,wXY/2,hz);
      auto* dl=new G4LogicalVolume(ds,mat,lvName.c_str());
      fScoringLVs.push_back(dl);
      curZ+=hz;
      new G4PVPlacement(nullptr,G4ThreeVector(0,0,curZ-wZ/2),
                        dl,pvName.c_str(),wLog,false,(G4int)i,true);
      curZ+=hz;
    }
    return wPhys;
  }

private:
  Options fOpts;
  G4Box* fWorldBox=nullptr;
  G4LogicalVolume* fShieldLV=nullptr;
  std::vector<G4LogicalVolume*> fScoringLVs;
};


// =============================================================================
// PrimaryGeneratorAction
// =============================================================================
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction(const Options& opts, RunAction* ra, DetectorConstruction* det)
  : fOpts(opts), fRunAction(ra), fDetector(det)
  {
    fGun = new G4ParticleGun(1);
    fGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    if(fOpts.spectrumFile.empty()) BuildGCR();
    else                           ReadFile(fOpts.spectrumFile);
  }
  ~PrimaryGeneratorAction() override {
    delete fGun; delete fPRand; delete fARand;
  }

  void GeneratePrimaries(G4Event* ev) override {
    // Place gun just inside the upstream world face. The world box can change
    // during sweep-mode geometry reinitialization, so always query the current one.
    auto* worldBox = GetCurrentWorldBox();
    if(worldBox){
      const G4double eps = 1.0*um;
      fGun->SetParticlePosition(G4ThreeVector(0,0,-worldBox->GetZHalfLength()+eps));
    }

    G4double tot=fTotP+fTotA;
    bool useP=(tot>0) ? (G4UniformRand()<fTotP/tot) : true;

    G4double E; G4ParticleDefinition* pd;
    if(useP){
      E=Sample(fPRand,fEP);
      pd=G4ParticleTable::GetParticleTable()->FindParticle("proton");
      RecordInput(true,E);
    } else {
      E=Sample(fARand,fEA);
      pd=G4ParticleTable::GetParticleTable()->FindParticle("alpha");
      RecordInput(false,E);
    }
    fGun->SetParticleDefinition(pd);
    fGun->SetParticleEnergy(E*MeV);
    fGun->GeneratePrimaryVertex(ev);
  }

private:
  void BuildGCR(){
    const G4int Ns=500;
    // Clamp user limits to [1 MeV, 1e5 MeV]
    G4double epLo=std::max(1.0,fOpts.eMinProton);
    G4double epHi=std::min(1.0e5,fOpts.eMaxProton);
    G4double eaLo=std::max(1.0,fOpts.eMinAlpha);
    G4double eaHi=std::min(1.0e5,fOpts.eMaxAlpha);

    fEP.resize(Ns); fEA.resize(Ns);
    std::vector<G4double> pF(Ns,0), aF(Ns,0);
    for(G4int i=0;i<Ns;++i){
      G4double fp=G4double(i)/(Ns-1);
      G4double Ep=epLo*std::pow(epHi/epLo,fp);
      G4double Ea=eaLo*std::pow(eaHi/eaLo,fp);
      fEP[i]=Ep; fEA[i]=Ea;
      pF[i]=GCR::ProtonWeight(Ep);
      aF[i]=GCR::AlphaWeight(Ea);
    }
    fPRand=new CLHEP::RandGeneral(pF.data(),Ns,1);
    fARand=new CLHEP::RandGeneral(aF.data(),Ns,1);
    fTotP=0; for(auto v:pF) fTotP+=v;
    fTotA=0; for(auto v:aF) fTotA+=v;
    G4cout<<"[GCR] Badhwar-O'Neill phi=550MV  "
          <<"proton ["<<epLo<<","<<epHi<<"] MeV  "
          <<"alpha ["<<eaLo<<","<<eaHi<<"] MeV"<<G4endl;
  }

  void ReadFile(const std::string& fn){
    std::ifstream f(fn);
    if(!f) G4Exception("PrimaryGeneratorAction","SpectrumFile",FatalException,("Cannot open spectrum '"+fn+"'").c_str());
    std::vector<G4double> pF,aF;
    std::string line;
    while(std::getline(f,line)){
      if(line.empty()||line[0]=='#') continue;
      std::istringstream ss(line); G4double E,p,a;
      if(!(ss>>E>>p>>a)) continue;
      if(E<fOpts.eMinProton||E>fOpts.eMaxProton) p=0;
      if(E<fOpts.eMinAlpha ||E>fOpts.eMaxAlpha)  a=0;
      fEP.push_back(E); fEA.push_back(E);
      pF.push_back(p);  aF.push_back(a);
    }
    if(fEP.empty()) G4Exception("PrimaryGeneratorAction","SpectrumFile",FatalException,"Spectrum file is empty or invalid.");
    fPRand=new CLHEP::RandGeneral(pF.data(),pF.size(),1);
    fARand=new CLHEP::RandGeneral(aF.data(),aF.size(),1);
    fTotP=0; for(auto v:pF) fTotP+=v;
    fTotA=0; for(auto v:aF) fTotA+=v;
    if(fTotP<=0 && fTotA<=0)
      G4Exception("PrimaryGeneratorAction","SpectrumFile",FatalException,"Spectrum weights are zero for both species after filtering.");
  }

  G4double Sample(CLHEP::RandGeneral* rng,const std::vector<G4double>& E){
    if(!rng||E.empty()) return 0;
    G4int i=static_cast<G4int>(rng->shoot()*E.size());
    i=std::max(0,std::min(i,(G4int)E.size()-1));
    return E[i];
  }

  inline void RecordInput(bool p,G4double E);   // defined after RunAction

  G4Box* GetCurrentWorldBox() const;

  Options fOpts;
  RunAction* fRunAction;
  DetectorConstruction* fDetector;
  G4ParticleGun* fGun;
  std::vector<G4double> fEP,fEA;
  CLHEP::RandGeneral* fPRand=nullptr;
  CLHEP::RandGeneral* fARand=nullptr;
  G4double fTotP=0,fTotA=0;
};

// =============================================================================
// RunAction
// =============================================================================
// Sweep result record: (thickness_mm, areal_density_g_cm2, doses_Gy_per_primary)
struct SweepPoint {
  G4double thickMM;
  G4double arealDensity;  // g/cm²
  std::vector<G4double> dose_Gy; // one entry per scoring volume
};

class RunAction : public G4UserRunAction {
public:
  RunAction(const Options& opts,
            DetectorConstruction* detector,
            const std::vector<std::string>& scoringNames,
            const std::vector<G4double>&    scoringThick)
  : fOpts(opts), fDetector(detector), fScoringNames(scoringNames), fScoringThick(scoringThick)
  {
    G4int N=SpecBins::N;
    fEdep.assign(scoringNames.size(),0);
    fLastDose.assign(scoringNames.size(),0);
    fInP.assign(N,0); fInA.assign(N,0);
    fOutP.assign(N,0); fOutA.assign(N,0); fOutN.assign(N,0);

    auto* am=G4AnalysisManager::Instance();
    am->SetVerboseLevel(0);
    am->SetFirstHistoId(0);
    am->SetDefaultFileType("csv");
    G4int nb=500; G4double eH=1e5;
    am->CreateH1("ProtonExitE",  "Proton KE at exit",  nb,0,eH,"MeV");
    am->CreateH1("AlphaExitE",   "Alpha KE at exit",   nb,0,eH,"MeV");
    am->CreateH1("NeutronExitE", "Neutron KE at exit", nb,0,eH,"MeV");
  }
  ~RunAction() override {}

  // ---- accumulators ----
  void AddInP  (G4double E){ fInP [SpecBins::Bin(E)]+=1; }
  void AddInA  (G4double E){ fInA [SpecBins::Bin(E)]+=1; }
  void AddOutP (G4double E){ fOutP[SpecBins::Bin(E)]+=1; }
  void AddOutA (G4double E){ fOutA[SpecBins::Bin(E)]+=1; }
  void AddOutN (G4double E){ fOutN[SpecBins::Bin(E)]+=1; }
  void AddEdep (std::size_t i,G4double edep){
    if(i<fEdep.size()) fEdep[i]+=edep/MeV;
  }

  // ---- sweep mode control ----
  void SetSweepMode(bool on,const std::string& mat,G4double tMM){
    fSweepMode=on; fCurrentMat=mat; fCurrentTmm=tMM;
  }
  void AppendSweepPoint(const SweepPoint& sp){ fSweepData.push_back(sp); }
  const std::vector<SweepPoint>& GetSweepData() const { return fSweepData; }
  const std::vector<G4double>&   GetLastDose()  const { return fLastDose; }

  // ---- LV pointer refresh (call from BeginOfRunAction and from main after reinit) ----
  void RefreshLVPointers(){
    fShieldLV = fDetector ? fDetector->GetShieldLV() : nullptr;
    fScoringLVs.clear();
    if(fDetector) fScoringLVs = fDetector->GetScoringLVs();
  }

  G4LogicalVolume*                     GetShieldLV()   const{ return fShieldLV; }
  const std::vector<G4LogicalVolume*>& GetScoringLVs() const{ return fScoringLVs; }

  // ---- G4 hooks ----
  void BeginOfRunAction(const G4Run*) override {
    RefreshLVPointers();

    // Open CSV output (per-run filename in sweep mode)
    auto* am=G4AnalysisManager::Instance();
    am->Reset();
    std::string csvName="shieldSimOutput";
    if(fSweepMode)
      csvName+="_"+SanitiseName(fCurrentMat)
               +"_"+FormatMM(fCurrentTmm)+"mm";
    csvName+=".csv";
    am->OpenFile(csvName);

    // Reset accumulators
    std::fill(fEdep.begin(),fEdep.end(),0);
    std::fill(fInP .begin(),fInP .end(),0);
    std::fill(fInA .begin(),fInA .end(),0);
    std::fill(fOutP.begin(),fOutP.end(),0);
    std::fill(fOutA.begin(),fOutA.end(),0);
    std::fill(fOutN.begin(),fOutN.end(),0);
  }

  void EndOfRunAction(const G4Run* run) override {
    G4int nEv=run->GetNumberOfEvent();
    if(nEv==0) return;

    // Compute dose per primary for each scoring volume.
    // Geant4 stores deposited energy and mass in its own consistent internal
    // units, so the safest approach is to let Geant4 handle the unit
    // conversion implicitly here and only convert to Gy when reporting.
    for(std::size_t i=0;i<fScoringLVs.size();++i){
      auto* lv=fScoringLVs[i];
      G4double dose=0.0;
      if(lv){
        // Recompute the mass from the current geometry after any sweep-driven
        // reinitialization.  The result is returned in Geant4 internal mass
        // units; fEdep is in Geant4 internal energy units, therefore fEdep/mass
        // is already in Geant4 dose units.
        G4double mass = lv->GetMass(true,false);
        if(mass>0.0) dose = fEdep[i]/mass/nEv;
      }
      fLastDose[i]=dose;
    }

    // Console report
    G4cout<<"\nDose per primary ["
          <<(fSweepMode ? fCurrentMat+" "+FormatMM(fCurrentTmm)+"mm" : "single run")
          <<"]:"<<G4endl;
    for(std::size_t i=0;i<fScoringNames.size();++i)
      G4cout<<"  "<<fScoringNames[i]<<": "
            <<std::setprecision(5)<<std::scientific
            <<fLastDose[i]/gray<<" Gy/primary"<<G4endl;

    // CSV histograms
    auto* am=G4AnalysisManager::Instance();
    am->Write(); am->CloseFile();

    // Tecplot spectrum file — append a new ZONE in sweep mode
    WriteSpectraTecplot(nEv);
  }

private:
  // ---- helpers ----
  static std::string SanitiseName(const std::string& s){
    std::string r=s;
    for(auto& c:r) if(c=='/'||c=='\\'||c==' ') c='_';
    return r;
  }
  static std::string FormatMM(G4double t){
    std::ostringstream ss; ss<<std::fixed<<std::setprecision(2)<<t;
    return ss.str();
  }

  void WriteSpectraTecplot(G4int nEv){
    std::string fname="shieldSim_spectra.dat";
    // In sweep mode: open for append after first point, create for first
    std::ios::openmode mode= (fSweepMode && !fFirstSpectraWrite)
                             ? (std::ios::out|std::ios::app)
                             :  std::ios::out;
    std::ofstream out(fname,mode);
    if(!out){ G4cerr<<"Cannot write "<<fname<<G4endl; return; }

    if(!fSweepMode || fFirstSpectraWrite){
      // Write file header once
      out<<"TITLE = \"GCR Shielding Simulation Spectra – Geant4/FTFP_BERT\"\n";
      out<<"# Differential spectra: counts / (MeV * event)\n";
      out<<"# Shield: "<<(fSweepMode?fCurrentMat:fOpts.shieldMaterial)<<"\n";
      out<<"# phi = 550 MV (Badhwar-O'Neill, solar minimum)\n";
      out<<"VARIABLES = \"Energy [MeV]\""
         <<" \"Input_Proton [1/(MeV evt)]\""
         <<" \"Input_Alpha [1/(MeV evt)]\""
         <<" \"Output_Proton [1/(MeV evt)]\""
         <<" \"Output_Alpha [1/(MeV evt)]\""
         <<" \"Output_Neutron [1/(MeV evt)]\"\n";
      if(fSweepMode) fFirstSpectraWrite=false;
    }

    // Zone title
    std::ostringstream zt;
    if(fSweepMode)
      zt<<fCurrentMat<<"_"<<FormatMM(fCurrentTmm)<<"mm";
    else
      zt<<fOpts.shieldMaterial<<"_"
        <<FormatMM(fOpts.shieldThickness/mm)<<"mm";
    zt<<" | GCR phi=550MV | "<<nEv<<" events";

    out<<"ZONE T=\""<<zt.str()<<"\","
       <<" I="<<SpecBins::N<<", DATAPACKING=POINT\n";
    out<<std::scientific<<std::setprecision(5);
    for(G4int i=0;i<SpecBins::N;++i){
      G4double dE=SpecBins::Width(i);
      G4double nr=dE*nEv;
      out<<std::setw(13)<<SpecBins::Center(i)
         <<std::setw(13)<<fInP [i]/nr
         <<std::setw(13)<<fInA [i]/nr
         <<std::setw(13)<<fOutP[i]/nr
         <<std::setw(13)<<fOutA[i]/nr
         <<std::setw(13)<<fOutN[i]/nr<<"\n";
    }
    out.close();
    G4cout<<"Spectra written: "<<fname
          <<(fSweepMode?" (appended zone)":"")<<G4endl;
  }

  // ---- data members ----
  Options                     fOpts;
  DetectorConstruction*       fDetector=nullptr;
  std::vector<std::string>    fScoringNames;
  std::vector<G4double>       fScoringThick;
  std::vector<G4LogicalVolume*> fScoringLVs;
  G4LogicalVolume*            fShieldLV=nullptr;
  std::vector<G4double>       fEdep;
  std::vector<G4double>       fLastDose;
  std::vector<G4double>       fInP,fInA,fOutP,fOutA,fOutN;
  // sweep state
  bool        fSweepMode=false;
  std::string fCurrentMat;
  G4double    fCurrentTmm=0;
  bool        fFirstSpectraWrite=true;
  std::vector<SweepPoint> fSweepData;
};

// =============================================================================
// PrimaryGeneratorAction::RecordInput  (RunAction now complete)
// =============================================================================
void PrimaryGeneratorAction::RecordInput(bool isProton,G4double E){
  if(isProton) fRunAction->AddInP(E);
  else         fRunAction->AddInA(E);
}

G4Box* PrimaryGeneratorAction::GetCurrentWorldBox() const {
  return fDetector ? fDetector->GetWorldBox() : nullptr;
}

// =============================================================================
// EventAction
// =============================================================================
class EventAction : public G4UserEventAction {
public:
  EventAction(){}
  ~EventAction() override {}
  void BeginOfEventAction(const G4Event*) override {}
  void EndOfEventAction  (const G4Event*) override {}
};

// =============================================================================
// SteppingAction  — uses RunAction's refreshed LV pointers
// =============================================================================
class SteppingAction : public G4UserSteppingAction {
public:
  explicit SteppingAction(RunAction* ra) : fRA(ra) {}
  ~SteppingAction() override {}

  void UserSteppingAction(const G4Step* step) override {
    auto* pre  = step->GetPreStepPoint();
    auto* post = step->GetPostStepPoint();
    auto* preVol  = pre ->GetPhysicalVolume()
                    ? pre ->GetPhysicalVolume()->GetLogicalVolume() : nullptr;
    auto* postVol = post->GetPhysicalVolume()
                    ? post->GetPhysicalVolume()->GetLogicalVolume() : nullptr;

    // Transmitted particles at the downstream shield face only. Reject side exits
    // and particles that turn around and leave through the entrance face.
    if(preVol==fRA->GetShieldLV() && postVol && preVol!=postVol){
      auto* shieldSolid = dynamic_cast<G4Box*>(fRA->GetShieldLV()->GetSolid());
      const auto& dir = post->GetMomentumDirection();
      const auto& pos = post->GetPosition();
      if(shieldSolid && dir.z() > 0.0){
        const G4double rearZ = shieldSolid->GetZHalfLength();
        const G4double tol = 1.0*um;
        if(std::abs(pos.z() - rearZ) <= tol){
          G4double ke=post->GetKineticEnergy()/MeV;
          const G4String& nm=step->GetTrack()->GetParticleDefinition()->GetParticleName();
          auto* am=G4AnalysisManager::Instance();
          if(nm=="proton") { am->FillH1(0,ke); fRA->AddOutP(ke); }
          else if(nm=="alpha")  { am->FillH1(1,ke); fRA->AddOutA(ke); }
          else if(nm=="neutron"){ am->FillH1(2,ke); fRA->AddOutN(ke); }
        }
      }
    }

    // Dose in scoring volumes (bug-2 fix: preVol)
    if(preVol){
      const auto& scLVs=fRA->GetScoringLVs();
      for(std::size_t i=0;i<scLVs.size();++i){
        if(preVol==scLVs[i]){
          G4double ed=step->GetTotalEnergyDeposit();
          if(ed>0) fRA->AddEdep(i,ed);
          break;
        }
      }
    }
  }
private:
  RunAction* fRA;
};

// =============================================================================
// Sweep Tecplot writer  (dose vs thickness)
// =============================================================================
static void WriteDoseSweepTecplot(
    const std::vector<SweepPoint>& data,
    const std::vector<std::string>& matNames,
    const std::string& shieldMat,
    const Options& opts)
{
  const std::string fname="shieldSim_dose_sweep.dat";
  std::ofstream out(fname);
  if(!out){ G4cerr<<"Cannot write "<<fname<<G4endl; return; }

  out<<"TITLE = \"GCR Dose vs Shield Thickness – Geant4/FTFP_BERT\"\n";
  out<<"# Dose normalised per primary GCR particle\n";
  out<<"# Shield material: "<<shieldMat<<"\n";
  out<<"# Events per point: "<<opts.nEvents<<"\n";
  out<<"# Proton energy range: ["<<opts.eMinProton<<", "<<opts.eMaxProton<<"] MeV\n";
  out<<"# Alpha  energy range: ["<<opts.eMinAlpha <<", "<<opts.eMaxAlpha <<"] MeV\n";
  out<<"# Scoring volumes:";
  for(const auto& n:matNames) out<<" "<<n;
  out<<"\n";

  // Variable list
  out<<"VARIABLES = \"Thickness [mm]\" \"Areal_Density [g/cm2]\"";
  for(const auto& n:matNames)
    out<<" \"Dose_"<<n<<"_perPrimary [Gy]\"";
  out<<"\n";

  out<<"ZONE T=\""<<shieldMat<<" shield | GCR phi=550MV | "
     <<opts.nEvents<<" evt/pt\","
     <<" I="<<data.size()<<", DATAPACKING=POINT\n";

  out<<std::scientific<<std::setprecision(6);
  for(const auto& pt:data){
    out<<std::setw(14)<<pt.thickMM
       <<std::setw(14)<<pt.arealDensity;
    for(G4double d:pt.dose_Gy)
      out<<std::setw(14)<<d/gray;
    out<<"\n";
  }
  out.close();
  G4cout<<"\nDose sweep written: "<<fname
        <<"  ("<<data.size()<<" thickness points)"<<G4endl;
}

// =============================================================================
// PrintHelp
// =============================================================================
static void PrintHelp(){
  std::cout<<"\nUsage: shieldSim [options]\n\n";
  std::cout<<"SPECTRUM OPTIONS\n";
  std::cout<<"  --spectrum=<file>        Three-column tabulated spectrum file\n";
  std::cout<<"                           (E MeV, protonFlux, alphaFlux).\n";
  std::cout<<"                           Omit for built-in GCR Badhwar-O'Neill\n";
  std::cout<<"                           (phi=550 MV, solar minimum).\n";
  std::cout<<"  --emin=<MeV>             Global energy minimum for both species\n";
  std::cout<<"                           (default: 10 MeV).\n";
  std::cout<<"  --emax=<MeV>             Global energy maximum for both species\n";
  std::cout<<"                           (default: 100000 MeV = 100 GeV).\n";
  std::cout<<"  --emin-p=<MeV>           Proton energy minimum (overrides --emin).\n";
  std::cout<<"  --emax-p=<MeV>           Proton energy maximum (overrides --emax).\n";
  std::cout<<"  --emin-a=<MeV>           Alpha energy minimum  (overrides --emin).\n";
  std::cout<<"  --emax-a=<MeV>           Alpha energy maximum  (overrides --emax).\n";
  std::cout<<"\nSINGLE-RUN GEOMETRY\n";
  std::cout<<"  --shield=<M>:<t>         Shield NIST material and thickness in mm.\n";
  std::cout<<"                           Default: G4_Al:2\n";
  std::cout<<"  --scoring=<M>:<t>,...    Comma-separated scoring slabs (material:mm).\n";
  std::cout<<"                           Default: G4_WATER:1,G4_Si:1\n";
  std::cout<<"  --events=<n>             Primary particles per run. Default: 10000.\n";
  std::cout<<"\nDOSE-VS-THICKNESS SWEEP  (--sweep enables this mode)\n";
  std::cout<<"  --sweep                  Enable sweep; runs BeamOn for each thickness\n";
  std::cout<<"                           and writes shieldSim_dose_sweep.dat.\n";
  std::cout<<"  --sweep-material=<M>     Shield material to sweep over\n";
  std::cout<<"                           (default: same as --shield material).\n";
  std::cout<<"  --sweep-tmin=<mm>        Minimum sweep thickness. Default: 0.5 mm.\n";
  std::cout<<"  --sweep-tmax=<mm>        Maximum sweep thickness. Default: 50 mm.\n";
  std::cout<<"  --sweep-n=<n>            Number of thickness points. Default: 10.\n";
  std::cout<<"  --sweep-log              Use log-spaced thicknesses (default: linear).\n";
  std::cout<<"\nOUTPUT FILES\n";
  std::cout<<"  shieldSim_spectra.dat    Tecplot: input/output spectra per run\n";
  std::cout<<"                           (multi-zone in sweep mode).\n";
  std::cout<<"  shieldSim_dose_sweep.dat Tecplot: dose vs thickness (sweep mode).\n";
  std::cout<<"  shieldSimOutput[*].csv   G4AnalysisManager exit histograms.\n";
  std::cout<<"\nEXAMPLES\n";
  std::cout<<"  # Single run, 2 mm Al, full GCR spectrum\n";
  std::cout<<"  ./shieldSim --shield=G4_Al:2 --events=50000\n\n";
  std::cout<<"  # SEP-range protons only (1-300 MeV)\n";
  std::cout<<"  ./shieldSim --emin-p=1 --emax-p=300 --emin-a=10 --emax-a=300\n\n";
  std::cout<<"  # Dose vs thickness sweep, Al 0.5-30 mm, 15 log-spaced points\n";
  std::cout<<"  ./shieldSim --sweep --sweep-material=G4_Al\\\n";
  std::cout<<"              --sweep-tmin=0.5 --sweep-tmax=30\\\n";
  std::cout<<"              --sweep-n=15 --sweep-log --events=20000\n\n";
  std::cout<<"  # Sweep + restricted energy range + custom scoring\n";
  std::cout<<"  ./shieldSim --sweep --sweep-material=G4_POLYETHYLENE\\\n";
  std::cout<<"              --sweep-tmin=1 --sweep-tmax=100 --sweep-n=20\\\n";
  std::cout<<"              --emin=10 --emax=10000\\\n";
  std::cout<<"              --scoring=G4_WATER:1,G4_Si:1 --events=50000\n\n";
}

// =============================================================================
// ParseArguments
// =============================================================================
Options ParseArguments(int argc, char** argv){
  Options o;
  o.scoringMaterials.push_back({"G4_WATER",1.*mm});
  o.scoringMaterials.push_back({"G4_Si",   1.*mm});

  // Helper for "key=value" parsing
  auto strVal=[](const std::string& arg,const std::string& prefix)->std::string{
    return arg.substr(prefix.size());
  };

  for(int i=1;i<argc;++i){
    std::string a=argv[i];
    if(a=="-h"||a=="-help"||a=="--help"){ o.showHelp=true; }
    else if(a.find("--spectrum=")==0){ o.spectrumFile=strVal(a,"--spectrum="); }
    // energy limits — global
    else if(a.find("--emin=")==0){ o.eMinProton=o.eMinAlpha=std::stod(strVal(a,"--emin=")); }
    else if(a.find("--emax=")==0){ o.eMaxProton=o.eMaxAlpha=std::stod(strVal(a,"--emax=")); }
    // per species
    else if(a.find("--emin-p=")==0){ o.eMinProton=std::stod(strVal(a,"--emin-p=")); }
    else if(a.find("--emax-p=")==0){ o.eMaxProton=std::stod(strVal(a,"--emax-p=")); }
    else if(a.find("--emin-a=")==0){ o.eMinAlpha =std::stod(strVal(a,"--emin-a=")); }
    else if(a.find("--emax-a=")==0){ o.eMaxAlpha =std::stod(strVal(a,"--emax-a=")); }
    // single-run geometry
    else if(a.find("--shield=")==0){
      std::string v=strVal(a,"--shield=");
      auto p=v.find(':');
      if(p!=std::string::npos){ o.shieldMaterial=v.substr(0,p);
                                 o.shieldThickness=std::stod(v.substr(p+1))*mm; }
      else o.shieldMaterial=v;
    }
    else if(a.find("--scoring=")==0){
      o.scoringMaterials.clear();
      std::stringstream ss(strVal(a,"--scoring=")); std::string item;
      while(std::getline(ss,item,',')){
        auto p=item.find(':');
        std::string mat=item; G4double t=1.;
        if(p!=std::string::npos){ mat=item.substr(0,p);
                                   t=std::stod(item.substr(p+1)); }
        o.scoringMaterials.push_back({mat,t*mm});
      }
    }
    else if(a.find("--events=")==0){ o.nEvents=std::stoi(strVal(a,"--events=")); }
    // sweep
    else if(a=="--sweep"){ o.doSweep=true; }
    else if(a.find("--sweep-material=")==0){ o.sweepMaterial=strVal(a,"--sweep-material="); }
    else if(a.find("--sweep-tmin=")==0)    { o.sweepTmin=std::stod(strVal(a,"--sweep-tmin=")); }
    else if(a.find("--sweep-tmax=")==0)    { o.sweepTmax=std::stod(strVal(a,"--sweep-tmax=")); }
    else if(a.find("--sweep-n=")==0)       { o.sweepN=std::stoi(strVal(a,"--sweep-n=")); }
    else if(a=="--sweep-log")              { o.sweepLog=true; }
    else { G4cout<<"Unknown option ignored: "<<a<<G4endl; }
  }

  if(o.sweepMaterial.empty()) o.sweepMaterial=o.shieldMaterial;

  if(o.nEvents<=0)
    G4Exception("ParseArguments","BadInput",FatalException,"--events must be > 0");
  if(o.eMinProton<=0 || o.eMinAlpha<=0)
    G4Exception("ParseArguments","BadInput",FatalException,"Energy minima must be > 0 MeV");
  if(o.eMinProton>o.eMaxProton)
    G4Exception("ParseArguments","BadInput",FatalException,"Require --emin-p <= --emax-p");
  if(o.eMinAlpha>o.eMaxAlpha)
    G4Exception("ParseArguments","BadInput",FatalException,"Require --emin-a <= --emax-a");
  if(o.shieldThickness<0)
    G4Exception("ParseArguments","BadInput",FatalException,"Shield thickness must be >= 0");
  if(o.scoringMaterials.empty())
    G4Exception("ParseArguments","BadInput",FatalException,"At least one scoring slab is required");
  for(const auto& sm : o.scoringMaterials){
    if(sm.second<=0)
      G4Exception("ParseArguments","BadInput",FatalException,"Scoring slab thickness must be > 0");
  }
  if(o.doSweep){
    if(o.sweepN<=0)
      G4Exception("ParseArguments","BadInput",FatalException,"--sweep-n must be > 0");
    if(o.sweepTmin<=0 || o.sweepTmax<=0)
      G4Exception("ParseArguments","BadInput",FatalException,"Sweep thicknesses must be > 0 mm");
    if(o.sweepTmin>o.sweepTmax)
      G4Exception("ParseArguments","BadInput",FatalException,"Require --sweep-tmin <= --sweep-tmax");
  }
  return o;
}

// =============================================================================
// main
// =============================================================================
int main(int argc,char** argv){
  Options opts=ParseArguments(argc,argv);
  if(opts.showHelp){ PrintHelp(); return 0; }

  // Collect scoring names and thicknesses (for RunAction constructor)
  std::vector<std::string> scNames;
  std::vector<G4double>    scThick;
  for(const auto& s:opts.scoringMaterials){
    scNames.push_back(s.first);
    scThick.push_back(s.second);
  }

  // Config summary
  G4cout<<"============ shieldSim configuration ============"<<G4endl;
  G4cout<<"Spectrum  : "
        <<(opts.spectrumFile.empty()?"GCR Badhwar-O'Neill phi=550MV"
                                    :opts.spectrumFile)<<G4endl;
  G4cout<<"E proton  : ["<<opts.eMinProton<<", "<<opts.eMaxProton<<"] MeV"<<G4endl;
  G4cout<<"E alpha   : ["<<opts.eMinAlpha <<", "<<opts.eMaxAlpha <<"] MeV"<<G4endl;
  G4cout<<"Events    : "<<opts.nEvents<<" per run"<<G4endl;
  G4cout<<"Scoring   :";
  for(const auto& s:opts.scoringMaterials)
    G4cout<<"  "<<s.first<<":"<<s.second/mm<<"mm";
  G4cout<<G4endl;
  if(opts.doSweep){
    G4cout<<"Mode      : SWEEP"<<G4endl;
    G4cout<<"Material  : "<<opts.sweepMaterial<<G4endl;
    G4cout<<"Thickness : "<<opts.sweepTmin<<" – "<<opts.sweepTmax<<" mm  ("
          <<opts.sweepN<<" points, "<<(opts.sweepLog?"log":"linear")<<")"<<G4endl;
  } else {
    G4cout<<"Mode      : single run"<<G4endl;
    G4cout<<"Shield    : "<<opts.shieldMaterial
          <<", "<<opts.shieldThickness/mm<<" mm"<<G4endl;
  }
  G4cout<<"================================================="<<G4endl;

  // Build thickness list for sweep
  std::vector<G4double> thicknesses;
  if(opts.doSweep){
    if(opts.sweepN==1){
      thicknesses.push_back(opts.sweepTmin);
    } else {
      for(G4int i=0;i<opts.sweepN;++i){
        G4double frac=G4double(i)/(opts.sweepN-1);
      G4double t;
      if(opts.sweepLog)
        t=opts.sweepTmin*std::pow(opts.sweepTmax/opts.sweepTmin,frac);
      else
        t=opts.sweepTmin+(opts.sweepTmax-opts.sweepTmin)*frac;
      thicknesses.push_back(t);
      }
    }
  } else {
    thicknesses.push_back(opts.shieldThickness/mm);
  }

  // Set initial shield thickness
  if(opts.doSweep){
    opts.shieldMaterial  = opts.sweepMaterial;
    opts.shieldThickness = thicknesses[0]*mm;
  }

  // Construct run manager and all user classes
  auto* runManager = new G4RunManager();

  auto* detector = new DetectorConstruction(opts);
  runManager->SetUserInitialization(detector);

  G4PhysListFactory factory;
  auto* physList=factory.GetReferencePhysList("FTFP_BERT");
  if(!physList) G4Exception("main","PhysList",FatalException,
                            "Cannot create FTFP_BERT.");
  physList->SetVerboseLevel(0);
  runManager->SetUserInitialization(physList);

  auto* runAction=new RunAction(opts,detector,scNames,scThick);
  runManager->SetUserAction(new PrimaryGeneratorAction(opts,runAction,detector));
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(new EventAction());
  runManager->SetUserAction(new SteppingAction(runAction));

  runManager->Initialize();

  // ---- Sweep / single-run loop ----
  for(G4int step=0;step<(G4int)thicknesses.size();++step){
    G4double tMM=thicknesses[step];

    if(step>0){
      // Update geometry for next thickness
      detector->SetShieldThickness(tMM*mm);
      runManager->ReinitializeGeometry();
    }

    if(opts.doSweep){
      runAction->SetSweepMode(true,opts.sweepMaterial,tMM);
      G4cout<<"\n>>> Sweep point "<<step+1<<"/"<<thicknesses.size()
            <<"  t = "<<std::fixed<<std::setprecision(2)<<tMM<<" mm"<<G4endl;
    }

    runManager->BeamOn(opts.nEvents);

    // Collect result for sweep Tecplot
    if(opts.doSweep){
      // Get areal density from the shield material
      auto* nist=G4NistManager::Instance();
      auto* mat=nist->FindOrBuildMaterial(opts.sweepMaterial);
      G4double rho_gcc = mat ? mat->GetDensity()/(g/cm3) : 0.;
      G4double areal   = rho_gcc * tMM * 0.1;  // g/cm² (1 mm = 0.1 cm)

      SweepPoint sp;
      sp.thickMM      = tMM;
      sp.arealDensity = areal;
      sp.dose_Gy      = runAction->GetLastDose();
      runAction->AppendSweepPoint(sp);
    }
  }

  // Write dose sweep Tecplot
  if(opts.doSweep)
    WriteDoseSweepTecplot(runAction->GetSweepData(),scNames,
                          opts.sweepMaterial,opts);

  delete runManager;
  G4cout<<"\nSimulation completed."<<G4endl;
  return 0;
}
