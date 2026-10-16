#include "DetectorConstruction.hh"

#include "AbsorberGeometry.hh"
#include "CaseDefinition.hh"
#include "MaterialCatalog.hh"
#include "MaterialFactory.hh"

#include <G4Box.hh>
#include <G4Exception.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4UserLimits.hh>
#include <G4VPhysicalVolume.hh>

DetectorConstruction::DetectorConstruction(const Options& opts) : fOpts(opts) {}
DetectorConstruction::~DetectorConstruction() {}

void DetectorConstruction::SetShieldThickness(G4double t){ fOpts.shieldThickness=t; }
void DetectorConstruction::SetActiveAbsorber(std::size_t index){ fOpts.activeAbsorberIndex=index; }
const Options& DetectorConstruction::GetOptions() const  { return fOpts; }
G4LogicalVolume* DetectorConstruction::GetShieldLV() const { return fShieldLV; }
const std::vector<G4LogicalVolume*>& DetectorConstruction::GetShieldLVs() const { return fShieldLVs; }
const std::vector<G4LogicalVolume*>& DetectorConstruction::GetScoringLVs() const { return fScoringLVs; }
G4Box* DetectorConstruction::GetWorldBox() const { return fWorldBox; }
G4double DetectorConstruction::GetSourcePlaneHalfX() const {
  if(fOpts.caseDefinition) return fOpts.caseDefinition->sourceSizeXmm*mm/2.0;
  return fWorldBox ? fWorldBox->GetXHalfLength() : 0.0;
}
G4double DetectorConstruction::GetSourcePlaneHalfY() const {
  if(fOpts.caseDefinition) return fOpts.caseDefinition->sourceSizeYmm*mm/2.0;
  return fWorldBox ? fWorldBox->GetYHalfLength() : 0.0;
}
G4double DetectorConstruction::GetSourcePlaneCenterX() const {
  return fOpts.caseDefinition ? fOpts.caseDefinition->sourceCenterXmm*mm : 0.0;
}
G4double DetectorConstruction::GetSourcePlaneCenterY() const {
  return fOpts.caseDefinition ? fOpts.caseDefinition->sourceCenterYmm*mm : 0.0;
}
G4double DetectorConstruction::GetSourcePlaneZ() const {
  if(fOpts.caseDefinition) return fOpts.caseDefinition->sourceZmm*mm;
  return fWorldBox ? -fWorldBox->GetZHalfLength()+1.0*um : 0.0;
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  // Construct() may be called multiple times in sweep mode.  Clear all cached
  // geometry pointers before rebuilding so that RunAction can refresh them for
  // the current geometry.
  fWorldBox = nullptr;
  fShieldLV = nullptr;
  fShieldLVs.clear();
  fScoringLVs.clear();

  // The versioned case path is deliberately separate from the compatibility
  // CLI path below.  Every dimension and material comes from the resolved case;
  // no legacy 5-cm world or 1-mm spacing default can leak into this branch.
  if(fOpts.caseDefinition){
    const auto& c=*fOpts.caseDefinition;
    if(fOpts.activeAbsorberIndex>=c.absorbers.size())
      G4Exception("Detector","AbsorberIndex",FatalException,"Active absorber index is out of range");

    auto* vac=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    auto* wSol=new G4Box("World",c.worldSizeXmm*mm/2.0,c.worldSizeYmm*mm/2.0,c.worldSizeZmm*mm/2.0);
    auto* wLog=new G4LogicalVolume(wSol,vac,"World");
    fWorldBox=wSol;
    auto* wPhys=new G4PVPlacement(nullptr,G4ThreeVector(),wLog,"World",nullptr,false,0);

    if(!fMaterialFactory)
      fMaterialFactory=std::make_unique<shield_sim::MaterialFactory>(*c.materials);
    G4double currentFront=c.shieldFrontZmm*mm;
    for(std::size_t i=0;i<c.shieldLayers.size();++i){
      const auto& layer=c.shieldLayers[i];
      auto* material=fMaterialFactory->Resolve(layer.materialId,"shield");
      const G4double thickness=layer.thicknessMm*mm;
      const std::string stableId="shield/"+layer.id;
      auto* solid=new G4Box((stableId+"/solid").c_str(),c.shieldSizeXmm*mm/2.0,
                            c.shieldSizeYmm*mm/2.0,thickness/2.0);
      auto* logical=new G4LogicalVolume(solid,material,(stableId+"/logical").c_str());
      if(c.maxStepMm>0.0) logical->SetUserLimits(new G4UserLimits(c.maxStepMm*mm));
      new G4PVPlacement(nullptr,G4ThreeVector(c.shieldCenterXmm*mm,c.shieldCenterYmm*mm,
                                               currentFront+thickness/2.0),logical,
                        (stableId+"/physical").c_str(),wLog,false,static_cast<G4int>(i),true);
      fShieldLVs.push_back(logical);
      currentFront+=thickness;
    }
    // The last ordered layer owns the downstream stack face used for exit
    // scoring.  Internal layer interfaces are never counted as shield exits.
    fShieldLV=fShieldLVs.back();

    fAbsorberGeometry=std::make_unique<shield_sim::SlabAbsorberGeometry>(
        c.absorbers[fOpts.activeAbsorberIndex]);
    const auto built=fAbsorberGeometry->Build(wLog,*fMaterialFactory,c.maxStepMm);
    for(const auto& region:built.sensitiveRegions) fScoringLVs.push_back(region.logicalVolume);
    return wPhys;
  }

  auto* nist=G4NistManager::Instance();

  // Finite transverse detector size.  In beam mode, dose per primary depends on
  // this detector area because the full scoring slab mass is used in Edep/mass.
  G4double wXY=5.*cm;
  G4double scSum=0;
  for(const auto& s:fOpts.scoringMaterials) scSum+=s.second;

  // Include a 1 mm upstream gap, the shield, scoring stack, and a 1 mm
  // downstream gap.  The source starts just inside the upstream world face.
  G4double wZ=2.*mm+fOpts.shieldThickness+scSum+2.*mm;

  auto* vac  =nist->FindOrBuildMaterial("G4_Galactic");
  auto* wSol =new G4Box("World",wXY/2,wXY/2,wZ/2);
  auto* wLog =new G4LogicalVolume(wSol,vac,"World");
  fWorldBox = wSol;
  auto* wPhys=new G4PVPlacement(nullptr,G4ThreeVector(),wLog,"World",nullptr,false,0);

  auto* shMat=FindOrBuildShieldMaterial(fOpts.shieldMaterial);
  if(!shMat) G4Exception("Detector","ShieldMat",FatalException,
      ("Shield material "+fOpts.shieldMaterial+" not found").c_str());

  // The shield local origin is its center.  Its downstream face is therefore at
  // local z = +shieldThickness/2.  This is important for SteppingAction, which
  // transforms global step-point coordinates into this local shield frame before
  // scoring transmitted particles.
  G4double shHz=fOpts.shieldThickness/2;
  auto* shS=new G4Box("Shield",wXY/2,wXY/2,shHz);
  auto* shL=new G4LogicalVolume(shS,shMat,"Shield");

  // Optional numerical step limiter.  The user limit itself is attached here,
  // while shieldSim.cc registers G4StepLimiterPhysics when --max-step is used.
  // This separation follows the Geant4 design: geometry volumes carry the
  // maximum permitted step, and the step-limiter process enforces it.  The
  // option is intended for Layer-3 numerical convergence tests, especially
  // thin-target TID and LET-spectrum checks where a very long charged-particle
  // step could smear local response estimates.
  if(fOpts.maxStepLength>0.0)
    shL->SetUserLimits(new G4UserLimits(fOpts.maxStepLength));

  fShieldLV = shL;
  fShieldLVs.push_back(shL);
  G4double zSh=1.*mm+shHz;
  new G4PVPlacement(nullptr,G4ThreeVector(0,0,zSh-wZ/2),shL,"Shield",wLog,false,0,true);

  // Downstream scoring stack.  Each slab is immediately adjacent to the
  // preceding volume.  Energy deposited in these logical volumes is accumulated
  // by SteppingAction and normalized by RunAction at the end of each run.
  //
  // Scoring slabs use FindOrBuildDetectorMaterial() rather than the shielding
  // resolver.  The detector resolver knows about tissue targets and electronics
  // materials such as BFO, CNS, SoftTissue, SiO2, GaAs, and InGaAs.  It still
  // falls back to the shielding catalog and to raw G4_* names, so a user can
  // deliberately score dose in a shielding material as well.
  G4double curZ=1.*mm+fOpts.shieldThickness;
  for(std::size_t i=0;i<fOpts.scoringMaterials.size();++i){
    const auto& sm=fOpts.scoringMaterials[i];
    auto* mat=FindOrBuildDetectorMaterial(sm.first);
    if(!mat) G4Exception("Detector","ScoreMat",FatalException,
        ("Scoring material "+sm.first+" not found").c_str());
    G4double hz=sm.second/2;
    std::string idx=std::to_string(i);
    std::string solidName="ScoringSolid_"+idx+"_"+sm.first;
    std::string lvName   ="ScoringLV_"+idx+"_"+sm.first;
    std::string pvName   ="ScoringPV_"+idx+"_"+sm.first;
    auto* ds=new G4Box(solidName.c_str(),wXY/2,wXY/2,hz);
    auto* dl=new G4LogicalVolume(ds,mat,lvName.c_str());

    // Apply the same optional step limiter to each scoring slab.  The scoring
    // volumes are where TID is accumulated, so this is the most important place
    // to test step-size convergence.  Leaving --max-step unset preserves the
    // normal physics-list stepping behavior.
    if(fOpts.maxStepLength>0.0)
      dl->SetUserLimits(new G4UserLimits(fOpts.maxStepLength));

    fScoringLVs.push_back(dl);
    curZ+=hz;
    new G4PVPlacement(nullptr,G4ThreeVector(0,0,curZ-wZ/2),
                      dl,pvName.c_str(),wLog,false,(G4int)i,true);
    curZ+=hz;
  }
  return wPhys;
}
