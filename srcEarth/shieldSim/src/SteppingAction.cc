#include "SteppingAction.hh"

#include "AncestryTrackingAction.hh"
#include "RunAction.hh"

#include <G4AffineTransform.hh>
#include <G4AnalysisManager.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4ParticleDefinition.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4SystemOfUnits.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4ThreeVector.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VProcess.hh>

#include <cmath>

SteppingAction::SteppingAction(RunAction* ra) : fRA(ra) {}
SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
  auto* pre  = step->GetPreStepPoint();
  auto* post = step->GetPostStepPoint();
  auto* preVol  = pre ->GetPhysicalVolume()
                  ? pre ->GetPhysicalVolume()->GetLogicalVolume() : nullptr;
  auto* postVol = post->GetPhysicalVolume()
                  ? post->GetPhysicalVolume()->GetLogicalVolume() : nullptr;

  // Score particles that leave the shield through the downstream face.  Side
  // exits and particles leaving through the upstream face are rejected using the
  // local crossing position and local momentum direction.
  if(preVol==fRA->GetShieldLV() && postVol && preVol!=postVol){
    auto* shieldSolid = dynamic_cast<G4Box*>(fRA->GetShieldLV()->GetSolid());
    if(shieldSolid){
      auto touchable = pre->GetTouchableHandle();
      auto* history = touchable->GetHistory();

      if(history){
        const G4AffineTransform& globalToLocal = history->GetTopTransform();
        const G4ThreeVector localPos =
            globalToLocal.TransformPoint(post->GetPosition());
        const G4ThreeVector localDir =
            globalToLocal.TransformAxis(post->GetMomentumDirection());

        if(localDir.z() > 0.0){
          const G4double rearZ = shieldSolid->GetZHalfLength();
          const G4double tol = 1.0*um;
          if(std::abs(localPos.z() - rearZ) <= tol){
            G4double ke=post->GetKineticEnergy()/MeV;
            const G4String& nm=step->GetTrack()->GetParticleDefinition()->GetParticleName();
            const auto* definition=step->GetTrack()->GetParticleDefinition();
            const auto* ancestry=dynamic_cast<const TrackAncestryInfo*>(
                step->GetTrack()->GetUserInformation());
            fRA->RecordShieldExit(nm,definition->GetPDGEncoding(),definition->GetAtomicNumber(),
                static_cast<int>(definition->GetAtomicMass()),step->GetTrack()->GetParentID()==0,
                ancestry?ancestry->rootParticle:nm,
                ancestry?ancestry->generation:(step->GetTrack()->GetParentID()==0?0:1),ke);
            auto* am=G4AnalysisManager::Instance();
            if(nm=="proton") {
              am->FillH1(0,ke); fRA->AddOutP(ke);
              if(step->GetTrack()->GetParentID()==0) fRA->AddPrimaryProtonExit(ke);
            }
            else if(nm=="alpha")  { am->FillH1(1,ke); fRA->AddOutA(ke); }
            else if(nm=="neutron"){ am->FillH1(2,ke); fRA->AddOutN(ke); }

            // Optional diagnostic used by Layer-2 tests.  The record is made
            // only after the same downstream-face criteria used for physics
            // scoring are satisfied, so the file is a direct audit trail of
            // the accepted exit crossings.
            const G4ThreeVector globalPos = post->GetPosition();
            fRA->RecordExitParticle(nm,ke,
                                    globalPos.x()/mm, globalPos.y()/mm, globalPos.z()/mm,
                                    localPos.x()/mm,  localPos.y()/mm,  localPos.z()/mm,
                                    localDir.x(),     localDir.y(),     localDir.z());
          }
        }
      }
    }
  }

  // Stage-5 response is defined on the absorber's upstream entry surface, not
  // at the shield rear face and not from energy deposited later in the slab.
  // Boundary classification is delegated to the absorber geometry provider.
  // slab.v1 accepts local z=-halfLength with local direction z>0; a future
  // arbitrary-shape provider can define an explicit entrance patch without
  // changing this transport hook or the damage-response calculation.
  if(postVol && preVol!=postVol){
    const auto& absorberLVs=fRA->GetScoringLVs();
    for(std::size_t absorber=0;absorber<absorberLVs.size();++absorber){
      if(postVol!=absorberLVs[absorber]) continue;
      const G4double fluenceWeight=
          fRA->AbsorberEntryFluenceWeightPerCm2(absorber,*post);
      if(fluenceWeight>0.0){
        const auto* definition=step->GetTrack()->GetParticleDefinition();
        fRA->RecordAbsorberEntry(absorber,definition->GetParticleName(),
                                 post->GetKineticEnergy()/MeV,fluenceWeight);
      }
      break;
    }
  }

  // Dose scoring.  Use the pre-step volume because the energy deposition for a
  // step belongs to the material where the step occurred, even if the post-step
  // point is already in the next volume at a boundary.
  if(preVol){
    const auto& shieldLayers=fRA->GetShieldLVs();
    for(std::size_t layer=0;layer<shieldLayers.size();++layer){
      if(preVol!=shieldLayers[layer]) continue;
      fRA->RecordLayerStep(layer,step->GetTotalEnergyDeposit(),
                          step->GetNonIonizingEnergyDeposit(),step->GetStepLength());
      const auto* track=step->GetTrack();
      if(track->GetParentID()!=0 && track->GetCurrentStepNumber()==1){
        const auto* definition=track->GetParticleDefinition();
        const auto* ancestry=dynamic_cast<const TrackAncestryInfo*>(track->GetUserInformation());
        const auto* creator=track->GetCreatorProcess();
        fRA->RecordSecondaryBirth(layer,definition->GetParticleName(),
            definition->GetPDGEncoding(),definition->GetAtomicNumber(),
            static_cast<int>(definition->GetAtomicMass()),
            ancestry?ancestry->rootParticle:"unavailable",
            ancestry?ancestry->generation:1,
            creator?creator->GetProcessName():"unavailable",
            track->GetVertexKineticEnergy()/MeV);
      }
      break;
    }
    const auto& scLVs=fRA->GetScoringLVs();
    for(std::size_t i=0;i<scLVs.size();++i){
      if(preVol==scLVs[i]){
        const G4double total=step->GetTotalEnergyDeposit();
        const G4double nonIon=step->GetNonIonizingEnergyDeposit();
        if(total>0 || nonIon>0){
          const auto* track=step->GetTrack();
          const auto* definition=track->GetParticleDefinition();
          const auto* ancestry=dynamic_cast<const TrackAncestryInfo*>(
              track->GetUserInformation());
          const auto* creator=track->GetCreatorProcess();
          const auto* defining=post->GetProcessDefinedStep();
          fRA->RecordAbsorberStep(i,definition->GetParticleName(),
              track->GetParentID()==0,
              ancestry?ancestry->rootParticle:definition->GetParticleName(),
              ancestry?ancestry->generation:(track->GetParentID()==0?0:1),
              creator?creator->GetProcessName():"primary_no_creator_process",
              defining?defining->GetProcessName():"no_step_defining_process",
              total,nonIon);
        }
        break;
      }
    }
  }
}
