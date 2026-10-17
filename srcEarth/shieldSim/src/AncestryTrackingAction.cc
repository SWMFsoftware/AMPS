#include "AncestryTrackingAction.hh"

#include <G4ParticleDefinition.hh>
#include <G4Track.hh>
#include <G4TrackingManager.hh>
#include <G4ios.hh>

void TrackAncestryInfo::Print() const {
  G4cout<<"root="<<rootParticle<<" generation="<<generation
        <<" parent="<<parentParticle<<G4endl;
}

void AncestryTrackingAction::PreUserTrackingAction(const G4Track* track){
  if(track->GetUserInformation()) return;
  auto* information=new TrackAncestryInfo();
  const auto* definition=track->GetParticleDefinition();
  information->rootParticle=definition->GetParticleName();
  information->rootPdg=definition->GetPDGEncoding();
  information->rootZ=definition->GetAtomicNumber();
  information->rootA=definition->GetAtomicMass();
  information->generation=track->GetParentID()==0 ? 0 : 1;
  information->parentParticle=track->GetParentID()==0 ? "none" : "unavailable";
  // The Geant4 callback uses a const track even though this is the supported
  // hook for attaching user tracking information.
  const_cast<G4Track*>(track)->SetUserInformation(information);
}

void AncestryTrackingAction::PostUserTrackingAction(const G4Track* track){
  const auto* parentInfo=dynamic_cast<const TrackAncestryInfo*>(track->GetUserInformation());
  auto* secondaries=fpTrackingManager->GimmeSecondaries();
  if(!secondaries) return;
  const auto* parentDefinition=track->GetParticleDefinition();
  for(G4Track* child:*secondaries){
    if(child->GetUserInformation()) continue;
    auto* information=new TrackAncestryInfo();
    if(parentInfo){
      information->rootParticle=parentInfo->rootParticle;
      information->rootPdg=parentInfo->rootPdg;
      information->rootZ=parentInfo->rootZ;
      information->rootA=parentInfo->rootA;
      information->generation=parentInfo->generation+1;
    } else {
      information->rootParticle=parentDefinition->GetParticleName();
      information->rootPdg=parentDefinition->GetPDGEncoding();
      information->rootZ=parentDefinition->GetAtomicNumber();
      information->rootA=parentDefinition->GetAtomicMass();
      information->generation=1;
    }
    information->parentParticle=parentDefinition->GetParticleName();
    information->parentPdg=parentDefinition->GetPDGEncoding();
    child->SetUserInformation(information);
  }
}
