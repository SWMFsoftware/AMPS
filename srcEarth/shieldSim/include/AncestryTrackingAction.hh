#ifndef SHIELDSIM_ANCESTRY_TRACKING_ACTION_HH
#define SHIELDSIM_ANCESTRY_TRACKING_ACTION_HH

/*
 * Propagated particle ancestry for Stage-3 response accounting.
 *
 * Geant4 exposes a track's immediate parent ID, but the parent track may no
 * longer exist when a downstream scorer is called. This user information is
 * copied from a parent to every secondary after tracking, preserving the root
 * primary identity and generation without relying on a global track-ID map.
 */

#include <G4UserTrackingAction.hh>
#include <G4VUserTrackInformation.hh>

#include <string>

class G4Track;

class TrackAncestryInfo final : public G4VUserTrackInformation {
public:
  std::string rootParticle;
  int rootPdg = 0;
  int rootZ = 0;
  int rootA = 0;
  std::string parentParticle;
  int parentPdg = 0;
  int generation = 0;
  void Print() const override;
};

class AncestryTrackingAction final : public G4UserTrackingAction {
public:
  void PreUserTrackingAction(const G4Track* track) override;
  void PostUserTrackingAction(const G4Track* track) override;
};

#endif
