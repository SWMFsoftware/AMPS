#ifndef SRCSEP_REDUCED_SHOCK_GEOMETRY_H
#define SRCSEP_REDUCED_SHOCK_GEOMETRY_H
#include "util/sep_reduced_geometry.h"
#include <list>
namespace SEP {class cFieldLine;
namespace ReducedShock {namespace Native {
// Retained immutable material realization for source/boundary queries. Live
// topology is published only after all particle readers join. Arc coordinates
// refer to departure chords; a segment fraction is a material label, while a
// time-dependent physical position must advect both bounding vertices.
struct MaterialView {
  IntervalView interval;
  std::vector<Position> vertices;
  std::vector<double> arc;
  std::vector<Connection> departureConnections;
  std::uint64_t geometryGeneration=0;
};
Core::Result<std::shared_ptr<const MaterialView>> CaptureMaterialView();
// Reconstruct derived connectivity without publishing geometry or altering a
// live reader. The checkpoint supplies decoded material positions/generation;
// the configured core event supplies the same physical functions and epochs.
Core::Result<std::shared_ptr<const MaterialView>> PrepareCheckpointView(
    const IntervalView&,std::uint64_t generation,const std::vector<Position>&);
void CommitCheckpointView(std::shared_ptr<const MaterialView>*) noexcept;
// Native ownership begins before source histories/sample buffers are opened.
// These methods consume the existing active runtime geometry configuration;
// no origin, wind, magnetic or thermal model is chosen by this adapter.
Core::Status BuildInitialLine(std::list<SEP::cFieldLine>*);
Core::Status RegisterLines(double epochS);
// Stage and collectively validate all initial native state and connectivity
// before publishing any datum or metadata. No physical timestep is advanced.
Core::Status InstallBackground(double epochS);
Core::Status Advance(double departureS,double arrivalS);
GeometryControls Controls(double intervalS,double segmentLengthM);
Core::Status RefreshConnections(double epochS);
Core::Result<std::vector<Connection>> ConnectionsAt(double epochS);
Core::Status UpdateShockDistances(bool report);
// First-order proximity estimate for an internal deterministic mover shell.
// This is not a stochastic first-passage or a shock acceleration operator.
Core::Result<double> CrossingTimeEstimate(double coordinate,double relativeStreamingMPerS);
}}}
#endif
