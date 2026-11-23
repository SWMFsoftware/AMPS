#ifndef AMPS_FIELDLINE_FRAME_H
#define AMPS_FIELDLINE_FRAME_H
#include <cmath>
#include <stdexcept>

namespace AMPS { namespace FieldLine {
// Wire values are explicit and shared with PT_wrapper.f90. A geometry owner
// declares the coordinate convention; particle mover selection never chooses
// it. Unknown is intentionally unusable, including when U happens to be zero.
enum class CoordinateFrame { Unknown=0, Fixed=1, PlasmaComoving=2 };
inline const char* CoordinateFrameName(CoordinateFrame frame) {
  switch(frame) {
    case CoordinateFrame::Fixed:return "fixed-geometry-v1";
    case CoordinateFrame::PlasmaComoving:return "plasma-comoving-v1";
    default:throw std::invalid_argument("field-line coordinate frame is undeclared");
  }
}
// physicalParallelVelocity is U dot tangent [m/s], not a particle velocity.
// Moving material labels already carry the bulk translation through geometry.
// Return exact zero in that frame without erasing physical U or grad(U), which
// remain needed for momentum evolution, source fluxes and wave transforms.
inline double CoordinateAdvection(CoordinateFrame frame,double physicalParallelVelocity) {
  CoordinateFrameName(frame);
  if(!std::isfinite(physicalParallelVelocity))
    throw std::invalid_argument("physical parallel plasma velocity is nonfinite");
  return frame==CoordinateFrame::PlasmaComoving?0.0:physicalParallelVelocity;
}
// Native S is segment index + fraction; producer node labels are contiguous
// integers. A prefix insertion/removal changes its index by oldFirst-newFirst,
// preserving the same material pair and fraction. Domain checks belong to the
// importing owner, which knows which segments survive and which are absorbed.
inline double RemapCoordinate(double coordinate,int oldFirst,int newFirst) {
  if(!std::isfinite(coordinate)||coordinate<0)
    throw std::invalid_argument("cannot remap an invalid material particle coordinate");
  const double result=coordinate+(static_cast<double>(oldFirst)-newFirst);
  if(!std::isfinite(result))throw std::invalid_argument("material coordinate remap overflow");
  return result;
}
} }
#endif
