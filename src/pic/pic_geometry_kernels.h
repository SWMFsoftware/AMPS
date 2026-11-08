#ifndef AMPS_PIC_GEOMETRY_KERNELS_H
#define AMPS_PIC_GEOMETRY_KERNELS_H

/*
 * Dependency-light geometric kernels shared by ray tracing and verification.
 *
 * Coordinates are in the caller's length unit and the ray is
 *     x(t) = origin + t * direction.
 * Therefore t carries the inverse unit of ``direction``; when ``direction``
 * is a unit vector, t is path length.  All bounds are closed within the
 * supplied non-negative tolerance.  Degenerate boxes, non-finite input and a
 * zero direction are typed failures instead of divisions by zero.
 */

#include <algorithm>
#include <cmath>
#include <limits>

namespace PIC {
namespace GeometryKernels {

enum class RayBoxStatus {
  Success,
  NoIntersection,
  InvalidInput
};

struct RayBoxInterval {
  RayBoxStatus status;
  double entry_parameter;
  double exit_parameter;
  int entry_face;
  int exit_face;

  RayBoxInterval()
      : status(RayBoxStatus::InvalidInput), entry_parameter(0.0),
        exit_parameter(0.0), entry_face(-1), exit_face(-1) {}
};

inline bool FiniteVector3(const double* value) {
  return value != 0 && std::isfinite(value[0]) && std::isfinite(value[1]) &&
         std::isfinite(value[2]);
}

/* Slab clipping with explicit face identities (0=-x,1=+x,...).  ``ray_only``
 * rejects an intersection lying wholly behind the origin. */
inline RayBoxInterval ClipRayToAxisAlignedBox(
    const double* box_min, const double* box_max, const double* origin,
    const double* direction, double tolerance, bool ray_only=true) {
  RayBoxInterval result;
  if (!FiniteVector3(box_min) || !FiniteVector3(box_max) ||
      !FiniteVector3(origin) || !FiniteVector3(direction) ||
      !std::isfinite(tolerance) || tolerance < 0.0) return result;

  double t_min=-std::numeric_limits<double>::infinity();
  double t_max= std::numeric_limits<double>::infinity();
  int face_min=-1,face_max=-1;
  bool nonzero_direction=false;
  for (int axis=0; axis<3; ++axis) {
    if (!(box_max[axis] > box_min[axis])) return result;
    const double lower=box_min[axis]-tolerance;
    const double upper=box_max[axis]+tolerance;
    const double d=direction[axis];
    if (d == 0.0) {
      if (origin[axis] < lower || origin[axis] > upper) {
        result.status=RayBoxStatus::NoIntersection;
        return result;
      }
      continue;
    }
    nonzero_direction=true;
    double t0=(box_min[axis]-origin[axis])/d;
    double t1=(box_max[axis]-origin[axis])/d;
    int f0=2*axis, f1=2*axis+1;
    if (t0>t1) {
      std::swap(t0,t1);
      std::swap(f0,f1);
    }
    if (t0>t_min) { t_min=t0; face_min=f0; }
    if (t1<t_max) { t_max=t1; face_max=f1; }
    if (t_min>t_max+tolerance) {
      result.status=RayBoxStatus::NoIntersection;
      return result;
    }
  }
  if (!nonzero_direction) return result;
  if (ray_only && t_max < -tolerance) {
    result.status=RayBoxStatus::NoIntersection;
    return result;
  }
  result.status=RayBoxStatus::Success;
  result.entry_parameter=t_min;
  result.exit_parameter=t_max;
  result.entry_face=face_min;
  result.exit_face=face_max;
  return result;
}

/* Return the first forward exit from a point inside/on the box.  Local face
 * coordinates use the historical AMPS ordering and lie in [0,1] to tolerance. */
inline RayBoxStatus GetRayBoxExit(
    const double* box_min, const double* box_max, const double* origin,
    const double* direction, double tolerance, double* exit_point,
    double* face_local, int& exit_face) {
  exit_face=-1;
  const RayBoxInterval interval=ClipRayToAxisAlignedBox(
      box_min,box_max,origin,direction,tolerance,true);
  if (interval.status!=RayBoxStatus::Success) return interval.status;
  for (int axis=0; axis<3; ++axis) {
    if (origin[axis] < box_min[axis]-tolerance ||
        origin[axis] > box_max[axis]+tolerance)
      return RayBoxStatus::NoIntersection;
  }
  if (!(interval.exit_parameter>=-tolerance) || interval.exit_face<0)
    return RayBoxStatus::NoIntersection;
  exit_face=interval.exit_face;
  for (int axis=0; axis<3; ++axis)
    exit_point[axis]=origin[axis]+interval.exit_parameter*direction[axis];

  const int normal_axis=exit_face/2;
  int tangent[2], n=0;
  for (int axis=0; axis<3; ++axis)
    if (axis!=normal_axis) tangent[n++]=axis;
  for (int i=0;i<2;++i)
    face_local[i]=(exit_point[tangent[i]]-box_min[tangent[i]])/
                  (box_max[tangent[i]]-box_min[tangent[i]]);
  return RayBoxStatus::Success;
}

}  // namespace GeometryKernels
}  // namespace PIC

#endif
