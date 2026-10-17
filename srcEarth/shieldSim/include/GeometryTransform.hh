#ifndef SHIELDSIM_GEOMETRY_TRANSFORM_HH
#define SHIELDSIM_GEOMETRY_TRANSFORM_HH

/*
 * GeometryTransform.hh
 * --------------------
 * Small, dependency-free rigid-transform mathematics used by the case parser,
 * geometry hashes, source generator, and analytic validation tests.
 *
 * Keeping this representation independent of Geant4 is intentional.  Input
 * validation and canonical hashing must not depend on Geant4 navigation state,
 * and exact transform tests should run without loading Geant4 libraries.
 * Quaternions represent an active rotation from an object's local right-handed
 * frame into the world frame.  A point is transformed as R(q)*p + translation.
 */

#include <array>

namespace shield_sim {

struct Vector3d {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct Quaterniond {
  // Scalar-first convention.  q and -q encode the same rotation; Canonical()
  // removes that ambiguity before serialization and content hashing.
  double w = 1.0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  double Norm() const;
  Quaterniond Canonical(double tolerance = 1.0e-12) const;
};

struct RigidTransform {
  Vector3d translationMm;
  Quaterniond rotation;

  Vector3d Rotate(const Vector3d& local) const;
  Vector3d InverseRotate(const Vector3d& worldVector) const;
  Vector3d TransformPoint(const Vector3d& local) const;
  Vector3d InverseTransformPoint(const Vector3d& world) const;

  // Columns are the local x/y/z unit axes expressed in world coordinates.
  std::array<Vector3d,3> Axes() const;
};

double Dot(const Vector3d& a, const Vector3d& b);
Vector3d operator+(const Vector3d& a, const Vector3d& b);
Vector3d operator-(const Vector3d& a, const Vector3d& b);
Vector3d operator*(double scale, const Vector3d& vector);

} // namespace shield_sim

#endif
