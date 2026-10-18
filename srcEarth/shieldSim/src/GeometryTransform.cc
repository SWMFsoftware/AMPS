#include "GeometryTransform.hh"

#include <cmath>
#include <stdexcept>

namespace shield_sim {

double Quaterniond::Norm() const {
  return std::sqrt(w*w + x*x + y*y + z*z);
}

Quaterniond Quaterniond::Canonical(double tolerance) const {
  const double norm=Norm();
  if(!std::isfinite(norm) || std::abs(norm-1.0)>tolerance)
    throw std::invalid_argument("rotation quaternion must be finite and unit-normalized");
  Quaterniond q{w/norm,x/norm,y/norm,z/norm};
  // The first materially nonzero component is made positive.  Thus q and -q
  // have identical canonical bytes and cannot create two cache keys for the
  // same physical placement.
  const double components[4]={q.w,q.x,q.y,q.z};
  for(double component:components){
    if(std::abs(component)<=tolerance) continue;
    if(component<0.0){ q.w=-q.w; q.x=-q.x; q.y=-q.y; q.z=-q.z; }
    break;
  }
  return q;
}

Vector3d RigidTransform::Rotate(const Vector3d& v) const {
  const Quaterniond q=rotation.Canonical();
  // Expanded q*v*q^-1 avoids a matrix dependency and is stable for unit q.
  const double xx=q.x*q.x, yy=q.y*q.y, zz=q.z*q.z;
  const double xy=q.x*q.y, xz=q.x*q.z, yz=q.y*q.z;
  const double wx=q.w*q.x, wy=q.w*q.y, wz=q.w*q.z;
  return {
    (1.0-2.0*(yy+zz))*v.x + 2.0*(xy-wz)*v.y + 2.0*(xz+wy)*v.z,
    2.0*(xy+wz)*v.x + (1.0-2.0*(xx+zz))*v.y + 2.0*(yz-wx)*v.z,
    2.0*(xz-wy)*v.x + 2.0*(yz+wx)*v.y + (1.0-2.0*(xx+yy))*v.z
  };
}

Vector3d RigidTransform::InverseRotate(const Vector3d& v) const {
  const Quaterniond q=rotation.Canonical();
  RigidTransform inverse;
  inverse.rotation={q.w,-q.x,-q.y,-q.z};
  return inverse.Rotate(v);
}

Vector3d RigidTransform::TransformPoint(const Vector3d& local) const {
  return Rotate(local)+translationMm;
}

Vector3d RigidTransform::InverseTransformPoint(const Vector3d& world) const {
  return InverseRotate(world-translationMm);
}

std::array<Vector3d,3> RigidTransform::Axes() const {
  return {Rotate({1.0,0.0,0.0}),Rotate({0.0,1.0,0.0}),Rotate({0.0,0.0,1.0})};
}

RigidTransform ComposeTransforms(const RigidTransform& parent,
                                 const RigidTransform& child) {
  const Quaterniond a=parent.rotation.Canonical();
  const Quaterniond b=child.rotation.Canonical();
  // Hamilton product for active rotations: R(a*b) == R(a)R(b).
  Quaterniond product{
    a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z,
    a.w*b.x+a.x*b.w+a.y*b.z-a.z*b.y,
    a.w*b.y-a.x*b.z+a.y*b.w+a.z*b.x,
    a.w*b.z+a.x*b.y-a.y*b.x+a.z*b.w};
  return {parent.TransformPoint(child.translationMm),product.Canonical()};
}

double Dot(const Vector3d& a,const Vector3d& b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
Vector3d operator+(const Vector3d& a,const Vector3d& b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
Vector3d operator-(const Vector3d& a,const Vector3d& b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
Vector3d operator*(double s,const Vector3d& v){ return {s*v.x,s*v.y,s*v.z}; }

} // namespace shield_sim
