#ifndef _SRC_GEOSPACE_GPU_ENGINE_MATH_H_
#define _SRC_GEOSPACE_GPU_ENGINE_MATH_H_

//======================================================================================
// EngineMath.h -- 3-vector helpers of the batched trajectory engine.
//
// Every operation reproduces the operand order of the legacy helpers in
// gridless/GridlessParticleMovers.h (add/sub/mul/dot/cross/norm/unit).  The engine relies
// on this to be bitwise identical to the legacy per-trajectory tracer in double/SI.
//======================================================================================

#include "EngineConfig.h"

namespace Geospace {
namespace GPU {

template <class R> struct Vec3 { R x,y,z; };

template <class R> EG_HD Vec3<R> vadd(const Vec3<R>& a,const Vec3<R>& b) {
  return {a.x+b.x,a.y+b.y,a.z+b.z};
}
template <class R> EG_HD Vec3<R> vsub(const Vec3<R>& a,const Vec3<R>& b) {
  return {a.x-b.x,a.y-b.y,a.z-b.z};
}
template <class R> EG_HD Vec3<R> vmul(R s,const Vec3<R>& a) {
  return {s*a.x,s*a.y,s*a.z};
}
template <class R> EG_HD R vdot(const Vec3<R>& a,const Vec3<R>& b) {
  return a.x*b.x+a.y*b.y+a.z*b.z;
}
template <class R> EG_HD Vec3<R> vcross(const Vec3<R>& a,const Vec3<R>& b) {
  return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};
}
template <class R> EG_HD R vnorm(const Vec3<R>& a) { return std::sqrt(vdot(a,a)); }
template <class R> EG_HD Vec3<R> vunit(const Vec3<R>& a) {
  const R n=vnorm(a);
  return (n>R(0)) ? vmul(R(1)/n,a) : Vec3<R>{R(0),R(0),R(0)};
}
template <class R> EG_HD bool vfinite(const Vec3<R>& a) {
  return efinite(a.x) && efinite(a.y) && efinite(a.z);
}

} // namespace GPU
} // namespace Geospace

#endif
