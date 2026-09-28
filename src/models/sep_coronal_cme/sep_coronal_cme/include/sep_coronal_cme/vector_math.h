#ifndef SEP_CORONAL_CME_VECTOR_MATH_H
#define SEP_CORONAL_CME_VECTOR_MATH_H

#include <cmath>

namespace SEP { namespace CoronalCME {

// A deliberately small neutral vector type keeps all shared kernels free of
// an AMPS mesh/vector ABI. Coordinates and components are SI unless the API
// explicitly describes a unit vector.
struct Vec3 {
  double x = 0.0, y = 0.0, z = 0.0;
};

inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
inline Vec3 operator*(double s, Vec3 a) { return {s*a.x,s*a.y,s*a.z}; }
inline Vec3 operator*(Vec3 a, double s) { return s*a; }
inline Vec3 operator/(Vec3 a, double s) { return {a.x/s,a.y/s,a.z/s}; }
inline double Dot(Vec3 a, Vec3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline Vec3 Cross(Vec3 a, Vec3 b) {
  return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};
}
inline double Norm(Vec3 a) { return std::sqrt(Dot(a,a)); }
inline Vec3 Unit(Vec3 a) { const double n=Norm(a); return n>0.0 ? a/n : Vec3{}; }

} }  // namespace SEP::CoronalCME
#endif
