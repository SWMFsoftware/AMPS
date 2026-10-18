//3D vector algebra used throughout the code
#ifndef VEC3_H
#define VEC3_H

#include <cmath>

struct Vec3 {
  double x,y,z;

  Vec3() : x(0.0),y(0.0),z(0.0) {}
  Vec3(double a,double b,double c) : x(a),y(b),z(c) {}

  double operator[](int i) const {return (i==0) ? x : ((i==1) ? y : z);}
  double& operator[](int i) {return (i==0) ? x : ((i==1) ? y : z);}

  Vec3 operator+(const Vec3& v) const {return Vec3(x+v.x,y+v.y,z+v.z);}
  Vec3 operator-(const Vec3& v) const {return Vec3(x-v.x,y-v.y,z-v.z);}
  Vec3 operator-() const {return Vec3(-x,-y,-z);}
  Vec3 operator*(double a) const {return Vec3(a*x,a*y,a*z);}
  Vec3 operator/(double a) const {return Vec3(x/a,y/a,z/a);}
  Vec3& operator+=(const Vec3& v) {x+=v.x,y+=v.y,z+=v.z; return *this;}
  Vec3& operator-=(const Vec3& v) {x-=v.x,y-=v.y,z-=v.z; return *this;}
  Vec3& operator*=(double a) {x*=a,y*=a,z*=a; return *this;}
};

inline Vec3 operator*(double a,const Vec3& v) {return v*a;}
inline double Dot(const Vec3& a,const Vec3& b) {return a.x*b.x+a.y*b.y+a.z*b.z;}
inline Vec3 Cross(const Vec3& a,const Vec3& b) {return Vec3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);}
inline double Length(const Vec3& a) {return std::sqrt(Dot(a,a));}
inline Vec3 Normalized(const Vec3& a) {return a/Length(a);}

inline Vec3 Min(const Vec3& a,const Vec3& b) {return Vec3(std::fmin(a.x,b.x),std::fmin(a.y,b.y),std::fmin(a.z,b.z));}
inline Vec3 Max(const Vec3& a,const Vec3& b) {return Vec3(std::fmax(a.x,b.x),std::fmax(a.y,b.y),std::fmax(a.z,b.z));}

//build an orthonormal frame (e0,e1) perpendicular to the unit vector e2
inline void OrthonormalFrame(const Vec3& e2,Vec3& e0,Vec3& e1) {
  Vec3 a=(std::fabs(e2.x)<0.9) ? Vec3(1.0,0.0,0.0) : Vec3(0.0,1.0,0.0);

  e0=Normalized(Cross(e2,a));
  e1=Cross(e2,e0);
}

#endif
