// Bitwise check of the remainder()-based fmod used by the OpenACC build
// (Geospace::GPU::efmod in app/gpu/EngineConfig.h) against std::fmod on the host.
//
//   g++ -O2 -std=c++17 -ffp-contract=off efmod_exactness.cpp -o efmod_exactness
//   ./efmod_exactness [samples]
//
// The engine uses fmod(|driftAccumulatedAngle|, 2*pi); the test covers that use and
// generic signed arguments, zeros, infinities, NaN and exact halfway cases.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>

template <class T> T efmod_remainder(T x,T y) {   // the EG_HAVE_OPENACC branch
  const T r=std::remainder(x,y);
  if (x>=T(0)) return (r<T(0)) ? r+std::fabs(y) : r;
  return (r>T(0)) ? r-std::fabs(y) : r;
}

template <class T> bool same(T a,T b) {
  if (std::isnan(a) && std::isnan(b)) return true;
  return std::memcmp(&a,&b,sizeof(T))==0;
}

template <class T> long long check(long long n,unsigned seed) {
  std::mt19937_64 g(seed);
  const T twoPi=T(2.0*3.141592653589793238462643383279502884);
  long long bad=0;
  auto test=[&](T x,T y) {
    const T a=std::fmod(x,y), b=efmod_remainder(x,y);
    if (!same(a,b)) {
      if (++bad<=10) std::printf("  MISMATCH x=%.17g y=%.17g fmod=%.17g efmod=%.17g\n",
                                 double(x),double(y),double(a),double(b));
    }
  };
  // engine use: |accumulated drift angle| in [0, 1e6] rad, y = 2*pi
  std::uniform_real_distribution<double> ang(0.0,1.0e6);
  for (long long i=0;i<n;++i) test(T(ang(g)),twoPi);
  // exact multiples / halfway points of 2*pi and neighbours
  for (int k=0;k<200000;++k) {
    const T m=T(k)*twoPi, h=(T(k)+T(0.5))*twoPi;
    test(m,twoPi); test(std::nextafter(m,T(0)),twoPi); test(std::nextafter(m,T(1e30)),twoPi);
    test(h,twoPi); test(std::nextafter(h,T(0)),twoPi); test(std::nextafter(h,T(1e30)),twoPi);
  }
  // random bit patterns (all exponents, both signs)
  for (long long i=0;i<n;++i) {
    T x,y;
    if (sizeof(T)==8) { uint64_t a=g(),b=g(); std::memcpy(&x,&a,8); std::memcpy(&y,&b,8); }
    else { uint32_t a=uint32_t(g()),b=uint32_t(g()); std::memcpy(&x,&a,4); std::memcpy(&y,&b,4); }
    test(x,y);
  }
  // special values
  const T inf=std::numeric_limits<T>::infinity(), nan=std::numeric_limits<T>::quiet_NaN();
  const T sp[]={T(0),-T(0),T(1),-T(1),twoPi,-twoPi,inf,-inf,nan,
                std::numeric_limits<T>::denorm_min(),std::numeric_limits<T>::min(),
                std::numeric_limits<T>::max()};
  for (T x:sp) for (T y:sp) test(x,y);
  return bad;
}

int main(int argc,char** argv) {
  const long long n=(argc>1) ? std::atoll(argv[1]) : 20000000LL;
  const long long bd=check<double>(n,12345u), bf=check<float>(n,54321u);
  std::printf("efmod (remainder-based) vs std::fmod: %lld samples per type + edge cases: "
              "double mismatches=%lld float mismatches=%lld\n",n,bd,bf);
  return (bd==0 && bf==0) ? 0 : 1;
}
