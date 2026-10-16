// Check frozen T05 + IGRF field evaluation independently of trajectory chaos.
// Link against the exact production Fortran objects after the mandated AMPS build:
// c++ -std=c++17 -O2 -pthread check_field_repeatability.cpp \
//   build/interface/T05.o build/interface/geopack2008.o -lgfortran -o CHECK
// This is a concurrency diagnostic, not an observational/field-accuracy test.
#include <array>
#include <atomic>
#include <cmath>
#include <cstring>
#include <iostream>
#include <thread>
#include <vector>

extern "C" {
void recalc_08_(int*,int*,int*,int*,int*,double*,double*,double*);
void igrf_gsw_08_(double*,double*,double*,double*,double*,double*);
void t04_s_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
}

using Vector=std::array<double,3>;
Vector evaluate(int index) {
  int option=0;
  double parmod[10]={2.0,-20.0,1.0,-3.0,0.2,0.2,0.2,0.2,0.2,0.2};
  double tilt=0.15;
  Vector point={{2.0+std::fmod(0.031*index,8.0),
                 -4.0+std::fmod(0.013*index,8.0),
                 -3.0+std::fmod(0.019*index,6.0)}};
  Vector external,internal;
  t04_s_(&option,parmod,&tilt,&point[0],&point[1],&point[2],
         &external[0],&external[1],&external[2]);
  igrf_gsw_08_(&point[0],&point[1],&point[2],&internal[0],&internal[1],&internal[2]);
  for (int k=0;k<3;++k) external[k]+=internal[k];
  return external;
}

int main() {
  // Use an epoch inside the frozen SA04 interval.  Geopack receives day-of-year,
  // so 249 is 6 September in non-leap year 2017; this is not a calendar day field.
  int year=2017,day=249,hour=11,minute=0,second=0;
  double vx=-400,vy=0,vz=0;
  recalc_08_(&year,&day,&hour,&minute,&second,&vx,&vy,&vz);
  const int count=20000,workers=3,repeats=3;
  std::vector<Vector> reference(count);
  for (int i=0;i<count;++i) reference[i]=evaluate(i);
  std::atomic<int> differences(0);
  std::vector<std::thread> threads;
  for (int worker=0;worker<workers;++worker) threads.emplace_back([&,worker]() {
    for (int pass=0;pass<repeats;++pass) for (int i=worker;i<count;i+=workers) {
      const Vector value=evaluate(i);
      if (std::memcmp(value.data(),reference[i].data(),sizeof(Vector))!=0 ||
          !std::isfinite(value[0]) || !std::isfinite(value[1]) || !std::isfinite(value[2]))
        ++differences;
    }
  });
  for (auto& thread:threads) thread.join();
  std::cout << "Frozen T05/IGRF serial-versus-threaded exact comparisons="
            << count*repeats << " differences=" << differences.load() << '\n';
  return differences.load()==0 ? 0 : 1;
}
