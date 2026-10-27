//======================================================================================
// test_mover_equivalence.cpp -- bitwise equivalence of the engine movers
// (gpu/EngineMovers.h) with the legacy StepParticleChecked()
// (gridless/GridlessParticleMovers.cpp) for every mover, including the
// stateful HYBRID RK4/GC4 selector, in a tilted dipole field.
//
// For each mover and seed a trajectory of up to 400 checked steps is integrated by both
// implementations from the same state with the same dt sequence; after every step the
// position, momentum and the "step accepted" flag must be bitwise identical.
//
// Build: see run_unit_tests.sh in this directory.
//======================================================================================

#include "GridlessParticleMovers.h"
#include "gpu/EngineMovers.h"
#include "gpu/EngineField.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>

static bool Same(double a,double b) { return std::memcmp(&a,&b,sizeof(double))==0; }
// Both NaN (payload/sign may differ): the tracer classifies either as NUMERICAL_FAILURE.
static bool SameOrBothNaN(double a,double b,int& nanPairs) {
  if (Same(a,b)) return true;
  if (std::isnan(a) && std::isnan(b)) { ++nanPairs; return true; }
  return false;
}

// Field shared by both implementations: identical arithmetic, so differences can only
// come from the movers.
static Geospace::GPU::EgDipole gDip;

class RefField : public IGridlessFieldEvaluator {
 public:
  void GetB_T(const V3& x,V3& B) const override {
    double xs[3]={x.x,x.y,x.z},b[3];
    Geospace::GPU::DipoleFieldSI(gDip,xs,b,false);
    B.x=b[0]; B.y=b[1]; B.z=b[2];
  }
};

struct EngField {
  void GetB(const Geospace::GPU::Vec3<double>& x,Geospace::GPU::Vec3<double>& B) {
    double xs[3]={x.x,x.y,x.z},b[3];
    Geospace::GPU::DipoleFieldSI(gDip,xs,b,false);
    B.x=b[0]; B.y=b[1]; B.z=b[2];
  }
};

int main() {
  const double Re=6371.0e3;
  const double M=1.0*Geospace::GPU::EgDipoleConstants::kM_E_Am2;
  const double tilt=11.0*Geospace::GPU::EgDipoleConstants::kPI/180.0;
  gDip.mx=M*std::sin(tilt); gDip.my=M*0.0; gDip.mz=M*std::cos(tilt);
  RefField ref;
  EngField eng;
  typedef Geospace::GPU::Movers<double,Geospace::GPU::UnitsSI> EM;
  const double q=1.602176565E-19,m0=1.007276466621*1.66053906660E-27,c=299792458.0;

  int failures=0,steps=0,rejected=0,gcSteps=0,nanPairs=0;
  for (int mover=0;mover<=8;++mover) {
    for (unsigned seed=1;seed<=40;++seed) {
      std::mt19937 rng(seed*97u+mover);
      std::uniform_real_distribution<double> u(0.0,1.0);
      const double r=(1.3+8.0*u(rng))*Re,th=std::acos(2*u(rng)-1),ph=6.283185307179586*u(rng);
      V3 xr{r*std::sin(th)*std::cos(ph),r*std::sin(th)*std::sin(ph),r*std::cos(th)};
      // HYBRID needs adiabatic orbits (rho/L<0.03) to enter its GC branch: low rigidity.
      const double R_GV=(mover==8) ? 0.0005+0.01*u(rng) : 0.05+20.0*u(rng)*u(rng);
      const double pm=(R_GV*1.0e9*q)/c;
      const double a=(mover==8) ? 0.6*u(rng) : std::acos(2*u(rng)-1);
      const double b=6.283185307179586*u(rng);
      V3 pr{pm*std::sin(a)*std::cos(b),pm*std::sin(a)*std::sin(b),pm*std::cos(a)};
      const double rInner=Re;
      Geospace::GPU::Vec3<double> xe{xr.x,xr.y,xr.z},pe{pr.x,pr.y,pr.z};
      Geospace::GPU::EgHybrid<double> hyb;
      Geospace::GPU::ResetHybrid(hyb,xe,rInner);
      ResetHybridTrajectoryContext(xr,rInner);
      const MoverType mt=static_cast<MoverType>(mover);
      for (int k=0;k<400;++k) {
        // dt from the local gyrofrequency (same formula in both), occasionally large to
        // exercise inner-sphere segment checks.
        V3 B; ref.GetB_T(xr,B);
        const double Bm=std::sqrt(B.x*B.x+B.y*B.y+B.z*B.z);
        const double gam=std::sqrt(1.0+(pm*pm)/((m0*c)*(m0*c)));
        double dt=0.15*gam*m0/(q*Bm)*(0.2+1.8*u(rng));
        if (k%37==5) dt*=50.0;
        bool useGCr=false,useGCe=false;
        if (mover==8) {
          useGCr=HybridPrepareStepUseGuidingCenter(xr,pr,q,ref);
          useGCe=EM::HybridPrepare(hyb,xe,pe,q,eng);
          if (useGCr!=useGCe) { ++failures; std::printf("hybrid decision mismatch m=%d s=%u k=%d\n",mover,seed,k); break; }
          if (useGCr) ++gcSteps;
        }
        const bool okr=StepParticleChecked(mt,xr,pr,q,m0,dt,ref,rInner);
        const bool oke=EM::StepChecked(mover,hyb,xe,pe,q,m0,dt,eng,rInner);
        ++steps;
        if (okr!=oke || !SameOrBothNaN(xr.x,xe.x,nanPairs) || !SameOrBothNaN(xr.y,xe.y,nanPairs) ||
            !SameOrBothNaN(xr.z,xe.z,nanPairs) || !SameOrBothNaN(pr.x,pe.x,nanPairs) ||
            !SameOrBothNaN(pr.y,pe.y,nanPairs) || !SameOrBothNaN(pr.z,pe.z,nanPairs)) {
          if (failures<10)
            std::printf("mismatch mover=%d seed=%u step=%d ok=%d/%d x=%.17g/%.17g\n",
                        mover,seed,k,okr,oke,xr.x,xe.x);
          ++failures;
          break;
        }
        if (!okr) { ++rejected; break; }
        // Stop at the first non-finite state: both tracers terminate there with
        // NUMERICAL_FAILURE.  Integrating further from NaN would only compare NaN
        // payload/sign bits, which carry no information.
        const double rr=std::sqrt(xr.x*xr.x+xr.y*xr.y+xr.z*xr.z);
        if (rr>40*Re || !std::isfinite(rr)) break;
      }
    }
  }
  std::printf("mover equivalence: %d checked steps (9 movers x 40 trajectories), %d inner-sphere "
              "rejections, %d HYBRID GC steps, %d NaN-in-both components, %d mismatches\n",
              steps,rejected,gcSteps,nanPairs,failures);
  return failures==0 ? 0 : 1;
}
