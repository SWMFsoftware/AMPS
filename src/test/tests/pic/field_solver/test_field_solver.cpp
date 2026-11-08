#include "test_harness.h"
#include "pic_field_solver_kernels.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace K=PIC::FieldSolverKernels;
static constexpr double pi=3.141592653589793238462643383279502884;

static double poisson_error(int intervals,double* final_residual=nullptr) {
  const double dx=1.0/intervals;
  std::vector<double> rhs(intervals-1),solution;
  for (int i=1;i<intervals;++i) rhs[i-1]=pi*pi*std::sin(pi*i*dx);
  const K::IterativeResult result=K::SolvePoissonDirichlet1D(
      rhs,dx,1.0e-13,4*intervals,solution);
  AMPS_EXPECT_TRUE(result.status==K::Status::Success);
  AMPS_EXPECT_TRUE(result.residual_history.size()>=2);
  if (final_residual) *final_residual=result.residual_history.back();
  double sum=0.0;
  for (int i=1;i<intervals;++i) {
    const double difference=solution[i-1]-std::sin(pi*i*dx);
    sum+=difference*difference;
  }
  return std::sqrt(sum/(intervals-1));
}

AMPS_TEST_CASE("PIC-FIELD-POISSON-RESIDUAL-001",
               "CG solves the independently assembled Dirichlet Poisson system") {
  double residual=0.0;
  const double error=poisson_error(64,&residual);
  AMPS_EXPECT_TRUE(error<2.1e-4);
  AMPS_EXPECT_TRUE(residual<1.0e-9);
}

AMPS_TEST_CASE("PIC-FIELD-POISSON-CONVERGENCE-001",
               "manufactured sine solution converges at the advertised second order") {
  const double e32=poisson_error(32),e64=poisson_error(64),e128=poisson_error(128);
  const double p1=std::log(e32/e64)/std::log(2.0);
  const double p2=std::log(e64/e128)/std::log(2.0);
  // The analytic sine is not a discrete eigenfunction with the same
  // eigenvalue; the leading truncation term is O(dx^2), with small O(dx^4)
  // corrections on the coarser pair.  A two-percent band is deliberately
  // tighter than confusing first- and second-order behavior while not using a
  // production-generated expected value.
  AMPS_EXPECT_TRUE(p1>1.98 && p1<2.02);
  AMPS_EXPECT_TRUE(p2>1.98 && p2<2.02);
}

AMPS_TEST_CASE("PIC-FIELD-WAVE-DISPERSION-001",
               "periodic Yee eigenmode follows its independent discrete dispersion law") {
  const int n=96,mode=5,steps=137;
  const double courant=0.7,kdx=2.0*pi*mode/n;
  const double omega_dt=2.0*std::asin(courant*std::sin(0.5*kdx));
  std::vector<double> electric(n),magnetic(n);
  for (int i=0;i<n;++i) {
    electric[i]=std::sin(kdx*i);
    magnetic[i]=std::sin(kdx*(i+0.5)+0.5*omega_dt); // B at n=-1/2
  }
  for (int step=0;step<steps;++step)
    AMPS_EXPECT_TRUE(K::AdvanceVacuumYee1D(electric,magnetic,courant)==K::Status::Success);
  double error2=0.0,reference2=0.0;
  for (int i=0;i<n;++i) {
    const double exact=std::sin(kdx*i-steps*omega_dt);
    error2+=(electric[i]-exact)*(electric[i]-exact);
    reference2+=exact*exact;
  }
  AMPS_EXPECT_TRUE(std::sqrt(error2/reference2)<2.0e-13);
}

AMPS_TEST_CASE("PIC-FIELD-WAVE-ENERGY-001",
               "Yee staggered invariant and particle-plus-field ledger close independently") {
  const int n=80,mode=3,steps=200;
  const double courant=0.55,kdx=2.0*pi*mode/n;
  const double omega_dt=2.0*std::asin(courant*std::sin(0.5*kdx));
  std::vector<double> electric(n),magnetic_old(n),magnetic_new;
  for (int i=0;i<n;++i) {
    electric[i]=std::sin(kdx*i);
    magnetic_old[i]=std::sin(kdx*(i+0.5)+0.5*omega_dt);
  }
  magnetic_new=magnetic_old;
  AMPS_EXPECT_TRUE(K::AdvanceVacuumYee1D(electric,magnetic_new,courant)==K::Status::Success);
  double invariant0=K::YeeInvariant(electric,magnetic_new,magnetic_old);
  for (int step=1;step<steps;++step) {
    magnetic_old=magnetic_new;
    AMPS_EXPECT_TRUE(K::AdvanceVacuumYee1D(electric,magnetic_new,courant)==K::Status::Success);
    AMPS_EXPECT_NEAR(K::YeeInvariant(electric,magnetic_new,magnetic_old),
                     invariant0,3.0e-12);
  }
  const std::vector<double> mass={2.0,4.0},speed2={9.0,1.0};
  const double before=K::ParticleFieldEnergy(10.0,7.0,mass,speed2);
  const std::vector<double> speed2_after={10.0,1.0}; // +1 J particle energy
  const double after=K::ParticleFieldEnergy(9.0,7.0,mass,speed2_after);
  AMPS_EXPECT_NEAR(after,before,0.0);
}

AMPS_TEST_CASE("PIC-FIELD-SOLVER-INVALID-001",
               "field kernels reject unstable, nonfinite, and impossible requests") {
  std::vector<double> solution;
  const std::vector<double> rhs={1.0,2.0,3.0};
  AMPS_EXPECT_TRUE(K::SolvePoissonDirichlet1D(rhs,0.0,1.0e-8,10,solution).status==
                   K::Status::InvalidArgument);
  std::vector<double> e(8,0.0),b(8,0.0);
  AMPS_EXPECT_TRUE(K::AdvanceVacuumYee1D(e,b,1.01)==K::Status::InvalidArgument);
  e[2]=std::numeric_limits<double>::quiet_NaN();
  AMPS_EXPECT_TRUE(K::AdvanceVacuumYee1D(e,b,0.5)==K::Status::InvalidArgument);
}
