#include "test_framework.hpp"
#include "../reference/shk05_oblique_v1.hpp"
#include "../reference/shk12_near_mach_v1.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_shock.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {
using Vec3=swcme::shock::Vec3;
using Primitive=swcme::shock::PrimitiveState;

constexpr double GAMMA=5.0/3.0;
constexpr double MP=swcme::constants::PROTON_MASS_KG;
constexpr double KB=swcme::constants::BOLTZMANN_J_K;
constexpr double MU0=swcme::constants::VACUUM_PERMEABILITY_N_A2;

inline double dot(const Vec3& a,const Vec3& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline Vec3 sub(const Vec3& a,const Vec3& b){return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};}
inline Vec3 scale(const Vec3& a,double s){return {{s*a[0],s*a[1],s*a[2]}};}
inline double norm(const Vec3& a){return std::sqrt(std::max(0.0,dot(a,a)));}

void print_check(const std::string& label,double actual,double expected,double tol,bool pass){
  std::cout<<"  "<<std::left<<std::setw(38)<<label<<std::right
           <<" actual="<<std::scientific<<std::setprecision(12)<<actual
           <<" expected="<<expected<<" tol="<<tol<<' '
           <<(pass?"PASS":"FAIL")<<'\n';
}

void expect_rel(swcme_test::Context& c,const std::string& label,
                double actual,double expected,double rel_tol,double floor=1e-30){
  const double scalev=std::max({std::abs(actual),std::abs(expected),floor});
  const double tol=rel_tol*scalev;
  const bool pass=std::isfinite(actual)&&std::isfinite(expected)&&std::abs(actual-expected)<=tol;
  print_check(label,actual,expected,tol,pass); c.expect_true(pass,label);
}

void expect_abs(swcme_test::Context& c,const std::string& label,
                double actual,double expected,double tol){
  const bool pass=std::isfinite(actual)&&std::isfinite(expected)&&std::abs(actual-expected)<=tol;
  print_check(label,actual,expected,tol,pass); c.expect_true(pass,label);
}

Primitive fixture(double theta_deg=40.0,double n_cm3=5.0,double B_nT=5.0,double T=1.2e5){
  Primitive s;
  const double th=theta_deg*swcme::constants::PI/180.0;
  const double n=n_cm3*1e6;
  s.rho_kg_m3=n*MP;
  s.pressure_Pa=n*KB*T;
  s.velocity_m_s={{4.0e5,0.0,0.0}};
  s.magnetic_T={{B_nT*1e-9*std::cos(th),B_nT*1e-9*std::sin(th),0.0}};
  return s;
}

// Independent fast-mode speed used only to construct threshold fixtures.  It
// intentionally does not call the production helper under test.
double reference_fast_speed(const Primitive& s,const Vec3& n,double gamma=GAMMA){
  const double B2=dot(s.magnetic_T,s.magnetic_T);
  const double vA2=B2/(MU0*s.rho_kg_m3);
  const double cs2=gamma*s.pressure_Pa/s.rho_kg_m3;
  const double cosBn=(B2>0.0)?dot(s.magnetic_T,n)/std::sqrt(B2):0.0;
  const double a=vA2+cs2;
  const double disc=std::max(0.0,a*a-4.0*vA2*cs2*cosBn*cosBn);
  return std::sqrt(0.5*(a+std::sqrt(disc)));
}

// Independent perpendicular-shock reference.  For B_n=0 and zero upstream
// tangential flow, ideal-MHD induction gives B2_t=r B1_t and mass continuity
// gives u2_n=u1_n/r.  Normal momentum then determines p2; the scalar energy
// equation below is solved by bisection.  Keeping this test-specific reduction
// independent of swcme_shock.hpp makes SHK04 a real validation rather than a
// second call to the production solver.
double perpendicular_reference_compression(const Primitive& up,double Vsh){
  const double u=Vsh-up.velocity_m_s[0]; // positive inflow magnitude
  const double B=std::abs(up.magnetic_T[1]);
  const double rmax=(GAMMA+1.0)/(GAMMA-1.0);
  auto q=[&](double r){
    const double p2=up.rho_kg_m3*u*u+up.pressure_Pa+B*B/(2*MU0)
                   -up.rho_kg_m3*u*u/r-r*r*B*B/(2*MU0);
    const double E1=0.5*up.rho_kg_m3*u*u
                   +GAMMA/(GAMMA-1.0)*up.pressure_Pa+B*B/MU0;
    const double E2=(1.0/r)*(0.5*up.rho_kg_m3*u*u/r
                   +GAMMA/(GAMMA-1.0)*p2+r*r*B*B/MU0);
    return (E2-E1)/(r-1.0);
  };
  double lo=1.0+1e-8, flo=q(lo), hi=rmax-1e-8, fhi=q(hi);
  if ((flo<0)==(fhi<0)) return 0.0;
  for(int i=0;i<160;++i){
    const double mid=0.5*(lo+hi), fm=q(mid);
    if ((flo<0)!=(fm<0)){hi=mid;fhi=fm;} else {lo=mid;flo=fm;}
  }
  (void)fhi;
  return 0.5*(lo+hi);
}

enum class StressStratum { NoShock, WeakLimit, Resolved, Conditioned };

struct StressCase {
  std::size_t index=0;
  StressStratum stratum=StressStratum::Resolved;
  Primitive upstream;
  Vec3 normal{{1.0,0.0,0.0}};
  double shock_speed_m_s=0.0;
  double gamma=GAMMA;
  double requested_fast_mach=0.0;
  swcme::shock::JumpResult result;
};

// Small fixed algorithm used instead of the standard-library distributions.
// The latter are permitted to map engine bits differently between library
// implementations, which would make an alleged fixed-seed campaign vary by
// platform.  SplitMix64 and the explicit 53-bit conversion below make every
// SHK15 primitive reproducible from the seed printed in the README.
class ShockStressRandom {
public:
  explicit ShockStressRandom(std::uint64_t seed): state_(seed) {}
  std::uint64_t next_u64(){
    state_+=0x9e3779b97f4a7c15ULL;
    std::uint64_t z=state_;
    z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
    z=(z^(z>>27))*0x94d049bb133111ebULL;
    return z^(z>>31);
  }
  double unit(){
    return static_cast<double>(next_u64()>>11)*0x1.0p-53;
  }
  double uniform(double lower,double upper){
    return lower+(upper-lower)*unit();
  }
  double log_uniform(double lower,double upper){
    return std::exp(uniform(std::log(lower),std::log(upper)));
  }
private:
  std::uint64_t state_;
};

Vec3 add(const Vec3& a,const Vec3& b){
  return {{a[0]+b[0],a[1]+b[1],a[2]+b[2]}};
}

Vec3 cross(const Vec3& a,const Vec3& b){
  return {{a[1]*b[2]-a[2]*b[1],
           a[2]*b[0]-a[0]*b[2],
           a[0]*b[1]-a[1]*b[0]}};
}

Vec3 unit_vector(const Vec3& value){
  const double magnitude=norm(value);
  return magnitude>0.0?scale(value,1.0/magnitude):Vec3{{1.0,0.0,0.0}};
}

// Form a deterministic orthonormal frame about an arbitrary random normal.
// Choosing the least-aligned Cartesian axis avoids a nearly zero cross product
// and also supports the arbitrary-orientation checks reused by SHK07/SHK08.
std::array<Vec3,2> tangential_basis(const Vec3& normal){
  Vec3 axis={{1.0,0.0,0.0}};
  if(std::abs(normal[1])<=std::abs(normal[0]) &&
     std::abs(normal[1])<=std::abs(normal[2])) axis={{0.0,1.0,0.0}};
  else if(std::abs(normal[2])<=std::abs(normal[0]) &&
          std::abs(normal[2])<=std::abs(normal[1])) axis={{0.0,0.0,1.0}};
  const Vec3 t1=unit_vector(cross(normal,axis));
  return {{t1,unit_vector(cross(normal,t1))}};
}

Vec3 random_unit_vector(ShockStressRandom& random){
  const double z=random.uniform(-1.0,1.0);
  const double phi=random.uniform(0.0,2.0*swcme::constants::PI);
  const double radius=std::sqrt(std::max(0.0,1.0-z*z));
  return {{radius*std::cos(phi),radius*std::sin(phi),z}};
}

bool finite_primitive(const Primitive& state){
  if(!std::isfinite(state.rho_kg_m3) || !std::isfinite(state.pressure_Pa))
    return false;
  for(double value:state.velocity_m_s) if(!std::isfinite(value)) return false;
  for(double value:state.magnetic_T) if(!std::isfinite(value)) return false;
  return true;
}

bool finite_jump(const swcme::shock::JumpResult& jump){
  const std::array<double,20> scalars={{
      jump.compression,jump.theta_Bn_rad,jump.fast_speed_m_s,jump.fast_mach,
      jump.shock_normal_speed_m_s,jump.upstream_inflow_normal_m_s,
      jump.root_bracket_lower_compression,jump.root_bracket_upper_compression,
      jump.root_bracket_width,jump.minimum_tangential_determinant_relative,
      jump.closest_tangential_determinant_relative,
      jump.selected_tangential_determinant_relative,jump.downstream_fast_mach,
      jump.downstream_normal_alfven_mach,jump.mass_residual,
      jump.normal_B_residual,jump.electric_residual,jump.momentum_residual,
      jump.energy_residual,jump.entropy_ratio}};
  return finite_primitive(jump.upstream) && finite_primitive(jump.downstream) &&
      std::all_of(scalars.begin(),scalars.end(),
                  [](double value){return std::isfinite(value);});
}

const char* stress_stratum_name(StressStratum stratum){
  switch(stratum){
    case StressStratum::NoShock: return "no_shock";
    case StressStratum::WeakLimit: return "weak_limit";
    case StressStratum::Resolved: return "resolved";
    case StressStratum::Conditioned: return "conditioned";
  }
  return "unknown";
}

// Build the high-count campaign once per validation process.  Later SHK06-08
// tests consume these exact serialized results, so the full suite pays the
// 100,000 nonlinear solves only once while each focused invocation remains
// independently executable.  The four strata guarantee substantial no-shock,
// sub-resolution weak, ordinary resolved, and determinant-conditioned cover.
const std::vector<StressCase>& shock_stress_cases(){
  static const std::vector<StressCase> cases=[] {
    constexpr std::size_t count=100000;
    constexpr std::uint64_t seed=0x53484b31355f7631ULL;
    ShockStressRandom random(seed);
    std::vector<StressCase> generated;
    generated.reserve(count);
    for(std::size_t index=0;index<count;++index){
      StressCase item;
      item.index=index;
      item.gamma=random.uniform(1.2,5.0/3.0);
      item.normal=random_unit_vector(random);
      const auto basis=tangential_basis(item.normal);

      // Density, field strength, and beta span four to five orders of
      // magnitude logarithmically.  Pressure is derived from beta, while the
      // independent density range makes the implied temperature another broad
      // log-distributed physical scale rather than a fixed hidden constant.
      item.upstream.rho_kg_m3=random.log_uniform(0.01,100.0)*1.0e6*MP;
      const double field_magnitude=random.log_uniform(0.1,100.0)*1.0e-9;
      double beta=random.log_uniform(1.0e-3,100.0);
      double cosine_theta=random.uniform(0.02,0.999);
      if(random.next_u64()&1ULL) cosine_theta=-cosine_theta;
      const double sine_theta=std::sqrt(std::max(0.0,1.0-cosine_theta*cosine_theta));
      const double field_azimuth=random.uniform(0.0,2.0*swcme::constants::PI);
      const Vec3 field_direction=add(scale(item.normal,cosine_theta),
          add(scale(basis[0],sine_theta*std::cos(field_azimuth)),
              scale(basis[1],sine_theta*std::sin(field_azimuth))));
      item.upstream.magnetic_T=scale(field_direction,field_magnitude);
      item.upstream.pressure_Pa=beta*field_magnitude*field_magnitude/(2.0*MU0);

      const double normal_flow=random.uniform(2.0e5,8.0e5);
      const double tangent_scale=random.uniform(-2.0e5,2.0e5);
      const double tangent_angle=random.uniform(0.0,2.0*swcme::constants::PI);
      item.upstream.velocity_m_s=add(scale(item.normal,normal_flow),
          add(scale(basis[0],tangent_scale*std::cos(tangent_angle)),
              scale(basis[1],tangent_scale*std::sin(tangent_angle))));

      if(index<25000){
        item.stratum=StressStratum::NoShock;
        item.requested_fast_mach=random.uniform(0.2,0.999999);
      } else if(index<50000){
        item.stratum=StressStratum::WeakLimit;
        item.requested_fast_mach=1.0+random.log_uniform(1.0e-10,1.0e-6);
      } else if(index<99900){
        item.stratum=StressStratum::Resolved;
        // The resolved stratum deliberately excludes the weak/parallel,
        // low-beta corner assigned to the adjacent explicit-limit strata.
        // This makes a generic WRONG_BRANCH outcome a genuine regression
        // rather than an undocumented expectation for an ill-conditioned
        // random draw, while the campaign as a whole retains the full ranges.
        beta=random.log_uniform(0.1,10.0);
        item.upstream.pressure_Pa=
            beta*field_magnitude*field_magnitude/(2.0*MU0);
        double resolved_cosine=random.uniform(0.05,0.95);
        if(random.next_u64()&1ULL) resolved_cosine=-resolved_cosine;
        const double resolved_sine=std::sqrt(
            std::max(0.0,1.0-resolved_cosine*resolved_cosine));
        const double resolved_azimuth=
            random.uniform(0.0,2.0*swcme::constants::PI);
        item.upstream.magnetic_T=scale(add(scale(item.normal,resolved_cosine),
            add(scale(basis[0],resolved_sine*std::cos(resolved_azimuth)),
                scale(basis[1],resolved_sine*std::sin(resolved_azimuth)))),
            field_magnitude);
        item.requested_fast_mach=1.0+random.log_uniform(1.0,5.0);
      } else {
        item.stratum=StressStratum::Conditioned;
        // Place the determinant zero on a scan node for a low-beta, nearly
        // parallel state.  Varying the node and polarity across 100 examples
        // stresses singular segmentation without relying on chance sampling.
        beta=random.log_uniform(1.0e-3,1.0e-2);
        const double angle=random.uniform(0.5,3.0)*swcme::constants::PI/180.0;
        const double polarity=(random.next_u64()&1ULL)?1.0:-1.0;
        item.upstream.magnetic_T=add(
            scale(item.normal,polarity*field_magnitude*std::cos(angle)),
            scale(basis[0],polarity*field_magnitude*std::sin(angle)));
        item.upstream.pressure_Pa=beta*field_magnitude*field_magnitude/(2.0*MU0);
        const double scan_index=250.0+static_cast<double>(index%301);
        const double fraction=scan_index/800.0;
        const double rmax=(item.gamma+1.0)/(item.gamma-1.0);
        const double singular_r=1.0+1.0e-9+
            (rmax-1.0-1.0e-9)*fraction*fraction;
        const double Bn=dot(item.upstream.magnetic_T,item.normal);
        const double inflow=std::sqrt(singular_r*Bn*Bn/
                                      (MU0*item.upstream.rho_kg_m3));
        item.requested_fast_mach=inflow/
            reference_fast_speed(item.upstream,item.normal,item.gamma);
      }

      const double fast=reference_fast_speed(item.upstream,item.normal,item.gamma);
      item.shock_speed_m_s=dot(item.upstream.velocity_m_s,item.normal)+
                            item.requested_fast_mach*fast;
      item.result=swcme::shock::solve_ideal_mhd_fast_shock(
          item.upstream,item.normal,item.shock_speed_m_s,item.gamma);
      generated.push_back(item);
    }
    return generated;
  }();
  return cases;
}

std::string stress_reproducer(const StressCase& item){
  std::ostringstream stream;
  stream<<"index="<<item.index<<" stratum="<<stress_stratum_name(item.stratum)
        <<" gamma="<<std::setprecision(17)<<item.gamma
        <<" requested_Mfast="<<item.requested_fast_mach
        <<" status="<<swcme::shock::solve_status_name(item.result.status)
        <<" normal=["<<item.normal[0]<<','<<item.normal[1]<<','<<item.normal[2]<<']'
        <<" rho="<<item.upstream.rho_kg_m3
        <<" pressure="<<item.upstream.pressure_Pa
        <<" velocity=["<<item.upstream.velocity_m_s[0]<<','
        <<item.upstream.velocity_m_s[1]<<','<<item.upstream.velocity_m_s[2]<<']'
        <<" magnetic=["<<item.upstream.magnetic_T[0]<<','
        <<item.upstream.magnetic_T[1]<<','<<item.upstream.magnetic_T[2]<<']'
        <<" shock_speed="<<item.shock_speed_m_s;
  return stream.str();
}

} // namespace

void test_shk01(swcme_test::Context& context){
  std::cout<<"SHK01 fast-shock existence and no-shock threshold\n";
  const Primitive up=fixture(45.0);
  const Vec3 n={{1,0,0}};
  const double cf=reference_fast_speed(up,n);
  const double V1n=up.velocity_m_s[0];

  for(double mf : {0.95,1.0,1.05}){
    const auto r=swcme::shock::solve_ideal_mhd_fast_shock(up,n,V1n+mf*cf,GAMMA);
    const bool expect_shock=mf>1.0;
    context.expect_true(r.has_shock==expect_shock,"fast-shock threshold classification");
    if(!expect_shock){
      expect_abs(context,"no-shock compression",r.compression,1.0,0.0);
      expect_rel(context,"no-shock downstream rho",r.downstream.rho_kg_m3,up.rho_kg_m3,1e-15);
      expect_rel(context,"no-shock downstream pressure",r.downstream.pressure_Pa,up.pressure_Pa,1e-15);
    } else {
      context.expect_true(r.solver_converged,"super-fast fixture must converge");
      context.expect_true(r.compression>1.0,"super-fast fixture must compress");
    }
  }

  // Permanent regression guard for the legacy 1-D bug: an empirical sheath
  // compression floor must not create a shock when V_sh=V_sw.
  swcme1d::Params p1;
  p1.V0_sh_kms=p1.V_sw_kms;
  p1.sheath_comp_floor=1.8;
  const swcme1d::StepState s1=swcme1d::Model(p1).prepare_step(0.0);
  context.expect_true(!s1.has_shock,"1-D floor must not manufacture a shock");
  expect_abs(context,"1-D no-shock rc despite floor",s1.rc,1.0,0.0);

  swcme3d::Params p3;
  p3.shape=swcme3d::ShockShape::Sphere;
  p3.V0_sh_kms=p3.V_sw_kms;
  p3.sheath_comp_floor=1.8;
  const swcme3d::Model m3(p3);
  const auto s3=m3.prepare_step(0.0);
  swcme3d::LocalShockState js;
  const double u[3]={1,0,0};
  context.expect_true(m3.shock_state_direction(s3,u,js),"3-D sphere surface exists");
  context.expect_true(!js.has_shock,"3-D floor must not manufacture a shock");
  expect_abs(context,"3-D no-shock rc despite floor",js.compression,1.0,0.0);
}

void test_shk02(swcme_test::Context& context){
  std::cout<<"SHK02 shock obliquity and polarity invariance\n";
  const Vec3 n={{1,0,0}};
  for(double theta : {0.0,30.0,60.0,90.0}){
    Primitive up=fixture(theta);
    auto a=swcme::shock::solve_ideal_mhd_fast_shock(up,n,1.2e6,GAMMA);
    for(double& b:up.magnetic_T)b=-b;
    auto b=swcme::shock::solve_ideal_mhd_fast_shock(up,n,1.2e6,GAMMA);
    expect_abs(context,"theta_Bn polarity invariant",a.theta_Bn_rad,b.theta_Bn_rad,2e-14);
    expect_abs(context,"theta_Bn known fixture",a.theta_Bn_rad,
               theta*swcme::constants::PI/180.0,2e-14);
  }
}

void test_shk03(swcme_test::Context& context){
  std::cout<<"SHK03 parallel fast-shock analytical limit\n";
  Primitive up=fixture(0.0);
  const Vec3 n={{1,0,0}};
  const double Vsh=9.0e5;
  const double U=Vsh-up.velocity_m_s[0];
  const double cs=std::sqrt(GAMMA*up.pressure_Pa/up.rho_kg_m3);
  const double Ms=U/cs;
  const double rref=((GAMMA+1.0)*Ms*Ms)/((GAMMA-1.0)*Ms*Ms+2.0);
  const double p_ratio=(2.0*GAMMA*Ms*Ms-(GAMMA-1.0))/(GAMMA+1.0);
  const double p2ref=p_ratio*up.pressure_Pa;
  const double V2ref=Vsh-(Vsh-up.velocity_m_s[0])/rref;

  const auto r=swcme::shock::solve_ideal_mhd_fast_shock(up,n,Vsh,GAMMA);
  context.expect_true(r.has_shock&&r.solver_converged,"parallel fast shock must converge");
  expect_rel(context,"parallel compression",r.compression,rref,2e-10);
  expect_rel(context,"parallel downstream pressure",r.downstream.pressure_Pa,p2ref,2e-10);
  expect_rel(context,"parallel downstream normal speed",r.downstream.velocity_m_s[0],V2ref,2e-10);
}

void test_shk04(swcme_test::Context& context){
  std::cout<<"SHK04 perpendicular MHD-shock independent benchmark\n";
  Primitive up=fixture(90.0);
  // Eliminate tiny cos(pi/2) roundoff so this is exactly B_n=0.
  up.magnetic_T={{0.0,5.0e-9,0.0}};
  const double Vsh=8.0e5;
  const double rref=perpendicular_reference_compression(up,Vsh);
  context.expect_true(rref>1.0,"independent perpendicular reference root exists");
  const auto r=swcme::shock::solve_ideal_mhd_fast_shock(up,{{1,0,0}},Vsh,GAMMA);
  context.expect_true(r.has_shock&&r.solver_converged,"perpendicular shock must converge");
  expect_rel(context,"perpendicular compression",r.compression,rref,2e-9);
  expect_rel(context,"perpendicular B_t compression",r.downstream.magnetic_T[1],
             rref*up.magnetic_T[1],2e-9);
  expect_rel(context,"perpendicular density compression",r.downstream.rho_kg_m3,
             rref*up.rho_kg_m3,2e-9);
}

void test_shk05(swcme_test::Context& context){
  std::cout<<"SHK05 independent oblique-MHD benchmark\n";
  namespace reference=swcme_test::shk05_reference_v1;
  constexpr std::size_t case_count=sizeof(reference::CASES)/
                                   sizeof(reference::CASES[0]);

  // These metadata assertions prevent a future edit from quietly replacing
  // the promised high-precision, versioned campaign with a few hand-entered
  // doubles.  The generator solves all eight RH equations with Decimal
  // precision 80 and records how many independent Newton seeds reached the
  // unique evolutionary fast branch.
  context.expect_true(reference::FIXTURE_VERSION==1,
                      "SHK05 fixture schema/version is pinned to v1");
  context.expect_true(reference::DECIMAL_PRECISION_DIGITS>=50,
                      "SHK05 reference solve uses at least 50 decimal digits");
  context.expect_true(case_count>=12,
                      "SHK05 matrix spans at least twelve oblique states");

  for(const reference::Case& expected : reference::CASES){
    Primitive upstream;
    upstream.rho_kg_m3=expected.upstream_rho_kg_m3;
    upstream.pressure_Pa=expected.upstream_pressure_Pa;
    for(int component=0;component<3;++component){
      upstream.velocity_m_s[component]=expected.upstream_velocity_m_s[component];
      upstream.magnetic_T[component]=expected.upstream_magnetic_T[component];
    }
    const Vec3 normal={{1.0,0.0,0.0}};
    const auto actual=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,normal,expected.shock_speed_m_s,expected.gamma);
    const std::string prefix=std::string("SHK05 ")+expected.id+" ";

    // Fixture-generation diagnostics are release data, not merely comments.
    // Every case must have a unique admissible root reached from multiple
    // initial guesses, an 80-digit equation residual below 1e-50, and a
    // nonsingular final Newton path.  A fixture that fails these checks is too
    // poorly conditioned to distinguish a production defect from reference
    // uncertainty and must be regenerated as a new version.
    context.expect_true(std::string(expected.branch)=="EVOLUTIONARY_FAST",
                        prefix+"reference branch is evolutionary fast");
    context.expect_true(expected.physical_root_count==1 &&
                            expected.physical_seed_count>=2 &&
                            expected.converged_seed_count>=
                                expected.physical_seed_count,
                        prefix+"independent seeds select one physical root");
    context.expect_true(expected.reference_max_residual<1.0e-50 &&
                            expected.minimum_newton_pivot>1.0e-10,
                        prefix+"reference residual and conditioning are qualified");

    context.expect_true(actual.has_shock && actual.solver_converged &&
                            actual.status==swcme::shock::SolveStatus::Solved,
                        prefix+"production solver selects a solved fast shock");
    if(!actual.solver_converged) continue;

    // Compare every downstream primitive, not only compression.  The 1e-7
    // validation-plan requirement is tightened to 2e-9 for these deliberately
    // well-conditioned fixtures.  Component-specific floors keep a true zero
    // from being assigned a meaningless relative error while still detecting
    // spurious cross-plane velocity or magnetic field.
    expect_rel(context,prefix+"upstream fast Mach",actual.fast_mach,
               expected.fast_mach,2.0e-10,1.0);
    expect_rel(context,prefix+"compression",actual.compression,
               expected.compression,2.0e-9,1.0);
    expect_rel(context,prefix+"downstream density",
               actual.downstream.rho_kg_m3,expected.downstream_rho_kg_m3,
               2.0e-9,1.0e-30);
    expect_rel(context,prefix+"downstream pressure",
               actual.downstream.pressure_Pa,expected.downstream_pressure_Pa,
               2.0e-9,1.0e-30);
    for(int component=0;component<3;++component){
      expect_rel(context,prefix+"downstream velocity["+
                     std::to_string(component)+"]",
                 actual.downstream.velocity_m_s[component],
                 expected.downstream_velocity_m_s[component],2.0e-9,1.0);
      expect_rel(context,prefix+"downstream magnetic["+
                     std::to_string(component)+"]",
                 actual.downstream.magnetic_T[component],
                 expected.downstream_magnetic_T[component],2.0e-9,1.0e-20);
    }
    expect_rel(context,prefix+"entropy ratio",actual.entropy_ratio,
               expected.entropy_ratio,2.0e-9,1.0);

    // Independently classify the returned branch from downstream
    // characteristics.  A production result is evolutionary fast only if
    // the upstream normal flow is super-fast, the downstream normal flow is
    // sub-fast but remains super-Alfvenic, and entropy increases.
    const double downstream_fast_speed=reference_fast_speed(
        actual.downstream,normal,expected.gamma);
    const double downstream_normal_flow=std::abs(
        actual.downstream.velocity_m_s[0]-expected.shock_speed_m_s);
    const double downstream_fast_mach=
        downstream_normal_flow/downstream_fast_speed;
    const double downstream_normal_alfven=
        std::abs(actual.downstream.magnetic_T[0])/
        std::sqrt(MU0*actual.downstream.rho_kg_m3);
    expect_rel(context,prefix+"downstream fast Mach",downstream_fast_mach,
               expected.downstream_fast_mach,2.0e-9,1.0);
    context.expect_true(actual.fast_mach>1.0 && downstream_fast_mach<1.0 &&
                            downstream_normal_flow>downstream_normal_alfven &&
                            actual.entropy_ratio>1.0,
                        prefix+"production branch is evolutionary fast");

    // Production conservation diagnostics are recomputed from the accepted
    // primitives by swcme_shock.hpp and are deliberately separate from the
    // frozen reference.  Requiring all five here makes a full-state mismatch
    // retain useful evidence about which invariant failed.
    context.expect_true(actual.mass_residual<=1.0e-9 &&
                            actual.normal_B_residual<=1.0e-10 &&
                            actual.electric_residual<=1.0e-8 &&
                            actual.momentum_residual<=1.0e-8 &&
                            actual.energy_residual<=1.0e-8,
                        prefix+"all production RH residuals pass");
  }
}

void test_shk06(swcme_test::Context& context){
  std::cout<<"SHK06 Rankine-Hugoniot mass-flux conservation\n";
  constexpr long double acceptance=1.0e-9L;
  long double maximum_residual=0.0L;
  std::size_t solved_count=0;
  std::size_t weak_count=0;
  std::size_t high_compression_count=0;
  std::size_t failures_reported=0;

  // Compute the invariant solely from the serialized upstream/downstream
  // primitives and public normal/speed inputs.  Long-double accumulation and
  // a symmetric larger-flux normalization keep this oracle independent of the
  // production residual field and meaningful for very small physical fluxes.
  auto check_mass_flux=[&](const Primitive& upstream,
                           const swcme::shock::JumpResult& result,
                           const Vec3& normal,double shock_speed,
                           const std::string& reproducer){
    if(result.status!=swcme::shock::SolveStatus::Solved) return;
    long double u1n=0.0L,u2n=0.0L;
    for(int component=0;component<3;++component){
      const long double shock_component=
          static_cast<long double>(shock_speed)*normal[component];
      u1n+=(static_cast<long double>(upstream.velocity_m_s[component])-
            shock_component)*normal[component];
      u2n+=(static_cast<long double>(result.downstream.velocity_m_s[component])-
            shock_component)*normal[component];
    }
    const long double flux1=static_cast<long double>(upstream.rho_kg_m3)*u1n;
    const long double flux2=
        static_cast<long double>(result.downstream.rho_kg_m3)*u2n;
    const long double scale_flux=std::max(
        {std::abs(flux1),std::abs(flux2),
         std::numeric_limits<long double>::min()});
    const long double residual=std::abs(flux1-flux2)/scale_flux;
    maximum_residual=std::max(maximum_residual,residual);
    ++solved_count;
    if(result.compression<1.01) ++weak_count;
    if(result.compression>4.0) ++high_compression_count;
    if((!std::isfinite(residual) || residual>acceptance) &&
       failures_reported<8){
      std::cerr<<"    SHK06 mass-flux failure residual="
               <<static_cast<double>(residual)<<' '<<reproducer<<'\n';
      ++failures_reported;
    }
  };

  for(const StressCase& item:shock_stress_cases()){
    check_mass_flux(item.upstream,item.result,item.normal,
                    item.shock_speed_m_s,stress_reproducer(item));
  }

  // Add the independently frozen SHK12 weak branches because the resolved
  // random stratum intentionally starts at Mach two.  These cases make the
  // mass-flux gate sensitive near compression one without relaxing SHK15's
  // requirement that all ordinary random states solve cleanly.
  for(const auto& expected:swcme_test::shk12_reference_v1::CASES){
    Primitive upstream;
    upstream.rho_kg_m3=expected.upstream_rho_kg_m3;
    upstream.pressure_Pa=expected.upstream_pressure_Pa;
    for(int component=0;component<3;++component){
      upstream.velocity_m_s[component]=expected.upstream_velocity_m_s[component];
      upstream.magnetic_T[component]=expected.upstream_magnetic_T[component];
    }
    const Vec3 normal={{1.0,0.0,0.0}};
    const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,normal,expected.shock_speed_m_s,expected.gamma);
    check_mass_flux(upstream,result,normal,expected.shock_speed_m_s,
                    std::string("fixture=")+expected.id);
  }

  std::cout<<"  independently checked="<<solved_count
           <<" weak="<<weak_count
           <<" high_compression="<<high_compression_count
           <<" max_normalized_residual="<<std::scientific
           <<static_cast<double>(maximum_residual)<<'\n';
  context.expect_true(solved_count>=49000,
                      "SHK06 independently checks the complete solved stress population");
  context.expect_true(weak_count>=8,
                      "SHK06 includes weak shocks below compression 1.01");
  context.expect_true(high_compression_count>=8,
                      "SHK06 includes high-compression shocks above four");
  context.expect_true(std::isfinite(maximum_residual) &&
                          maximum_residual<=acceptance,
                      "SHK06 independent normalized mass flux meets 1e-9 threshold");

  // Integration check: RESOLVED_COMPRESSION reaches the exact RH state at
  // the inner edge of its finite C1 shock layer.  The mathematical shock state
  // itself remains exact through shock_state_direction().
  swcme3d::Params p;
  p.shape=swcme3d::ShockShape::Sphere; p.r0_Rs=20.0; p.V0_sh_kms=1200.0;
  p.V_sw_kms=400.0; p.Gamma_kmInv=1e-8; p.sin_theta=1.0;
  // The canonical default is SOURCE/SHOCK_ONLY after Fix 14. This test is
  // specifically a resolved-compression integration check, so both members of
  // the alternative FULL_ICME/RESOLVED_COMPRESSION pair are explicit.
  p.region_mode=swcme::regions::Mode::FullICME;
  p.shock_acceleration_mode=swcme::acceleration::Mode::ResolvedCompression;
  const swcme3d::Model m(p); const auto s=m.prepare_step(0.0);
  const double u[3]={1,0,0}; swcme3d::LocalShockState js;
  context.expect_true(m.shock_state_direction(s,u,js)&&js.has_shock&&js.solver_converged,
                      "3-D apex shock state exists");
  const auto boundaries=swcme::regions::make_boundaries(js.Rdir_m,s.region_config);
  const double rq=js.Rdir_m-0.5*boundaries.smooth_shock_width_m;
  const double x=rq, y=0.0,z=0.0;
  double nq=0,vx=0,vy=0,vz=0,bx=0,by=0,bz=0;
  m.evaluate_cartesian_with_B(s,&x,&y,&z,&nq,&vx,&vy,&vz,&bx,&by,&bz,1);
  expect_rel(context,"3-D resolved-shock inner density",nq,js.downstream_n_m3,2e-7);
  expect_rel(context,"3-D resolved-shock inner Vx",vx,js.downstream.velocity_m_s[0],2e-7,1.0);
  expect_rel(context,"3-D resolved-shock inner Vy",vy,js.downstream.velocity_m_s[1],2e-7,1.0);
  expect_rel(context,"3-D resolved-shock inner Bx",bx,js.downstream.magnetic_T[0],2e-7,1e-30);
  expect_rel(context,"3-D resolved-shock inner By",by,js.downstream.magnetic_T[1],2e-7,1e-30);
}

void test_shk07(swcme_test::Context& context){
  std::cout<<"SHK07 Rankine-Hugoniot normal magnetic-field continuity\n";
  constexpr long double acceptance=1.0e-10L;
  long double maximum_residual=0.0L;
  std::size_t checked=0;
  std::size_t positive_polarity=0;
  std::size_t negative_polarity=0;
  std::size_t failures_reported=0;

  // Evaluate B.n from the serialized Cartesian vectors using long-double
  // products.  The scale contains the physical normal field and a tiny total-
  // field guard, so nearly perpendicular samples remain finite without making
  // an absolute tesla tolerance part of this dimensionless invariant.
  auto check_normal_field=[&](const Primitive& upstream,
                              const swcme::shock::JumpResult& result,
                              const Vec3& normal,
                              const std::string& reproducer){
    if(result.status!=swcme::shock::SolveStatus::Solved) return;
    long double B1n=0.0L,B2n=0.0L,B1sq=0.0L,B2sq=0.0L;
    for(int component=0;component<3;++component){
      const long double B1=upstream.magnetic_T[component];
      const long double B2=result.downstream.magnetic_T[component];
      B1n+=B1*normal[component];
      B2n+=B2*normal[component];
      B1sq+=B1*B1;
      B2sq+=B2*B2;
    }
    const long double field_scale=std::max(std::sqrt(B1sq),std::sqrt(B2sq));
    const long double scale_normal=std::max(
        {std::abs(B1n),std::abs(B2n),1.0e-12L*field_scale,1.0e-300L});
    const long double residual=std::abs(B1n-B2n)/scale_normal;
    maximum_residual=std::max(maximum_residual,residual);
    ++checked;
    if(B1n>=0.0L) ++positive_polarity; else ++negative_polarity;
    if((!std::isfinite(residual) || residual>acceptance) &&
       failures_reported<8){
      std::cerr<<"    SHK07 normal-B failure residual="
               <<static_cast<double>(residual)<<' '<<reproducer<<'\n';
      ++failures_reported;
    }
  };

  for(const StressCase& item:shock_stress_cases()){
    check_normal_field(item.upstream,item.result,item.normal,
                       stress_reproducer(item));
  }

  std::size_t covariance_pairs=0;
  for(const StressCase& item:shock_stress_cases()){
    if(item.result.status!=swcme::shock::SolveStatus::Solved ||
       covariance_pairs>=256) continue;

    // Cyclic coordinate permutation is a proper three-dimensional rotation
    // (determinant +1).  Applying it to normal, velocity, and magnetic field
    // supplies a direct rotational-covariance pair without transforming any
    // scalar shock input or expected result through production helpers.
    auto rotate=[](const Vec3& value){
      return Vec3{{value[2],value[0],value[1]}};
    };
    Primitive rotated=item.upstream;
    rotated.velocity_m_s=rotate(item.upstream.velocity_m_s);
    rotated.magnetic_T=rotate(item.upstream.magnetic_T);
    const Vec3 rotated_normal=rotate(item.normal);
    const auto rotated_result=swcme::shock::solve_ideal_mhd_fast_shock(
        rotated,rotated_normal,item.shock_speed_m_s,item.gamma);
    check_normal_field(rotated,rotated_result,rotated_normal,
                       "rotated "+stress_reproducer(item));
    context.expect_true(
        rotated_result.status==swcme::shock::SolveStatus::Solved &&
            std::abs(rotated_result.compression-item.result.compression)<=
                2.0e-10*std::max(1.0,item.result.compression),
        "SHK07 proper rotation preserves solved compression");

    Primitive reversed=item.upstream;
    reversed.magnetic_T=scale(reversed.magnetic_T,-1.0);
    const auto reversed_result=swcme::shock::solve_ideal_mhd_fast_shock(
        reversed,item.normal,item.shock_speed_m_s,item.gamma);
    check_normal_field(reversed,reversed_result,item.normal,
                       "polarity-reversed "+stress_reproducer(item));
    context.expect_true(
        reversed_result.status==swcme::shock::SolveStatus::Solved &&
            std::abs(reversed_result.compression-item.result.compression)<=
                2.0e-10*std::max(1.0,item.result.compression),
        "SHK07 magnetic-polarity reversal preserves solved compression");
    ++covariance_pairs;
  }

  std::cout<<"  independently checked="<<checked
           <<" rotation/polarity pairs="<<covariance_pairs
           <<" positive_Bn="<<positive_polarity
           <<" negative_Bn="<<negative_polarity
           <<" max_normalized_residual="<<std::scientific
           <<static_cast<double>(maximum_residual)<<'\n';
  context.expect_true(checked>=49000,
                      "SHK07 independently checks the solved stress population");
  context.expect_true(covariance_pairs==256,
                      "SHK07 completes all proper-rotation and polarity pairs");
  context.expect_true(positive_polarity>10000 && negative_polarity>10000,
                      "SHK07 covers both normal-field polarities broadly");
  context.expect_true(std::isfinite(maximum_residual) &&
                          maximum_residual<=acceptance,
                      "SHK07 independent normal-B residual meets 1e-10 threshold");
}

void test_shk08(swcme_test::Context& context){
  std::cout<<"SHK08 Rankine-Hugoniot tangential electric-field conservation\n";
  constexpr long double acceptance=1.0e-8L;
  long double maximum_component_residual=0.0L;
  std::size_t checked_states=0;
  std::size_t checked_components=0;
  std::size_t deterministic_states=0;
  std::size_t failures_reported=0;

  // Reconstruct E=-u x B in Cartesian coordinates from public primitives.
  // Neither the production cross-product helper nor electric_residual enters
  // this oracle.  Each of two test-owned tangential axes is checked separately
  // so cancellation between components cannot hide a sign or ordering error.
  auto check_electric_field=[&](const Primitive& upstream,
                                const swcme::shock::JumpResult& result,
                                const Vec3& normal,double shock_speed,
                                const std::string& reproducer){
    if(result.status!=swcme::shock::SolveStatus::Solved) return;
    const auto basis=tangential_basis(normal);
    std::array<long double,3> u1{{0.0L,0.0L,0.0L}};
    std::array<long double,3> u2{{0.0L,0.0L,0.0L}};
    for(int component=0;component<3;++component){
      const long double shock_component=
          static_cast<long double>(shock_speed)*normal[component];
      u1[component]=upstream.velocity_m_s[component]-shock_component;
      u2[component]=result.downstream.velocity_m_s[component]-shock_component;
    }
    auto electric=[](const std::array<long double,3>& velocity,
                     const Vec3& magnetic){
      return std::array<long double,3>{{
          -(velocity[1]*magnetic[2]-velocity[2]*magnetic[1]),
          -(velocity[2]*magnetic[0]-velocity[0]*magnetic[2]),
          -(velocity[0]*magnetic[1]-velocity[1]*magnetic[0])}};
    };
    const auto E1=electric(u1,upstream.magnetic_T);
    const auto E2=electric(u2,result.downstream.magnetic_T);
    long double E1sq=0.0L,E2sq=0.0L;
    for(int component=0;component<3;++component){
      E1sq+=E1[component]*E1[component];
      E2sq+=E2[component]*E2[component];
    }
    const long double vector_scale=std::max(std::sqrt(E1sq),std::sqrt(E2sq));
    for(int tangent=0;tangent<2;++tangent){
      long double component1=0.0L,component2=0.0L;
      for(int component=0;component<3;++component){
        component1+=E1[component]*basis[tangent][component];
        component2+=E2[component]*basis[tangent][component];
      }
      const long double component_scale=std::max(
          {std::abs(component1),std::abs(component2),
           1.0e-12L*vector_scale,1.0e-300L});
      const long double residual=
          std::abs(component1-component2)/component_scale;
      maximum_component_residual=std::max(maximum_component_residual,residual);
      ++checked_components;
      if((!std::isfinite(residual) || residual>acceptance) &&
         failures_reported<8){
        std::cerr<<"    SHK08 tangential component "<<tangent
                 <<" failure residual="<<static_cast<double>(residual)
                 <<' '<<reproducer<<'\n';
        ++failures_reported;
      }
    }
    ++checked_states;
  };

  for(const StressCase& item:shock_stress_cases()){
    check_electric_field(item.upstream,item.result,item.normal,
                         item.shock_speed_m_s,stress_reproducer(item));
  }

  // The frozen SHK05 matrix supplies deterministic oblique states whose full
  // downstream primitives originate in an 80-digit independent eight-equation
  // solve.  Re-solving those inputs here complements random stress with named,
  // reviewable regression fixtures.
  for(const auto& expected:swcme_test::shk05_reference_v1::CASES){
    Primitive upstream;
    upstream.rho_kg_m3=expected.upstream_rho_kg_m3;
    upstream.pressure_Pa=expected.upstream_pressure_Pa;
    for(int component=0;component<3;++component){
      upstream.velocity_m_s[component]=expected.upstream_velocity_m_s[component];
      upstream.magnetic_T[component]=expected.upstream_magnetic_T[component];
    }
    const Vec3 normal={{1.0,0.0,0.0}};
    const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,normal,expected.shock_speed_m_s,expected.gamma);
    check_electric_field(upstream,result,normal,expected.shock_speed_m_s,
                         std::string("fixture=")+expected.id);
    if(result.status==swcme::shock::SolveStatus::Solved) ++deterministic_states;
  }

  std::size_t covariance_pairs=0;
  for(const StressCase& item:shock_stress_cases()){
    if(item.result.status!=swcme::shock::SolveStatus::Solved ||
       covariance_pairs>=256) continue;
    auto rotate=[](const Vec3& value){
      return Vec3{{value[2],value[0],value[1]}};
    };
    Primitive rotated=item.upstream;
    rotated.velocity_m_s=rotate(item.upstream.velocity_m_s);
    rotated.magnetic_T=rotate(item.upstream.magnetic_T);
    const Vec3 rotated_normal=rotate(item.normal);
    const auto rotated_result=swcme::shock::solve_ideal_mhd_fast_shock(
        rotated,rotated_normal,item.shock_speed_m_s,item.gamma);
    check_electric_field(rotated,rotated_result,rotated_normal,
                         item.shock_speed_m_s,
                         "rotated "+stress_reproducer(item));
    context.expect_true(
        rotated_result.status==swcme::shock::SolveStatus::Solved &&
            std::abs(rotated_result.compression-item.result.compression)<=
                2.0e-10*std::max(1.0,item.result.compression),
        "SHK08 proper rotation preserves solved compression");

    Primitive reversed=item.upstream;
    reversed.magnetic_T=scale(reversed.magnetic_T,-1.0);
    const auto reversed_result=swcme::shock::solve_ideal_mhd_fast_shock(
        reversed,item.normal,item.shock_speed_m_s,item.gamma);
    check_electric_field(reversed,reversed_result,item.normal,
                         item.shock_speed_m_s,
                         "polarity-reversed "+stress_reproducer(item));
    context.expect_true(
        reversed_result.status==swcme::shock::SolveStatus::Solved &&
            std::abs(reversed_result.compression-item.result.compression)<=
                2.0e-10*std::max(1.0,item.result.compression),
        "SHK08 magnetic-polarity reversal preserves solved compression");
    ++covariance_pairs;
  }

  std::cout<<"  independently checked states="<<checked_states
           <<" components="<<checked_components
           <<" deterministic="<<deterministic_states
           <<" rotation/polarity pairs="<<covariance_pairs
           <<" max_component_residual="<<std::scientific
           <<static_cast<double>(maximum_component_residual)<<'\n';
  context.expect_true(checked_states>=49000 &&
                          checked_components==2*checked_states,
                      "SHK08 checks two tangential components for every solved state");
  context.expect_true(deterministic_states>=12,
                      "SHK08 includes all named independent oblique fixtures");
  context.expect_true(covariance_pairs==256,
                      "SHK08 completes all proper-rotation and polarity pairs");
  context.expect_true(std::isfinite(maximum_component_residual) &&
                          maximum_component_residual<=acceptance,
                      "SHK08 independent tangential-E components meet 1e-8 threshold");
}

void test_shk09(swcme_test::Context& context){
  std::cout<<"SHK09 Rankine-Hugoniot momentum-flux conservation\n";
  constexpr long double acceptance=1.0e-8L;
  std::array<long double,3> maximum_residual{{0.0L,0.0L,0.0L}};
  std::size_t checked_states=0;
  std::size_t deterministic_states=0;
  std::size_t failures_reported=0;

  // Rebuild the complete ideal-MHD momentum flux from serialized primitive
  // records.  The implementation below deliberately avoids production vector
  // helpers and the stored momentum_residual.  Long-double accumulation makes
  // the comparison sensitive to missing pressure or magnetic-stress terms
  // rather than to ordinary binary64 cancellation in Cartesian projection.
  auto check_momentum=[&](const Primitive& upstream,
                          const swcme::shock::JumpResult& result,
                          const Vec3& normal,double shock_speed,
                          const std::string& reproducer){
    if(result.status!=swcme::shock::SolveStatus::Solved) return;
    const auto tangents=tangential_basis(normal);
    const std::array<Vec3,3> axes={{normal,tangents[0],tangents[1]}};

    struct FluxTerms {
      std::array<long double,3> total{{0.0L,0.0L,0.0L}};
      std::array<long double,3> dynamic_scale{{0.0L,0.0L,0.0L}};
      std::array<long double,3> thermal_scale{{0.0L,0.0L,0.0L}};
      std::array<long double,3> magnetic_scale{{0.0L,0.0L,0.0L}};
    };
    auto flux=[&](const Primitive& state){
      FluxTerms value;
      std::array<long double,3> velocity{{0.0L,0.0L,0.0L}};
      long double un=0.0L,Bn=0.0L,B2=0.0L;
      for(int component=0;component<3;++component){
        velocity[component]=state.velocity_m_s[component]-
            static_cast<long double>(shock_speed)*normal[component];
        un+=velocity[component]*normal[component];
        Bn+=static_cast<long double>(state.magnetic_T[component])*
            normal[component];
        B2+=static_cast<long double>(state.magnetic_T[component])*
            state.magnetic_T[component];
      }
      for(int axis_index=0;axis_index<3;++axis_index){
        long double ua=0.0L,Ba=0.0L,na=0.0L;
        for(int component=0;component<3;++component){
          ua+=velocity[component]*axes[axis_index][component];
          Ba+=static_cast<long double>(state.magnetic_T[component])*
              axes[axis_index][component];
          na+=static_cast<long double>(normal[component])*
              axes[axis_index][component];
        }
        const long double dynamic=state.rho_kg_m3*un*ua;
        const long double thermal=na*state.pressure_Pa;
        const long double magnetic=na*B2/(2.0L*MU0)-Bn*Ba/MU0;
        value.total[axis_index]=dynamic+thermal+magnetic;
        value.dynamic_scale[axis_index]=std::abs(dynamic);
        value.thermal_scale[axis_index]=std::abs(thermal);
        value.magnetic_scale[axis_index]=
            std::abs(na*B2/(2.0L*MU0))+std::abs(Bn*Ba/MU0);
      }
      return value;
    };

    const FluxTerms before=flux(upstream);
    const FluxTerms after=flux(result.downstream);
    for(int component=0;component<3;++component){
      // Sum the dynamic, thermal, and magnetic magnitudes on the more energetic
      // side.  This documented scale cannot become spuriously small through
      // cancellation of physical terms in the conserved total itself.
      const long double before_scale=before.dynamic_scale[component]+
          before.thermal_scale[component]+before.magnetic_scale[component];
      const long double after_scale=after.dynamic_scale[component]+
          after.thermal_scale[component]+after.magnetic_scale[component];
      const long double physical_scale=std::max(
          {before_scale,after_scale,std::numeric_limits<long double>::min()});
      const long double residual=
          std::abs(before.total[component]-after.total[component])/physical_scale;
      maximum_residual[component]=std::max(maximum_residual[component],residual);
      if((!std::isfinite(residual) || residual>acceptance) &&
         failures_reported<8){
        std::cerr<<"    SHK09 component="<<component
                 <<" residual="<<static_cast<double>(residual)
                 <<" dynamic_scale="
                 <<static_cast<double>(std::max(before.dynamic_scale[component],
                                                 after.dynamic_scale[component]))
                 <<" thermal_scale="
                 <<static_cast<double>(std::max(before.thermal_scale[component],
                                                 after.thermal_scale[component]))
                 <<" magnetic_scale="
                 <<static_cast<double>(std::max(before.magnetic_scale[component],
                                                 after.magnetic_scale[component]))
                 <<' '<<reproducer<<'\n';
        ++failures_reported;
      }
    }
    ++checked_states;
  };

  for(const StressCase& item:shock_stress_cases()){
    check_momentum(item.upstream,item.result,item.normal,item.shock_speed_m_s,
                   stress_reproducer(item));
  }
  for(const auto& expected:swcme_test::shk05_reference_v1::CASES){
    Primitive upstream;
    upstream.rho_kg_m3=expected.upstream_rho_kg_m3;
    upstream.pressure_Pa=expected.upstream_pressure_Pa;
    for(int component=0;component<3;++component){
      upstream.velocity_m_s[component]=expected.upstream_velocity_m_s[component];
      upstream.magnetic_T[component]=expected.upstream_magnetic_T[component];
    }
    const Vec3 normal={{1.0,0.0,0.0}};
    const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,normal,expected.shock_speed_m_s,expected.gamma);
    check_momentum(upstream,result,normal,expected.shock_speed_m_s,
                   std::string("fixture=")+expected.id);
    if(result.status==swcme::shock::SolveStatus::Solved) ++deterministic_states;
  }

  std::cout<<"  independently checked="<<checked_states
           <<" deterministic="<<deterministic_states
           <<" max_normal="<<std::scientific
           <<static_cast<double>(maximum_residual[0])
           <<" max_t1="<<static_cast<double>(maximum_residual[1])
           <<" max_t2="<<static_cast<double>(maximum_residual[2])<<'\n';
  context.expect_true(checked_states>=49000,
                      "SHK09 checks the complete solved stress population");
  context.expect_true(deterministic_states>=12,
                      "SHK09 checks all named high-precision oblique fixtures");
  context.expect_true(std::all_of(maximum_residual.begin(),maximum_residual.end(),
                      [](long double value){
                        return std::isfinite(value) && value<=acceptance;
                      }),
                      "SHK09 all three independent momentum components meet 1e-8");
}

void test_shk10(swcme_test::Context& context){
  std::cout<<"SHK10 Rankine-Hugoniot total-energy-flux conservation\n";
  constexpr long double acceptance=1.0e-8L;
  long double maximum_residual=0.0L;
  double maximum_production_residual=0.0;
  std::size_t checked_states=0;
  std::size_t deterministic_states=0;
  std::size_t failures_reported=0;

  // Keep the four physical contributions separate.  This independent
  // expression is algebraically equivalent to the ideal-MHD total-energy
  // flux, but it neither calls total_energy_flux_normal() nor consumes the
  // stored energy residual.  Long-double products reduce evaluation
  // cancellation after the binary64 primitive state has been serialized.
  auto check_energy=[&](const Primitive& upstream,
                        const swcme::shock::JumpResult& result,
                        const Vec3& normal,double shock_speed,double gamma,
                        const std::string& reproducer){
    if(result.status!=swcme::shock::SolveStatus::Solved) return;
    struct EnergyTerms {
      long double kinetic=0.0L;
      long double enthalpy=0.0L;
      long double magnetic_advection=0.0L;
      long double magnetic_work=0.0L;
      long double total() const {
        return kinetic+enthalpy+magnetic_advection+magnetic_work;
      }
      long double scale() const {
        return std::abs(kinetic)+std::abs(enthalpy)+
               std::abs(magnetic_advection)+std::abs(magnetic_work);
      }
    };
    auto terms=[&](const Primitive& state){
      std::array<long double,3> velocity{{0.0L,0.0L,0.0L}};
      long double un=0.0L,Bn=0.0L,u2=0.0L,B2=0.0L,u_dot_B=0.0L;
      for(int component=0;component<3;++component){
        velocity[component]=state.velocity_m_s[component]-
            static_cast<long double>(shock_speed)*normal[component];
        un+=velocity[component]*normal[component];
        Bn+=static_cast<long double>(state.magnetic_T[component])*
            normal[component];
        u2+=velocity[component]*velocity[component];
        B2+=static_cast<long double>(state.magnetic_T[component])*
            state.magnetic_T[component];
        u_dot_B+=velocity[component]*state.magnetic_T[component];
      }
      EnergyTerms value;
      value.kinetic=un*0.5L*state.rho_kg_m3*u2;
      value.enthalpy=un*(static_cast<long double>(gamma)/(gamma-1.0L))*
                     state.pressure_Pa;
      value.magnetic_advection=un*B2/MU0;
      value.magnetic_work=-Bn*u_dot_B/MU0;
      return value;
    };

    const EnergyTerms before=terms(upstream);
    const EnergyTerms after=terms(result.downstream);
    const long double physical_scale=std::max(
        {before.scale(),after.scale(),std::numeric_limits<long double>::min()});
    const long double residual=
        std::abs(before.total()-after.total())/physical_scale;
    maximum_residual=std::max(maximum_residual,residual);
    maximum_production_residual=
        std::max(maximum_production_residual,result.energy_residual);
    ++checked_states;
    if((!std::isfinite(residual) || residual>acceptance) &&
       failures_reported<8){
      std::cerr<<"    SHK10 residual="<<static_cast<double>(residual)
               <<" kinetic="
               <<static_cast<double>(std::max(std::abs(before.kinetic),
                                               std::abs(after.kinetic)))
               <<" enthalpy="
               <<static_cast<double>(std::max(std::abs(before.enthalpy),
                                               std::abs(after.enthalpy)))
               <<" magnetic_advection="
               <<static_cast<double>(std::max(
                    std::abs(before.magnetic_advection),
                    std::abs(after.magnetic_advection)))
               <<" magnetic_work="
               <<static_cast<double>(std::max(std::abs(before.magnetic_work),
                                               std::abs(after.magnetic_work)))
               <<' '<<reproducer<<'\n';
      ++failures_reported;
    }
  };

  for(const StressCase& item:shock_stress_cases()){
    check_energy(item.upstream,item.result,item.normal,item.shock_speed_m_s,
                 item.gamma,stress_reproducer(item));
  }
  for(const auto& expected:swcme_test::shk05_reference_v1::CASES){
    Primitive upstream;
    upstream.rho_kg_m3=expected.upstream_rho_kg_m3;
    upstream.pressure_Pa=expected.upstream_pressure_Pa;
    for(int component=0;component<3;++component){
      upstream.velocity_m_s[component]=expected.upstream_velocity_m_s[component];
      upstream.magnetic_T[component]=expected.upstream_magnetic_T[component];
    }
    const Vec3 normal={{1.0,0.0,0.0}};
    const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,normal,expected.shock_speed_m_s,expected.gamma);
    check_energy(upstream,result,normal,expected.shock_speed_m_s,
                 expected.gamma,std::string("fixture=")+expected.id);
    context.expect_true(expected.reference_max_residual<1.0e-50,
                        std::string("SHK10 high-precision residual ")+expected.id);
    if(result.status==swcme::shock::SolveStatus::Solved) ++deterministic_states;
  }

  std::cout<<"  independently checked="<<checked_states
           <<" deterministic="<<deterministic_states
           <<" max_long_double_residual="<<std::scientific
           <<static_cast<double>(maximum_residual)
           <<" max_production_residual="<<maximum_production_residual<<'\n';
  context.expect_true(checked_states>=49000,
                      "SHK10 checks the complete solved stress population");
  context.expect_true(deterministic_states>=12,
                      "SHK10 checks all named high-precision oblique fixtures");
  context.expect_true(std::isfinite(maximum_residual) &&
                          maximum_residual<=acceptance,
                      "SHK10 independent energy flux meets 1e-8 threshold");
}

void test_shk11(swcme_test::Context& context){
  std::cout<<"SHK11 physical admissibility and entropy increase\n";
  std::size_t solved=0,no_shock=0,numerical_limit=0;
  std::size_t rejected_generic=0,failures_reported=0;
  double minimum_entropy_ratio=std::numeric_limits<double>::max();
  double maximum_compression_fraction=0.0;

  // Apply a test-owned branch classification to every campaign outcome.  The
  // characteristic speeds come from the independent reference_fast_speed()
  // formula, while entropy and normal Alfven speed are reconstructed directly
  // from serialized primitives.  This prevents the production status flag
  // from serving as its own admissibility oracle.
  auto check_case=[&](const Primitive& upstream,
                      const swcme::shock::JumpResult& result,
                      const Vec3& normal,double shock_speed,double gamma,
                      const std::string& reproducer){
    bool pass=true;
    if(result.status==swcme::shock::SolveStatus::Solved){
      ++solved;
      const double gamma_bound=(gamma+1.0)/(gamma-1.0);
      const double upstream_fast=reference_fast_speed(upstream,normal,gamma);
      const double downstream_fast=
          reference_fast_speed(result.downstream,normal,gamma);
      double upstream_un=0.0,downstream_un=0.0,downstream_Bn=0.0;
      for(int component=0;component<3;++component){
        const double shock_component=shock_speed*normal[component];
        upstream_un+=(upstream.velocity_m_s[component]-shock_component)*
                     normal[component];
        downstream_un+=(result.downstream.velocity_m_s[component]-shock_component)*
                       normal[component];
        downstream_Bn+=result.downstream.magnetic_T[component]*normal[component];
      }
      const double downstream_normal_alfven=std::abs(downstream_Bn)/
          std::sqrt(MU0*result.downstream.rho_kg_m3);
      const double entropy1=upstream.pressure_Pa/
          std::pow(upstream.rho_kg_m3,gamma);
      const double entropy2=result.downstream.pressure_Pa/
          std::pow(result.downstream.rho_kg_m3,gamma);
      const double entropy_ratio=entropy2/entropy1;
      minimum_entropy_ratio=std::min(minimum_entropy_ratio,entropy_ratio);
      maximum_compression_fraction=std::max(
          maximum_compression_fraction,result.compression/gamma_bound);
      pass=result.has_shock && result.solver_converged &&
          upstream.rho_kg_m3>0.0 && upstream.pressure_Pa>0.0 &&
          result.downstream.rho_kg_m3>upstream.rho_kg_m3 &&
          result.downstream.pressure_Pa>0.0 &&
          result.compression>1.0 &&
          result.compression<=gamma_bound*(1.0+1.0e-12) &&
          entropy_ratio>=1.0-1.0e-10 &&
          std::abs(upstream_un)>upstream_fast &&
          std::abs(downstream_un)<=downstream_fast*(1.0+2.0e-10) &&
          (downstream_normal_alfven==0.0 ||
           std::abs(downstream_un)>=downstream_normal_alfven*(1.0-2.0e-10)) &&
          result.evolutionary_fast_branch;
    } else if(result.status==swcme::shock::SolveStatus::NoShock){
      ++no_shock;
      pass=!result.has_shock && result.solver_converged &&
           result.compression==1.0 &&
           result.downstream.rho_kg_m3==upstream.rho_kg_m3 &&
           result.downstream.pressure_Pa==upstream.pressure_Pa &&
           result.downstream.velocity_m_s==upstream.velocity_m_s &&
           result.downstream.magnetic_T==upstream.magnetic_T;
    } else if(result.status==
                  swcme::shock::SolveStatus::NumericallyUnresolvedWeakShock ||
              result.status==swcme::shock::SolveStatus::NumericallySingular){
      ++numerical_limit;
      pass=result.has_shock && !result.solver_converged && finite_jump(result);
    } else {
      ++rejected_generic;
      pass=false;
    }
    if(!pass && failures_reported<8){
      std::cerr<<"    SHK11 admissibility failure "<<reproducer<<'\n';
      ++failures_reported;
    }
    return pass;
  };

  std::size_t campaign_passed=0;
  for(const StressCase& item:shock_stress_cases()){
    if(check_case(item.upstream,item.result,item.normal,item.shock_speed_m_s,
                  item.gamma,stress_reproducer(item))) ++campaign_passed;
  }

  std::size_t reference_passed=0;
  for(const auto& expected:swcme_test::shk05_reference_v1::CASES){
    Primitive upstream;
    upstream.rho_kg_m3=expected.upstream_rho_kg_m3;
    upstream.pressure_Pa=expected.upstream_pressure_Pa;
    for(int component=0;component<3;++component){
      upstream.velocity_m_s[component]=expected.upstream_velocity_m_s[component];
      upstream.magnetic_T[component]=expected.upstream_magnetic_T[component];
    }
    const Vec3 normal={{1.0,0.0,0.0}};
    const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,normal,expected.shock_speed_m_s,expected.gamma);
    const bool metadata=expected.physical_root_count==1 &&
        expected.physical_seed_count>0 && expected.converged_seed_count>1 &&
        std::string(expected.branch)=="EVOLUTIONARY_FAST";
    if(metadata && check_case(upstream,result,normal,expected.shock_speed_m_s,
                              expected.gamma,
                              std::string("fixture=")+expected.id)){
      ++reference_passed;
    }
  }

  std::cout<<"  campaign_passed="<<campaign_passed
           <<" solved="<<solved<<" no_shock="<<no_shock
           <<" numerical_limit="<<numerical_limit
           <<" generic_rejection="<<rejected_generic
           <<" min_entropy_ratio="<<std::scientific<<minimum_entropy_ratio
           <<" max_compression/bound="<<maximum_compression_fraction<<'\n';
  context.expect_true(campaign_passed==shock_stress_cases().size(),
                      "SHK11 every stress outcome satisfies its independent contract");
  context.expect_true(reference_passed>=12,
                      "SHK11 all named reference roots are uniquely evolutionary fast");
  context.expect_true(rejected_generic==0,
                      "SHK11 has no generic or silently substituted rejected root");
  context.expect_true(solved>=49000 && no_shock==25000 && numerical_limit>=25000,
                      "SHK11 exercises solved, no-shock, and explicit numerical limits");
}

void test_shk12(swcme_test::Context& context){
  std::cout<<"SHK12 two-sided near-Mach-one shock limit\n";
  namespace reference=swcme_test::shk12_reference_v1;
  constexpr std::size_t case_count=sizeof(reference::CASES)/
                                   sizeof(reference::CASES[0]);
  const Vec3 n={{1,0,0}};

  // The frozen values are generated with an 80-digit, eight-variable Newton
  // solve and logarithmic continuation.  Pinning the schema and matrix size
  // prevents a future fixture edit from silently reducing the required beta,
  // obliquity, gamma, or Mach-excess coverage.
  context.expect_true(reference::FIXTURE_VERSION==1,
                      "SHK12 fixture schema/version is pinned to v1");
  context.expect_true(reference::DECIMAL_PRECISION_DIGITS>=50,
                      "SHK12 reference uses at least 50 decimal digits");
  context.expect_true(case_count>=48,
                      "SHK12 reference contains eight six-point families");

  double previous_compression=0.0;
  double previous_excess=0.0;
  for(std::size_t index=0;index<case_count;++index){
    const reference::Case& expected=reference::CASES[index];
    Primitive upstream;
    upstream.rho_kg_m3=expected.upstream_rho_kg_m3;
    upstream.pressure_Pa=expected.upstream_pressure_Pa;
    for(int component=0;component<3;++component){
      upstream.velocity_m_s[component]=expected.upstream_velocity_m_s[component];
      upstream.magnetic_T[component]=expected.upstream_magnetic_T[component];
    }
    const auto actual=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,n,expected.shock_speed_m_s,expected.gamma);
    const std::string prefix=std::string("SHK12 ")+expected.id+" ";

    // Reference conditioning and residual metadata are executable acceptance
    // criteria.  A changed fixture therefore cannot hide a poorly converged
    // high-precision solve behind rounded binary64 values.
    context.expect_true(expected.reference_max_residual<1.0e-50,
                        prefix+"reference residual is below 1e-50");
    context.expect_true(expected.minimum_newton_pivot>1.0e-12,
                        prefix+"reference Jacobian remains conditioned");
    context.expect_true(actual.status==swcme::shock::SolveStatus::Solved &&
                            actual.has_shock && actual.solver_converged,
                        prefix+"well-conditioned weak branch is solved");

    // The plan requires the full primitive state, not compression alone.  The
    // 2e-7 component tolerance leaves a small binary64 bracket margin while
    // remaining independent of and much tighter than visually smooth trends.
    expect_rel(context,prefix+"Mach excess",actual.fast_mach-1.0,
               expected.mach_excess,2.0e-7,1.0e-15);
    expect_rel(context,prefix+"compression",actual.compression,
               expected.compression,2.0e-7);
    expect_rel(context,prefix+"downstream density",actual.downstream.rho_kg_m3,
               expected.downstream_rho_kg_m3,2.0e-7,1.0e-30);
    expect_rel(context,prefix+"downstream pressure",actual.downstream.pressure_Pa,
               expected.downstream_pressure_Pa,2.0e-7,1.0e-30);
    for(int component=0;component<3;++component){
      expect_rel(context,prefix+"downstream velocity["+
                     std::to_string(component)+"]",
                 actual.downstream.velocity_m_s[component],
                 expected.downstream_velocity_m_s[component],2.0e-7,1.0);
      expect_rel(context,prefix+"downstream magnetic["+
                     std::to_string(component)+"]",
                 actual.downstream.magnetic_T[component],
                 expected.downstream_magnetic_T[component],2.0e-7,1.0e-15);
    }

    // Published iteration/bracket diagnostics prove scalar convergence and
    // make any future weak-limit regression directly reproducible.
    context.expect_true(actual.root_iterations>0 && actual.root_iterations<=120,
                        prefix+"iteration count is bounded");
    context.expect_true(actual.root_bracket_lower_compression<=actual.compression &&
                            actual.compression<=actual.root_bracket_upper_compression &&
                            actual.root_bracket_width>=0.0 &&
                            actual.root_bracket_width<=3.0e-13*
                                std::max(1.0,actual.compression),
                        prefix+"final bracket contains the reported root");

    // Each family is stored from larger to smaller Mach excess.  Compression
    // must decrease toward one without an empirical floor; the final 1e-5
    // reference is close enough to expose even a small imposed compression.
    if(index%6!=0){
      context.expect_true(expected.mach_excess<previous_excess &&
                              actual.compression<previous_compression,
                          prefix+"compression approaches one monotonically");
    }
    if(index%6==5){
      context.expect_true(actual.compression-1.0<2.0e-5,
                          prefix+"smallest resolved shock approaches identity");
    }
    previous_excess=expected.mach_excess;
    previous_compression=actual.compression;
  }

  auto make_sweep_state=[](double beta,double theta_deg){
    Primitive state;
    const double B=5.0e-9;
    const double theta=theta_deg*swcme::constants::PI/180.0;
    const double phi=25.0*swcme::constants::PI/180.0;
    state.rho_kg_m3=5.0e6*MP;
    state.pressure_Pa=beta*B*B/(2.0*MU0);
    // Zero normal flow makes the exactly critical floating-point fixture
    // unambiguous; nonzero tangential flow still exercises the vector jump.
    state.velocity_m_s={{0.0,2.0e4,-1.0e4}};
    state.magnetic_T={{B*std::cos(theta),
                       B*std::sin(theta)*std::cos(phi),
                       B*std::sin(theta)*std::sin(phi)}};
    return state;
  };

  // Sweep both sides of Mfast=1 down to 1e-12 for every beta/angle/gamma
  // family.  Subcritical and exactly critical states are physical NoShock;
  // small positive excesses are physical shocks that are explicitly marked
  // numerically unresolved, never silently folded into NoShock.
  for(double beta : {0.01,0.1,1.0,10.0}){
    for(double theta : {1.0,30.0,60.0,89.0}){
      for(double gamma : {1.4,5.0/3.0}){
        const Primitive upstream=make_sweep_state(beta,theta);
        const double fast=reference_fast_speed(upstream,n,gamma);
        for(double excess : {-0.5,-0.1,-1.0e-2,-1.0e-4,-1.0e-8,-1.0e-12,0.0}){
          const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
              upstream,n,(1.0+excess)*fast,gamma);
          context.expect_true(result.status==swcme::shock::SolveStatus::NoShock &&
                                  !result.has_shock && result.solver_converged &&
                                  result.compression==1.0,
                              "SHK12 subcritical/critical state is NoShock");
        }
        for(double excess : {1.0e-12,1.0e-10,1.0e-8,1.0e-7,1.0e-6}){
          const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
              upstream,n,(1.0+excess)*fast,gamma);
          context.expect_true(
              result.status==swcme::shock::SolveStatus::NumericallyUnresolvedWeakShock &&
                  result.has_shock && !result.solver_converged &&
                  std::isfinite(result.downstream.rho_kg_m3) &&
                  std::isfinite(result.downstream.pressure_Pa),
              "SHK12 sub-resolution supercritical state is explicitly unresolved");
        }
      }
    }
  }

  // These low-beta, nearly parallel points expose a near-singular tangential
  // system.  The old outermost-root rule returned a discontinuous r~4--6
  // state.  They must retain a finite diagnostic candidate but carry the
  // unresolved status so downstream SEP physics cannot consume that branch.
  for(const auto parameters :
      {std::array<double,3>{{0.01,1.4,1.0e-2}},
       std::array<double,3>{{0.01,5.0/3.0,1.0e-1}}}){
    const Primitive upstream=make_sweep_state(parameters[0],1.0);
    const double fast=reference_fast_speed(upstream,n,parameters[1]);
    const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,n,(1.0+parameters[2])*fast,parameters[1]);
    context.expect_true(
        result.status==swcme::shock::SolveStatus::NumericallyUnresolvedWeakShock &&
            result.has_shock && !result.solver_converged &&
            std::isfinite(result.compression) &&
            std::isfinite(result.downstream.rho_kg_m3) &&
            std::isfinite(result.downstream.pressure_Pa),
        std::string("SHK12 discontinuous near-singular branch is explicitly unresolved; status=")+
            swcme::shock::solve_status_name(result.status)+
            " compression="+std::to_string(result.compression));
  }
}

void test_shk16(swcme_test::Context& context){
  std::cout<<"SHK16 near-singular tangential system\n";
  const Vec3 normal={{1.0,0.0,0.0}};
  Primitive upstream=fixture(5.0,5.0,5.0,2.0e4);

  // Choose the shock-frame inflow from an independently specified singular
  // compression.  For the tangential 2x2 system, det=0 when
  // rho*u1n^2/r=Bn^2/mu0.  Selecting one of the production scan abscissae is
  // deliberate: it proves that a sampled pole splits the residual sequence
  // instead of allowing a false sign change to bridge the discontinuity.
  const double rmax=(GAMMA+1.0)/(GAMMA-1.0);
  const double scan_fraction=0.5;
  const double singular_compression=1.0+1.0e-9+
      (rmax-1.0-1.0e-9)*scan_fraction*scan_fraction;
  const double Bn=dot(upstream.magnetic_T,normal);
  const double inflow=std::sqrt(singular_compression*Bn*Bn/
                                (MU0*upstream.rho_kg_m3));
  const double shock_speed=upstream.velocity_m_s[0]+inflow;

  // Exercise both sides of the pole at logarithmically decreasing offsets.
  // The expected signed conditioning is calculated here from the analytic
  // determinant rather than copied from candidate_for_compression().  The two
  // closest offsets also probe the explicit 1e-12 binary64 safety boundary.
  for(double offset : {1.0e-3,1.0e-6,1.0e-9,1.0e-13}){
    for(double side : {-1.0,1.0}){
      const double compression=singular_compression*(1.0+side*offset);
      const auto candidate=swcme::shock::detail::candidate_for_compression(
          upstream,normal,shock_speed,GAMMA,compression);
      const double inertial=upstream.rho_kg_m3*inflow*inflow/compression;
      const double magnetic=Bn*Bn/MU0;
      const double reference_relative=(inertial-magnetic)/
          std::max({std::abs(inertial),std::abs(magnetic),1.0e-300});
      const std::string label="SHK16 offset="+std::to_string(offset)+
                              (side<0.0?" below":" above");
      expect_abs(context,label+" signed determinant",
                 candidate.tangential_determinant_relative,
                 reference_relative,2.0e-15);
      context.expect_true(
          (candidate.tangential_determinant_relative>0.0)==(side<0.0),
          label+" preserves determinant side");
      context.expect_true(candidate.singular==(offset<=1.0e-12),
                          label+" has stable singular classification");
    }
  }

  // One-ULP perturbations are the most hostile representable changes to the
  // constructed family.  Every outcome must remain explicitly classified,
  // finite, and, if a physical root is returned, on the evolutionary fast
  // branch with a bracket confined to one side of the determinant pole.
  const std::array<double,3> speeds={{
      std::nextafter(shock_speed,-std::numeric_limits<double>::infinity()),
      shock_speed,
      std::nextafter(shock_speed,std::numeric_limits<double>::infinity())}};
  swcme::shock::SolveStatus baseline_status=swcme::shock::SolveStatus::InvalidInput;
  for(std::size_t index=0;index<speeds.size();++index){
    const auto result=swcme::shock::solve_ideal_mhd_fast_shock(
        upstream,normal,speeds[index],GAMMA);
    if(index==0) baseline_status=result.status;
    context.expect_true(result.status==baseline_status,
                        "SHK16 one-ULP perturbations retain status");
    context.expect_true(result.status==swcme::shock::SolveStatus::Solved ||
                            result.status==swcme::shock::SolveStatus::NumericallySingular,
                        "SHK16 result is solved or explicitly singular");
    context.expect_true(result.encountered_tangential_singularity,
                        "SHK16 sampled pole is retained in diagnostics");
    context.expect_true(std::isfinite(result.minimum_tangential_determinant_relative) &&
                            std::isfinite(result.closest_tangential_determinant_relative) &&
                            std::isfinite(result.root_bracket_lower_compression) &&
                            std::isfinite(result.root_bracket_upper_compression),
                        "SHK16 diagnostic payload remains finite");
    if(result.status==swcme::shock::SolveStatus::Solved){
      context.expect_true(result.evolutionary_fast_branch,
                          "SHK16 accepted result is evolutionary fast");
      context.expect_true(result.root_bracket_upper_compression<singular_compression ||
                              result.root_bracket_lower_compression>singular_compression,
                          "SHK16 root bracket never crosses singular pole");
    }
  }
}

void test_shk15(swcme_test::Context& context){
  std::cout<<"SHK15 high-count deterministic random shock stress\n";
  const auto& cases=shock_stress_cases();
  std::array<std::size_t,9> status_counts{{0,0,0,0,0,0,0,0,0}};
  std::size_t finite_count=0;
  std::size_t expected_classification_count=0;
  std::size_t conditioned_count=0;
  std::size_t reported_failures=0;

  for(const StressCase& item:cases){
    const auto status_index=static_cast<std::size_t>(item.result.status);
    if(status_index<status_counts.size()) ++status_counts[status_index];
    if(finite_jump(item.result)) ++finite_count;
    if(item.result.encountered_tangential_singularity) ++conditioned_count;

    bool expected=false;
    if(item.stratum==StressStratum::NoShock){
      expected=item.result.status==swcme::shock::SolveStatus::NoShock &&
          !item.result.has_shock && item.result.solver_converged;
    } else if(item.stratum==StressStratum::WeakLimit){
      expected=item.result.status==
          swcme::shock::SolveStatus::NumericallyUnresolvedWeakShock &&
          item.result.has_shock && !item.result.solver_converged;
    } else {
      // Resolved and determinant-conditioned inputs may either produce the
      // verified evolutionary fast branch or one of the two documented
      // supported-limit rejections.  Generic bracket, reconstruction,
      // conservation, and wrong-branch failures remain forbidden.
      expected=(item.result.status==swcme::shock::SolveStatus::Solved &&
                    item.result.has_shock && item.result.solver_converged &&
                    item.result.evolutionary_fast_branch) ||
          ((item.result.status==swcme::shock::SolveStatus::NumericallySingular ||
            item.result.status==
                swcme::shock::SolveStatus::NumericallyUnresolvedWeakShock) &&
                    item.result.has_shock && !item.result.solver_converged);
    }
    if(expected) ++expected_classification_count;
    else if(reported_failures<8){
      // A complete, fixed-seed primitive record is a directly reusable
      // reproducer.  Limiting diagnostics prevents a systemic regression from
      // flooding CI while preserving the first examples for minimization and
      // independent high-precision analysis.
      std::cerr<<"    SHK15 unexpected case: "<<stress_reproducer(item)<<'\n';
      ++reported_failures;
    }
  }

  std::cout<<"  seed=0x53484b31355f7631 cases="<<cases.size()
           <<" solved="
           <<status_counts[static_cast<std::size_t>(swcme::shock::SolveStatus::Solved)]
           <<" no_shock="
           <<status_counts[static_cast<std::size_t>(swcme::shock::SolveStatus::NoShock)]
           <<" weak_limit="
           <<status_counts[static_cast<std::size_t>(
                  swcme::shock::SolveStatus::NumericallyUnresolvedWeakShock)]
           <<" singular="
           <<status_counts[static_cast<std::size_t>(
                  swcme::shock::SolveStatus::NumericallySingular)]
           <<" conditioned_trials="<<conditioned_count<<'\n';
  context.expect_true(cases.size()>=100000,
                      "SHK15 executes at least 100,000 physical inputs");
  context.expect_true(finite_count==cases.size(),
                      "SHK15 emits no NaN or infinite diagnostic/primitive values");
  context.expect_true(expected_classification_count==cases.size(),
                      "SHK15 every input has an expected physical or supported-limit status");
  context.expect_true(conditioned_count>=90,
                      "SHK15 includes at least 90 sampled singular systems");
  context.expect_true(
      status_counts[static_cast<std::size_t>(swcme::shock::SolveStatus::NoPhysicalBracket)]==0 &&
      status_counts[static_cast<std::size_t>(swcme::shock::SolveStatus::InvalidAcceptedState)]==0 &&
      status_counts[static_cast<std::size_t>(swcme::shock::SolveStatus::ConservationFailure)]==0 &&
      status_counts[static_cast<std::size_t>(swcme::shock::SolveStatus::WrongBranch)]==0,
      "SHK15 has zero unclassified bracket, reconstruction, conservation, or branch failures");
}
