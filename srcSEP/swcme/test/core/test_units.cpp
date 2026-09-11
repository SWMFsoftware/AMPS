#include "test_framework.hpp"

#include <swcme1d.hpp>
#include <swcme3d.hpp>
#include <swcme_constants.hpp>
#include <swcme_units.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

constexpr double REF_AU_M = 149597870700.0;
constexpr double REF_RS_M = 695700000.0;
constexpr double REF_PI = 3.141592653589793238462643383279502884;
constexpr double ROUNDOFF_REL_TOL =
    64.0 * std::numeric_limits<double>::epsilon();

struct Metrics {
  double max_forward_error = 0.0;
  double max_roundtrip_error = 0.0;
  double max_model_difference = 0.0;
  double alfven_error = 0.0;
};

double relerr(double a, double b) {
  const double e=std::abs(a-b);
  return b==0.0 ? e : e/std::abs(b);
}

void check(swcme_test::Context& context, Metrics& metrics,
           const std::string& name, double value, double reference,
           double tolerance=ROUNDOFF_REL_TOL) {
  const double error=relerr(value,reference);
  const bool pass=std::isfinite(value) && error<=tolerance;
  metrics.max_forward_error=std::max(metrics.max_forward_error,error);
  std::cout << std::left << std::setw(35) << name << std::right
            << " value=" << std::scientific << std::setprecision(12) << value
            << " ref=" << reference << " rel_err=" << error
            << " tol=" << tolerance << ' ' << (pass?"PASS":"FAIL") << '\n';
  context.record_result(pass);
}

void check_roundtrip(swcme_test::Context& context, Metrics& metrics,
                     const std::string& name, double original,
                     double recovered) {
  const double error=relerr(recovered,original);
  const bool pass=std::isfinite(recovered) && error<=ROUNDOFF_REL_TOL;
  metrics.max_roundtrip_error=std::max(metrics.max_roundtrip_error,error);
  std::cout << std::left << std::setw(35) << name << std::right
            << " original=" << std::scientific << std::setprecision(12)
            << original << " recovered=" << recovered << " rel_err=" << error
            << ' ' << (pass?"PASS":"FAIL") << '\n';
  context.record_result(pass);
}

void check_model_pair(swcme_test::Context& context, Metrics& metrics,
                      const std::string& name, double a, double b) {
  const double error=relerr(a,b);
  const bool pass=std::isfinite(a) && std::isfinite(b) &&
                  error<=ROUNDOFF_REL_TOL;
  metrics.max_model_difference=std::max(metrics.max_model_difference,error);
  std::cout << std::left << std::setw(35) << name << std::right
            << " 1D=" << std::scientific << std::setprecision(12) << a
            << " 3D=" << b << " rel_diff=" << error << ' '
            << (pass?"PASS":"FAIL") << '\n';
  context.record_result(pass);
}

struct Ambient {
  double density_m3=0.0;
  double B_T=0.0;
};

Ambient ambient_1d(double n_cm3, double B_nT) {
  swcme1d::Params p; p.n1AU_cm3=n_cm3; p.B1AU_nT=B_nT;
  const swcme1d::Model model(p);
  const auto S=model.prepare_step(0.0);
  const double r=swcme::constants::AU_M;
  double n=0,V=0,Br=0,Bp=0,Bm=0;
  model.evaluate_radii_with_B_div(S,&r,&n,&V,&Br,&Bp,&Bm,nullptr,1);
  return {n,Bm};
}

Ambient ambient_3d(double n_cm3, double B_nT) {
  swcme3d::Params p; p.n1AU_cm3=n_cm3; p.B1AU_nT=B_nT;
  const swcme3d::Model model(p);
  const auto S=model.prepare_step(0.0);
  const double x=swcme::constants::AU_M,y=0,z=0;
  double n=0,vx=0,vy=0,vz=0,bx=0,by=0,bz=0;
  model.evaluate_cartesian_with_B(S,&x,&y,&z,&n,&vx,&vy,&vz,&bx,&by,&bz,1);
  return {n,std::sqrt(bx*bx+by*by+bz*bz)};
}

}  // namespace

void test_cfg02(swcme_test::Context& context) {
  Metrics metrics;
  std::cout << "CFG02 centralized forward conversions\n";

  const double velocities[]={0.0,1.0,321.0,400.0,1500.0,3000.0};
  for (double x:velocities) {
    check(context,metrics,"km/s -> m/s",
          swcme::units::km_per_s_to_m_per_s(x),x*1000.0,0.0);
  }
  const double fields[]={0.0,0.1,1.0,5.0,1000.0};
  for (double x:fields) {
    check(context,metrics,"nT -> T",swcme::units::nT_to_T(x),x*1.0e-9);
  }
  const double densities[]={0.0,0.01,1.0,5.0,1.0e4};
  for (double x:densities) {
    check(context,metrics,"cm^-3 -> m^-3",swcme::units::cm3_to_m3(x),x*1.0e6);
  }
  const double au_values[]={0.0,0.1,0.5,1.0,2.0,5.0};
  for (double x:au_values) {
    check(context,metrics,"AU -> m",swcme::units::au_to_m(x),x*REF_AU_M);
  }
  const double rs_values[]={0.0,1.0,2.0,20.0,REF_AU_M/REF_RS_M};
  for (double x:rs_values) {
    check(context,metrics,"Rs -> m",swcme::units::solar_radii_to_m(x),x*REF_RS_M);
  }
  const double hours[]={0.0,0.5,1.0,24.0,72.0};
  for (double x:hours) {
    check(context,metrics,"hour -> s",swcme::units::hours_to_seconds(x),x*3600.0,0.0);
  }
  const double degrees[]={0.0,40.0,90.0,180.0};
  for (double x:degrees) {
    check(context,metrics,"degree -> rad",swcme::units::degrees_to_radians(x),
          x*REF_PI/180.0);
  }
  const double inverse_km[]={0.0,1.0e-8,2.5e-5};
  for (double x:inverse_km) {
    check(context,metrics,"km^-1 -> m^-1",
          swcme::units::km_inverse_to_m_inverse(x),x*1.0e-3);
  }

  std::cout << "CFG02 centralized round-trip conversions\n";
  check_roundtrip(context,metrics,"velocity round trip",321.0,
      swcme::units::m_per_s_to_km_per_s(swcme::units::km_per_s_to_m_per_s(321.0)));
  check_roundtrip(context,metrics,"field round trip",5.0,
      swcme::units::T_to_nT(swcme::units::nT_to_T(5.0)));
  check_roundtrip(context,metrics,"density round trip",8.0,
      swcme::units::m3_to_cm3(swcme::units::cm3_to_m3(8.0)));
  check_roundtrip(context,metrics,"AU round trip",0.5,
      swcme::units::m_to_au(swcme::units::au_to_m(0.5)));
  check_roundtrip(context,metrics,"Rs round trip",20.0,
      swcme::units::m_to_solar_radii(swcme::units::solar_radii_to_m(20.0)));
  check_roundtrip(context,metrics,"time round trip",24.0,
      swcme::units::seconds_to_hours(swcme::units::hours_to_seconds(24.0)));
  check_roundtrip(context,metrics,"angle round trip",40.0,
      swcme::units::radians_to_degrees(swcme::units::degrees_to_radians(40.0)));
  check_roundtrip(context,metrics,"inverse length round trip",2.5e-5,
      swcme::units::m_inverse_to_km_inverse(
          swcme::units::km_inverse_to_m_inverse(2.5e-5)));

  std::cout << "CFG02 1-D/3-D prepared-state consistency\n";
  // Only physically valid values enter prepare_step().  CFG01 separately
  // verifies that zero/negative invalid model inputs are rejected rather than
  // silently changed by unit conversion.
  for (double v : {1.0,321.0,400.0,1500.0,3000.0}) {
    swcme1d::Params p1; p1.V_sw_kms=v; p1.V0_sh_kms=std::max(1500.0,v+1.0);
    swcme3d::Params p3; p3.V_sw_kms=v; p3.V0_sh_kms=std::max(1500.0,v+1.0);
    check_model_pair(context,metrics,"prepared V_sw",
                     swcme1d::Model(p1).prepare_step(0.0).V_up_ms,
                     swcme3d::Model(p3).prepare_step(0.0).V_sw_ms);
  }
  for (double n : {0.01,1.0,5.0,8.0,100.0,1.0e4}) {
    check_model_pair(context,metrics,"ambient density",ambient_1d(n,5.0).density_m3,
                     ambient_3d(n,5.0).density_m3);
  }
  for (double b : {0.1,1.0,5.0,7.0,100.0,1000.0}) {
    check_model_pair(context,metrics,"ambient |B|",ambient_1d(5.0,b).B_T,
                     ambient_3d(5.0,b).B_T);
  }

  std::cout << "CFG02 Alfven-speed dimensional smoke test\n";
  const Ambient a1=ambient_1d(5.0,5.0);
  const Ambient a3=ambient_3d(5.0,5.0);
  const double rho_ref=swcme::units::cm3_to_m3(5.0)*
                       swcme::constants::PROTON_MASS_KG;
  const double B_ref=swcme::units::nT_to_T(5.0);
  const double va_ref=B_ref/std::sqrt(
      swcme::constants::VACUUM_PERMEABILITY_N_A2*rho_ref);
  const double va1=swcme::physics::alfven_speed_m_s(
      a1.B_T,a1.density_m3*swcme::constants::PROTON_MASS_KG);
  const double va3=swcme::physics::alfven_speed_m_s(
      a3.B_T,a3.density_m3*swcme::constants::PROTON_MASS_KG);
  check(context,metrics,"1D Alfven speed",va1,va_ref);
  check(context,metrics,"3D Alfven speed",va3,va_ref);
  check_model_pair(context,metrics,"1D/3D Alfven speed",va1,va3);
  metrics.alfven_error=std::max(relerr(va1,va_ref),relerr(va3,va_ref));

  std::cout << "CFG02 metrics"
            << " max_forward=" << std::scientific << metrics.max_forward_error
            << " max_roundtrip=" << metrics.max_roundtrip_error
            << " max_1d3d=" << metrics.max_model_difference
            << " alfven=" << metrics.alfven_error << '\n';
}
