// Stage12 core verification of the explicit D10-1 upstream material profile.
// Every actual value comes from the new canonical Profile/Event public API.
// Expected values use the separate long-double, test-only oracle header. No
// production derivative/map/force/residual helper computes its own reference.
#include "material_upstream_oracle.hpp"
#include "swcme_material_upstream.hpp"
#include "swcme1d.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace P=swcme::material1d;
namespace O=swcme_material_oracle;
using L=O::Scalar;
struct Metric {L value=0,limit=0;std::string relation="<=",units="dimensionless";bool initialized=false;};
struct Checks {
  unsigned checks=0,failures=0;
  std::map<std::string,Metric> metrics;
  void require(bool condition,const std::string& reason) {
    ++checks;if(!condition){++failures;std::cerr<<"ASSERT "<<reason<<'\n';}
  }
  bool status(const swcme::ModelStatus& status,const std::string& reason) {
    require(status.ok(),reason+": "+status.summary());return status.ok();
  }
  void upper(const std::string& name,L value,L limit,const std::string& units="dimensionless") {
    auto& m=metrics[name];m.value=m.initialized?std::max(m.value,value):value;
    m.initialized=true;m.limit=limit;m.units=units;
    require(std::isfinite(value)&&value<=limit,name);
  }
  void lower(const std::string& name,L value,L limit) {
    auto& m=metrics[name];m.value=m.initialized?std::min(m.value,value):value;
    m.initialized=true;m.limit=limit;m.relation=">";
    require(std::isfinite(value)&&value>limit,name);
  }
  void strict_upper(const std::string& name,L value,L limit) {
    auto& m=metrics[name];m.value=m.initialized?std::max(m.value,value):value;
    m.initialized=true;m.limit=limit;m.relation="<";
    require(std::isfinite(value)&&value<limit,name);
  }
  void near(const std::string& name,L actual,L reference,L physicalScale) {
    // Never normalize a tiny SI quantity with max(1,...). Zero tensor entries
    // use the WHOLE physical field/r scale declared before the run, so a T/m
    // derivative cannot pass merely because its dimensionful value is small.
    const L error=physicalScale>0?std::abs(actual-reference)/physicalScale:
        (actual==reference?0:std::numeric_limits<L>::infinity());
    upper(name,error,O::ExactTolerance);
  }
  void report() const {
    for(const auto& row:metrics)
      std::cout<<"MATERIAL_METRIC "<<row.first<<' '<<std::scientific<<std::setprecision(18)
               <<row.second.value<<' '<<row.second.relation<<' '<<row.second.limit<<' '<<row.second.units<<'\n';
  }
};

P::Config config(int polarity=1,P::Vector axis={{0,0,1}}) {
  P::Config c;c.inner_radius_m=0.30*static_cast<double>(O::AU);c.outer_radius_m=1.20*static_cast<double>(O::AU);
  c.reference_radius_m=static_cast<double>(O::AU);c.wind_m_s=static_cast<double>(O::Speed);
  c.reference_electron_density_m3=static_cast<double>(O::ElectronDensity);
  c.reference_electron_temperature_k=static_cast<double>(O::ElectronTemperature);
  c.reference_proton_temperature_k=static_cast<double>(O::ProtonTemperature);
  c.reference_field_t=static_cast<double>(O::FieldMagnitude);c.reference_sin_colatitude=1;
  c.polarity=polarity;c.rotation_axis=axis;c.rotation_rad_s=static_cast<double>(O::Rotation);
  c.solar_gravity_m3_s2=static_cast<double>(O::Gravity);c.gamma=static_cast<double>(O::Gamma);return c;
}
O::Vector wide(const P::Vector& v) {return {{v[0],v[1],v[2]}};}
P::Vector narrow(const O::Vector& v) {return {{static_cast<double>(v[0]),static_cast<double>(v[1]),static_cast<double>(v[2])}};}
O::Vector unit(const P::Vector& axis) {const auto v=wide(axis);return O::scaled(1/O::norm(v),v);}
P::Vector position(L radius,const O::Vector& direction) {return narrow(O::scaled(radius,direction));}
const std::vector<O::Vector> directions={{{1,0,0}},{{0,1,0}},{{0,0,1}},{{-1,0,0}},
    {{0.36L,-0.48L,0.80L}},{{0.36L,-0.48L,-0.80L}},{{-0.60L,0,0.80L}}};
void vectors(Checks& q,const std::string& name,const P::Vector& actual,const O::Vector& expected,L scale) {
  for(unsigned i=0;i<3;++i)q.near(name,actual[i],expected[i],scale);
}
void primitives(Checks& q,const P::Sample& a,const O::State& e) {
  q.near("electron_density",a.electron_density_m3,e.electrons,e.electrons);
  q.near("proton_density",a.proton_density_m3,e.protons,e.protons);
  q.near("mass_density",a.mass_density_kg_m3,e.density,e.density);
  q.near("electron_temperature",a.electron_temperature_k,e.electronTemperature,e.electronTemperature);
  q.near("proton_temperature",a.proton_temperature_k,e.protonTemperature,e.protonTemperature);
  q.near("pressure",a.pressure_pa,e.pressure,e.pressure);
  q.near("sound_speed",a.sound_speed_m_s,e.soundSpeed,e.soundSpeed);
  vectors(q,"position",a.position_m,e.position,e.radius);
  vectors(q,"velocity",a.velocity_m_s,e.velocity,O::Speed);
  vectors(q,"magnetic_field",a.magnetic_field_t,e.magnetic,O::norm(e.magnetic));
}

template<std::size_t N> void append(std::vector<double>& out,const std::array<double,N>& value) {
  out.insert(out.end(),value.begin(),value.end());
}
std::vector<double> image(const P::Sample& s) {
  std::vector<double> a;for(const auto* v:{&s.position_m,&s.velocity_m_s,&s.magnetic_field_t,
      &s.pressure_gradient_pa_m_inv,&s.mass_density_gradient_kg_m4,&s.electron_temperature_gradient_k_m,
      &s.proton_temperature_gradient_k_m,&s.external_force_n_m3})append(a,*v);
  append(a,s.velocity_gradient_s_inv);append(a,s.magnetic_gradient_t_m_inv);
  for(double x:{s.electron_density_m3,s.proton_density_m3,s.electron_temperature_k,s.proton_temperature_k,
      s.mass_density_kg_m3,s.pressure_pa,s.sound_speed_m_s,s.external_work_w_m3,s.heat_w_m3})a.push_back(x);
  return a;
}
std::vector<double> image(const P::Motion& s) {
  std::vector<double> a;append(a,s.position_m);append(a,s.velocity_m_s);append(a,s.deformation);a.push_back(s.jacobian);return a;
}
P::Sample sentinel() {
  P::Sample s;for(auto* v:{&s.position_m,&s.velocity_m_s,&s.magnetic_field_t,&s.pressure_gradient_pa_m_inv,
      &s.mass_density_gradient_kg_m4,&s.electron_temperature_gradient_k_m,&s.proton_temperature_gradient_k_m,
      &s.external_force_n_m3})v->fill(-17);
  s.velocity_gradient_s_inv.fill(-19);s.magnetic_gradient_t_m_inv.fill(-23);
  s.electron_density_m3=s.proton_density_m3=s.electron_temperature_k=s.proton_temperature_k=
      s.mass_density_kg_m3=s.pressure_pa=s.sound_speed_m_s=s.external_work_w_m3=s.heat_w_m3=-29;return s;
}
P::Motion motion_sentinel() {P::Motion m;m.position_m.fill(-17);m.velocity_m_s.fill(-19);m.deformation.fill(-23);m.jacobian=-29;return m;}

// MAT1D01 — Check the approved SI anchors through actual public samples. Latitude changes only
// the query geometry: the declared equatorial field normalization remains fixed, and both radial
// polarities retain their signed components.
void test01(Checks& q) {
  for(int polarity:{1,-1})for(const auto axis:{P::Vector{{0,0,1}},P::Vector{{2,-1,2}}}) {
    const auto c=config(polarity,axis);const P::Profile p(c);const auto k=unit(axis);
    for(L r:{0.30L,0.50L,1.00L,1.20L})for(std::size_t i=0;i<directions.size();++i) {
      // Closed endpoints use exactly representable Cartesian axial radii.
      // Interior nonaxial points exercise every component without pretending
      // a rounded radius outside the closed interval is still covered.
      if((r==0.30L||r==1.20L)&&i>=4)continue;
      const auto x=position(r*O::AU,directions[i]);P::Sample actual;
      if(q.status(p.sample(x,actual),"sample approved primitive"))primitives(q,actual,O::evaluate(wide(x),polarity,k));
    }
  }
}

// MAT1D02 — A radial shell has rho*r^2 constant. The independently integrated inventory and
// signed radial mass flux must close before and after a complete300s material translation; no
// particle-domain clipping is allowed.
void test02(Checks& q) {
  const P::Profile p(config());const L inner=.35L*O::AU,outer=1.10L*O::AU;
  constexpr unsigned intervals=64;L inventory=0;
  for(unsigned i=0;i<intervals;++i) {
    const L r=inner+(outer-inner)*(static_cast<L>(i)+.5L)/intervals;
    P::Sample s;if(!q.status(p.sample(position(r,{{1,0,0}}),s),"mass inventory sample"))continue;
    inventory+=s.mass_density_kg_m3*r*r*(outer-inner)/intervals;
    const L flux=s.mass_density_kg_m3*O::dot(wide(s.velocity_m_s),{{1,0,0}})*r*r;
    q.upper("mass_flux",std::abs(flux/O::massFluxPerSolidAngle()-1),O::BudgetTolerance);
  }
  q.upper("shell_mass_inventory",std::abs(inventory/O::massPerSolidAngle(inner,outer)-1),O::BudgetTolerance);
  // Constant U moves both material shell boundaries by the same distance.
  // Its independent inventory therefore remains constant; there is no hidden
  // clipping at the provider boundary or replacement by an Eulerian cell mass.
  const L dt=300,shift=O::Speed*dt;L moved=0;
  for(unsigned i=0;i<intervals;++i) {
    const L r=inner+shift+(outer-inner)*(static_cast<L>(i)+.5L)/intervals;
    P::Sample s;if(q.status(p.sample(position(r,{{1,0,0}}),s),"moved shell mass"))
      moved+=s.mass_density_kg_m3*r*r*(outer-inner)/intervals;
  }
  q.upper("material_shell_mass",std::abs((moved-inventory)/inventory),O::BudgetTolerance);
}

// MAT1D03 — Check all nine entries of the actual deformation, its determinant and inverse, then
// independent mass and signed Cauchy pushes. A scalar radial displacement alone cannot qualify a
// curved material tube.
void test03(Checks& q) {
  for(int sign:{1,-1}) {
    const P::Profile p(config(sign));
    for(const auto& e:directions)for(L r:{.40L,.80L})for(double elapsed:{0.0,300.0,7200.0}) {
      const auto x=position(r*O::AU,e);const auto expected=O::move(wide(x),elapsed);
      P::Motion moved;if(!q.status(p.move(x,elapsed,moved),"exact material move"))continue;
      vectors(q,"material_position",moved.position_m,expected.position,O::norm(expected.position));
      vectors(q,"material_velocity",moved.velocity_m_s,O::scaled(O::Speed,O::scaled(1/O::norm(wide(x)),wide(x))),O::Speed);
      q.near("material_jacobian",moved.jacobian,expected.jacobian,expected.jacobian);
      for(unsigned i=0;i<9;++i)q.near("material_deformation",moved.deformation[i],expected.deformation[i],std::max(1.0L,std::sqrt(expected.jacobian)));
      P::Sample a,b;if(!q.status(p.sample(x,a),"departure sample")||!q.status(p.sample(moved.position_m,b),"target sample"))continue;
      q.upper("material_mass_push",std::abs(b.mass_density_kg_m3*expected.jacobian/a.mass_density_kg_m3-1),O::BudgetTolerance);
      const auto cauchy=O::scaled(1/expected.jacobian,O::product(expected.deformation,wide(a.magnetic_field_t)));
      q.upper("signed_cauchy_field",O::norm(O::subtract(wide(b.magnetic_field_t),cauchy))/O::norm(cauchy),O::BudgetTolerance);
      P::Motion inverse;if(q.status(p.move(moved.position_m,-elapsed,inverse),"inverse material map")) {
        vectors(q,"material_inverse",inverse.position_m,wide(x),O::norm(wide(x)));
        for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j) {
          L product=0;for(unsigned k=0;k<3;++k)product+=inverse.deformation[3*i+k]*moved.deformation[3*k+j];
          q.near("deformation_inverse",product,i==j?1:0,1);
        }
      }
    }
  }
}

// MAT1D04 — The stationary Parker field obeys divB=0 and U.gradB=(B.grad)U-B divU. Rebuild every
// component from exported tensors, and separately preserve signed radial/winding flux and
// polarity-even force.
void test04(Checks& q) {
  for(const auto& e:directions)for(L r:{.35L,.60L,1.10L}) {
    P::Sample a,b;const P::Profile plus(config()),minus(config(-1));const auto x=position(r*O::AU,e);
    if(!q.status(plus.sample(x,a),"positive field")||!q.status(minus.sample(x,b),"negative field"))continue;
    const auto reference=O::evaluate(wide(x),1,{{0,0,1}});const auto er=reference.radial;
    const L scale=O::norm(reference.magnetic)/reference.radius;
    L divB=0,divU=0;for(unsigned i=0;i<3;++i){divB+=a.magnetic_gradient_t_m_inv[3*i+i];divU+=a.velocity_gradient_s_inv[3*i+i];}
    q.upper("div_B",std::abs(divB)/scale,O::BudgetTolerance);
    q.near("div_U",divU,reference.divergence,O::Speed/reference.radius);
    for(unsigned i=0;i<3;++i) {
      L induction=divU*a.magnetic_field_t[i];
      for(unsigned j=0;j<3;++j)induction+=a.velocity_m_s[j]*a.magnetic_gradient_t_m_inv[3*i+j]-a.velocity_gradient_s_inv[3*i+j]*a.magnetic_field_t[j];
      q.upper("ideal_induction",std::abs(induction)/(O::Speed*scale),O::BudgetTolerance);
      q.near("signed_polarity",a.magnetic_field_t[i],-b.magnetic_field_t[i],O::norm(reference.magnetic));
      q.near("polarity_even_force",a.external_force_n_m3[i],b.external_force_n_m3[i],O::norm(reference.force));
    }
    const L br=O::dot(wide(a.magnetic_field_t),er);
    const L brRef=O::FieldMagnitude/std::sqrt(1+std::pow(O::Rotation*O::AU/O::Speed,2));
    q.near("signed_radial_flux",reference.radius*reference.radius*br,brRef*O::AU*O::AU,std::abs(brRef*O::AU*O::AU));
    const auto transverse=O::subtract(wide(a.magnetic_field_t),O::scaled(br,er));
    const auto expected=O::scaled(-brRef*O::AU*O::AU*O::Rotation,O::cross({{0,0,1}},er));
    const auto wound=O::scaled(reference.radius*O::Speed,transverse);
    q.upper("signed_winding_flux",O::norm(O::subtract(wound,expected))/(O::FieldMagnitude*O::AU*O::Speed),O::BudgetTolerance);
  }
}

// MAT1D05 — Verify all Cartesian tensor/scalar-gradient entries, including exact polar limits.
// Independent centered samples must retain the original DIV03 second-order refinement gate;
// neither production gradient helper supplies the oracle.
void test05(Checks& q) {
  for(const auto axis:{P::Vector{{0,0,1}},P::Vector{{2,-1,2}}}) {
    const P::Profile p(config(1,axis));
    auto gradientDirections=directions;
    gradientDirections.push_back(unit(axis));gradientDirections.push_back(O::scaled(-1,unit(axis)));
    // At the exact rotation axis B_phi itself vanishes, but its Cartesian
    // transverse derivatives do not. Both poles therefore need full-tensor
    // checks; checking field magnitude alone would miss that regular limit.
    for(const auto& e:gradientDirections) {
      const auto x=position(.70L*O::AU,e);P::Sample s;if(!q.status(p.sample(x,s),"gradient sample"))continue;
      const auto reference=O::evaluate(wide(x),1,unit(axis));
      const L uScale=O::Speed/reference.radius,bScale=O::norm(reference.magnetic)/reference.radius;
      for(unsigned i=0;i<9;++i) {
        q.near("velocity_tensor",s.velocity_gradient_s_inv[i],reference.velocityGradient[i],uScale);
        q.near("magnetic_tensor",s.magnetic_gradient_t_m_inv[i],reference.magneticGradient[i],bScale);
      }
      vectors(q,"pressure_gradient",s.pressure_gradient_pa_m_inv,reference.pressureGradient,reference.pressure/reference.radius);
      vectors(q,"density_gradient",s.mass_density_gradient_kg_m4,O::scaled(-2*reference.density/reference.radius,reference.radial),reference.density/reference.radius);
      vectors(q,"electron_temperature_gradient",s.electron_temperature_gradient_k_m,O::scaled(-4*reference.electronTemperature/(3*reference.radius),reference.radial),reference.electronTemperature/reference.radius);
      vectors(q,"proton_temperature_gradient",s.proton_temperature_gradient_k_m,O::scaled(-4*reference.protonTemperature/(3*reference.radius),reference.radial),reference.protonTemperature/reference.radius);
    }
  }
  // Inherit DIV03's 1.90<order<2.10 gate. Fixed h/r powers are chosen before
  // candidate output, well above roundoff and strictly inside the coverage.
  // Derivatives use actual exported samples; expectations are the independent
  // analytical tensor, never the production gradient under test.
  const P::Vector axis{{2,-1,2}};const P::Profile p(config(1,axis));
  const auto x=position(.70L*O::AU,{{.36L,-.48L,.80L}});const auto reference=O::evaluate(wide(x),1,unit(axis));
  L previousU=0,previousB=0,previousP=0;
  for(int level:{8,9,10}) {
    const double h=std::ldexp(static_cast<double>(reference.radius),-level);
    L errorU=0,errorB=0,errorP=0;
    for(unsigned j=0;j<3;++j) {
      auto left=x,right=x;left[j]-=h;right[j]+=h;P::Sample a,b;
      if(!q.status(p.sample(left,a),"left finite difference")||!q.status(p.sample(right,b),"right finite difference"))continue;
      for(unsigned i=0;i<3;++i) {
        const L du=(static_cast<L>(b.velocity_m_s[i])-a.velocity_m_s[i])/(2*h)-reference.velocityGradient[3*i+j];
        const L db=(static_cast<L>(b.magnetic_field_t[i])-a.magnetic_field_t[i])/(2*h)-reference.magneticGradient[3*i+j];
        errorU+=du*du;errorB+=db*db;
      }
      const L dp=(static_cast<L>(b.pressure_pa)-a.pressure_pa)/(2*h)-reference.pressureGradient[j];errorP+=dp*dp;
    }
    errorU=std::sqrt(errorU);errorB=std::sqrt(errorB);errorP=std::sqrt(errorP);
    if(level>8) {
      for(const auto& row:{std::pair<const char*,std::pair<L,L>>{"velocity",{previousU,errorU}},
          {"magnetic",{previousB,errorB}},{"pressure",{previousP,errorP}}}) {
        const L order=std::log(row.second.first/row.second.second)/std::log(2.0L);
        q.lower(std::string(row.first)+"_derivative_order_lower",order,O::MinimumDerivativeOrder);
        q.strict_upper(std::string(row.first)+"_derivative_order_upper",order,O::MaximumDerivativeOrder);
      }
    }
    previousU=errorU;previousB=errorB;previousP=errorP;
  }
}

// MAT1D06 — Entropy is advected with exactly zero heat while external mechanical work remains
// separate. Independently integrate exported work and gravity and reconstruct total energy flux,
// including Poynting transport.
void test06(Checks& q) {
  const P::Profile p(config());const L a=.35L*O::AU,b=1.10L*O::AU;
  const auto entropyReference=O::entropy(O::evaluate({{O::AU,0,0}},1,{{0,0,1}}));
  for(L radius:{.30L,.50L,1.00L,1.20L}) {
    P::Sample s;if(!q.status(p.sample(position(radius*O::AU,{{1,0,0}}),s),"entropy sample"))continue;
    q.upper("entropy",std::abs(s.pressure_pa/std::pow(static_cast<L>(s.mass_density_kg_m3),O::Gamma)/entropyReference-1),O::BudgetTolerance);
    q.require(s.heat_w_m3==0.0,"thermal heat input is exactly zero, distinct from mechanical work");
  }
  // Composite Simpson with a fixed8192-panel grid. For the two known smooth
  // integrands r^-2 and r^-7/3, its a-priori fourth-derivative remainder is
  // below the frozen1e-10 budget over this interval. No candidate residual
  // selects the panel count. Integrate ACTUAL exported work and gravity power.
  constexpr unsigned panels=8192;const L h=(b-a)/panels;L external=0,gravity=0;
  for(unsigned i=0;i<=panels;++i) {
    const L r=a+i*h;P::Sample s;if(!q.status(p.sample(position(r,{{1,0,0}}),s),"shell work sample"))continue;
    const L weight=(i==0||i==panels)?1:(i%2?4:2);
    external+=weight*r*r*s.external_work_w_m3;
    gravity+=weight*(-s.mass_density_kg_m3*O::Speed*O::Gravity);
  }
  external*=h/3;gravity*=h/3;
  P::Sample first,last;q.status(p.sample(position(a,{{1,0,0}}),first),"inner energy flux");q.status(p.sample(position(b,{{1,0,0}}),last),"outer energy flux");
  auto flux=[](const P::Sample& s) {
    const auto x=wide(s.position_m);const L r=O::norm(x);const auto er=O::scaled(1/r,x);const auto u=wide(s.velocity_m_s),mag=wide(s.magnetic_field_t);
    return r*r*(O::dot(u,er)*(0.5L*s.mass_density_kg_m3*O::dot(u,u)+O::Gamma/(O::Gamma-1)*s.pressure_pa+O::dot(mag,mag)/O::Mu0)-O::dot(mag,er)*O::dot(u,mag)/O::Mu0);
  };
  const L fluxDifference=flux(last)-flux(first),referenceGravity=O::shellGravityPowerPerSolidAngle(a,b);
  const L referenceThermal=O::shellPressureWorkPerSolidAngle(a,b);
  const L scale=std::abs(referenceGravity)+std::abs(referenceThermal);
  q.upper("integrated_gravity_work",std::abs(gravity-referenceGravity)/scale,O::BudgetTolerance);
  q.upper("integrated_external_work",std::abs(external-(referenceThermal-referenceGravity))/scale,O::BudgetTolerance);
  q.upper("total_energy_budget",std::abs(fluxDifference-external-gravity)/scale,O::BudgetTolerance);
}

// MAT1D07 — Constant radial U has zero material acceleration but generally requires nonzero
// force. Off-equatorial Lorentz force is transverse; scalar equatorial checks alone would miss
// its omission.
void test07(Checks& q) {
  for(int sign:{1,-1})for(const auto axis:{P::Vector{{0,0,1}},P::Vector{{2,-1,2}}}) {
    const P::Profile p(config(sign,axis));
    for(L r:{.35L,.70L,1.10L})for(const auto& e:directions) {
      const auto x=position(r*O::AU,e);const auto ref=O::evaluate(wide(x),sign,unit(axis));P::Sample s;
      if(!q.status(p.sample(x,s),"force/work sample"))continue;
      const L gravityScale=ref.density*O::Gravity/(ref.radius*ref.radius);
      const L forceScale=O::norm(ref.pressureGradient)+gravityScale+O::norm(ref.lorentz);
      vectors(q,"external_force",s.external_force_n_m3,ref.force,forceScale);
      q.near("external_work",s.external_work_w_m3,ref.work,O::Speed*forceScale);
      const O::Vector curl{{s.magnetic_gradient_t_m_inv[7]-s.magnetic_gradient_t_m_inv[5],
          s.magnetic_gradient_t_m_inv[2]-s.magnetic_gradient_t_m_inv[6],s.magnetic_gradient_t_m_inv[3]-s.magnetic_gradient_t_m_inv[1]}};
      const auto lorentz=O::scaled(1/O::Mu0,O::cross(curl,wide(s.magnetic_field_t)));
      q.upper("lorentz_independent",O::norm(O::subtract(lorentz,ref.lorentz))/forceScale,O::BudgetTolerance);
      O::Vector momentum{};
      for(unsigned i=0;i<3;++i) {
        momentum[i]=-s.pressure_gradient_pa_m_inv[i]-s.mass_density_kg_m3*O::Gravity/(ref.radius*ref.radius)*ref.radial[i]+lorentz[i]+s.external_force_n_m3[i];
        for(unsigned j=0;j<3;++j)momentum[i]-=s.mass_density_kg_m3*s.velocity_m_s[j]*s.velocity_gradient_s_inv[3*i+j];
      }
      q.upper("full_vector_force_balance",O::norm(momentum)/forceScale,O::BudgetTolerance);
    }
  }
}

P::EventConfig event_config(swcme::kinematics::Mode mode=swcme::kinematics::Mode::Ballistic) {
  P::EventConfig e;e.identity="D10-1-controlled-upstream";e.begin_s=0;e.end_s=1000;
  e.front.mode=mode;e.front.r0_m=.50*static_cast<double>(O::AU);e.front.V0_m_s=1.0e6;
  e.front.Vsw_m_s=static_cast<double>(O::Speed);e.front.Gamma_m_inv=8.0e-11;
  e.front.extrapolation=swcme::kinematics::ExtrapolationPolicy::OutsideTime;
  if(mode==swcme::kinematics::Mode::DataDriven)for(double t:{0.0,300.0,600.0,1000.0}) {
    e.front.data_time_s.push_back(t);e.front.data_radius_m.push_back(e.front.r0_m+e.front.V0_m_s*t);
  }
  return e;
}

// MAT1D08 — Spatial provider support and event time support are closed and distinct. Require
// complete caller-supplied physical steps, all vertices and optional shock queries. Large epochs
// must consume explicitdt, never rounded end-begin.
void test08(Checks& q) {
  for(const auto mode:{swcme::kinematics::Mode::Ballistic,swcme::kinematics::Mode::DBM,swcme::kinematics::Mode::DataDriven}) {
    const P::Event e(config(),event_config(mode));const auto fingerprint=e.fingerprint();
    const P::Vector x{{.80*static_cast<double>(O::AU),0,0}};P::Sample sample;
    q.status(e.sample(x,0,sample),"closed time begin");q.status(e.sample(x,1000,sample),"closed time end");
    q.status(e.require_step(0,300,300,{x},true),"whole fixed300-second step with front");
    for(const auto times:{std::array<double,3>{{0,299,300}},{{300,0,300}},{{0,300,0}},{{0,300,150}}})
      q.require(!e.require_step(times[0],times[1],times[2],{x},true).ok(),"partial, reversed or nonpositive fixed step rejected");
    auto untouched=motion_sentinel();const auto before=image(untouched);
    const P::Vector nearEdge{{static_cast<double>(1.20L*O::AU-150*O::Speed),0,0}};
    q.require(!e.advance_fixed(nearEdge,0,300,300,untouched).ok()&&image(untouched)==before,
              "whole trajectory outside coverage rejects without shortened motion or output mutation");
    q.require(!e.require_step(0,300,300,{x,nearEdge},false).ok(),"all vertices covered for full step, not only the first");
    auto out=sentinel();const auto saved=image(out);
    for(double epoch:{std::nextafter(0.0,-1.0),std::nextafter(1000.0,2000.0),std::numeric_limits<double>::infinity()})
      q.require(!e.sample(x,epoch,out).ok()&&image(out)==saved,"outside/nonfinite event query is atomic");
    q.require(e.fingerprint()==fingerprint,"checked rejection does not mutate provider identity");
  }
  // Provider support covers the shock too only when explicitly requested.
  // An initially covered front which leaves during the full300-second step
  // must fail atomically; shortening the physical dt would hide the defect.
  auto edge=event_config();edge.front.r0_m=config().outer_radius_m-150*edge.front.V0_m_s;
  const P::Event edgeEvent(config(),edge);const P::Vector innerVertex{{.80*static_cast<double>(O::AU),0,0}};
  q.status(edgeEvent.require_step(0,300,300,{innerVertex},false),"wind-only whole step does not require absent shock query coverage");
  q.require(!edgeEvent.require_step(0,300,300,{innerVertex},true).ok(),"front leaving provider support rejects full required shock step");
  P::Shock unchanged;unchanged.epoch_s=-17;unchanged.position_m.fill(-19);unchanged.outward_normal.fill(-23);unchanged.jump.compression=-29;
  const auto savedPosition=unchanged.position_m,savedNormal=unchanged.outward_normal;
  q.require(!edgeEvent.shock_at({{1,0,0}},300,unchanged).ok()&&unchanged.epoch_s==-17&&
      unchanged.position_m==savedPosition&&unchanged.outward_normal==savedNormal&&unchanged.jump.compression==-29,
      "unsupported shock query does not mutate output or extrapolate its endpoint");
  auto outside=sentinel();const auto savedSample=image(outside);
  const P::Vector releaseOutside{{std::nextafter(config().outer_radius_m,std::numeric_limits<double>::infinity()),0,0}};
  q.require(!edgeEvent.sample(releaseOutside,0,outside).ok()&&image(outside)==savedSample,
      "independently requested release/insert position must itself be covered");

  // Use a genuine covered data-driven front at large absolute epochs; no
  // fictitious giant ballistic radius is introduced to test clock precision.
  auto c=event_config(swcme::kinematics::Mode::DataDriven);const double begin=1.0e12,dt=.01,end=begin+dt;
  c.begin_s=begin;c.end_s=begin+1;c.front.data_time_s={begin,begin+1};c.front.data_radius_m={.50*static_cast<double>(O::AU),.50*static_cast<double>(O::AU)+1.0e6};
  const P::Event e(config(),c);const P::Vector x{{.80*static_cast<double>(O::AU),0,0}};P::Motion m;
  q.require(end-begin!=dt,"large epoch fixture actually distinguishes stored endpoints from explicit physical duration");
  q.status(e.require_step(begin,end,dt,{x},true),"large epoch fixed step is covered");
  if(q.status(e.advance_fixed(x,begin,end,dt,m),"advance uses explicit physical duration"))
    vectors(q,"large_epoch_fixed_duration",m.position_m,O::move(wide(x),dt).position,O::norm(wide(x)));
  const double wrongEnd=std::nextafter(end,begin);
  q.require(!e.require_step(begin,wrongEnd,dt,{x},true).ok(),"wrong endpoint is not silently synchronized by shortening dt");
}

struct Flux {O::Vector electric,momentum;L mass=0,bn=0,energy=0,energyScale=0,momentumScale=0,electricScale=0;};
Flux independent_flux(const swcme::shock::PrimitiveState& s,const O::Vector& n,L speed) {
  Flux f;const auto u=O::subtract(wide(s.velocity_m_s),O::scaled(speed,n));const auto b=wide(s.magnetic_T);
  const L un=O::dot(u,n),bn=O::dot(b,n),u2=O::dot(u,u),b2=O::dot(b,b);
  f.mass=s.rho_kg_m3*un;f.bn=bn;
  f.electric=O::cross(u,b);f.electric=O::subtract(f.electric,O::scaled(O::dot(f.electric,n),n));
  f.electricScale=O::norm(u)*O::norm(b);
  const auto advected=O::scaled(s.rho_kg_m3*un,u),normal=O::scaled(s.pressure_Pa+b2/(2*O::Mu0),n),tension=O::scaled(bn/O::Mu0,b);
  for(unsigned i=0;i<3;++i)f.momentum[i]=advected[i]+normal[i]-tension[i];
  f.momentumScale=O::norm(advected)+O::norm(normal)+O::norm(tension);
  const L kinetic=.5L*s.rho_kg_m3*u2*un,enthalpy=O::Gamma/(O::Gamma-1)*s.pressure_Pa*un,
      magnetic=b2*un/O::Mu0,work=-bn*O::dot(u,b)/O::Mu0;
  f.energy=kinetic+enthalpy+magnetic+work;
  f.energyScale=std::abs(kinetic)+std::abs(enthalpy)+std::abs(magnetic)+std::abs(work);return f;
}

// MAT1D09 — A mathematical surface is insufficient: require a genuine canonical fast shock and
// independently rebuild signed mass, normalB, tangentialE, vector momentum and total-energy
// fluxes in the shock-relative frame.
void test09(Checks& q) {
  for(int sign:{1,-1})for(const auto mode:{swcme::kinematics::Mode::Ballistic,swcme::kinematics::Mode::DBM,swcme::kinematics::Mode::DataDriven}) {
    const auto ec=event_config(mode);const P::Event e(config(sign),ec);
    for(double t:{0.0,600.0})for(const auto& ray:{directions[0],directions[4],directions[5]}) {
      P::Shock shock;if(!q.status(e.shock_at(narrow(ray),t,shock),"genuine shock query"))continue;
      L radius=ec.front.r0_m+ec.front.V0_m_s*t,speed=ec.front.V0_m_s;
      if(mode==swcme::kinematics::Mode::DBM) {
        const L excess=ec.front.V0_m_s-O::Speed,gamma=ec.front.Gamma_m_inv;
        radius=ec.front.r0_m+O::Speed*t+std::log1p(gamma*excess*t)/gamma;
        speed=O::Speed+excess/(1+gamma*excess*t);
      }
      vectors(q,"front_inertial_position",shock.position_m,O::scaled(radius,ray),radius);
      q.near("front_inertial_speed",shock.jump.shock_normal_speed_m_s,speed,speed);
      const auto ambient=O::evaluate(wide(shock.position_m),sign,{{0,0,1}});
      const auto& j=shock.jump;q.require(j.has_shock&&j.solver_converged&&j.status==swcme::shock::SolveStatus::Solved,
          "actual superfast front remains genuine solved shock, not ambient/no-shock fallback");
      q.require(j.compression>1&&j.entropy_ratio>1&&j.evolutionary_fast_branch,"compression, entropy and fast evolutionary branch retained");
      q.near("shock_upstream_density",j.upstream.rho_kg_m3,ambient.density,ambient.density);
      q.near("shock_upstream_pressure",j.upstream.pressure_Pa,ambient.pressure,ambient.pressure);
      vectors(q,"shock_upstream_velocity",j.upstream.velocity_m_s,ambient.velocity,O::Speed);
      vectors(q,"shock_upstream_field",j.upstream.magnetic_T,ambient.magnetic,O::norm(ambient.magnetic));
      const L sound2=O::Gamma*ambient.pressure/ambient.density,alfven2=O::dot(ambient.magnetic,ambient.magnetic)/(O::Mu0*ambient.density),normalAlfven2=std::pow(O::dot(ambient.magnetic,ray),2)/(O::Mu0*ambient.density);
      const L fast=std::sqrt(.5L*(sound2+alfven2+std::sqrt(std::pow(sound2+alfven2,2)-4*sound2*normalAlfven2)));
      q.near("independent_fast_speed",j.fast_speed_m_s,fast,fast);
      const auto a=independent_flux(j.upstream,wide(shock.outward_normal),speed),b=independent_flux(j.downstream,wide(shock.outward_normal),speed);
      q.upper("RH_mass",std::abs(a.mass-b.mass)/std::max(std::abs(a.mass),std::abs(b.mass)),O::ShockMassTolerance);
      q.upper("RH_normal_B",std::abs(a.bn-b.bn)/std::max(O::norm(wide(j.upstream.magnetic_T)),O::norm(wide(j.downstream.magnetic_T))),O::ShockNormalBTolerance);
      q.upper("RH_tangential_E",O::norm(O::subtract(a.electric,b.electric))/(a.electricScale+b.electricScale),O::ShockElectricTolerance);
      q.upper("RH_vector_momentum",O::norm(O::subtract(a.momentum,b.momentum))/(a.momentumScale+b.momentumScale),O::ShockMomentumTolerance);
      q.upper("RH_total_energy",std::abs(a.energy-b.energy)/(a.energyScale+b.energyScale),O::ShockEnergyTolerance);
      q.upper("canonical_mass_residual",j.mass_residual,O::ShockMassTolerance);q.upper("canonical_normal_B_residual",j.normal_B_residual,O::ShockNormalBTolerance);
      q.upper("canonical_electric_residual",j.electric_residual,O::ShockElectricTolerance);q.upper("canonical_momentum_residual",j.momentum_residual,O::ShockMomentumTolerance);q.upper("canonical_energy_residual",j.energy_residual,O::ShockEnergyTolerance);
    }
  }
  auto slow=event_config();slow.front.V0_m_s=static_cast<double>(O::Speed);const P::Event noShock(config(),slow);P::Shock state;
  if(q.status(noShock.shock_at({{1,0,0}},300,state),"subfast surface query")) {
    q.require(!state.jump.has_shock&&state.jump.solver_converged&&state.jump.status==swcme::shock::SolveStatus::NoShock&&state.jump.compression==1,
              "subfast mathematical front cannot manufacture a shock/compression floor");
    const auto& a=state.jump.upstream;const auto& b=state.jump.downstream;
    q.require(a.rho_kg_m3==b.rho_kg_m3&&a.pressure_Pa==b.pressure_Pa&&a.velocity_m_s==b.velocity_m_s&&a.magnetic_T==b.magnetic_T,
              "NoShock preserves complete upstream one-sided state exactly");
  }
}

// MAT1D10 — Malformed or unsupported contracts fail before publication. Preserve every sample
// field on checked rejection, and verify the old Leblanc diagnostic keeps both its values and
// visible nonconserved radial mass flux.
void test10(Checks& q) {
  auto rejects=[&](const P::Config& c){bool caught=false;try{const P::Profile p(c);(void)p;}catch(const std::invalid_argument&){caught=true;}q.require(caught,"invalid/missing explicit physical contract rejected");};
  rejects(P::Config{});
  for(const auto mutation:{0,1,2,3,4,5,6,7,8}) {
    auto c=config();
    if(mutation==0)c.wind_m_s=0;
    if(mutation==1)c.polarity=0;
    if(mutation==2)c.rotation_axis={{0,0,0}};
    if(mutation==3)c.gamma=1;
    if(mutation==4)c.reference_electron_temperature_k=0;
    if(mutation==5)c.reference_proton_temperature_k=0;
    if(mutation==6)c.reference_sin_colatitude=-1;
    if(mutation==7)c.outer_radius_m=c.inner_radius_m;
    if(mutation==8)c.solar_gravity_m3_s2=0;
    rejects(c);
  }
  const P::Profile p(config());auto s=sentinel();const auto before=image(s);
  for(const auto x:{P::Vector{{std::nextafter(config().inner_radius_m,0.0),0,0}},P::Vector{{std::nextafter(config().outer_radius_m,std::numeric_limits<double>::infinity()),0,0}},P::Vector{{0,0,0}},P::Vector{{std::numeric_limits<double>::quiet_NaN(),1,2}}})
    q.require(!p.sample(x,s).ok()&&image(s)==before,"outside/invalid position cannot partially overwrite sample");
  for(int mutation:{0,1,2,3,4}) {
    auto e=event_config();
    if(mutation==0)e.identity.clear();
    if(mutation==1)e.end_s=-1;
    if(mutation==2)e.front.mode=static_cast<swcme::kinematics::Mode>(-1);
    if(mutation==3)e.front.extrapolation=swcme::kinematics::ExtrapolationPolicy::Ballistic;
    if(mutation==4)e.front.Vsw_m_s+=1;
    bool caught=false;try{const P::Event rejected(config(),e);(void)rejected;}catch(const std::invalid_argument&){caught=true;}
    q.require(caught,"unsupported event mode/extrapolation/identity/time/wind is rejected");
  }
  // The old diagnostic is still Leblanc + constant radial wind. Verify both
  // its independent value and its NONZERO radial mass-flux defect. Its current
  // output must not be silently redefined as the newly admitted r^-2 profile.
  swcme1d::Params parameters;parameters.region_mode=swcme::regions::Mode::ShockOnly;
  const swcme1d::Model legacy(parameters);const auto prepared=legacy.prepare_step(0);
  L firstFlux=0,lastFlux=0;
  const auto leblanc=[](L r){const L rs=r/695700000.0L;return 3.3e5L/(rs*rs)+4.1e6L/std::pow(rs,4)+8e7L/std::pow(rs,6);};
  for(L radius:{.30L,1.20L}) {
    const double r=static_cast<double>(radius*O::AU);double n=0,u=0;
    q.status(legacy.evaluate_radii_fast_checked(prepared,&r,&n,&u,1),"legacy diagnostic retained");
    const L expected=static_cast<L>(parameters.n1AU_cm3)*1e6L*leblanc(r)/leblanc(O::AU);
    q.near("legacy_leblanc_unchanged",n,expected,expected);
    const L flux=n*r*r*u;if(firstFlux==0)firstFlux=flux;else lastFlux=flux;
  }
  q.lower("legacy_mass_flux_defect_visible",std::abs(lastFlux/firstFlux-1),O::BudgetTolerance);
}

struct Test {const char* id;const char* name;void(*run)(Checks&);};
const Test tests[]={
 {"MAT1D01","Approved primitive/composition/EOS anchors",test01},
 {"MAT1D02","Independent radial mass inventory and flux",test02},
 {"MAT1D03","Full material deformation and signed Cauchy transport",test03},
 {"MAT1D04","Signed magnetic flux and ideal induction",test04},
 {"MAT1D05","All Cartesian gradients and second-order differentiation",test05},
 {"MAT1D06","Adiabatic entropy and total-energy/work budget",test06},
 {"MAT1D07","Full off-equator vector force and mechanical work",test07},
 {"MAT1D08","Closed whole-step temporal/spatial coverage and large epoch",test08},
 {"MAT1D09","Genuine canonical shock with independent RH fluxes",test09},
 {"MAT1D10","Invalid/unsupported and preserved legacy diagnostics",test10}};
int run(const Test& t) {
  Checks q;try{t.run(q);}catch(const std::exception& error){std::cerr<<"ERROR "<<error.what()<<'\n';return 2;}
  q.report();std::cout<<'['<<t.id<<"] "<<(q.failures?"FAIL":"PASS")<<" checks="<<q.checks<<" failures="<<q.failures<<" name="<<t.name<<'\n';return q.failures?1:0;
}
}
int main(int argc,char** argv) {
  if(argc==2&&std::string(argv[1])=="--list") {for(const auto& t:tests)std::cout<<t.id<<" MATERIAL1D "<<t.name<<'\n';return 0;}
  if(argc==3&&std::string(argv[1])=="--test") {for(const auto& t:tests)if(t.id==std::string(argv[2]))return run(t);std::cerr<<"unknown material test ID\n";return 2;}
  if(argc==1||(argc==2&&std::string(argv[1])=="--all")) {int code=0;for(const auto& t:tests)code=std::max(code,run(t));return code;}
  std::cerr<<"usage: test_material_upstream [--list|--all|--test MAT1Dxx]\n";return 2;
}
