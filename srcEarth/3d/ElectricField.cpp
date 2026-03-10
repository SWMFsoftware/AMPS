#include "ElectricField.h"

#include <cmath>
#include <algorithm>

#include "pic.h"
#include "../../interface/T96Interface.h"
#include "../../interface/T05Interface.h"
#include "../gridless/DipoleInterface.h"

namespace {
inline double sqr(double x) { return x*x; }
const double kOmegaEarth = 7.2921159e-5;

void Cross(const double a[3],const double b[3],double c[3]) {
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

void EvalDipoleSI(double B[3],const double x[3],const EarthUtil::AmpsParam& prm) {
  // Reuse the analytic dipole helper that already exists in srcEarth/gridless.
  // That interface owns the canonical dipole constants and the tilt semantics
  // used elsewhere in srcEarth, so using it here keeps the 3D initialization
  // path consistent with the dipole-only cutoff/density workflows.
  //
  // IMPORTANT:
  // The gridless dipole interface does not expose a free function named
  // `GetMagneticFieldDipole(...)`. Instead, it keeps a small shared state in
  // `Earth::GridlessMode::Dipole::gParams` and provides the public helpers
  // `SetMomentScale`, `SetTiltDeg`, and `GetB_Tesla`.
  //
  // Therefore, for the 3D field initialization path we temporarily configure
  // the dipole helper from the parsed input parameters and then evaluate the
  // field at the requested SI position.
  Earth::GridlessMode::Dipole::SetMomentScale(prm.field.dipoleMoment_Me);
  Earth::GridlessMode::Dipole::SetTiltDeg(prm.field.dipoleTilt_deg);
  Earth::GridlessMode::Dipole::GetB_Tesla(x,B);
}
} // namespace

namespace Earth {
namespace Mode3D {

void EvaluateBackgroundMagneticFieldSI(double B[3],const double xGSM_SI[3],const EarthUtil::AmpsParam& prm) {
  const std::string model=EarthUtil::ToUpper(prm.field.model);
  if (model=="DIPOLE") { EvalDipoleSI(B,xGSM_SI,prm); return; }
  if (model=="T96") { ::T96::GetMagneticField(B,const_cast<double*>(xGSM_SI)); return; }
  if (model=="T05") { ::T05::GetMagneticField(B,const_cast<double*>(xGSM_SI)); return; }
  B[0]=B[1]=B[2]=0.0;
}

// Corotation: E = -(Omega x r) x B. This is a clean SI-space implementation
// of Earth corotation suitable for initializing a mesh field. The literature
// often discusses corotation together with the Volland-Stern family of inner-
// magnetosphere convection models: Volland (1973), Stern (1975), Maynard &
// Chen (1975), and later model comparisons such as Pierrard et al. (2008).
static void EvalCorotationSI(double E[3],const double x[3],const EarthUtil::AmpsParam& prm) {
  double B[3],vx[3],e[3],omega[3]={0.0,0.0,kOmegaEarth};
  EvaluateBackgroundMagneticFieldSI(B,x,prm);
  Cross(omega,x,vx);
  Cross(vx,B,e);
  E[0]=-prm.efield.corotationScale*e[0];
  E[1]=-prm.efield.corotationScale*e[1];
  E[2]=-prm.efield.corotationScale*e[2];
}
static double EvalCorotationPotentialV(const double x[3],const EarthUtil::AmpsParam& prm) {(void)x;(void)prm; return 0.0;}

// Volland-Stern-type volumetric potential. Classical published forms are most
// commonly given in equatorial/L-shell variables. For AMPS mesh initialization
// we need an everywhere-defined 3D recipe, so we use a simple cylindrical
// extension with the standard shielded two-cell local-time dependence.
static double EvalVSPotentialV(const double x[3],const EarthUtil::AmpsParam& prm) {
  const double rho_m=std::sqrt(sqr(x[0])+sqr(x[1]));
  const double rho_re=std::max(rho_m/_EARTH__RADIUS_,prm.efield.lMin);
  const double phi=std::atan2(x[1],x[0]);
  const double ampV=1000.0*prm.efield.vsPotential_kV*prm.efield.vsScale;
  return ampV*std::pow(prm.efield.vsReferenceL/rho_re,prm.efield.vsGamma)*std::sin(phi);
}

static void EvalVollandSternSI(double E[3],const double x[3],const EarthUtil::AmpsParam& prm) {
  const double r=std::sqrt(sqr(x[0])+sqr(x[1])+sqr(x[2]));
  const double h=std::max(1.0e3,1.0e-5*std::max(r,_EARTH__RADIUS_));
  E[0]=E[1]=E[2]=0.0;
  for (int idim=0;idim<3;idim++) {
    double xp[3]={x[0],x[1],x[2]},xm[3]={x[0],x[1],x[2]};
    xp[idim]+=h; xm[idim]-=h;
    E[idim]=-(EvalVSPotentialV(xp,prm)-EvalVSPotentialV(xm,prm))/(2.0*h);
  }
}

void EvaluateElectricFieldSI(double E[3],const double xGSM_SI[3],const EarthUtil::AmpsParam& prm) {
  E[0]=E[1]=E[2]=0.0;
  const std::string model=EarthUtil::ToUpper(prm.efield.model);
  if (model=="NONE" || model.empty()) return;
  const double r=std::sqrt(sqr(xGSM_SI[0])+sqr(xGSM_SI[1])+sqr(xGSM_SI[2]));
  if (r<std::max(1.0,prm.efield.rMin_km*1000.0)) return;
  if (model=="COROTATION") { EvalCorotationSI(E,xGSM_SI,prm); return; }
  if (model=="VOLLAND_STERN") { EvalVollandSternSI(E,xGSM_SI,prm); return; }
  if (model=="COROTATION_VOLLAND_STERN") {
    double a[3],b[3];
    EvalCorotationSI(a,xGSM_SI,prm);
    EvalVollandSternSI(b,xGSM_SI,prm);
    for (int i=0;i<3;i++) E[i]=a[i]+b[i];
  }
}

double EvaluateElectricPotential_V(const double xGSM_SI[3],const EarthUtil::AmpsParam& prm) {
  const std::string model=EarthUtil::ToUpper(prm.efield.model);
  if (model=="COROTATION") return EvalCorotationPotentialV(xGSM_SI,prm);
  if (model=="VOLLAND_STERN") return EvalVSPotentialV(xGSM_SI,prm);
  if (model=="COROTATION_VOLLAND_STERN") return EvalCorotationPotentialV(xGSM_SI,prm)+EvalVSPotentialV(xGSM_SI,prm);
  return 0.0;
}

} // namespace Mode3D
} // namespace Earth
