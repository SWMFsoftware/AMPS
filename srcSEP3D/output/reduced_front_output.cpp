#include "reduced_front_output.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace SEP3D { namespace Output {

std::string SerializeReducedFrontTecplot(
    const SEP::CoronaSwcme::ShockFront::Epoch& epoch,
    const SEP::CoronaSwcme::ShockFront::Configuration& configuration) {
  namespace SF=SEP::CoronaSwcme::ShockFront;
  const int nPolar=configuration.polarCells;
  const int nAzimuth=configuration.azimuthCells;
  const std::size_t expected=static_cast<std::size_t>(nPolar)*nAzimuth;
  if(nPolar<2||nAzimuth<3||epoch.records.size()!=expected)return {};

  std::ostringstream out;
  out<<std::setprecision(17);
  out<<"TITLE=\"srcSEP3D reduced prescribed front; RH limits valid only when shock_accepted=1\"\n"
     <<"VARIABLES=\"x_m\",\"y_m\",\"z_m\",\"normal_x\",\"normal_y\","
       "\"normal_z\",\"normal_speed_m_s\",\"quadrature_area_m2\","
       "\"time_s\",\"generation\",\"status_code\",\"shock_accepted\","
       "\"upstream_valid\",\"downstream_valid\",\"fast_mach\","
       "\"theta_Bn_rad\",\"theta_Bn_valid\",\"density_compression\","
       "\"magnetic_compression\",\"magnetic_compression_valid\","
       "\"rho1_kg_m3\",\"p1_Pa\",\"Tproton1_K\","
       "\"u1x_m_s\",\"u1y_m_s\",\"u1z_m_s\","
       "\"b1x_T\",\"b1y_T\",\"b1z_T\","
       "\"rho2_kg_m3\",\"p2_Pa\","
       "\"u2x_m_s\",\"u2y_m_s\",\"u2z_m_s\","
       "\"b2x_T\",\"b2y_T\",\"b2z_T\"\n"
     <<"ZONE T=\"front-t"<<epoch.trajectory.timeS<<"-g"<<epoch.generation
     <<"\", N="<<expected<<", E="
     <<static_cast<std::size_t>(nPolar-1)*nAzimuth
     <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n"
     <<"AUXDATA event_identity=\""<<epoch.eventIdentity<<"\"\n"
     <<"AUXDATA phase=\""<<SF::Name(epoch.trajectory.phase)<<"\"\n"
     <<"AUXDATA volume_role=\"ambient-reference-only\"\n"
     <<"AUXDATA no_shock_fill=\"upstream-ambient-visualization-placeholder\"\n";

  for(const auto& record:epoch.records) {
    const bool upstreamValid=record.status!=SF::FrontStatus::OutsideFrontSupport&&
        record.status!=SF::FrontStatus::BelowPhysicalInnerBoundary&&
        record.status!=SF::FrontStatus::AmbientUnavailable;
    // A status label alone is insufficient: an accepted record must also own
    // the complete RH state.  The provider establishes this invariant, while
    // the conjunction here prevents a malformed record from publishing a
    // shock_accepted=1 flag next to visualization fallback values.
    const bool accepted=record.status==SF::FrontStatus::SolvedFastShock&&
        record.downstreamValid;
    const bool thetaValid=upstreamValid&&record.magneticDirectionValid;
    const bool magneticCompressionValid=accepted&&
        record.diagnostics.magneticCompressionValid;
    const double theta=thetaValid?std::acos(std::max(0.0,std::min(1.0,
        std::abs(record.signedMagneticNormalCosine)))):0.0;
    // These finite neutral values are deliberately an output convention, not
    // an alternate physical closure.  In particular, Mf=0 prevents a plotting
    // package from drawing a non-shock patch as super-fast, while X=CB=1 says
    // that the display fallback is the unchanged upstream ambient.  Consumers
    // must inspect shock_accepted/downstream_valid before interpreting any
    // downstream column as a shock limit.
    const double outputMach=accepted?record.fastMach:0.0;
    const double densityCompression=accepted?
        record.jump.compressionRatio:1.0;
    const double magneticCompression=magneticCompressionValid?
        record.diagnostics.magneticCompression:1.0;
    const auto& x=record.geometry.positionM;
    const auto& n=record.geometry.outwardNormal;
    out<<x.x<<' '<<x.y<<' '<<x.z<<' '<<n.x<<' '<<n.y<<' '<<n.z<<' '
       <<record.geometry.normalSpeedMPerS<<' '<<record.geometry.areaM2<<' '
       <<epoch.trajectory.timeS<<' '<<epoch.generation<<' '
       <<static_cast<int>(record.status)<<' '<<(accepted?1:0)<<' '
       <<(upstreamValid?1:0)<<' '<<(accepted?1:0)<<' '
       <<outputMach<<' '<<theta<<' '<<(thetaValid?1:0)<<' '
       <<densityCompression<<' '<<magneticCompression<<' '
       <<(magneticCompressionValid?1:0)<<' ';
    if(upstreamValid) {
      const auto& u=record.upstream.velocityMPerS;
      const auto& b=record.upstream.magneticFieldT;
      out<<record.upstream.plasma.massDensityKgM3<<' '
         <<record.upstream.plasma.pressurePa<<' '
         <<record.upstream.protonTemperatureK<<' '
         <<u.x<<' '<<u.y<<' '<<u.z<<' '
         <<b.x<<' '<<b.y<<' '<<b.z<<' ';
    } else {
      // No ambient state exists to copy.  Zero is only a finite Tecplot
      // placeholder and upstream_valid=0 makes it unusable as plasma data.
      for(int i=0;i<9;++i)out<<0.0<<' ';
    }
    if(accepted) {
      const auto& u=record.jump.downstream.velocityMPerS;
      const auto& b=record.jump.downstream.magneticFieldT;
      out<<record.jump.downstream.massDensityKgM3<<' '
         <<record.jump.downstream.pressurePa<<' '
         <<u.x<<' '<<u.y<<' '<<u.z<<' '
         <<b.x<<' '<<b.y<<' '<<b.z;
    } else if(upstreamValid) {
      // Finite no-shock display encoding requested for Tecplot.  Copying the
      // upstream primitive makes X=CB=1 algebraically consistent.  The two
      // validity flags above stay zero, so this must never be called a
      // downstream state, sheath, or ejecta value.
      const auto& u=record.upstream.velocityMPerS;
      const auto& b=record.upstream.magneticFieldT;
      out<<record.upstream.plasma.massDensityKgM3<<' '
         <<record.upstream.plasma.pressurePa<<' '
         <<u.x<<' '<<u.y<<' '<<u.z<<' '
         <<b.x<<' '<<b.y<<' '<<b.z;
    } else {
      for(int i=0;i<8;++i)out<<0.0<<(i==7?'\n':' ');
      continue;
    }
    out<<'\n';
  }

  // Records are ordered [polar][azimuth].  Close only the periodic azimuthal
  // seam; the two finite-width support rims remain open.  This avoids drawing
  // a fictitious cap across the unsupported rear of the generating sphere.
  for(int i=0;i<nPolar-1;++i)for(int j=0;j<nAzimuth;++j) {
    const int next=(j+1)%nAzimuth;
    const int a=i*nAzimuth+j+1;
    const int b=i*nAzimuth+next+1;
    const int c=(i+1)*nAzimuth+next+1;
    const int d=(i+1)*nAzimuth+j+1;
    out<<a<<' '<<b<<' '<<c<<' '<<d<<'\n';
  }
  return out.str();
}

} } // namespace SEP3D::Output
