#include "reduced_front_output.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

namespace SEP3D { namespace Output { namespace {

// Tecplot requires BLOCK packing when an FE zone mixes nodal coordinates and
// cell-centred variables. Six values per text line keeps the product readable
// without assigning any physical meaning to line breaks inside a data block.
template<class Getter>
void WriteBlock(std::ostringstream* out,std::size_t count,Getter value) {
  for(std::size_t i=0;i<count;++i) {
    *out<<value(i);
    *out<<(i+1==count||i%6==5?'\n':' ');
  }
}

} // namespace

std::string SerializeReducedFrontTecplot(
    const SEP::CoronaSwcme::ShockFront::Epoch& epoch,
    const SEP::CoronaSwcme::ShockFront::Configuration& configuration) {
  namespace SF=SEP::CoronaSwcme::ShockFront;
  constexpr std::size_t kCellVariables=36;
  const std::size_t nodeCount=epoch.vertices.size();
  const std::size_t faceCount=epoch.triangles.size();
  if(configuration.surfaceTopology!="triangular-sse-cap-v1"||nodeCount<4||
      faceCount<3||epoch.records.size()!=faceCount)return {};

  // Assemble face values before writing any bytes. The output convention is
  // deliberately finite for Tecplot/VisIt, but validity remains explicit:
  // only SolvedFastShock plus a complete RH state sets shock_accepted=1.
  // Every other physical/numerical status receives Mach=0, compression=1 and
  // (when available) an unchanged upstream primitive in the display-only
  // downstream columns. Those placeholders are not downstream CME plasma.
  std::vector<std::array<double,kCellVariables>> cells(faceCount);
  for(std::size_t face=0;face<faceCount;++face) {
    const auto& triangle=epoch.triangles[face];
    const auto& record=epoch.records[face];
    if(triangle.stableId!=record.geometry.stableId||
        triangle.curvedAreaM2!=record.geometry.areaM2||
        triangle.vertex[0]>=nodeCount||triangle.vertex[1]>=nodeCount||
        triangle.vertex[2]>=nodeCount)return {};
    const bool upstreamValid=record.status!=SF::FrontStatus::OutsideFrontSupport&&
        record.status!=SF::FrontStatus::BelowPhysicalInnerBoundary&&
        record.status!=SF::FrontStatus::AmbientUnavailable;
    const bool accepted=record.status==SF::FrontStatus::SolvedFastShock&&
        record.downstreamValid;
    const bool thetaValid=upstreamValid&&record.magneticDirectionValid;
    const bool magneticCompressionValid=accepted&&
        record.diagnostics.magneticCompressionValid;
    const double theta=thetaValid?std::acos(std::max(0.0,std::min(1.0,
        std::abs(record.signedMagneticNormalCosine)))):0.0;
    const auto& n=record.geometry.outwardNormal;
    auto& row=cells[face];
    row[0]=n.x;row[1]=n.y;row[2]=n.z;
    row[3]=record.geometry.normalSpeedMPerS;
    row[4]=triangle.curvedAreaM2;
    row[5]=triangle.planarAreaM2;
    row[6]=epoch.trajectory.timeS;
    row[7]=static_cast<double>(epoch.generation);
    row[8]=static_cast<double>(triangle.stableId);
    row[9]=static_cast<double>(static_cast<int>(record.status));
    row[10]=accepted?1.0:0.0;
    row[11]=upstreamValid?1.0:0.0;
    row[12]=accepted?1.0:0.0;
    row[13]=accepted?record.fastMach:0.0;
    row[14]=theta;
    row[15]=thetaValid?1.0:0.0;
    row[16]=accepted?record.jump.compressionRatio:1.0;
    row[17]=magneticCompressionValid?
        record.diagnostics.magneticCompression:1.0;
    row[18]=magneticCompressionValid?1.0:0.0;
    if(upstreamValid) {
      const auto& u=record.upstream.velocityMPerS;
      const auto& b=record.upstream.magneticFieldT;
      row[19]=record.upstream.plasma.massDensityKgM3;
      row[20]=record.upstream.plasma.pressurePa;
      row[21]=record.upstream.protonTemperatureK;
      row[22]=u.x;row[23]=u.y;row[24]=u.z;
      row[25]=b.x;row[26]=b.y;row[27]=b.z;
    }
    if(accepted) {
      const auto& u=record.jump.downstream.velocityMPerS;
      const auto& b=record.jump.downstream.magneticFieldT;
      row[28]=record.jump.downstream.massDensityKgM3;
      row[29]=record.jump.downstream.pressurePa;
      row[30]=u.x;row[31]=u.y;row[32]=u.z;
      row[33]=b.x;row[34]=b.y;row[35]=b.z;
    } else if(upstreamValid) {
      const auto& u=record.upstream.velocityMPerS;
      const auto& b=record.upstream.magneticFieldT;
      row[28]=record.upstream.plasma.massDensityKgM3;
      row[29]=record.upstream.plasma.pressurePa;
      row[30]=u.x;row[31]=u.y;row[32]=u.z;
      row[33]=b.x;row[34]=b.y;row[35]=b.z;
    }
    for(double value:row)if(!std::isfinite(value))return {};
  }

  std::ostringstream out;
  out<<std::setprecision(17);
  out<<"TITLE=\"srcSEP3D reduced prescribed front; RH limits valid only when shock_accepted=1\"\n"
     <<"VARIABLES=\"x_m\",\"y_m\",\"z_m\",\"normal_x\",\"normal_y\","
       "\"normal_z\",\"normal_speed_m_s\",\"quadrature_area_m2\","
       "\"planar_chord_area_m2\",\"time_s\",\"generation\","
       "\"triangle_stable_id\",\"status_code\",\"shock_accepted\","
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
     <<"\", N="<<nodeCount<<", E="<<faceCount
     <<", DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE, "
       "VARLOCATION=([4-39]=CELLCENTERED)\n"
     <<"AUXDATA event_identity=\""<<epoch.eventIdentity<<"\"\n"
     <<"AUXDATA phase=\""<<SF::Name(epoch.trajectory.phase)<<"\"\n"
     <<"AUXDATA surface_topology=\""<<configuration.surfaceTopology<<"\"\n"
     <<"AUXDATA volume_role=\"ambient-reference-only\"\n"
     <<"AUXDATA no_shock_fill=\"upstream-ambient-visualization-placeholder\"\n";

  // Only Cartesian coordinates are nodal. All normals, areas, status flags
  // and plasma/RH values are face-centred so no visualization interpolation
  // can turn a shock/no-shock boundary into a fictitious intermediate state.
  WriteBlock(&out,nodeCount,[&](std::size_t i){return epoch.vertices[i].positionM.x;});
  WriteBlock(&out,nodeCount,[&](std::size_t i){return epoch.vertices[i].positionM.y;});
  WriteBlock(&out,nodeCount,[&](std::size_t i){return epoch.vertices[i].positionM.z;});
  for(std::size_t variable=0;variable<kCellVariables;++variable)
    WriteBlock(&out,faceCount,[&](std::size_t i){return cells[i][variable];});
  for(const auto& triangle:epoch.triangles)out
      <<triangle.vertex[0]+1<<' '<<triangle.vertex[1]+1<<' '
      <<triangle.vertex[2]+1<<'\n';
  return out.str();
}

} } // namespace SEP3D::Output
