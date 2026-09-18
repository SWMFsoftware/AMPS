#ifndef _SEP_FLUX_TUBE_GEOMETRY_CORE_H_
#define _SEP_FLUX_TUBE_GEOMETRY_CORE_H_

#include "sep_physical_units.h"
#include "sep_common_header_path.h"

#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace FieldLine {
namespace FluxTubeGeometryCore {

// Every field line owns a conserved magnetic flux in webers.  Generation and
// provenance make remap/restart decisions auditable instead of reconstructing
// an implicit area normalization independently in each consumer.
struct MagneticFluxRecord {
  double magneticFluxWb = 0.0;
  std::uint64_t generation = 0;
  std::string provenance;
};

Transport::Status ValidateMagneticFluxRecord(const MagneticFluxRecord& record);
Units::AreaM2 AreaFromMagneticFlux(const MagneticFluxRecord& record,
                                  Units::MagneticFieldT local_abs_B_T);
Transport::Status RedistributeMagneticFlux(
    const MagneticFluxRecord& parent,
    const std::vector<double>& fractions,
    std::vector<MagneticFluxRecord>* children);
Transport::Status SerializeMagneticFluxTable(
    const std::vector<MagneticFluxRecord>& records, std::string* text);
Transport::Status DeserializeMagneticFluxTable(
    const std::string& text, std::vector<MagneticFluxRecord>* records);

// Return the local tube area from conservation of magnetic flux,
// A_local |B_local| = A_reference |B_reference|.  Invalid or zero fields are
// rejected because inventing an area profile would silently corrupt every
// downstream volume normalization.
Units::AreaM2 AreaFromMagneticFlux(
    Units::AreaM2 reference_area_m2,
    Units::MagneticFieldT reference_abs_B_T,
    Units::MagneticFieldT local_abs_B_T);

// Integrate a linearly varying cross-sectional area over a segment.  Field-line
// vertices store the endpoint state, so the trapezoid is the conservative and
// grid-convergent finite-volume representation of that state.
Units::VolumeM3 IntegrateLinearArea(
    Units::AreaM2 begin_area_m2,
    Units::AreaM2 end_area_m2,
    Units::LengthM segment_length_m,
    double begin_fraction = 0.0,
    double end_fraction = 1.0);

// Integrate an arbitrary smooth area profile from its begin/mid/end samples.
// The PIC adapter evaluates those samples from A|B|=constant (or the configured
// explicit profile), making this the production segment-volume quadrature.
Units::VolumeM3 IntegrateSimpsonArea(
    Units::AreaM2 begin_area_m2,
    Units::AreaM2 middle_area_m2,
    Units::AreaM2 end_area_m2,
    Units::LengthM interval_length_m);

// A shock moving normally through a flux tube sweeps A v dt cubic metres.
Units::VolumeM3 SweptVolume(
    Units::AreaM2 area_m2,
    Units::SpeedMPerS normal_speed_m_s,
    Units::TimeS time_s);

// Shared physical-source helpers keep injection efficiency and statistical
// particle weight independent of the background provider used to obtain the
// same density and shock state.
double InjectedPhysicalParticleCount(
    Units::NumberDensityPerM3 density_per_m3,
    Units::VolumeM3 swept_volume_m3,
    double injection_efficiency);
double MacroparticleWeight(double physical_particle_count,
                           long macroparticle_count);
double SpectralNumberDensityPerJ(Units::NumberDensityPerM3 density_per_m3,
                                 double normalized_probability_per_J);

}  // namespace FluxTubeGeometryCore
}  // namespace FieldLine
}  // namespace SEP

#endif
