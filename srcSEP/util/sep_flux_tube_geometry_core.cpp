#include "sep_flux_tube_geometry_core.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace {

void RequireFinitePositive(double value, const char* name) {
  if (!std::isfinite(value) || value <= 0.0) {
    throw std::invalid_argument(name);
  }
}

}  // namespace

SEP::Transport::Status
SEP::FieldLine::FluxTubeGeometryCore::ValidateMagneticFluxRecord(
    const MagneticFluxRecord& record) {
  if (!std::isfinite(record.magneticFluxWb) || record.magneticFluxWb<=0.0 ||
      record.provenance.empty())
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
        "field-line magnetic flux must be positive, finite, and provenance tagged");
  return Transport::Status::Ok();
}

SEP::Units::AreaM2
SEP::FieldLine::FluxTubeGeometryCore::AreaFromMagneticFlux(
    const MagneticFluxRecord& record, Units::MagneticFieldT local_abs_B_T) {
  const Transport::Status status=ValidateMagneticFluxRecord(record);
  if (!status.ok()) throw std::invalid_argument(status.message);
  RequireFinitePositive(local_abs_B_T.Value(),
                        "local |B| must be finite and positive");
  return Units::AreaM2(record.magneticFluxWb/local_abs_B_T.Value());
}

SEP::Transport::Status
SEP::FieldLine::FluxTubeGeometryCore::RedistributeMagneticFlux(
    const MagneticFluxRecord& parent, const std::vector<double>& fractions,
    std::vector<MagneticFluxRecord>* children) {
  const Transport::Status status=ValidateMagneticFluxRecord(parent);
  if (!status.ok() || !children || fractions.empty())
    return status.ok() ? Transport::Status::Error(
        Transport::StatusCode::InvalidArgument,"flux redistribution output is invalid")
        : status;
  long double sum=0.0L;
  for (std::size_t i=0;i<fractions.size();++i) {
    if (!std::isfinite(fractions[i]) || fractions[i]<=0.0)
      return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                      "flux fractions must be positive and finite");
    sum+=fractions[i];
  }
  if (std::fabs(static_cast<double>(sum)-1.0)>1.0e-12)
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                    "flux redistribution fractions must sum to one");
  children->assign(fractions.size(),MagneticFluxRecord());
  double assigned=0.0;
  for (std::size_t i=0;i<fractions.size();++i) {
    MagneticFluxRecord& child=(*children)[i];
    child.magneticFluxWb=i+1==fractions.size()
        ? parent.magneticFluxWb-assigned
        : parent.magneticFluxWb*fractions[i];
    assigned+=child.magneticFluxWb;
    child.generation=parent.generation+1;
    child.provenance=parent.provenance+";conservative-refinement";
  }
  return Transport::Status::Ok();
}

SEP::Transport::Status
SEP::FieldLine::FluxTubeGeometryCore::SerializeMagneticFluxTable(
    const std::vector<MagneticFluxRecord>& records, std::string* text) {
  if (!text) return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                             "flux checkpoint output is null");
  std::ostringstream out;
  out << "SEP_FIELD_LINE_FLUX 1\n" << records.size() << '\n'
      << std::setprecision(17);
  for (std::size_t i=0;i<records.size();++i) {
    const Transport::Status valid=ValidateMagneticFluxRecord(records[i]);
    if (!valid.ok() || records[i].provenance.find('\n')!=std::string::npos)
      return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                      "flux record cannot be serialized");
    out << records[i].magneticFluxWb << ' ' << records[i].generation << ' '
        << records[i].provenance.size() << ' ' << records[i].provenance << '\n';
  }
  *text=out.str();
  return Transport::Status::Ok();
}

SEP::Transport::Status
SEP::FieldLine::FluxTubeGeometryCore::DeserializeMagneticFluxTable(
    const std::string& text, std::vector<MagneticFluxRecord>* records) {
  if (!records) return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                                "flux checkpoint destination is null");
  std::istringstream in(text); std::string magic; int version=0; std::size_t n=0;
  if (!(in>>magic>>version>>n) || magic!="SEP_FIELD_LINE_FLUX" || version!=1)
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                    "flux checkpoint header is invalid");
  std::vector<MagneticFluxRecord> parsed(n);
  for (std::size_t i=0;i<n;++i) {
    std::size_t length=0;
    if (!(in>>parsed[i].magneticFluxWb>>parsed[i].generation>>length))
      return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                      "flux checkpoint is truncated");
    if (length==0)
      return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                      "flux checkpoint provenance is empty");
    in.get();
    parsed[i].provenance.resize(length);
    in.read(&parsed[i].provenance[0],static_cast<std::streamsize>(length));
    if (!in || !ValidateMagneticFluxRecord(parsed[i]).ok())
      return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                      "flux checkpoint record is invalid");
    if (in.peek()=='\n') in.get();
  }
  *records=parsed;
  return Transport::Status::Ok();
}

SEP::Units::AreaM2
SEP::FieldLine::FluxTubeGeometryCore::AreaFromMagneticFlux(
    Units::AreaM2 reference_area_m2,
    Units::MagneticFieldT reference_abs_B_T,
    Units::MagneticFieldT local_abs_B_T) {
  RequireFinitePositive(reference_area_m2.Value(),
                        "reference area must be finite and positive");
  RequireFinitePositive(reference_abs_B_T.Value(),
                        "reference |B| must be finite and positive");
  RequireFinitePositive(local_abs_B_T.Value(),
                        "local |B| must be finite and positive");
  return Units::AreaM2(reference_area_m2.Value() *
                       reference_abs_B_T.Value() / local_abs_B_T.Value());
}

SEP::Units::VolumeM3
SEP::FieldLine::FluxTubeGeometryCore::IntegrateLinearArea(
    Units::AreaM2 begin_area_m2,
    Units::AreaM2 end_area_m2,
    Units::LengthM segment_length_m,
    double begin_fraction,
    double end_fraction) {
  RequireFinitePositive(begin_area_m2.Value(),
                        "begin area must be finite and positive");
  RequireFinitePositive(end_area_m2.Value(),
                        "end area must be finite and positive");
  RequireFinitePositive(segment_length_m.Value(),
                        "segment length must be finite and positive");
  if (!std::isfinite(begin_fraction) || !std::isfinite(end_fraction) ||
      begin_fraction < 0.0 || end_fraction > 1.0 ||
      end_fraction <= begin_fraction) {
    throw std::invalid_argument("segment fractions must satisfy 0 <= begin < end <= 1");
  }

  const double delta_area_m2 = end_area_m2.Value() - begin_area_m2.Value();
  const double area_at_begin_m2 =
      begin_area_m2.Value() + begin_fraction * delta_area_m2;
  const double area_at_end_m2 =
      begin_area_m2.Value() + end_fraction * delta_area_m2;
  const double partial_length_m =
      (end_fraction - begin_fraction) * segment_length_m.Value();
  return Units::VolumeM3(
      0.5 * (area_at_begin_m2 + area_at_end_m2) * partial_length_m);
}

SEP::Units::VolumeM3
SEP::FieldLine::FluxTubeGeometryCore::IntegrateSimpsonArea(
    Units::AreaM2 begin_area_m2,
    Units::AreaM2 middle_area_m2,
    Units::AreaM2 end_area_m2,
    Units::LengthM interval_length_m) {
  RequireFinitePositive(begin_area_m2.Value(),
                        "begin area must be finite and positive");
  RequireFinitePositive(middle_area_m2.Value(),
                        "middle area must be finite and positive");
  RequireFinitePositive(end_area_m2.Value(),
                        "end area must be finite and positive");
  RequireFinitePositive(interval_length_m.Value(),
                        "interval length must be finite and positive");
  return Units::VolumeM3(interval_length_m.Value() *
      (begin_area_m2.Value() + 4.0*middle_area_m2.Value() +
       end_area_m2.Value()) / 6.0);
}

SEP::Units::VolumeM3 SEP::FieldLine::FluxTubeGeometryCore::SweptVolume(
    Units::AreaM2 area_m2,
    Units::SpeedMPerS normal_speed_m_s,
    Units::TimeS time_s) {
  RequireFinitePositive(area_m2.Value(), "area must be finite and positive");
  RequireFinitePositive(normal_speed_m_s.Value(),
                        "speed must be finite and positive");
  RequireFinitePositive(time_s.Value(), "time must be finite and positive");
  return Units::VolumeM3(area_m2.Value() * normal_speed_m_s.Value() *
                         time_s.Value());
}

double SEP::FieldLine::FluxTubeGeometryCore::InjectedPhysicalParticleCount(
    Units::NumberDensityPerM3 density_per_m3,
    Units::VolumeM3 swept_volume_m3,
    double injection_efficiency) {
  if (!std::isfinite(density_per_m3.Value()) || density_per_m3.Value() < 0.0 ||
      !std::isfinite(swept_volume_m3.Value()) || swept_volume_m3.Value() < 0.0 ||
      !std::isfinite(injection_efficiency) || injection_efficiency < 0.0 ||
      injection_efficiency > 1.0) {
    throw std::invalid_argument("source density, volume, and efficiency are invalid");
  }
  return density_per_m3.Value() * swept_volume_m3.Value() *
         injection_efficiency;
}

double SEP::FieldLine::FluxTubeGeometryCore::MacroparticleWeight(
    double physical_particle_count,
    long macroparticle_count) {
  if (!std::isfinite(physical_particle_count) ||
      physical_particle_count < 0.0 || macroparticle_count <= 0) {
    throw std::invalid_argument("physical and macroparticle counts are invalid");
  }
  return physical_particle_count / static_cast<double>(macroparticle_count);
}

double SEP::FieldLine::FluxTubeGeometryCore::SpectralNumberDensityPerJ(
    Units::NumberDensityPerM3 density_per_m3,
    double normalized_probability_per_J) {
  if (!std::isfinite(density_per_m3.Value()) || density_per_m3.Value() < 0.0 ||
      !std::isfinite(normalized_probability_per_J) ||
      normalized_probability_per_J < 0.0) {
    throw std::invalid_argument("density and spectral probability are invalid");
  }
  return density_per_m3.Value() * normalized_probability_per_J;
}
