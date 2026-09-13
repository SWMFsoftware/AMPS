#include "sep.h"

#include <cmath>
#include <stdexcept>

namespace {

// The historical implementation implicitly assigned a one-metre radius to the
// first point.  Preserve that normalization explicitly as an SI reference area
// while replacing the dimensionally incorrect "radius grows as r squared"
// rule.  Applications with a measured tube area should set it before startup.
double g_reference_area_m2 = Pi;
SEP::FieldLine::FluxTubeGeometry::ExplicitAreaProfile g_explicit_area = NULL;

double MagneticMagnitudeT(PIC::FieldLine::cFieldLineVertex* vertex) {
  if (vertex == NULL) return 0.0;
  const double* B_T = vertex->GetMagneticField();
  if (B_T == NULL) return 0.0;
  const double magnitude_T =
      std::sqrt(B_T[0]*B_T[0] + B_T[1]*B_T[1] + B_T[2]*B_T[2]);
  return std::isfinite(magnitude_T) ? magnitude_T : 0.0;
}

double ExplicitAreaM2(const double* x_m, int iFieldLine) {
  if (g_explicit_area == NULL || x_m == NULL) return 0.0;
  const double area_m2 = g_explicit_area(x_m, iFieldLine);
  if (!std::isfinite(area_m2) || area_m2 <= 0.0) {
    throw std::runtime_error(
        "FluxTubeGeometry explicit area profile returned an invalid area");
  }
  return area_m2;
}

}  // namespace

void SEP::FieldLine::FluxTubeGeometry::SetReferenceAreaM2(
    double reference_area_m2) {
  if (!std::isfinite(reference_area_m2) || reference_area_m2 <= 0.0) {
    throw std::invalid_argument(
        "FluxTubeGeometry reference area must be finite and positive");
  }
  g_reference_area_m2 = reference_area_m2;
}

void SEP::FieldLine::FluxTubeGeometry::SetExplicitAreaProfile(
    ExplicitAreaProfile profile) {
  if (profile == NULL) {
    throw std::invalid_argument(
        "FluxTubeGeometry explicit area profile cannot be null");
  }
  g_explicit_area = profile;
}

void SEP::FieldLine::FluxTubeGeometry::ClearExplicitAreaProfile() {
  g_explicit_area = NULL;
}

double SEP::FieldLine::FluxTubeGeometry::AreaAtVertexM2(
    PIC::FieldLine::cFieldLineVertex* vertex,
    int iFieldLine) {
  if (vertex == NULL || iFieldLine < 0 ||
      iFieldLine >= PIC::FieldLine::nFieldLine) {
    throw std::invalid_argument("FluxTubeGeometry vertex or field-line id is invalid");
  }

  PIC::FieldLine::cFieldLineSegment* first_segment =
      PIC::FieldLine::FieldLinesAll[iFieldLine].GetFirstSegment();
  if (first_segment == NULL) {
    throw std::runtime_error("FluxTubeGeometry field line has no segments");
  }

  const double reference_abs_B_T = MagneticMagnitudeT(first_segment->GetBegin());
  const double local_abs_B_T = MagneticMagnitudeT(vertex);
  if (reference_abs_B_T > 0.0 && local_abs_B_T > 0.0) {
    return FluxTubeGeometryCore::AreaFromMagneticFlux(
        Units::AreaM2(g_reference_area_m2),
        Units::MagneticFieldT(reference_abs_B_T),
        Units::MagneticFieldT(local_abs_B_T)).Value();
  }

  // A missing magnetic magnitude does not justify silently fabricating an r^2
  // profile.  Requiring an explicit area callback makes the physical closure a
  // visible application configuration decision.
  const double explicit_area_m2 = ExplicitAreaM2(vertex->GetX(), iFieldLine);
  if (explicit_area_m2 > 0.0) return explicit_area_m2;
  throw std::runtime_error(
      "FluxTubeGeometry requires finite |B| or an explicit area profile");
}

double SEP::FieldLine::FluxTubeGeometry::AreaAtSegmentFractionM2(
    PIC::FieldLine::cFieldLineSegment* segment,
    int iFieldLine,
    double fraction) {
  if (segment == NULL || !std::isfinite(fraction) || fraction < 0.0 ||
      fraction > 1.0 || iFieldLine < 0 ||
      iFieldLine >= PIC::FieldLine::nFieldLine) {
    throw std::invalid_argument("FluxTubeGeometry segment fraction is invalid");
  }
  PIC::FieldLine::cFieldLineSegment* first_segment =
      PIC::FieldLine::FieldLinesAll[iFieldLine].GetFirstSegment();
  const double reference_abs_B_T =
      first_segment ? MagneticMagnitudeT(first_segment->GetBegin()) : 0.0;
  double local_B_T[3] = {0.0,0.0,0.0};
  segment->GetMagneticField(fraction, local_B_T);
  const double local_abs_B_T = std::sqrt(
      local_B_T[0]*local_B_T[0] + local_B_T[1]*local_B_T[1] +
      local_B_T[2]*local_B_T[2]);

  if (reference_abs_B_T > 0.0 && std::isfinite(local_abs_B_T) &&
      local_abs_B_T > 0.0) {
    return FluxTubeGeometryCore::AreaFromMagneticFlux(
        Units::AreaM2(g_reference_area_m2),
        Units::MagneticFieldT(reference_abs_B_T),
        Units::MagneticFieldT(local_abs_B_T)).Value();
  }

  double x_m[3] = {0.0,0.0,0.0};
  segment->GetCartesian(x_m, fraction);
  const double explicit_area_m2 = ExplicitAreaM2(x_m, iFieldLine);
  if (explicit_area_m2 > 0.0) return explicit_area_m2;
  throw std::runtime_error(
      "FluxTubeGeometry requires finite |B| or an explicit area profile");
}

double SEP::FieldLine::FluxTubeGeometry::PartialSegmentVolumeM3(
    PIC::FieldLine::cFieldLineSegment* segment,
    int iFieldLine,
    double begin_fraction,
    double end_fraction) {
  if (segment == NULL) {
    throw std::invalid_argument("FluxTubeGeometry segment cannot be null");
  }
  if (!std::isfinite(begin_fraction) || !std::isfinite(end_fraction) ||
      begin_fraction < 0.0 || end_fraction > 1.0 ||
      end_fraction <= begin_fraction) {
    throw std::invalid_argument(
        "FluxTubeGeometry fractions must satisfy 0 <= begin < end <= 1");
  }

  // Simpson integration samples the actual flux-conserving area at both ends
  // and the midpoint.  It therefore converges for curved B(s) without replacing
  // A|B|=constant by an unrelated conical-radius assumption.
  const double middle_fraction = 0.5*(begin_fraction+end_fraction);
  const double begin_area_m2 =
      AreaAtSegmentFractionM2(segment,iFieldLine,begin_fraction);
  const double middle_area_m2 =
      AreaAtSegmentFractionM2(segment,iFieldLine,middle_fraction);
  const double end_area_m2 =
      AreaAtSegmentFractionM2(segment,iFieldLine,end_fraction);
  const double partial_length_m =
      (end_fraction-begin_fraction)*segment->GetLength();
  return FluxTubeGeometryCore::IntegrateSimpsonArea(
      Units::AreaM2(begin_area_m2), Units::AreaM2(middle_area_m2),
      Units::AreaM2(end_area_m2),
      Units::LengthM(partial_length_m)).Value();
}

double SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(
    PIC::FieldLine::cFieldLineSegment* segment,
    int iFieldLine) {
  return PartialSegmentVolumeM3(segment, iFieldLine, 0.0, 1.0);
}

double SEP::FieldLine::FluxTubeGeometry::SweptVolumeM3(
    PIC::FieldLine::cFieldLineSegment* segment,
    int iFieldLine,
    double fraction,
    double normal_speed_m_s,
    double time_s) {
  return FluxTubeGeometryCore::SweptVolume(
      Units::AreaM2(AreaAtSegmentFractionM2(segment, iFieldLine, fraction)),
      Units::SpeedMPerS(normal_speed_m_s), Units::TimeS(time_s)).Value();
}
