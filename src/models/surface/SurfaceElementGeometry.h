/*
 * Dependency-light geometry for a latitude/longitude spherical surface.
 *
 * This header deliberately contains no PIC, MPI, or configured-application
 * dependency.  Surface source models need the same geometry in preprocessing,
 * normalization, and tests, where constructing the full AMR internal-boundary
 * object would introduce unrelated global state.
 */
#ifndef AMPS_MODELS_SURFACE_ELEMENT_GEOMETRY_H
#define AMPS_MODELS_SURFACE_ELEMENT_GEOMETRY_H

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace Surface {
namespace Geometry {

enum class ZenithSpacing {
  EqualCosine,
  EqualAngle
};

struct ElementIndex {
  std::size_t zenith;
  std::size_t azimuth;
};

/**
 * Analytic partition of a sphere into zenith-by-azimuth surface elements.
 *
 * The sphere radius is in metres and each returned area is in m^2.  For an
 * element bounded by cos(theta)_north and cos(theta)_south and by an azimuth
 * width d(phi), the exact curved area is
 *
 *   A = R^2 d(phi) [cos(theta)_north - cos(theta)_south].
 *
 * EqualCosine therefore produces equal-area elements.  EqualAngle uses
 * theta edges separated by pi/Nz and has smaller polar cells.  The flattened
 * ordering follows the existing AMPS spherical-boundary convention,
 * `flat = zenith + Nz * azimuth`; changing this order would silently attach
 * source strengths to the wrong physical cells.
 *
 * A zero/nonfinite radius or zero resolution is invalid geometry and raises
 * std::invalid_argument.  Invalid indices raise std::out_of_range.  Those
 * typed failures are preferable in dependency-light preprocessing to a NaN
 * area propagating into a source probability table.
 */
class SphericalElementGrid {
 public:
  SphericalElementGrid(std::size_t zenith_count,
                       std::size_t azimuth_count,
                       double radius_m,
                       ZenithSpacing spacing)
      : zenith_count_(zenith_count), azimuth_count_(azimuth_count),
        radius_m_(radius_m), spacing_(spacing) {
    if (zenith_count_ == 0 || azimuth_count_ == 0)
      throw std::invalid_argument("surface resolution must be positive");
    if (!(radius_m_ > 0.0) || !std::isfinite(radius_m_))
      throw std::invalid_argument("sphere radius must be finite and positive");
    if (spacing_ != ZenithSpacing::EqualCosine &&
        spacing_ != ZenithSpacing::EqualAngle)
      throw std::invalid_argument("unknown spherical zenith spacing");
    if (zenith_count_ > std::numeric_limits<std::size_t>::max() /
                             azimuth_count_)
      throw std::overflow_error("surface element count overflows size_t");
  }

  std::size_t zenith_count() const { return zenith_count_; }
  std::size_t azimuth_count() const { return azimuth_count_; }
  std::size_t size() const { return zenith_count_ * azimuth_count_; }
  double radius_m() const { return radius_m_; }

  std::size_t flatten(std::size_t zenith, std::size_t azimuth) const {
    validate_pair(zenith, azimuth);
    return zenith + zenith_count_ * azimuth;
  }

  ElementIndex unflatten(std::size_t element) const {
    if (element >= size())
      throw std::out_of_range("surface element index is out of range");
    ElementIndex result;
    result.azimuth = element / zenith_count_;
    result.zenith = element - zenith_count_ * result.azimuth;
    return result;
  }

  double area(std::size_t zenith, std::size_t azimuth) const {
    validate_pair(zenith, azimuth);
    (void)azimuth;  // Every azimuth bin has the same angular width.

    const double pi = 3.141592653589793238462643383279502884;
    const double dphi = 2.0 * pi / static_cast<double>(azimuth_count_);
    double cosine_north;
    double cosine_south;

    if (spacing_ == ZenithSpacing::EqualCosine) {
      const double dcos = 2.0 / static_cast<double>(zenith_count_);
      cosine_north = 1.0 - dcos * static_cast<double>(zenith);
      cosine_south = 1.0 - dcos * static_cast<double>(zenith + 1);
    }
    else {
      const double dtheta = pi / static_cast<double>(zenith_count_);
      cosine_north = std::cos(dtheta * static_cast<double>(zenith));
      cosine_south = std::cos(dtheta * static_cast<double>(zenith + 1));
    }

    return radius_m_ * radius_m_ * dphi *
           (cosine_north - cosine_south);
  }

  double area(std::size_t element) const {
    const ElementIndex index = unflatten(element);
    return area(index.zenith, index.azimuth);
  }

  /** Unit outward normal at the angular centre of an element. */
  std::array<double, 3> centre_normal(std::size_t element) const {
    const ElementIndex index = unflatten(element);
    const double pi = 3.141592653589793238462643383279502884;
    const double phi = 2.0 * pi *
        (static_cast<double>(index.azimuth) + 0.5) /
        static_cast<double>(azimuth_count_);
    double cosine_theta;

    if (spacing_ == ZenithSpacing::EqualCosine) {
      cosine_theta = 1.0 - 2.0 *
          (static_cast<double>(index.zenith) + 0.5) /
          static_cast<double>(zenith_count_);
    }
    else {
      const double theta = pi *
          (static_cast<double>(index.zenith) + 0.5) /
          static_cast<double>(zenith_count_);
      cosine_theta = std::cos(theta);
    }

    // Guard the square root against a one-ulp overshoot at extreme, but valid,
    // resolutions.  This is a roundoff correction, not a geometry clamp.
    const double sine_theta = std::sqrt(
        std::fmax(0.0, 1.0 - cosine_theta * cosine_theta));
    return {{sine_theta * std::cos(phi), sine_theta * std::sin(phi),
             cosine_theta}};
  }

 private:
  void validate_pair(std::size_t zenith, std::size_t azimuth) const {
    if (zenith >= zenith_count_ || azimuth >= azimuth_count_)
      throw std::out_of_range("surface zenith/azimuth index is out of range");
  }

  std::size_t zenith_count_;
  std::size_t azimuth_count_;
  double radius_m_;
  ZenithSpacing spacing_;
};

}  // namespace Geometry
}  // namespace Surface

#endif
