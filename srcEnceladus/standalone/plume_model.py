"""Steady, collisionless, gravity-free Enceladus plume kernels in SI units.

This module has no AMPS dependency. ``Q`` always means physical molecules/s;
it is not the unnormalised amplitude used by the historical fitting program.
"""
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erfc, gamma, hyperu

KB = 1.380649e-23
POINT_SOURCE_LAW = "tenishev2010_point"
SURFACE_FLUX_LAW = "surface_flux_maxwellian"
SOURCE_LAWS = (POINT_SOURCE_LAW, SURFACE_FLUX_LAW)


def unit(vector):
    x = np.asarray(vector, dtype=float)
    if x.shape != (3,) or not np.isfinite(x).all() or np.linalg.norm(x) == 0:
        raise ValueError("Expected a finite nonzero three-vector")
    return x / np.linalg.norm(x)


def spherical_source(radius_m, latitude_deg, west_longitude_deg,
                     tilt_deg, azimuth_deg):
    """Return position, radial normal, and jet axis in IAU_ENCELADUS.

    Latitude is signed planetocentric latitude, west longitude increases
    westward, tilt is from the radial normal, and azimuth is clockwise from
    local north. Thus the published positive south-latitude magnitudes must be
    supplied as negative values.
    """
    values = [radius_m, latitude_deg, west_longitude_deg, tilt_deg, azimuth_deg]
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite source geometry")
    if radius_m <= 0 or not -90 <= latitude_deg <= 90 or not 0 <= tilt_deg <= 90:
        raise ValueError("Invalid radius, latitude, or outward tilt")
    lat, lon, tilt, az = np.deg2rad(
        [latitude_deg, -west_longitude_deg, tilt_deg, azimuth_deg])
    normal = np.array([np.cos(lat)*np.cos(lon),
                       np.cos(lat)*np.sin(lon), np.sin(lat)])
    north = np.array([-np.sin(lat)*np.cos(lon),
                      -np.sin(lat)*np.sin(lon), np.cos(lat)])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    axis = (np.cos(tilt)*normal
            + np.sin(tilt)*(np.cos(az)*north + np.sin(az)*east))
    return radius_m*normal, normal, unit(axis)


def radial_moment(power, cosine, speed_ratio):
    """Return K_p(c,s)=integral_0^inf t^p exp(-t^2+2sct-s^2) dt."""
    if power not in (2, 3) or not np.isfinite(speed_ratio) or speed_ratio < 0:
        raise ValueError("Require p=2 or 3 and finite speed ratio >= 0")
    c = np.asarray(cosine, dtype=float)
    if not np.isfinite(c).all() or np.any(np.abs(c) > 1+1e-12):
        raise ValueError("Invalid direction cosine")
    c = np.clip(c, -1, 1)
    a = speed_ratio*c
    result = np.empty_like(a)

    mask = a >= -2
    z = a[mask]
    j0 = np.sqrt(np.pi)/2*erfc(-z)
    j1 = 0.5*np.exp(-z*z)+z*j0
    j2 = 0.5*j0+z*j1
    jp = j2 if power == 2 else j1+z*j2
    result[mask] = np.exp(-speed_ratio**2*(1-c[mask]**2))*jp

    z = a[~mask]
    result[~mask] = (np.exp(-speed_ratio**2)*gamma(power+1)/2**(power+1)
                     * hyperu((power+1)/2, 0.5, z*z))
    if not np.isfinite(result).all() or np.any(result < 0):
        raise ArithmeticError("Radial velocity moment failed")
    return result


@lru_cache(maxsize=16)
def angular_grid(order):
    x, w = leggauss(order)
    return (x+1)/2, w/2, np.pi*(x+1), np.pi*w


def _angular_integral(speed_ratio, axis_normal, source_law, order):
    mu, wmu, phi, wphi = angular_grid(order)
    c = (axis_normal*mu[:, None]
         + np.sqrt(max(0.0, 1-axis_normal**2))
         * np.sqrt(1-mu[:, None]**2)*np.cos(phi)[None, :])
    integrand = radial_moment(3, c, speed_ratio)
    if source_law == SURFACE_FLUX_LAW:
        integrand = mu[:, None]*integrand
    elif source_law != POINT_SOURCE_LAW:
        raise ValueError("Unknown source law")
    return float(wmu @ integrand @ wphi)


@lru_cache(maxsize=4096)
def normalization(speed_ratio, axis_normal, source_law=POINT_SOURCE_LAW):
    """Flux normalization over the outward body hemisphere."""
    if source_law not in SOURCE_LAWS:
        raise ValueError("Unknown source law: "+str(source_law))
    if not 0 <= speed_ratio <= 20 or not 0 <= axis_normal <= 1:
        raise ValueError("Supported speed ratio is [0,20]; axis must point out")
    previous = _angular_integral(speed_ratio, axis_normal, source_law, 24)
    for order in (48, 96, 192, 384):
        current = _angular_integral(speed_ratio, axis_normal, source_law, order)
        if current > 0 and abs(current-previous) <= 1e-9*current:
            return current
        previous = current
    raise ArithmeticError("Source normalization did not converge")


@dataclass(frozen=True)
class Plume:
    position_m: tuple
    normal: tuple
    axis: tuple
    temperature_K: float
    bulk_speed_m_s: float
    molecular_mass_kg: float
    source_law: str = POINT_SOURCE_LAW

    def __post_init__(self):
        p = np.asarray(self.position_m, dtype=float)
        if p.shape != (3,) or not np.isfinite(p).all() or np.linalg.norm(p) == 0:
            raise ValueError("Invalid source position")
        for value in (self.normal, self.axis):
            if abs(np.linalg.norm(value)-1) > 1e-10:
                raise ValueError("Normal and axis must be unit vectors")
            unit(value)
        values = [self.temperature_K, self.bulk_speed_m_s, self.molecular_mass_kg]
        if (not np.isfinite(values).all() or self.temperature_K <= 0
                or self.molecular_mass_kg <= 0 or self.bulk_speed_m_s < 0):
            raise ValueError("Require finite T,m>0 and u>=0")
        if np.dot(self.normal, self.axis) < -1e-12:
            raise ValueError("Source axis points into the body")
        if self.source_law not in SOURCE_LAWS:
            raise ValueError("Unknown source law")
        _ = self.norm

    @property
    def thermal_speed(self):
        return math.sqrt(2*KB*self.temperature_K/self.molecular_mass_kg)

    @property
    def speed_ratio(self):
        return self.bulk_speed_m_s/self.thermal_speed

    @property
    def norm(self):
        axis_normal = float(np.clip(np.dot(self.normal, self.axis), 0, 1))
        return normalization(self.speed_ratio, axis_normal, self.source_law)

    def _rays(self, points):
        x = np.asarray(points, dtype=float)
        if x.shape[-1:] != (3,) or not np.isfinite(x).all():
            raise ValueError("Finite positions with final dimension 3 required")
        displacement = x-np.asarray(self.position_m)
        radius = np.linalg.norm(displacement, axis=-1)
        if np.any(radius == 0):
            raise ValueError("Point-source singularity at a source location")
        ray = displacement/radius[..., None]
        mu = ray @ np.asarray(self.normal)
        return radius, ray @ np.asarray(self.axis), mu

    def _law_factor(self, mu):
        if self.source_law == POINT_SOURCE_LAW:
            return np.ones_like(mu)
        return np.maximum(mu, 0.0)

    def density_per_rate(self, points):
        """Number density [m^-3] per unit total production [s^-1]."""
        radius, cosine, mu = self._rays(points)
        outward = mu >= 0
        return np.where(
            outward,
            self._law_factor(mu)*radial_moment(2, cosine, self.speed_ratio)
            / (self.norm*self.thermal_speed*radius*radius),
            0.0)

    def radial_flux_per_rate(self, points):
        """Source-centered radial number flux per unit production."""
        radius, cosine, mu = self._rays(points)
        outward = mu >= 0
        return np.where(
            outward,
            self._law_factor(mu)*radial_moment(3, cosine, self.speed_ratio)
            / (self.norm*radius*radius),
            0.0)


@dataclass(frozen=True)
class CompositeSource:
    """One fitted production rate distributed over normalized point sources."""
    components: tuple
    weights: tuple

    def __post_init__(self):
        if not self.components or len(self.components) != len(self.weights):
            raise ValueError("Composite source needs equally sized components/weights")
        weights = np.asarray(self.weights, dtype=float)
        if (not np.isfinite(weights).all() or np.any(weights < 0)
                or not np.isclose(weights.sum(), 1, rtol=0, atol=1e-12)):
            raise ValueError("Composite weights must be nonnegative and sum to one")

    def density_per_rate(self, points):
        return sum(w*p.density_per_rate(points)
                   for p, w in zip(self.components, self.weights))

    def radial_flux_per_rate(self, points):
        return sum(w*p.radial_flux_per_rate(points)
                   for p, w in zip(self.components, self.weights))


def line_source(points_m, temperature_K, bulk_speed_m_s, molecular_mass_kg,
                source_law=POINT_SOURCE_LAW):
    """Construct a vertically directed, uniform-production polyline source."""
    points = np.asarray(points_m, dtype=float)
    if (points.ndim != 2 or points.shape[1] != 3 or len(points) < 2
            or not np.isfinite(points).all()):
        raise ValueError("A line source needs at least two finite 3-D nodes")
    segment_length = np.linalg.norm(np.diff(points, axis=0), axis=1)
    if np.any(segment_length <= 0):
        raise ValueError("Line source contains a zero-length segment")
    weights = np.zeros(len(points))
    weights[:-1] += segment_length/2
    weights[1:] += segment_length/2
    weights /= weights.sum()
    components = []
    for point in points:
        normal = unit(point)
        components.append(Plume(tuple(point), tuple(normal), tuple(normal),
                                temperature_K, bulk_speed_m_s,
                                molecular_mass_kg, source_law))
    return CompositeSource(tuple(components), tuple(weights))


def design_matrix(sources, points, reference_radius_m, background="both",
                  power_law_exponent=-2.0):
    """Build density columns for nonnegative linear inference.

    Source columns multiply Q/1e28. A constant background column multiplies
    c0 [m^-3]. A power-law column multiplies the density at the reference
    radius, so n_bg = a*(r/reference_radius)**power_law_exponent.
    """
    x = np.asarray(points, dtype=float)
    if x.ndim != 2 or x.shape[1] != 3 or not np.isfinite(x).all():
        raise ValueError("points must be a finite (N,3) array")
    if not np.isfinite(reference_radius_m) or reference_radius_m <= 0:
        raise ValueError("Positive reference radius required")
    radius = np.linalg.norm(x, axis=1)
    if np.any(radius <= reference_radius_m):
        raise ValueError("Observation is inside the adopted reference sphere")
    if background not in ("none", "constant", "power_law", "both"):
        raise ValueError("Unknown background selection")
    if not np.isfinite(power_law_exponent) or power_law_exponent >= 0:
        raise ValueError("Background power-law exponent must be finite and negative")
    columns = [source.density_per_rate(x)*1e28 for source in sources]
    if background in ("constant", "both"):
        columns.append(np.ones(len(x)))
    if background in ("power_law", "both"):
        columns.append((radius/reference_radius_m)**power_law_exponent)
    if not columns:
        raise ValueError("At least one source or background term is required")
    return np.column_stack(columns)


def line_of_sight_design_matrix(
        sources, tangent_points_m, directions, reference_radius_m,
        background="power_law", power_law_exponent=-2.0, order=64,
        batch_size=128, observer_coordinate_from_tangent_m=None):
    """Integrate density columns along occultation rays.

    ``tangent_points_m`` are closest-approach points to the body centre and
    ``directions`` are unit vectors along each ray. Source columns multiply
    ``Q/1e28`` and therefore have units m^-2. If observer distance is omitted,
    the integration is over an infinite chord. Otherwise the signed coordinate
    of the finite observer along the source-pointing ray sets the lower limit;
    it is negative when the tangent point lies ahead of the observer and
    positive when the infinite-line tangent lies behind it.
    A power-law background column
    multiplies density at the reference radius. A spatially constant density
    is deliberately unsupported because its infinite chord integral diverges.

    The substitution ``s=b*tan(theta)`` maps each infinite line to a finite
    Gauss-Legendre interval and preserves the 1/r^2 tails analytically well.
    """
    tangent = np.asarray(tangent_points_m, dtype=float)
    ray = np.asarray(directions, dtype=float)
    if (tangent.ndim != 2 or tangent.shape[1] != 3
            or ray.shape != tangent.shape or not np.isfinite(tangent).all()
            or not np.isfinite(ray).all()):
        raise ValueError("Finite tangent points and directions must be (N,3)")
    if (not np.isfinite(reference_radius_m) or reference_radius_m <= 0
            or not isinstance(order, int) or order < 8
            or not isinstance(batch_size, int) or batch_size < 1):
        raise ValueError("Invalid radius or quadrature controls")
    if background not in ("none", "power_law"):
        raise ValueError("Occultation chords support only none or power_law background")
    if not np.isfinite(power_law_exponent) or power_law_exponent >= -1:
        raise ValueError("Infinite chord requires a power-law exponent below -1")
    norm = np.linalg.norm(ray, axis=1)
    if np.any(norm == 0):
        raise ValueError("Ray direction cannot be zero")
    ray = ray/norm[:, None]
    impact = np.linalg.norm(tangent, axis=1)
    if np.any(impact <= reference_radius_m):
        raise ValueError("Occultation chord intersects the reference sphere")
    perpendicular_error = np.abs(np.sum(tangent*ray, axis=1))/impact
    if np.any(perpendicular_error > 1.0e-7):
        raise ValueError("Ray direction is not perpendicular to tangent radius")

    if observer_coordinate_from_tangent_m is None:
        observer_coordinate = np.full(len(tangent), -np.inf)
    else:
        observer_coordinate = np.asarray(
            observer_coordinate_from_tangent_m, dtype=float)
        if (observer_coordinate.shape != (len(tangent),)
                or not np.isfinite(observer_coordinate).all()):
            raise ValueError("Finite observer coordinates must match all rays")

    nodes, node_weights = leggauss(order)
    column_count = len(sources)+(1 if background == "power_law" else 0)
    if column_count == 0:
        raise ValueError("At least one source or background term is required")
    result = np.empty((len(tangent), column_count), dtype=float)

    for start in range(0, len(tangent), batch_size):
        stop = min(start+batch_size, len(tangent))
        b = impact[start:stop]
        finite = np.isfinite(observer_coordinate[start:stop])
        theta_low = np.full(stop-start, -0.5*np.pi)
        theta_low[finite] = np.arctan(
            observer_coordinate[start:stop][finite]/b[finite])
        theta_high = np.full(stop-start, 0.5*np.pi)
        theta_half = 0.5*(theta_high-theta_low)
        theta_mid = 0.5*(theta_high+theta_low)
        theta = theta_mid[:, None]+theta_half[:, None]*nodes[None, :]
        distance = b[:, None]*np.tan(theta)
        weights = (b[:, None]*theta_half[:, None]
                   /(np.cos(theta)**2)*node_weights[None, :])
        points = (tangent[start:stop, None, :]
                  + distance[:, :, None]*ray[start:stop, None, :])
        flat = points.reshape(-1, 3)
        for index, source in enumerate(sources):
            density = source.density_per_rate(flat).reshape(stop-start, order)
            result[start:stop, index] = np.sum(density*weights, axis=1)*1.0e28
        if background == "power_law":
            radius = np.linalg.norm(points, axis=2)
            density = (radius/reference_radius_m)**power_law_exponent
            result[start:stop, -1] = np.sum(density*weights, axis=1)
    if not np.isfinite(result).all() or np.any(result < 0):
        raise ArithmeticError("Line-of-sight integration failed")
    return result


def outward_ray_design_matrix(
        sources, origins_m, directions, reference_radius_m,
        background="power_law", power_law_exponent=-2.0, order=64,
        batch_size=128):
    """Integrate density from finite observers along outward UVIS rays."""
    origin = np.asarray(origins_m, dtype=float)
    ray = np.asarray(directions, dtype=float)
    if (origin.ndim != 2 or origin.shape[1] != 3 or ray.shape != origin.shape
            or not np.isfinite(origin).all() or not np.isfinite(ray).all()):
        raise ValueError("Finite origins and directions must be (N,3)")
    if (not np.isfinite(reference_radius_m) or reference_radius_m <= 0
            or not isinstance(order, int) or order < 8
            or not isinstance(batch_size, int) or batch_size < 1):
        raise ValueError("Invalid radius or quadrature controls")
    if background not in ("none", "power_law"):
        raise ValueError("Occultation rays support only none or power_law background")
    if not np.isfinite(power_law_exponent) or power_law_exponent >= -1:
        raise ValueError("Infinite ray requires a power-law exponent below -1")
    norm = np.linalg.norm(ray, axis=1)
    if np.any(norm == 0):
        raise ValueError("Ray direction cannot be zero")
    ray = ray/norm[:, None]
    radius0 = np.linalg.norm(origin, axis=1)
    if np.any(radius0 <= reference_radius_m):
        raise ValueError("UVIS observer is inside the reference sphere")
    if np.any(np.sum(origin*ray, axis=1) < -1.0e-8*radius0):
        raise ValueError("outward_ray_design_matrix received an inward ray")

    nodes, node_weights = leggauss(order)
    scaled_distance = (1.0+nodes)/(1.0-nodes)
    scaled_jacobian = 2.0/(1.0-nodes)**2
    column_count = len(sources)+(1 if background == "power_law" else 0)
    if column_count == 0:
        raise ValueError("At least one source or background term is required")
    result = np.empty((len(origin), column_count), dtype=float)
    for start in range(0, len(origin), batch_size):
        stop = min(start+batch_size, len(origin))
        scale = radius0[start:stop]
        distance = scale[:, None]*scaled_distance[None, :]
        weights = scale[:, None]*scaled_jacobian[None, :]*node_weights[None, :]
        points = (origin[start:stop, None, :]
                  + distance[:, :, None]*ray[start:stop, None, :])
        flat = points.reshape(-1, 3)
        for index, source in enumerate(sources):
            density = source.density_per_rate(flat).reshape(stop-start, order)
            result[start:stop, index] = np.sum(density*weights, axis=1)*1.0e28
        if background == "power_law":
            radius = np.linalg.norm(points, axis=2)
            density = (radius/reference_radius_m)**power_law_exponent
            result[start:stop, -1] = np.sum(density*weights, axis=1)
    if not np.isfinite(result).all() or np.any(result < 0):
        raise ArithmeticError("Outward-ray integration failed")
    return result
