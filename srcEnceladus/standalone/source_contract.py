"""Dependency-free source membership and spherical projection contracts.

Case names identify physical bases, not interchangeable lists of coefficients.
These canonical IDs match the published one-based Porco selection and the
explicit opt-in E17 candidate. Rates and fitted background terms never enter
this module; source positions are IAU_ENCELADUS Cartesian metres.
"""
import math

HYBRID = "tenishev2014_hybrid_8"
PORCO = "porco2014_98"
SOLAR = "portyankina2022_solar_41"
CANDIDATE = "porco2014_98_plus_e17_candidate"
HYBRID_IDS = (
    "jet:I_baghdad", "jet:II_damascus", "jet:III_damascus",
    "jet:IV_alexandria", "jet:V_cairo", "jet:VI_baghdad",
    "jet:VII_baghdad", "jet:VIII_cairo", "stripe:alexandria",
    "stripe:baghdad", "stripe:cairo", "stripe:damascus",
    "stripe:baghdad_extra", "stripe:cairo_extra", "stripe:damascus_extra")
SOLAR_IDS = (1, 4, 6, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 21, 24, 25,
             27, 28, 32, 35, 36, 37, 40, 54, 60, 61, 62, 63, 66, 67, 68,
             73, 84, 85, 86, 89, 90, 91, 92, 93, 94)
CASES = (HYBRID, PORCO, SOLAR, CANDIDATE)


def source_ids(case):
    """Return every required physical component, including zero-rate sources."""
    if case == HYBRID:
        return HYBRID_IDS
    if case not in CASES:
        raise ValueError("unsupported source case: " + str(case))
    ids = SOLAR_IDS if case == SOLAR else range(1, 99)
    result = tuple("jet:porco_{:03d}".format(i) for i in ids)
    return result + (("jet:e17_candidate_baghdad",) if case == CANDIDATE else ())


def project_vertices(points, radius):
    """Map x to R*x/|x| while preserving latitude/longitude and node order.

    R is a positive radius in metres. This is the user-selected local vertical
    on a sphere; no constant altitude translation or topographic normal is
    inferred. The same routine serves fitting and export so the arc-length
    weights are formed from the same projected vertices in both models.
    """
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError("reference radius must be finite and positive")
    projected = []
    for point in points:
        if len(point) != 3 or not all(math.isfinite(float(x)) for x in point):
            raise ValueError("stripe vertex must contain three finite coordinates")
        norm = math.sqrt(math.fsum(float(x) ** 2 for x in point))
        if norm == 0:
            raise ValueError("stripe vertex cannot be at the centre")
        vertex = [float(x) * (radius / norm) for x in point]
        if projected and vertex == projected[-1]:
            raise ValueError("projected stripe has a zero-length segment")
        projected.append(vertex)
    if len(projected) < 2:
        raise ValueError("stripe needs at least two vertices")
    return projected
