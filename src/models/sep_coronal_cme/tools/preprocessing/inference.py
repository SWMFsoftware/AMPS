"""Explicit offline inference operators with full covariance propagation.

These operators consume already decoded, metadata-owned observations. No image
projection, instrument coordinate system, or unit convention is guessed. The
ellipse is a projected contour fit; it is NOT claimed to reconstruct a 3-D CME
from one viewpoint. A 3-D ellipsoid fit uses independently reconstructed points
and a declared center/orientation. Those construction authorities are retained.
"""
from __future__ import annotations
import math
from .core import (require, finite, linear_fit, covariance, propagate, matvec,
                   transpose, matmul, load_json, sha256_bytes, validate_metadata,
                   publish_bundle, digest)


def real_harmonic(l, m, sine, theta, phi):
    """Same real orthonormal/Condon-Shortley convention as the C++ PFSS API."""
    x = math.cos(theta)
    p = 1.0
    for k in range(1, m+1):
        p *= -(2*k-1)*math.sqrt(max(0.0, 1-x*x))
    if l > m:
        previous, current = p, x*(2*m+1)*p
        for n in range(m+2, l+1):
            previous, current = current, ((2*n-1)*x*current-(n+m-1)*previous)/(n-m)
        p = current
    norm = math.sqrt((2*l+1)/(4*math.pi)*(1 if m == 0 else 2)*
                     math.factorial(l-m)/math.factorial(l+m))
    return norm*p*(math.sin(m*phi) if sine else math.cos(m*phi))


def fit_magnetogram(source, parameters):
    degree = parameters["maximum_degree"]
    require(isinstance(degree, int) and 1 <= degree <= 12, "unsupported harmonic fit degree")
    theta, phi, field = source["theta_rad"], source["phi_rad"], source["radial_field_T"]
    require(len(theta) == len(phi) == len(field), "magnetogram coordinate/value dimensions differ")
    finite(theta+phi+field)
    require(all(0 <= t <= math.pi for t in theta), "invalid magnetogram colatitude")
    modes = [(l, m, sine) for l in range(1, degree+1) for m in range(l+1)
             for sine in ([False] if m == 0 else [False, True])]
    # Fit a monopole simultaneously, then remove only that explicit mode. This
    # estimates its uncertainty rather than subtracting an unrecorded map mean.
    modes = [(0, 0, False)]+modes
    design = [[real_harmonic(l, m, sine, t, p) for l, m, sine in modes] for t, p in zip(theta, phi)]
    fit = linear_fit(design, field, source["covariance_T2"])
    require(abs(fit["parameters"][0]) <= parameters["maximum_removed_monopole_T"],
            "magnetogram flux-balance correction exceeds preregistered bound")
    return {"schema": "sep-real-harmonic-coefficients-v1", "normalization": "real-orthonormal-condon-shortley",
            "modes": [{"degree": l, "order": m, "sine": sine, "coefficient_T": value}
                      for (l, m, sine), value in zip(modes[1:], fit["parameters"][1:])],
            "covariance_T2": [row[1:] for row in fit["covariance"][1:]],
            "removed_monopole_T": fit["parameters"][0], "fit_chi2": fit["chi2"]}


def ellipse_parameters(coefficients):
    a, b, c, d, e = coefficients
    determinant = a*c-b*b/4
    require(a*c > 0 and determinant > 0, "conic is not a bounded ellipse")
    center = [(-c*d+b*e/2)/(2*determinant), (b*d/2-a*e)/(2*determinant)]
    level = 1+a*center[0]**2+b*center[0]*center[1]+c*center[1]**2
    if a < 0:
        a, b, c, level = -a, -b, -c, -level
    splitting = math.hypot(a-c, b)
    small, large = (a+c-splitting)/2, (a+c+splitting)/2
    require(small > 0 and level > 0, "degenerate ellipse")
    angle = (0.5*math.atan2(b, a-c)+math.pi/2) % math.pi
    return center+[math.sqrt(level/small), math.sqrt(level/large), angle]


def fit_image_ellipse(source, parameters):
    points = source["contour_xy"]
    scale = parameters["coordinate_scale_m"]
    require(scale > 0 and math.isfinite(scale), "image needs explicit metric/pixel scale")
    finite(v for point in points for v in point)
    require(all(len(point) == 2 for point in points), "ellipse contour must be 2-D")
    # Normalize coordinates to avoid ill-conditioned meter-scale conic columns.
    normalization = max(abs(v) for point in points for v in point)
    require(normalization > 0, "empty/degenerate image contour")
    design = [[(x/normalization)**2, x*y/normalization**2, (y/normalization)**2,
               x/normalization, y/normalization] for x, y in points]
    fit = linear_fit(design, [1.0]*len(points), source["implicit_equation_covariance"])
    physical = ellipse_parameters(fit["parameters"])
    jacobian = [[] for _ in physical]
    for i, value in enumerate(fit["parameters"]):
        step = 1e-6*max(abs(value), 1e-2)
        plus, minus = list(fit["parameters"]), list(fit["parameters"])
        plus[i] += step; minus[i] -= step
        fplus, fminus = ellipse_parameters(plus), ellipse_parameters(minus)
        for j in range(5):
            diff = fplus[j]-fminus[j]
            if j == 4:
                diff = (diff+math.pi/2) % math.pi-math.pi/2
            jacobian[j].append(diff/(2*step))
    factors = [normalization*scale]*4+[1.0]
    return {"schema": "sep-projected-ellipse-fit-v1",
            "parameters": [x*f for x, f in zip(physical, factors)],
            "parameter_names": ["center_x_m", "center_y_m", "major_axis_m", "minor_axis_m", "angle_rad"],
            "covariance": propagate([[x*factors[i] for x in row] for i, row in enumerate(jacobian)], fit["covariance"]),
            "projection": parameters["projection"], "fit_chi2": fit["chi2"]}


def validate_rotation(rotation):
    require(len(rotation) == 3 and all(len(row) == 3 for row in rotation), "3-D rotation matrix required")
    finite(v for row in rotation for v in row)
    gram = matmul(rotation, transpose(rotation))
    require(all(abs(gram[i][j]-float(i == j)) < 1e-10 for i in range(3) for j in range(3)), "frame transform is not orthonormal")
    a = rotation
    determinant = a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])
    require(abs(determinant-1) < 1e-10, "frame transform is a reflection")


def transform_ephemeris(source, parameters):
    rotation = parameters["rotation"]
    validate_rotation(rotation)
    require(parameters["source_frame"] != "" and parameters["target_frame"] != "", "missing frame names")
    times = source["time_s"]
    finite(times)
    require(len(times) >= 1 and all(a < b for a, b in zip(times, times[1:])), "ephemeris time order/duplicates")
    require(len(source["position_m"]) == len(source["velocity_m_per_s"]) == len(source["covariance"]) == len(times), "ephemeris array dimensions differ")
    jacobian = [row+[0.0]*3 for row in rotation]+[[0.0]*3+row for row in rotation]
    positions, velocities, covariances = [], [], []
    for x, u, c in zip(source["position_m"], source["velocity_m_per_s"], source["covariance"]):
        finite(x+u); require(len(x) == len(u) == 3, "ephemeris needs full position/velocity vectors")
        covariance(c, 6)
        positions.append(matvec(rotation, x)); velocities.append(matvec(rotation, u))
        covariances.append(propagate(jacobian, c))
    return {"schema": "sep-observer-ephemeris-v1", "time_s": times, "position_m": positions,
            "velocity_m_per_s": velocities, "covariance": covariances, "transform": parameters}


def fit_power_law(source, parameters):
    radius, values = source["radius_m"], source["values"]
    reference = parameters["reference_radius_m"]
    finite(radius+values)
    require(reference > 0 and len(radius) == len(values) and all(r > 0 for r in radius) and all(v > 0 for v in values), "power-law support/values invalid")
    c = covariance(source["covariance"], len(values))
    require(all(math.sqrt(c[i][i])/values[i] < 0.5 for i in range(len(values))), "log-linear first-order uncertainty is too large")
    log_cov = [[c[i][j]/(values[i]*values[j]) for j in range(len(values))] for i in range(len(values))]
    fit = linear_fit([[1.0, math.log(r/reference)] for r in radius], [math.log(v) for v in values], log_cov)
    return {"schema": "sep-positive-power-law-profile-v1", "reference_radius_m": reference,
            "log_reference_value": fit["parameters"][0], "exponent": fit["parameters"][1],
            "log_parameter_covariance": fit["covariance"], "support_m": [min(radius), max(radius)],
            "uncertainty_approximation": "first-order-log-linear", "fit_chi2": fit["chi2"]}


def evaluate_power_law(profile, radii):
    # Prepared candidate builders can export exact density at the requested
    # radio/front nodes instead of forcing a coronal density into a power law.
    # Covariance already includes their declared geometry/closure uncertainty.
    if profile["schema"] == "sep-front-sampled-density-v1":
        require(profile["front_radii_m"] == radii, "sampled density belongs to another front support")
        values = profile["density_m3"]
        finite(values)
        require(len(values) == len(radii) and all(n > 0 for n in values), "sampled density is incomplete/nonpositive")
        covariance(profile["covariance_m6"], len(values))
        return values, profile["covariance_m6"]
    require(profile["schema"] == "sep-positive-power-law-profile-v1", "unsupported density inference profile")
    require(all(profile["support_m"][0] <= r <= profile["support_m"][1] for r in radii), "profile extrapolation outside declared support")
    design = [[1.0, math.log(r/profile["reference_radius_m"])] for r in radii]
    values = [math.exp(profile["log_reference_value"]+profile["exponent"]*row[1]) for row in design]
    jacobian = [[v*x for x in row] for v, row in zip(values, design)]
    return values, propagate(jacobian, profile["log_parameter_covariance"])


def fit_ellipsoid(source, parameters):
    # The independent reconstruction owns center/attitude. We fit all three
    # positive squared-axis coefficients in that explicit body frame and
    # propagate their correlated uncertainty to physical semiaxes.
    rotation = parameters["body_to_frame_rotation"]
    validate_rotation(rotation)
    center = parameters["center_m"]
    finite(center)
    require(len(center) == 3, "ellipsoid center dimension")
    axes_scale = parameters["normalization_m"]
    require(axes_scale > 0, "ellipsoid normalization must be positive")
    body = [matvec(transpose(rotation), [(x-c)/axes_scale for x, c in zip(point, center)]) for point in source["surface_points_m"]]
    fit = linear_fit([[x*x for x in point] for point in body], [1.0]*len(body), source["implicit_equation_covariance"])
    require(all(x > 0 for x in fit["parameters"]), "reconstruction is not an ellipsoid")
    axes = [axes_scale/math.sqrt(x) for x in fit["parameters"]]
    jacobian = [[(-0.5*axes[i]/fit["parameters"][i] if i == j else 0.0) for j in range(3)] for i in range(3)]
    return {"schema": "sep-reconstructed-ellipsoid-fit-v1", "center_m": center,
            "body_to_frame_rotation": rotation, "axes_m": axes, "axis_covariance_m2": propagate(jacobian, fit["covariance"]),
            "center_attitude_authority": parameters["reconstruction_asset_sha256"],
            "uncertainty_conditioning": "fixed-declared-center-and-attitude", "fit_chi2": fit["chi2"]}


def fit_plasma_sheet(source, parameters):
    # Gaussian enhancement above an independently declared background.
    distance, measured = source["distance_m"], source["density_m3"]
    background = parameters["background_density_m3"]
    enhancement = [n-background for n in measured]
    require(all(n > 0 for n in enhancement), "plasma-sheet enhancement must be positive")
    c = covariance(source["covariance_m6"], len(measured))
    log_cov = [[c[i][j]/(enhancement[i]*enhancement[j]) for j in range(len(measured))] for i in range(len(measured))]
    scale = parameters["distance_scale_m"]
    require(scale > 0, "plasma-sheet distance scale required")
    fit = linear_fit([[1.0, (d/scale)**2] for d in distance], [math.log(n) for n in enhancement], log_cov)
    slope = fit["parameters"][1]
    require(slope < 0, "plasma sheet does not have a finite Gaussian width")
    width = scale/math.sqrt(-2*slope)
    amplitude = math.exp(fit["parameters"][0])
    return {"schema": "sep-plasma-sheet-profile-v1", "background_density_m3": background,
            "enhancement_density_m3": amplitude, "width_m": width,
            "parameter_covariance": propagate([[amplitude, 0], [0, -width/(2*slope)]], fit["covariance"]),
            "support_m": [min(distance), max(distance)]}


def shock_compression(source, parameters):
    n1, n2 = source["upstream_density_m3"], source["downstream_density_m3"]
    finite([n1, n2]); require(n1 > 0 and n2 > n1, "compression observation is not compressive")
    c = covariance(source["joint_covariance_m6"], 2)
    ratio = n2/n1
    return {"schema": "sep-compression-observation-v1", "compression": ratio,
            "variance": propagate([[-n2/n1**2, 1/n1]], c)[0][0], "definition": parameters["definition"]}


def fold_response(source, parameters):
    response = parameters["response_matrix"]
    require(parameters["response_sha256"] == digest(response) and parameters["input_definition"] and parameters["output_definition"], "missing/mismatched instrument-response authority")
    finite(v for row in response for v in row)
    require(all(v >= 0 for row in response for v in row), "negative response weight")
    covariance(source["covariance"], len(source["values"]))
    return {"schema": "sep-instrument-folded-product-v1", "values": matvec(response, source["values"]),
            "covariance": propagate(response, source["covariance"]), "response": parameters}


def validate_table(source, parameters):
    # Closure/critical-Mach/front/comparison tables remain source owned. This
    # operator validates declared columns, SI units, support and covariance;
    # it does not invent an inference rule or silently interpolate/extrapolate.
    required = parameters["column_units"]
    require(set(source["columns"]) == set(required), "table columns differ from registered definitions")
    lengths = {len(values) for values in source["columns"].values()}
    require(len(lengths) == 1 and next(iter(lengths)) > 0, "table dimensions/support invalid")
    finite(v for values in source["columns"].values() for v in values)
    require(all(unit not in {"", "unknown"} for unit in required.values()), "ambiguous table units")
    if "time_s" in source["columns"]:
        times = source["columns"]["time_s"]
        require(all(a < b for a, b in zip(times, times[1:])), "nonmonotone table cadence")
    for key in parameters.get("positive_columns", []):
        require(key in required and all(v > 0 for v in source["columns"][key]), "nonpositive physical table column")
    require("variable_definitions" in parameters and "interpolation" in parameters, "table needs definitions and interpolation authority")
    uncertainty = parameters["uncertainty_model"]
    if uncertainty == "covariance":
        columns = parameters["uncertainty_columns"]
        require(columns and set(columns) <= set(required) and len(columns) == len(set(columns)), "uncertainty columns are missing/duplicated")
        covariance(source["covariance"], next(iter(lengths))*len(columns))
        evidence = {"covariance": source["covariance"], "column_order": columns}
    else:
        require(uncertainty == "named-ensemble" and source["ensemble_members"], "missing uncertainty/ensemble model")
        members = source["ensemble_members"]
        require(len({m["id"] for m in members}) == len(members), "duplicate uncertainty ensemble member")
        for member in members:
            require(set(member["columns"]) == set(required) and all(len(v) == next(iter(lengths)) for v in member["columns"].values()), "ensemble support mismatch")
            finite(v for values in member["columns"].values() for v in values)
        evidence = {"ensemble_members": members}
    return {"schema": parameters["output_schema"], "columns": source["columns"], "declarations": parameters, "uncertainty": evidence}

def fit_kinematics(source, parameters):
    # Simultaneous GLS retains correlations between center and axes. The
    # independently reconstructed attitude is conditioned on, not inferred
    # from a single viewpoint or silently merged with this kinematic fit.
    degree, epoch, scale = parameters["degree"], parameters["reference_time_s"], parameters["time_scale_s"]
    require(degree in {1, 2} and scale > 0, "kinematic polynomial needs a linear/quadratic preregistered degree and time scale")
    times = source["time_s"]
    finite(times)
    require(times and all(a < b for a, b in zip(times, times[1:])), "kinematic cadence/order invalid")
    require(len(source["center_m"]) == len(source["axes_m"]) == len(times), "kinematic observation dimensions")
    values, design = [], []
    for time, center, axes in zip(times, source["center_m"], source["axes_m"]):
        require(len(center) == len(axes) == 3 and all(a > 0 for a in axes), "front center/axes are incomplete or nonphysical")
        finite(center+axes)
        powers = [((time-epoch)/scale)**k for k in range(degree+1)]
        for variable, value in enumerate(center+axes):
            row = [0.0]*(6*(degree+1))
            row[variable*(degree+1):(variable+1)*(degree+1)] = powers
            design.append(row); values.append(value)
    fit = linear_fit(design, values, source["joint_center_axes_covariance_m2"])
    coefficients = [fit["parameters"][i*(degree+1):(i+1)*(degree+1)] for i in range(6)]
    for polynomial in coefficients[3:]:
        support = [(times[0]-epoch)/scale, (times[-1]-epoch)/scale]
        if degree == 2 and polynomial[2] != 0:
            root = -polynomial[1]/(2*polynomial[2])
            if support[0] < root < support[1]: support.append(root)
        require(all(sum(c*x**k for k, c in enumerate(polynomial)) > 0 for x in support), "fitted axis becomes nonpositive inside time coverage")
    return {"schema": "sep-front-kinematic-inference-v1", "reference_time_s": epoch, "time_scale_s": scale,
            "degree": degree, "center_axes_coefficients_m": coefficients, "joint_coefficient_covariance": fit["covariance"],
            "attitude_authority": parameters["attitude_asset_sha256"], "support_time_s": [times[0], times[-1]],
            "uncertainty_conditioning": "fixed-declared-attitude", "fit_chi2": fit["chi2"]}


def fit_wsa_speeds(source, parameters):
    # The shape/hyperparameters are preregistered ensemble authorities. Only
    # slow speed and the positive fast-minus-slow amplitude are inferred here.
    # Retaining their joint covariance avoids falsely independent speed bounds.
    expansion, distance, speed = source["expansion_factor"], source["footpoint_distance_rad"], source["speed_m_per_s"]
    finite(expansion+distance+speed)
    require(len(expansion) == len(distance) == len(speed) and all(f > 0 for f in expansion), "WSA input support invalid")
    alpha, beta, width, power, exponent = [parameters[k] for k in ("alpha", "beta", "width_rad", "distance_power", "outer_exponent")]
    require(alpha > 0 and 0 <= beta < 1 and width > 0 and power > 0 and exponent > 0 and all(d >= 0 for d in distance), "WSA-like shape controls invalid")
    shape = [(1+f)**(-alpha)*(1-beta*math.exp(-(d/width)**power))**exponent for f, d in zip(expansion, distance)]
    fit = linear_fit([[1, value] for value in shape], speed, source["speed_covariance"])
    slow, amplitude = fit["parameters"]
    require(slow > 0 and amplitude > 0, "WSA inferred target speeds are not positive/ordered")
    return {"schema": "sep-wsa-target-speed-inference-v1", "slow_speed_m_per_s": slow,
            "fast_speed_m_per_s": slow+amplitude,
            "speed_covariance": propagate([[1, 0], [1, 1]], fit["covariance"]),
            "shape_authority": parameters, "fit_chi2": fit["chi2"]}


def critical_mach_table(source, parameters):
    beta, angle, values = source["beta"], source["obliquity_rad"], source["values"]
    finite(beta+angle+[v for row in values for v in row])
    require(beta and angle and all(x >= 0 for x in beta) and all(0 <= x <= math.pi/2 for x in angle), "critical-Mach beta/obliquity domain invalid")
    require(all(a < b for a, b in zip(beta, beta[1:])) and all(a < b for a, b in zip(angle, angle[1:])), "critical-Mach axes must be ordered")
    require(len(values) == len(beta) and all(len(row) == len(angle) for row in values) and all(v > 0 for row in values for v in row), "critical-Mach table values/shape invalid")
    require(parameters["mach_convention"] in {"fast", "total-alfven", "normal-alfven"} and parameters["gamma_ad"] > 1 and parameters["interpolation"] == "bounded-bilinear", "critical-Mach convention/EOS/interpolation missing")
    covariance(source["value_covariance"], len(beta)*len(angle))
    return dict(source, schema="sep-critical-mach-table-v1", declarations=parameters)


def formation_constraint(source, parameters):
    times = source["time_s"]
    finite(times)
    require(times and all(a < b for a, b in zip(times, times[1:])), "formation constraint time support invalid")
    if parameters["kind"] == "radio-frequency-time":
        require(len(source["frequency_hz"]) == len(times), "radio time dimensions differ")
        covariance(source["frequency_covariance_hz2"], len(times))
        require(all(f > 0 for f in source["frequency_hz"]), "nonpositive radio frequency")
        require("harmonic_hypotheses" in source, "radio harmonic ambiguity must be explicit")
    else:
        require(parameters["kind"] == "joint-preinferred-height" and source["independence"] == "joint-density-inference" and
                len(source["height_m"]) == len(times) and source["density_model_sha256"], "invalid/independent pre-inferred height")
        covariance(source["joint_height_log_density_covariance"], len(times)+1)
    return dict(source, schema="sep-formation-constraint-v1", declarations=parameters)


OPERATORS = {
    "magnetogram-harmonics": (fit_magnetogram, "T"),
    "image-ellipse": (fit_image_ellipse, "pixel-or-m-explicit-scale"),
    "ephemeris-transform": (transform_ephemeris, "SI-position-velocity"),
    "power-law-profile": (fit_power_law, None),
    "reconstructed-ellipsoid": (fit_ellipsoid, "m"),
    "plasma-sheet-width": (fit_plasma_sheet, "m^-3"),
    "shock-compression": (shock_compression, "m^-3"),
    "instrument-response": (fold_response, None),
    "declared-table": (validate_table, "explicit-column-SI"),
    "wsa-speed-fit": (fit_wsa_speeds, "m/s"),
    "critical-mach-table": (critical_mach_table, "dimensionless"),
    "formation-constraint": (formation_constraint, "explicit-radio-or-height-SI"),
    "front-kinematic-fit": (fit_kinematics, "m"),
}

def preprocess_job(job, base, destination):
    require(job["schema"] == "sep-observation-preprocess-job-v1", "unsupported preprocessing job")
    assets, sources = [], {}
    for request in job["assets"]:
        m = dict(request["metadata"])
        validate_metadata(m)
        require(m["processing_version"] == job["processing_version"] and m["epoch_utc"] == job["epoch_utc"], "mixed processing versions/epochs")
        source_path = Path(base)/request["source_file"]
        source_bytes = source_path.read_bytes()
        require(sha256_bytes(source_bytes) == m["source_sha256"], "source checksum does not match supplied bytes")
        source = load_json(source_path)
        require(request["operator"] in OPERATORS, "unknown inference operator")
        function, unit = OPERATORS[request["operator"]]
        if unit is not None:
            require(m["units"] == unit, "ambiguous/incompatible input units")
        else:
            require(m["units"] in {"m^-3", "K", "m/s", "Pa", "differential-intensity-SI", "counts/s"}, "unsupported scalar-product units")
        if request["operator"] == "ephemeris-transform":
            require(request["parameters"]["transform_epoch_utc"] == job["epoch_utc"], "frame transformation belongs to another epoch")
            require(m["frame"] == request["parameters"]["source_frame"], "ephemeris source frame mismatch")
            require(request["parameters"]["target_frame"] == job["output_frame"], "ephemeris target frame mismatch")
        else:
            require(m["frame"] == job["output_frame"], "asset frame differs from declared job frame")
        product = function(source, request["parameters"])
        assets.append({"metadata": m, "operator": request["operator"],
                       "parameters": request["parameters"], "product": product})
        sources[m["source_sha256"]] = source_bytes
    return publish_bundle(destination, assets, sources, job["processing_version"])

# Path is intentionally imported here: the job API resolves files relative to
# an explicit offline job location, not a process-dependent application cwd.
from pathlib import Path
