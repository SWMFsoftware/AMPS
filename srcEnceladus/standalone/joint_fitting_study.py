#!/usr/bin/env python3
"""Reproducible E3/E5 plus UVIS plume fitting and cross-prediction study."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import lsq_linear, nnls

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from fitting import (SOURCE_CASES, load_config, matrix_from_config,
                     solve_amplitudes, sources_from_config)  # noqa: E402
from plume_fit import load_density  # noqa: E402
from plume_model import (line_of_sight_design_matrix,
                         outward_ray_design_matrix)  # noqa: E402


def sha256(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            value.update(block)
    return value.hexdigest()


def json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def ray_directions(tangent, time):
    """Reproduce the SwRI occultation-plane construction.

    For a tangent track r(t), the line of sight is normal to both r and dr/dt.
    Its sign is irrelevant to an infinite chord integral; signs are made
    continuous only to make the exported geometry easier to inspect.
    """
    tangent = np.asarray(tangent, dtype=float)
    time = np.asarray(time, dtype=float)
    velocity = np.gradient(tangent, time, axis=0, edge_order=2)
    direction = np.cross(tangent, velocity)
    norm = np.linalg.norm(direction, axis=1)
    if np.any(norm == 0) or not np.isfinite(norm).all():
        raise ValueError("Degenerate UVIS tangent track")
    direction /= norm[:, None]
    for index in range(1, len(direction)):
        if np.dot(direction[index-1], direction[index]) < 0:
            direction[index] *= -1
    orthogonality = np.abs(np.sum(direction*tangent, axis=1))
    orthogonality /= np.linalg.norm(tangent, axis=1)
    if orthogonality.max() > 1.0e-10:
        raise ArithmeticError("UVIS line of sight is not tangent")
    return direction, velocity, float(orthogonality.max())


def bounded_weighted_fit(matrix, values, sigma, signed_last=False):
    matrix = np.asarray(matrix, dtype=float)
    values = np.asarray(values, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    weighted = matrix/sigma[:, None]
    target = values/sigma
    if signed_last:
        # Eliminate the unconstrained reference-count offset analytically.
        # Solving the projected source problem with the Lawson-Hanson active
        # set is robust when there are more, and nearly dependent, jet columns
        # than UVIS samples (the 98-jet/2007-occultation case).
        source = weighted[:, :-1]
        baseline = weighted[:, -1]
        baseline_norm2 = float(baseline@baseline)
        if not np.isfinite(baseline_norm2) or baseline_norm2 <= 0:
            raise ValueError("Signed nuisance column must be finite and nonzero")
        projected_source = source-np.outer(
            baseline, (baseline@source)/baseline_norm2)
        projected_target = target-baseline*(baseline@target)/baseline_norm2
        source_scale = np.linalg.norm(projected_source, axis=0)
        source_observable = source_scale > 0
        safe_source_scale = np.where(source_observable, source_scale, 1.0)
        scaled_source = (projected_source[:, source_observable]
                         / safe_source_scale[source_observable])
        source_coefficient = np.zeros(source.shape[1])
        if scaled_source.shape[1]:
            solution, _ = nnls(
                scaled_source, projected_target, maxiter=10000)
            source_coefficient[source_observable] = (
                solution/safe_source_scale[source_observable])
        baseline_coefficient = float(
            baseline@(target-source@source_coefficient)/baseline_norm2)
        coefficient = np.concatenate(
            (source_coefficient, [baseline_coefficient]))
        diagnostic = np.column_stack(
            (scaled_source, baseline/np.sqrt(baseline_norm2)))
        observable = np.concatenate((source_observable, [True]))
        solver = "projected_nnls_with_analytic_signed_nuisance"
    else:
        scale = np.linalg.norm(weighted, axis=0)
        observable = scale > 0
        safe = np.where(observable, scale, 1.0)
        diagnostic = weighted[:, observable]/safe[observable]
        fit = lsq_linear(diagnostic, target, bounds=(0, np.inf),
                         tol=1.0e-12, max_iter=10000)
        if not fit.success:
            raise RuntimeError(fit.message)
        coefficient = np.zeros(matrix.shape[1])
        coefficient[observable] = fit.x/safe[observable]
        solver = "scipy_lsq_linear_trf"
    prediction = matrix@coefficient
    singular = np.linalg.svd(diagnostic, compute_uv=False)
    relative = singular/singular[0] if len(singular) else singular
    return {
        "solver": solver,
        "coefficient": coefficient,
        "prediction": prediction,
        "residual": (prediction-values)/sigma,
        "chi2": float(np.sum(((prediction-values)/sigma)**2)),
        "rank": int(np.linalg.matrix_rank(diagnostic)),
        "condition_number_scaled": (
            float(singular[0]/singular[-1]) if len(singular) and singular[-1] > 0
            else float("inf")),
        "effective_rank": {
            "1e-2": int(np.count_nonzero(relative >= 1.0e-2)),
            "1e-3": int(np.count_nonzero(relative >= 1.0e-3)),
            "1e-4": int(np.count_nonzero(relative >= 1.0e-4)),
        },
    }


def load_uvis(path):
    rows = read_csv(path)
    tangent = np.asarray([[float(row[name]) for name in
                           ("tangent_x_m", "tangent_y_m", "tangent_z_m")]
                          for row in rows])
    time = np.asarray([float(row["time_source_s"]) for row in rows])
    direction, velocity, orthogonality = ray_directions(tangent, time)
    return {
        "rows": rows,
        "time": np.asarray([float(row["time_from_min_altitude_s"]) for row in rows]),
        "altitude_m": np.asarray([float(row["tangent_altitude_m"]) for row in rows]),
        "tangent": tangent,
        "direction": direction,
        "track_velocity": velocity,
        "tau": np.asarray([float(row["optical_depth_raw"]) for row in rows]),
        "sigma_tau": np.asarray([float(row["poisson_sigma_optical_depth"])
                                 for row in rows]),
        "orthogonality_max": orthogonality,
    }


def occultation_matrix(sources, data, reference_radius, exponent, order,
                       observer_geometry=None):
    if observer_geometry is None:
        matrix = line_of_sight_design_matrix(
            sources, data["tangent"], data["direction"], reference_radius,
            background="power_law", power_law_exponent=exponent,
            order=order, batch_size=64)
        return matrix, {
            "geometry": "infinite chord from SwRI tangent-track normal",
            "infinite_chord_rows": len(data["tangent"]),
            "finite_tangent_rows": 0,
            "outward_observer_rows": 0,
            "ray_tangent_orthogonality_max": data["orthogonality_max"],
        }

    coordinate = observer_geometry["coordinate"]
    sun_direction = observer_geometry["direction"]
    if len(coordinate) != len(data["tangent"]):
        raise ValueError("Observer geometry row mismatch")
    inward = coordinate < 0
    outward = ~inward
    matrix = np.empty((len(coordinate), len(sources)+1))
    orthogonality = 0.0
    if inward.any():
        tangent = data["tangent"][inward]
        direction = sun_direction[inward].copy()
        direction -= (np.sum(direction*tangent, axis=1)
                      / np.sum(tangent*tangent, axis=1))[:, None]*tangent
        direction /= np.linalg.norm(direction, axis=1)[:, None]
        orthogonality = float(np.max(np.abs(np.sum(direction*tangent, axis=1))
                                     / np.linalg.norm(tangent, axis=1)))
        matrix[inward] = line_of_sight_design_matrix(
            sources, tangent, direction, reference_radius,
            background="power_law", power_law_exponent=exponent,
            order=order, batch_size=64,
            observer_coordinate_from_tangent_m=coordinate[inward])
    if outward.any():
        origin = data["tangent"][outward]
        direction = sun_direction[outward].copy()
        radial_dot = np.sum(origin*direction, axis=1)
        direction[radial_dot < 0] *= -1
        matrix[outward] = outward_ray_design_matrix(
            sources, origin, direction, reference_radius,
            background="power_law", power_law_exponent=exponent,
            order=order, batch_size=64)
    return matrix, {
        "geometry": "finite Cassini-to-Sun ray from SPICE",
        "infinite_chord_rows": 0,
        "finite_tangent_rows": int(inward.sum()),
        "outward_observer_rows": int(outward.sum()),
        "ray_tangent_orthogonality_max": orthogonality,
    }


def fit_uvis(column_matrix, tau, sigma_tau, published_column,
             published_sigma, labels):
    # The broadband response is not calibrated. Fit only its shape and a
    # signed reference-count offset, then attach the published peak-column
    # scale once. This avoids pretending each tau sample is an independent
    # absolute column measurement.
    augmented = np.column_stack((column_matrix, np.ones(len(tau))))
    shape = bounded_weighted_fit(augmented, tau, sigma_tau, signed_last=True)
    shape_coeff = shape["coefficient"][:-1]
    signal = column_matrix@shape_coeff
    maximum = float(signal.max())
    if maximum <= 0:
        raise RuntimeError("UVIS fit has no positive plume signal")
    physical_scale = published_column/maximum
    physical_coeff = shape_coeff*physical_scale
    physical_column = column_matrix@physical_coeff
    shape_dof = len(tau)-shape["rank"]
    return {
        "solver": shape["solver"],
        "labels": labels,
        "shape_coefficients": shape_coeff,
        "physical_coefficients": physical_coeff,
        "physical_column_m2": physical_column,
        "tau_prediction": shape["prediction"],
        "tau_plume_prediction": signal,
        "tau_baseline": float(shape["coefficient"][-1]),
        "chi2_shape": shape["chi2"],
        "shape_dof": shape_dof,
        "shape_reduced_chi2": (
            shape["chi2"]/shape_dof if shape_dof > 0 else None),
        "rank": shape["rank"],
        "effective_rank": shape["effective_rank"],
        "condition_number_scaled": shape["condition_number_scaled"],
        "published_peak_column_m2": published_column,
        "published_peak_column_sigma_m2": published_sigma,
        "absolute_scale_fractional_sigma": published_sigma/published_column,
        "normalization_method": (
            "shape fit in optical-depth space; common coefficient scale chosen "
            "so the modeled sampled peak equals the published column"),
    }


def full_vector(labels, values, full_labels):
    table = dict(zip(labels, values))
    return np.asarray([table.get(label, 0.0) for label in full_labels])


def uvis_shape_metric(matrix, physical_vector, tau, sigma, published_column):
    column = matrix@physical_vector
    augmented = np.column_stack((column, np.ones(len(column))))
    fit = bounded_weighted_fit(augmented, tau, sigma, signed_last=True)
    peak = float(column.max())
    return {
        "chi2_per_sample": fit["chi2"]/len(tau),
        "chi2": fit["chi2"],
        "best_tau_per_column_m2": float(fit["coefficient"][0]),
        "tau_baseline": float(fit["coefficient"][1]),
        "predicted_peak_column_m2": peak,
        "peak_column_ratio_to_published": peak/published_column,
    }


def density_metric(matrix, vector, density, sigma):
    prediction = matrix@vector
    residual = (prediction-density)/sigma
    positive = (prediction > 0) & (density > 0)
    return {
        "chi2": float(residual@residual),
        "chi2_per_sample": float(residual@residual/len(density)),
        "median_log10_prediction_over_observation": (
            float(np.median(np.log10(prediction[positive]/density[positive])))
            if positive.any() else None),
    }


def solve_combined_blocks(blocks, balance_by_dataset=True):
    """Fit one nonnegative vector to several observation blocks.

    Each block is ``(name, matrix, value, sigma)``.  Dataset balancing
    multiplies every point sigma by sqrt(N_dataset), so the minimized objective
    is the sum of the mean weighted-square residual from each dataset.  This
    prevents the 1,133-row solar occultation from overwhelming the 97-row
    stellar occultation solely because it has more samples.
    """
    if not blocks:
        raise ValueError("At least one dataset block is required")
    columns = np.asarray(blocks[0][1]).shape[1]
    matrices, values, sigmas, prepared = [], [], [], []
    weighting = {}
    samples = {}
    for name, matrix, value, sigma in blocks:
        matrix = np.asarray(matrix, dtype=float)
        value = np.asarray(value, dtype=float)
        sigma = np.asarray(sigma, dtype=float)
        if (matrix.ndim != 2 or matrix.shape[1] != columns
                or value.shape != (len(matrix),)
                or sigma.shape != value.shape or len(value) == 0):
            raise ValueError("Invalid combined-fit block: "+str(name))
        multiplier = np.sqrt(len(value)) if balance_by_dataset else 1.0
        matrices.append(matrix)
        values.append(value)
        sigmas.append(sigma*multiplier)
        prepared.append((name, matrix, value, sigma))
        weighting[name] = {
            "samples": len(value),
            "sigma_multiplier": float(multiplier),
            "objective_contribution": (
                "mean weighted-square residual"
                if balance_by_dataset else "sum weighted-square residual"),
        }
        samples[name] = len(value)
    fit = solve_amplitudes(
        np.vstack(matrices), np.concatenate(values), np.concatenate(sigmas))
    fit["datasets"] = list(samples)
    fit["dataset_weighting"] = weighting
    fit["balance_by_dataset"] = bool(balance_by_dataset)
    fit["per_dataset"] = {}
    for name, matrix, value, sigma in prepared:
        residual = (matrix@fit["coefficients"]-value)/sigma
        chi2 = float(residual@residual)
        fit["per_dataset"][name] = {
            "samples": len(value),
            "chi2": chi2,
            "chi2_per_sample": chi2/len(value),
            "balanced_objective_contribution": (
                chi2/len(value) if balance_by_dataset else chi2),
        }
    contribution = sum(
        item["balanced_objective_contribution"]
        for item in fit["per_dataset"].values())
    if not np.isclose(contribution, fit["chi2"], rtol=2.0e-12, atol=1.0e-12):
        raise ArithmeticError("Combined-fit objective bookkeeping mismatch")
    return fit


def normalized_source_shape(vector, labels):
    """Return a full-length, source-only vector with unit source sum.

    A source coefficient is expressed in units of 1e28 H2O/s.  Removing the
    fitted background terms and normalizing the remaining coefficients makes
    the vector a spatial source shape.  Multiplying the returned vector by a
    single coefficient therefore gives an independently fitted total source
    rate for another observation.
    """
    vector = np.asarray(vector, dtype=float)
    if vector.shape != (len(labels),):
        raise ValueError("Source-shape vector and labels do not match")
    source = np.asarray([
        label.startswith(("jet:", "stripe:")) for label in labels])
    if not source.any():
        raise ValueError("No physical source components in labels")
    if (not np.isfinite(vector).all()) or np.any(vector[source] < 0):
        raise ValueError("Source-shape coefficients must be finite and nonnegative")
    total = float(vector[source].sum())
    if total <= 0:
        raise ValueError("Source-shape vector has zero total source rate")
    shape = np.zeros(len(labels))
    shape[source] = vector[source]/total
    return shape


def fit_rate_for_fixed_source_shape(matrix, values, sigma, labels, shape_vector):
    """Fit one total source rate while retaining free background amplitudes.

    The relative jet/stripe amplitudes are frozen to ``shape_vector``.  The
    first design column is their unit-sum response, so its fitted coefficient
    is the total physical source rate in units of 1e28 H2O/s.  Backgrounds
    remain observation-specific nuisance parameters and are never folded into
    the source-rate estimate.
    """
    matrix = np.asarray(matrix, dtype=float)
    values = np.asarray(values, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    if (matrix.ndim != 2 or matrix.shape[1] != len(labels)
            or values.shape != (len(matrix),) or sigma.shape != values.shape):
        raise ValueError("Invalid fixed-shape rate-fit inputs")
    shape = normalized_source_shape(shape_vector, labels)
    background_index = [
        index for index, label in enumerate(labels)
        if label.startswith("background:")
        and np.any(np.asarray(matrix[:, index]) != 0.0)
    ]
    collapsed = np.column_stack(
        [matrix@shape]+[matrix[:, index] for index in background_index])
    fit = solve_amplitudes(collapsed, values, sigma)
    coefficient = np.asarray(fit["coefficients"])
    full = shape*coefficient[0]
    for output_index, label_index in enumerate(background_index, start=1):
        full[label_index] = coefficient[output_index]
    return {
        "total_source_rate_s": float(coefficient[0]*1.0e28),
        "total_source_coefficient_1e28_s-1": float(coefficient[0]),
        "source_shape_fractions": {
            label: float(shape[index]) for index, label in enumerate(labels)
            if label.startswith(("jet:", "stripe:"))
        },
        "background_coefficients": {
            labels[label_index]: float(coefficient[output_index])
            for output_index, label_index in enumerate(background_index, start=1)
        },
        "full_coefficients": full,
        "chi2": fit["chi2"],
        "chi2_per_sample": fit["chi2"]/len(values),
        "rank": fit["rank"],
        "samples": len(values),
        "condition_number_scaled": fit["condition_number_scaled"],
    }


def uvis_column_block(matrix, data, fit, full_labels, uvis_labels):
    """Convert one UVIS shape fit to the documented physical-column proxy."""
    plume_tau = np.maximum(data["tau"]-fit["tau_baseline"], 0.0)
    tau_peak = float(plume_tau.max())
    if not np.isfinite(tau_peak) or tau_peak <= 0:
        raise RuntimeError("UVIS reference-subtracted profile has no positive peak")
    column_scale = fit["published_peak_column_m2"]/tau_peak
    pseudo_column = plume_tau*column_scale
    pseudo_sigma = data["sigma_tau"]*column_scale
    full = np.zeros((len(matrix), len(full_labels)))
    for index, label in enumerate(uvis_labels):
        full[:, full_labels.index(label)] = matrix[:, index]
    diagnostic = {
        "tau_baseline_from_independent_shape_fit": fit["tau_baseline"],
        "nonnegative_reference_subtracted_samples": int(
            np.count_nonzero(plume_tau)),
        "tau_to_column_scale_m-2": float(column_scale),
        "published_peak_column_m2": fit["published_peak_column_m2"],
    }
    return full, pseudo_column, pseudo_sigma, diagnostic


def selected_columns(labels, variant):
    if variant == "full":
        return np.ones(len(labels), dtype=bool)
    if variant == "no_stripes":
        return np.asarray([not label.startswith("stripe:") for label in labels])
    if variant == "no_jets":
        return np.asarray([not label.startswith("jet:") for label in labels])
    if variant == "background_only":
        return np.asarray([label.startswith("background:") for label in labels])
    raise ValueError(variant)


def write_csv(path, rows):
    if not rows:
        raise ValueError("Cannot write empty table")
    with Path(path).open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--source-case", choices=SOURCE_CASES)
    parser.add_argument("--jet-catalog", type=Path)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--quadrature-order", type=int, default=96)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error("output directory already exists")
    args.output_dir.mkdir(parents=True)

    config = load_config(args.config, args.source_case, args.jet_catalog)
    sources, full_labels = sources_from_config(config)
    source_labels = [label for label in full_labels
                     if label.startswith(("jet:", "stripe:"))]
    jet_count = sum(label.startswith("jet:") for label in source_labels)
    stripe_count = sum(label.startswith("stripe:") for label in source_labels)
    uvis_labels = source_labels+["background:density_at_reference_radius_m-3"]
    constraints = {row["event_id"]: row for row in read_csv(
        args.data_root/"uvis_published_column_constraints.csv")}

    density_data = {}
    density_fits = {}
    density_matrices = {}
    for flyby in ("e3", "e5"):
        observation = args.data_root/f"{flyby}_h2o_sticking_corrected.csv"
        metadata = args.data_root/f"{flyby}_h2o_sticking_corrected.metadata.json"
        x, y, sigma, meta = load_density(observation, metadata)
        matrix, labels = matrix_from_config(config, x)
        if labels != full_labels:
            raise AssertionError("density label mismatch")
        fit = solve_amplitudes(matrix, y, sigma)
        density_data[flyby] = {"x": x, "y": y, "sigma": sigma,
                               "rows": read_csv(observation), "metadata": meta}
        density_matrices[flyby] = matrix
        density_fits[flyby] = fit

    uvis_files = {
        "uvis_2007_zeta_orionis": "uvis_2007_zeta_orionis_optical_depth.csv",
        "uvis_2010_e10_solar": "uvis_2010_e10_solar_optical_depth.csv",
    }
    uvis_data, uvis_matrices, uvis_fits, convergence = {}, {}, {}, {}
    for event, name in uvis_files.items():
        data = load_uvis(args.data_root/name)
        observer_geometry = None
        if event == "uvis_2010_e10_solar":
            geometry_rows = read_csv(
                args.data_root/"uvis_2010_e10_solar_observer_range.csv")
            geometry_time = np.asarray([float(row["time_source_s"])
                                        for row in geometry_rows])
            source_time = np.asarray([float(row["time_source_s"])
                                      for row in data["rows"]])
            if not np.allclose(geometry_time, source_time, rtol=0, atol=1.0e-5):
                raise ValueError("SPICE observer times do not match UVIS rows")
            observer_geometry = {
                "coordinate": np.asarray([
                    float(row["observer_coordinate_from_tangent_m"])
                    for row in geometry_rows]),
                "direction": np.asarray([[float(row[name]) for name in
                    ("ray_ux_table", "ray_uy_table", "ray_uz_table")]
                    for row in geometry_rows]),
            }
        matrix, geometry_diagnostic = occultation_matrix(
            sources, data, config["reference_radius_m"],
            config["background"]["power_law_exponent"],
            args.quadrature_order, observer_geometry)
        constraint = constraints[event]
        fit = fit_uvis(
            matrix, data["tau"], data["sigma_tau"],
            float(constraint["published_h2o_column_m2"]),
            float(constraint["published_h2o_column_sigma_m2"]), uvis_labels)
        # Check the quadrature on a deterministic 12-ray subset at a higher order.
        index = np.linspace(0, len(data["tangent"])-1, 12).round().astype(int)
        subset = dict(data)
        for key in ("tangent", "direction"):
            subset[key] = data[key][index]
        subset["orthogonality_max"] = data["orthogonality_max"]
        observer_subset = None
        if observer_geometry is not None:
            observer_subset = {key: value[index]
                               for key, value in observer_geometry.items()}
        check, _ = occultation_matrix(
            sources, subset, config["reference_radius_m"],
            config["background"]["power_law_exponent"],
            args.quadrature_order+32, observer_subset)
        base = matrix[index]
        threshold = 1.0e-10*np.max(check, axis=0, keepdims=True)
        valid = check > threshold
        relative = np.abs(base-check)/np.maximum(check, threshold)
        convergence[event] = {
            "comparison_order": args.quadrature_order+32,
            "rays": index,
            "max_relative_difference_material_columns": float(relative[valid].max()),
            "median_relative_difference_material_columns": float(np.median(relative[valid])),
            **geometry_diagnostic,
        }
        base_combined = base@fit["physical_coefficients"]
        check_combined = check@fit["physical_coefficients"]
        material_combined = check_combined > 1.0e-6*check_combined.max()
        combined_relative = np.abs(base_combined-check_combined)/np.maximum(
            check_combined, 1.0e-6*check_combined.max())
        convergence[event].update({
            "max_relative_difference_fitted_total_column": float(
                combined_relative[material_combined].max()),
            "median_relative_difference_fitted_total_column": float(
                np.median(combined_relative[material_combined])),
        })
        uvis_data[event], uvis_matrices[event], uvis_fits[event] = data, matrix, fit

    vectors = {
        "E3_independent": density_fits["e3"]["coefficients"],
        "E5_independent": density_fits["e5"]["coefficients"],
    }
    for event, fit in uvis_fits.items():
        vectors[event+"_independent"] = full_vector(
            fit["labels"], fit["physical_coefficients"], full_labels)

    # Fit one shared source vector within each instrument family.  E3 and E5
    # have equal sample counts, so balancing does not change their relative
    # weighting.  UVIS shapes are converted to physical-column proxies using
    # the independently fitted reference offset and one published peak scale.
    inms_blocks = []
    for flyby in ("e3", "e5"):
        data = density_data[flyby]
        inms_blocks.append((flyby.upper(), density_matrices[flyby],
                            data["y"], data["sigma"]))
    inms_combined = solve_combined_blocks(inms_blocks)
    vectors["INMS_combined"] = inms_combined["coefficients"]

    uvis_blocks = []
    uvis_full_blocks = []
    uvis_column_diagnostics = {}
    for event in uvis_files:
        data, matrix, fit = uvis_data[event], uvis_matrices[event], uvis_fits[event]
        full, pseudo_column, pseudo_sigma, diagnostic = uvis_column_block(
            matrix, data, fit, full_labels, uvis_labels)
        uvis_blocks.append((event, matrix, pseudo_column, pseudo_sigma))
        uvis_full_blocks.append((event, full, pseudo_column, pseudo_sigma))
        uvis_column_diagnostics[event] = diagnostic
    uvis_combined = solve_combined_blocks(uvis_blocks)
    uvis_combined["component_labels"] = uvis_labels
    vectors["UVIS_combined"] = full_vector(
        uvis_labels, uvis_combined["coefficients"], full_labels)

    # Existing all-observation balanced diagnostic, now built from the same
    # checked blocks as the two instrument-family combined fits.
    steady = solve_combined_blocks(inms_blocks+uvis_full_blocks)
    vectors["steady_balanced_diagnostic"] = steady["coefficients"]

    # Separate spatial-shape inference from epoch-to-epoch normalization.
    # This is the quantity needed for an orbital/tidal source-rate comparison:
    # a shared normalized shape is learned from a declared training set, then
    # exactly one total source-rate coefficient is refitted to each observation.
    target_blocks = {
        "E3": (density_matrices["e3"], density_data["e3"]["y"],
               density_data["e3"]["sigma"]),
        "E5": (density_matrices["e5"], density_data["e5"]["y"],
               density_data["e5"]["sigma"]),
    }
    for event, full, pseudo_column, pseudo_sigma in uvis_full_blocks:
        target_blocks[event] = (full, pseudo_column, pseudo_sigma)

    independent_keys = {
        "E3": "E3_independent",
        "E5": "E5_independent",
        "uvis_2007_zeta_orionis": "uvis_2007_zeta_orionis_independent",
        "uvis_2010_e10_solar": "uvis_2010_e10_solar_independent",
    }
    rate_strategies = {
        "independent_each_observation": {
            "description": (
                "Each observation independently determines both its relative "
                "source shape and total source rate."),
            "shape_training_datasets": "target observation only",
            "target_fits": {},
        }
    }
    for target, vector_key in independent_keys.items():
        vector = np.asarray(vectors[vector_key])
        shape = normalized_source_shape(vector, full_labels)
        matrix, values, sigma = target_blocks[target]
        residual = (matrix@vector-values)/sigma
        rate_strategies["independent_each_observation"]["target_fits"][target] = {
            "total_source_rate_s": float(
                sum(value for label, value in zip(full_labels, vector)
                    if label.startswith(("jet:", "stripe:")))*1.0e28),
            "source_shape_fractions": {
                label: float(shape[index])
                for index, label in enumerate(full_labels)
                if label.startswith(("jet:", "stripe:"))
            },
            "chi2": float(residual@residual),
            "chi2_per_sample": float(residual@residual/len(values)),
            "samples": len(values),
        }

    shared_shape_definitions = {
        "shared_shape_all_observations": (
            "E3, E5, 2007 UVIS, and 2010 UVIS",
            vectors["steady_balanced_diagnostic"],
            "All four observations determine one relative source shape; one "
            "total source rate is then fitted independently to each observation."),
        "shared_shape_inms_only": (
            "E3 and E5 INMS",
            vectors["INMS_combined"],
            "E3 and E5 determine one relative source shape; one total source "
            "rate is then fitted independently to every INMS and UVIS observation."),
        "shared_shape_uvis_only": (
            "2007 and 2010 UVIS",
            vectors["UVIS_combined"],
            "The two UVIS occultations determine one relative source shape; "
            "one total source rate is then fitted independently to every INMS "
            "and UVIS observation."),
    }
    for strategy, (training, shape_vector, description) in shared_shape_definitions.items():
        rate_strategies[strategy] = {
            "description": description,
            "shape_training_datasets": training,
            "target_fits": {},
        }
        for target, (matrix, values, sigma) in target_blocks.items():
            rate_strategies[strategy]["target_fits"][target] = (
                fit_rate_for_fixed_source_shape(
                    matrix, values, sigma, full_labels, shape_vector))

    rate_strategy_rows = []
    for strategy, strategy_result in rate_strategies.items():
        for target, fit in strategy_result["target_fits"].items():
            rate_strategy_rows.append({
                "strategy": strategy,
                "shape_training_datasets": strategy_result["shape_training_datasets"],
                "target": target,
                "total_source_rate_s": fit["total_source_rate_s"],
                "chi2": fit["chi2"],
                "chi2_per_sample": fit["chi2_per_sample"],
                "samples": fit["samples"],
            })

    cross_rows = []
    cross = {}
    for model_name, vector in vectors.items():
        cross[model_name] = {}
        for flyby in ("e3", "e5"):
            metric = density_metric(density_matrices[flyby], vector,
                                    density_data[flyby]["y"],
                                    density_data[flyby]["sigma"])
            target = flyby.upper()
            cross[model_name][target] = metric
            cross_rows.append({"model": model_name, "target": target,
                               "metric": "density_chi2_per_sample",
                               "value": metric["chi2_per_sample"]})
        uv_vector = np.asarray([vector[full_labels.index(label)]
                                for label in uvis_labels])
        for event in uvis_files:
            metric = uvis_shape_metric(
                uvis_matrices[event], uv_vector, uvis_data[event]["tau"],
                uvis_data[event]["sigma_tau"],
                uvis_fits[event]["published_peak_column_m2"])
            cross[model_name][event] = metric
            cross_rows.append({"model": model_name, "target": event,
                               "metric": "uvis_shape_chi2_per_sample",
                               "value": metric["chi2_per_sample"]})
            cross_rows.append({"model": model_name, "target": event,
                               "metric": "peak_column_ratio_to_published",
                               "value": metric["peak_column_ratio_to_published"]})

    comparison_rows = []
    for variant in ("full", "no_stripes", "no_jets", "background_only"):
        keep = selected_columns(full_labels, variant)
        for flyby in ("e3", "e5"):
            fit = solve_amplitudes(density_matrices[flyby][:, keep],
                                   density_data[flyby]["y"],
                                   density_data[flyby]["sigma"])
            comparison_rows.append({
                "dataset": flyby.upper(), "variant": variant,
                "free_nonnegative_components": int(keep.sum()),
                "chi2": fit["chi2"], "dof": len(density_data[flyby]["y"])-fit["rank"],
                "reduced_chi2": fit["chi2"]/max(len(density_data[flyby]["y"])-fit["rank"], 1),
            })
        keep_uvis = np.asarray([label in set(np.asarray(full_labels)[keep])
                                for label in uvis_labels])
        for event in uvis_files:
            fit = fit_uvis(
                uvis_matrices[event][:, keep_uvis], uvis_data[event]["tau"],
                uvis_data[event]["sigma_tau"],
                uvis_fits[event]["published_peak_column_m2"],
                uvis_fits[event]["published_peak_column_sigma_m2"],
                list(np.asarray(uvis_labels)[keep_uvis]))
            comparison_rows.append({
                "dataset": event, "variant": variant,
                "free_nonnegative_components": int(keep_uvis.sum()),
                "chi2": fit["chi2_shape"], "dof": fit["shape_dof"],
                "reduced_chi2": (
                    fit["shape_reduced_chi2"]
                    if fit["shape_reduced_chi2"] is not None else ""),
            })

    amplitude_rows = []
    for model_name, vector in vectors.items():
        for label, value in zip(full_labels, vector):
            amplitude_rows.append({
                "model": model_name, "component": label,
                "coefficient": value,
                "source_rate_s": value*1.0e28 if label.startswith(("jet:", "stripe:")) else "",
            })
    write_csv(args.output_dir/"component_amplitudes.csv", amplitude_rows)
    write_csv(args.output_dir/"model_comparison.csv", comparison_rows)
    write_csv(args.output_dir/"cross_prediction_metrics.csv", cross_rows)
    write_csv(args.output_dir/"source_rate_shape_strategies.csv",
              rate_strategy_rows)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    for ax, flyby, color in zip(axes[0], ("e3", "e5"), ("#d73027", "#4575b4")):
        data = density_data[flyby]
        time = np.asarray([float(row["time_ca_s"]) for row in data["rows"]])
        ax.errorbar(time, data["y"], yerr=data["sigma"], fmt=".", ms=3,
                    color="black", elinewidth=.4, label="INMS qualified H2O")
        ax.semilogy(time, density_matrices[flyby]@vectors[flyby.upper()+"_independent"],
                    color=color, lw=2, label="independent fit")
        ax.semilogy(time, density_matrices[flyby]@vectors["INMS_combined"],
                    color="#756bb1", lw=1.5, label="combined INMS fit")
        ax.semilogy(time, density_matrices[flyby]@vectors["steady_balanced_diagnostic"],
                    color="#31a354", lw=1.2, label="all-data combined diagnostic")
        ax.set(title=flyby.upper(), xlabel="Time from flyby reference (s)",
               ylabel="H2O number density (m$^{-3}$)")
        ax.grid(True, which="both", alpha=.25); ax.legend(fontsize=8)
    for ax, event, color in zip(axes[1], uvis_files, ("#984ea3", "#ff7f00")):
        data, fit = uvis_data[event], uvis_fits[event]
        ax.errorbar(data["time"], data["tau"], yerr=data["sigma_tau"],
                    fmt=".", ms=2, color="black", elinewidth=.3,
                    label="UVIS optical-depth proxy")
        ax.plot(data["time"], fit["tau_prediction"], color=color, lw=1.8,
                label="independent shape fit")
        for model_name, model_label, model_color, width in (
                ("UVIS_combined", "combined UVIS fit", "#756bb1", 1.4),
                ("steady_balanced_diagnostic", "all-data combined diagnostic",
                 "#31a354", 1.1)):
            physical = np.asarray([
                vectors[model_name][full_labels.index(label)]
                for label in uvis_labels])
            metric = uvis_shape_metric(
                uvis_matrices[event], physical, data["tau"],
                data["sigma_tau"], fit["published_peak_column_m2"])
            column = uvis_matrices[event]@physical
            ax.plot(data["time"],
                    metric["best_tau_per_column_m2"]*column
                    + metric["tau_baseline"], color=model_color, lw=width,
                    label=model_label)
        ax.set(title=event.replace("uvis_", "UVIS ").replace("_", " "),
               xlabel="Time from minimum tangent altitude (s)",
               ylabel="Optical-depth proxy")
        ax.grid(True, alpha=.25); ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output_dir/"fit_profiles.png", dpi=180)
    fig.savefig(args.output_dir/"fit_profiles.pdf")
    plt.close(fig)

    totals = {}
    for name, vector in vectors.items():
        totals[name] = float(sum(value for label, value in zip(full_labels, vector)
                                 if label.startswith(("jet:", "stripe:")))*1.0e28)
    summary = {
        "study": "E3/E5 local-density and 2007/2010 UVIS occultation fitting",
        "config": str(args.config.resolve()),
        "config_sha256": sha256(args.config),
        "software_sha256": {
            str(path): sha256(path) for path in (
                Path(__file__).resolve(),
                ROOT/"fitting.py",
                ROOT/"plume_fit.py",
                ROOT/"plume_model.py")
        },
        "source_case": config["source_case"],
        "jet_catalog": config.get("_jet_catalog_path"),
        "jet_catalog_sha256": (
            sha256(config["_jet_catalog_path"])
            if config.get("_jet_catalog_path") else None),
        "data_root": str(args.data_root.resolve()),
        "input_sha256": {
            str(path.resolve()): sha256(path)
            for path in sorted(args.data_root.iterdir()) if path.is_file()
        },
        "quadrature_order": args.quadrature_order,
        "component_labels": full_labels,
        "total_source_rate_s": totals,
        "inms_fits": density_fits,
        "uvis_fits": uvis_fits,
        "inms_combined_fit": inms_combined,
        "uvis_combined_fit": uvis_combined,
        "uvis_column_proxy_diagnostics": uvis_column_diagnostics,
        "steady_balanced_diagnostic": steady,
        "source_rate_shape_strategies": rate_strategies,
        "cross_prediction": cross,
        "quadrature_convergence": convergence,
        "interpretation_limits": [
            "E3/E5 are digitized sticking-corrected profiles, not independent raw-count inversions.",
            "UVIS resolved profiles constrain shape in optical-depth space; one published peak/representative column supplies each absolute scale.",
            "The INMS-, UVIS-, and all-data combined fits minimize a sum of per-dataset mean weighted-square residuals; they are balanced diagnostics, not formal joint likelihoods.",
            "Shared-shape rate strategies normalize only jet/stripe coefficients; each observation refits one total source coefficient plus its own nonnegative background nuisance amplitudes.",
            "Individual source rates are not identifiable when the scaled design matrix is ill-conditioned.",
            "The model is steady, collisionless, gravity-free, and omits time-dependent inlet and UVIS spectral-response physics.",
        ],
    }
    with (args.output_dir/"study_results.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, default=json_default, allow_nan=False)
        stream.write("\n")

    report = [
        "# E3/E5 and UVIS fitting study", "",
        "This is a reproducible Stage-0 diagnostic using the steady, collisionless, "
        "gravity-free multi-source model. E3/E5 are fitted in local-density space. "
        "UVIS is fitted in optical-depth-shape space and receives one external "
        "published column normalization per event.", "",
        "## Source case", "",
        f"`{config['source_case']}` resolves to {jet_count} discrete jets, "
        f"{stripe_count} stripe components, and "
        f"{len(full_labels)-len(source_labels)} density-background components. "
        "All source amplitudes are refitted for each independent dataset; catalog "
        "membership does not import fixed production rates.", "",
        "## Independent fits", "",
        "| Dataset-derived model | Total discrete+stripe source rate (s^-1) | Common scale systematic | Fit statistic |",
        "|---|---:|---:|---:|",
    ]
    report.append(
        f"| E3 independent | {totals['E3_independent']:.6e} | 23% | "
        f"reduced chi-square {density_fits['e3']['chi2']/(len(density_data['e3']['y'])-density_fits['e3']['rank']):.3f} |")
    report.append(
        f"| E5 independent | {totals['E5_independent']:.6e} | 23% | "
        f"reduced chi-square {density_fits['e5']['chi2']/(len(density_data['e5']['y'])-density_fits['e5']['rank']):.3f} |")
    for event in uvis_files:
        key = event+"_independent"
        reduced = uvis_fits[event]["shape_reduced_chi2"]
        statistic = (
            f"shape reduced chi-square {reduced:.3f}"
            if reduced is not None else
            f"shape chi-square {uvis_fits[event]['chi2_shape']:.3f}; "
            "zero residual degrees of freedom")
        report.append(
            f"| {event} independent | {totals[key]:.6e} | "
            f"{100*uvis_fits[event]['absolute_scale_fractional_sigma']:.1f}% | "
            f"{statistic} |")
    report += [
        "", "## Combined fits", "",
        "Each combined fit uses one physical source/background vector. Dataset "
        "balancing makes each included observation contribute its mean weighted-"
        "square residual, rather than weight in direct proportion to its number "
        "of samples. E3 and E5 both contain 79 samples, so balancing does not "
        "change their relative weight. UVIS first uses each independent signed "
        "reference offset and published peak column to form the documented "
        "physical-column proxy; the two proxy profiles then share one vector.", "",
        "The table reports the objective in the space actually fitted: local "
        "density for INMS and the documented physical-column proxy for UVIS. "
        "A dash means that dataset was not included in that combined fit.", "",
        "| Combined model | Total source rate (s^-1) | E3 density chi-square/sample | E5 density chi-square/sample | 2007 UVIS column-proxy chi-square/sample | 2010 UVIS column-proxy chi-square/sample |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for model_name, display, combined_fit in (
            ("INMS_combined", "E3+E5 combined", inms_combined),
            ("UVIS_combined", "2007+2010 UVIS combined", uvis_combined),
            ("steady_balanced_diagnostic", "All four datasets combined", steady)):
        objective = combined_fit["per_dataset"]
        cells = []
        for dataset in ("E3", "E5", "uvis_2007_zeta_orionis",
                        "uvis_2010_e10_solar"):
            cells.append(
                f"{objective[dataset]['chi2_per_sample']:.3f}"
                if dataset in objective else "--")
        report.append(
            f"| {display} | {totals[model_name]:.6e} | "
            +" | ".join(cells)+" |")

    report += [
        "", "## Source-shape strategies for orbital-phase comparisons", "",
        "These rates separate spatial-shape inference from temporal "
        "normalization. For each shared-shape strategy, only the normalized "
        "relative jet/stripe amplitudes are transferred. Each observation "
        "then refits one nonnegative total source coefficient and its own "
        "background nuisance amplitudes. UVIS fits use the documented "
        "physical-column proxy, since the raw optical-depth scale is not an "
        "absolute column calibration.", "",
        "| Shape strategy | E3 rate (s^-1) | E5 rate (s^-1) | 2007 UVIS rate (s^-1) | 2010 UVIS rate (s^-1) |",
        "|---|---:|---:|---:|---:|",
    ]
    strategy_display = {
        "independent_each_observation": "Independent shape for every observation",
        "shared_shape_all_observations": "Shape fitted to all observations",
        "shared_shape_inms_only": "Shape fitted to all INMS observations",
        "shared_shape_uvis_only": "Shape fitted to all UVIS observations",
    }
    for strategy, display in strategy_display.items():
        target_fits = rate_strategies[strategy]["target_fits"]
        report.append(
            f"| {display} | "
            f"{target_fits['E3']['total_source_rate_s']:.6e} | "
            f"{target_fits['E5']['total_source_rate_s']:.6e} | "
            f"{target_fits['uvis_2007_zeta_orionis']['total_source_rate_s']:.6e} | "
            f"{target_fits['uvis_2010_e10_solar']['total_source_rate_s']:.6e} |")

    report += ["", "## Model-component ablation", "",
               "The table reports the direct chi-square increase when a source family is removed. "
               "Lower is better; the UVIS values refer only to optical-depth shape.", "",
               "| Dataset | Full chi-square | No stripes | No jets | Background only |",
               "|---|---:|---:|---:|---:|"]
    for dataset in ("E3", "E5", "uvis_2007_zeta_orionis", "uvis_2010_e10_solar"):
        values = {row["variant"]: row["chi2"] for row in comparison_rows
                  if row["dataset"] == dataset}
        report.append(
            f"| {dataset} | {values['full']:.3f} | {values['no_stripes']:.3f} | "
            f"{values['no_jets']:.3f} | {values['background_only']:.3f} |")

    report += ["", "## Cross-epoch test", "",
               "A single frozen source map does not transfer between these observations. "
               "The balanced shared-source diagnostic has chi-square per sample of "
               f"{cross['steady_balanced_diagnostic']['E3']['chi2_per_sample']:.2f} (E3), "
               f"{cross['steady_balanced_diagnostic']['E5']['chi2_per_sample']:.2f} (E5), "
               f"{cross['steady_balanced_diagnostic']['uvis_2007_zeta_orionis']['chi2_per_sample']:.2f} (2007 UVIS), and "
               f"{cross['steady_balanced_diagnostic']['uvis_2010_e10_solar']['chi2_per_sample']:.2f} (2010 UVIS).", "",
               "The E3+E5 combined vector has chi-square per sample of "
               f"{cross['INMS_combined']['E3']['chi2_per_sample']:.2f} and "
               f"{cross['INMS_combined']['E5']['chi2_per_sample']:.2f} on its two "
               "fitted INMS profiles. The combined UVIS vector has shape chi-square "
               "per sample of "
               f"{cross['UVIS_combined']['uvis_2007_zeta_orionis']['chi2_per_sample']:.2f} "
               "and "
               f"{cross['UVIS_combined']['uvis_2010_e10_solar']['chi2_per_sample']:.2f} "
               "on its two fitted occultations.", "",
               f"The E3-derived map overpredicts the two published UVIS peak columns by factors "
               f"{cross['E3_independent']['uvis_2007_zeta_orionis']['peak_column_ratio_to_published']:.2f} and "
               f"{cross['E3_independent']['uvis_2010_e10_solar']['peak_column_ratio_to_published']:.2f}; "
               f"the E5-derived map overpredicts them by factors "
               f"{cross['E5_independent']['uvis_2007_zeta_orionis']['peak_column_ratio_to_published']:.2f} and "
               f"{cross['E5_independent']['uvis_2010_e10_solar']['peak_column_ratio_to_published']:.2f}. "
               "This is evidence against a time-invariant source-amplitude map under the present "
               "model, but it does not by itself distinguish true temporal variability from missing "
               "collisions, calibration response, or source geometry.", "",
               "The two UVIS-derived maps transfer more closely in absolute column: the 2007 map "
               f"predicts {cross['uvis_2007_zeta_orionis_independent']['uvis_2010_e10_solar']['peak_column_ratio_to_published']:.2f} times the 2010 column, "
               "and the 2010 map predicts "
               f"{cross['uvis_2010_e10_solar_independent']['uvis_2007_zeta_orionis']['peak_column_ratio_to_published']:.2f} times the 2007 column. "
               "Their transferred shapes are still statistically poor.", "",
               "## Identifiability", "",
               f"The E3 and E5 scaled design matrices have condition numbers "
               f"{density_fits['e3']['condition_number_scaled']:.3e} and "
               f"{density_fits['e5']['condition_number_scaled']:.3e}. Their effective ranks at "
               f"the 1e-3 relative singular-value cutoff are "
               f"{density_fits['e3']['effective_rank_by_relative_singular_cutoff']['1e-3']} and "
               f"{density_fits['e5']['effective_rank_by_relative_singular_cutoff']['1e-3']} out of "
               f"{len(full_labels)} fitted columns. Individual source rates therefore must not be "
               "interpreted as uniquely measured.", ""]
    if stripe_count:
        report += [
            "The `no_stripes` ablation removes all line-stripe components, while "
            "`no_jets` removes all discrete point sources. The numerical changes "
            "in the table, rather than source labels alone, quantify which family "
            "is required by each observation.", ""]
    else:
        report += [
            "This case contains no separate line-stripe components, so `no_stripes` "
            "is identical to `full`. The `no_jets` and `background_only` rows both "
            "remove every physical source and therefore are also identical.", ""]
    report += ["## Numerical verification", ""]
    for event in uvis_files:
        item = convergence[event]
        report.append(
            f"- {event}: tangent orthogonality error {item['ray_tangent_orthogonality_max']:.3e}; "
            f"fitted total-column quadrature difference at order {args.quadrature_order} versus "
            f"{item['comparison_order']} is at most "
            f"{item['max_relative_difference_fitted_total_column']:.3e} on the sampled material rays.")
    report += ["", "## Numerical and inference cautions", ""]
    report.extend(f"- {item}" for item in summary["interpretation_limits"])
    report += ["", "See `source_rate_shape_strategies.csv`, "
               "`model_comparison.csv`, `cross_prediction_metrics.csv`, "
               "`component_amplitudes.csv`, and `study_results.json` for the "
               "full diagnostics.", ""]
    (args.output_dir/"fitting_study.md").write_text("\n".join(report), encoding="utf-8")
    print(json.dumps({"output": str(args.output_dir),
                      "total_source_rate_s": totals,
                      "quadrature_convergence": convergence},
                     indent=2, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
