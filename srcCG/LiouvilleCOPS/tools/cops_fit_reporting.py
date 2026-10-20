#!/usr/bin/env python3
"""Deterministic Stage-A5 tables and figures for one LiouvilleCOPS fit.

This module does not choose data, fit parameters, or regularization.  It turns
already selected observations and already evaluated model/activity arrays into
self-identifying report products.  Every table header and PNG metadata block
carries the same SHA-256 report identity, allowing readers to detect accidental
mixing of profiles, coefficients, degrees, or lambda values.
"""

import datetime as dt
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from cops_io import canonical_sha256
from cops_selection import parse_utc


REPORT_SCHEMA = "liouville_cops.fit_report"
REPORT_SCHEMA_VERSION = 1


def report_identity(payload):
    """Return the canonical identity of the scientific choices in a report."""
    return canonical_sha256({"schema": REPORT_SCHEMA,
                             "schema_version": REPORT_SCHEMA_VERSION,
                             "payload": payload})


def _metadata(identity, description):
    """Portable PNG text metadata accepted by Matplotlib's PNG backend."""
    return {"Title": description, "Description": "%s v%d; identity=%s" %
            (REPORT_SCHEMA, REPORT_SCHEMA_VERSION, identity)}


def _header(stream, identity, description):
    stream.write("# ReportSchema=%s ReportSchemaVersion=%d ReportIdentity=%s\n" %
                 (REPORT_SCHEMA, REPORT_SCHEMA_VERSION, identity))
    stream.write("# " + description + "\n")


def _utc_array(values):
    base = dt.datetime(2016, 1, 1)
    return np.array([np.datetime64(base + dt.timedelta(seconds=parse_utc(value)))
                     for value in values])


def _interval_times(intervals):
    return [(_utc_array([start])[0], _utc_array([stop])[0])
            for start, stop in intervals]


def _shade_intervals(axis, intervals, x_kind, utc, altitude, label=None):
    """Shade UTC intervals on either the time or trajectory-altitude axis."""
    for index, (start, stop) in enumerate(_interval_times(intervals)):
        if x_kind == "time":
            lo, hi = start, stop
        else:
            mask = (utc >= start) & (utc <= stop) & np.isfinite(altitude)
            if not np.any(mask):
                continue
            lo, hi = np.nanmin(altitude[mask]), np.nanmax(altitude[mask])
        axis.axvspan(lo, hi, color="#f2a900", alpha=0.10,
                    label=label if index == 0 else None, zorder=0)


def write_pressure_products(directory, identity, identity_label, point_ids, utc_text,
                            altitude, predictions, references, used, saturated,
                            reasons, transition_intervals=()):
    """Write per-point residuals and time/altitude comparison figures."""
    directory = Path(directory)
    utc = _utc_array(utc_text)
    altitude = np.asarray(altitude, dtype=float)
    saturated = np.asarray(saturated, dtype=bool)
    pred = {key: np.asarray(value, dtype=float) for key, value in predictions.items()}
    ref = {key: np.asarray(value, dtype=float) for key, value in references.items()}
    use = {key: np.asarray(value, dtype=bool) for key, value in used.items()}

    table_path = directory / "pressure_residuals.dat"
    with table_path.open("w") as stream:
        _header(stream, identity,
                "model/reference pressures and log residuals; non-selected residuals are nan")
        stream.write("#iPoint UTC altitude_m gauge selected saturated observed_Pa "
                     "model_Pa log_residual normalized_log_residual reason\n")
        for point in range(len(utc_text)):
            for gauge in ("NG", "RG"):
                valid = (use[gauge][point] and ref[gauge][point] > 0.0 and
                         pred[gauge][point] > 0.0 and
                         np.isfinite(pred[gauge][point]))
                residual = (np.log(pred[gauge][point] / ref[gauge][point])
                            if valid else np.nan)
                count = np.count_nonzero(use[gauge])
                normalized = residual / np.sqrt(count) if count else np.nan
                stream.write("%d %s %.9e %s %d %d %.15e %.15e %s %s %s\n" %
                             (point_ids[point], utc_text[point], altitude[point], gauge,
                              int(use[gauge][point]),
                              int(saturated[point] and gauge == "NG"),
                              ref[gauge][point], pred[gauge][point],
                              "%.15e" % residual if np.isfinite(residual) else "nan",
                              "%.15e" % normalized if np.isfinite(normalized) else "nan",
                              reasons[gauge][point]))

    colors = {"NG": "#1f5fbf", "RG": "#d62728"}
    markers = {"NG": "o", "RG": "s"}
    outputs = [table_path.name]
    for x_kind, x, xlabel, filename in (
            ("time", utc, "UTC on 30 September 2016",
             "pressure_comparison_time.png"),
            ("altitude", altitude, "altitude above shape model [m]",
             "pressure_comparison_altitude.png")):
        fig, axes = plt.subplots(2, 1, figsize=(13.5, 9), sharex=True)
        for axis, gauge in zip(axes, ("NG", "RG")):
            available = (ref[gauge] > 0.0) & np.isfinite(ref[gauge])
            selected = available & use[gauge]
            excluded = available & ~use[gauge]
            normal_selected = selected & (~saturated if gauge == "NG" else
                                          np.ones_like(selected, dtype=bool))
            axis.plot(x, pred[gauge], color=colors[gauge], lw=1.35,
                      label="fitted model", zorder=2)
            axis.plot(x[normal_selected], ref[gauge][normal_selected],
                      markers[gauge], ms=3.2, color="black",
                      label="selected observation (%d)" % selected.sum(), zorder=4)
            axis.plot(x[excluded], ref[gauge][excluded], markers[gauge], ms=4.0,
                      mfc="none", mec="0.55", mew=0.8,
                      label="excluded observation (%d)" % excluded.sum(), zorder=3)
            if gauge == "NG":
                sat_available = available & saturated
                axis.plot(x[sat_available], ref[gauge][sat_available], "x",
                          color="#e68a00", ms=5.2, mew=1.0,
                          label="saturated NG (%d; selected %d)" %
                                (sat_available.sum(), (sat_available & selected).sum()),
                          zorder=5)
            _shade_intervals(axis, transition_intervals, x_kind, utc, altitude,
                             "instrument-transition interval")
            axis.set_yscale("log")
            axis.set_ylabel("%s pressure [Pa]" % gauge)
            axis.grid(True, which="both", alpha=0.22)
            axis.legend(fontsize=8, loc="best", frameon=False)
        if x_kind == "time":
            axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        else:
            axes[-1].set_xscale("log")
        axes[-1].set_xlabel(xlabel)
        fig.suptitle("COPS observations and fitted pressures\n%s" % identity_label,
                     fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(directory / filename, dpi=130,
                    metadata=_metadata(identity, identity_label))
        plt.close(fig)
        outputs.append(filename)

    fig, axes = plt.subplots(2, 2, figsize=(15, 8), sharey="row")
    for row, gauge in enumerate(("NG", "RG")):
        residual = np.full(len(utc), np.nan)
        mask = use[gauge] & (ref[gauge] > 0.0) & (pred[gauge] > 0.0)
        residual[mask] = np.log(pred[gauge][mask] / ref[gauge][mask])
        for col, (x_kind, x, xlabel) in enumerate((
                ("time", utc, "UTC"), ("altitude", altitude, "altitude [m]"))):
            axis = axes[row, col]
            axis.axhline(0.0, color="0.45", lw=0.8)
            axis.plot(x[mask], residual[mask], ".", ms=3, color=colors[gauge])
            if gauge == "NG":
                sm = mask & saturated
                axis.plot(x[sm], residual[sm], "x", ms=5, color="#e68a00")
            _shade_intervals(axis, transition_intervals, x_kind, utc, altitude)
            axis.set_xlabel(xlabel)
            axis.set_ylabel("ln(model/observed) %s" % gauge)
            axis.grid(True, alpha=0.22)
            if x_kind == "time":
                axis.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            else:
                axis.set_xscale("log")
    fig.suptitle("Selected-data log residuals\n%s" % identity_label, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    residual_plot = "pressure_residuals.png"
    fig.savefig(directory / residual_plot, dpi=130,
                metadata=_metadata(identity, identity_label))
    plt.close(fig)
    outputs.append(residual_plot)
    return outputs


def _region_row(name, mask, area, baseline, residual, total_linear, clipped):
    region_area = np.sum(area[mask])
    all_area = np.sum(area)
    if region_area <= 0.0:
        return [name, 0.0] + [np.nan] * 8
    mean = lambda value: np.sum(area[mask] * value[mask]) / region_area
    return [name, region_area / all_area, mean(baseline), mean(residual),
            mean(total_linear), mean(clipped), mean(baseline ** 2),
            mean(residual ** 2), mean(total_linear ** 2),
            np.sum(area[mask & (total_linear <= 0.0)]) / region_area]


def write_activity_products(directory, identity, identity_label, centroid, area,
                            seen, baseline, residual, total_linear, clipped,
                            degree_components, template_components=()):
    """Write degree/region diagnostics and baseline-residual-total maps."""
    directory = Path(directory)
    centroid = np.asarray(centroid, dtype=float)
    area = np.asarray(area, dtype=float)
    seen = np.asarray(seen, dtype=bool)
    fields = [np.asarray(value, dtype=float)
              for value in (baseline, residual, total_linear, clipped)]
    baseline, residual, total_linear, clipped = fields
    baseline_power = np.sum(area * baseline ** 2)

    # This is the authoritative A5 facet table for both fitting modes.  The
    # older ActivityEnhancementFacets.dat remains an A3/A4 exchange diagnostic
    # and exists only when the optional enhancement is enabled.
    facet_path = directory / "activity_facets.dat"
    with facet_path.open("w") as stream:
        _header(stream, identity,
                "H2O activity at shape-facet centroids; seen means a facet "
                "contributes to at least one loaded selected-data kernel")
        stream.write("#iFacet x_m y_m z_m area_m2 seen legacy_linear "
                     "residual_linear total_linear clipped_activity\n")
        for index, (xyz, facet_area) in enumerate(zip(centroid, area)):
            stream.write("%d %.9e %.9e %.9e %.9e %d %.15e %.15e %.15e %.15e\n" %
                         (index, xyz[0], xyz[1], xyz[2], facet_area,
                          int(seen[index]), baseline[index], residual[index],
                          total_linear[index], clipped[index]))

    degree_path = directory / "activity_degree_power.dat"
    with degree_path.open("w") as stream:
        _header(stream, identity,
                "unclipped activity-residual power by harmonic degree")
        stream.write("#degree n_modes coefficient_power spherical_field_power mesh_area_mean_square relative_mesh_power\n")
        for item in degree_components:
            degree, coefficients, component = item[:3]
            spherical_power = item[3] if len(item) > 3 else np.nan
            coefficients = np.asarray(coefficients, dtype=float)
            component = np.asarray(component, dtype=float)
            mesh_power = np.sum(area * component ** 2) / np.sum(area)
            stream.write("%d %d %.15e %.15e %.15e %.15e\n" %
                         (degree, len(coefficients), np.sum(coefficients ** 2),
                          spherical_power, mesh_power,
                          np.sum(area * component ** 2) / baseline_power))

    template_path = directory / "activity_template_power.dat"
    with template_path.open("w") as stream:
        _header(stream, identity,
                "fitted localized-template amplitudes and unclipped activity power")
        stream.write("#template coefficient mesh_area_mean_square relative_mesh_power\n")
        for name, coefficient, component in template_components:
            component = np.asarray(component, dtype=float)
            mesh_power = np.sum(area * component ** 2) / np.sum(area)
            stream.write("%s %.15e %.15e %.15e\n" %
                         (name, coefficient, mesh_power,
                          np.sum(area * component ** 2) / baseline_power))

    region_path = directory / "activity_region_diagnostics.dat"
    rows = [_region_row("all", np.ones(len(area), dtype=bool), area,
                        baseline, residual, total_linear, clipped),
            _region_row("seen", seen, area, baseline, residual,
                        total_linear, clipped),
            _region_row("unseen", ~seen, area, baseline, residual,
                        total_linear, clipped)]
    with region_path.open("w") as stream:
        _header(stream, identity,
                "area-weighted H2O activity diagnostics by trajectory visibility")
        stream.write("#region area_fraction mean_baseline mean_residual mean_total_linear mean_clipped baseline_mean_square residual_mean_square total_mean_square clipped_area_fraction\n")
        for row in rows:
            stream.write("%s " % row[0] + " ".join("%.15e" % value
                                                    for value in row[1:]) + "\n")

    radius = np.linalg.norm(centroid, axis=1)
    longitude = np.degrees(np.arctan2(centroid[:, 1], centroid[:, 0]))
    latitude = np.degrees(np.arcsin(np.clip(centroid[:, 2] / radius, -1.0, 1.0)))
    lon_edges = np.linspace(-180.0, 180.0, 241)
    lat_edges = np.linspace(-90.0, 90.0, 121)
    area_grid = np.histogram2d(latitude, longitude,
                               bins=(lat_edges, lon_edges), weights=area)[0]

    def grid(value):
        weighted = np.histogram2d(latitude, longitude,
                                  bins=(lat_edges, lon_edges),
                                  weights=area * value)[0]
        return np.divide(weighted, area_grid,
                         out=np.full_like(weighted, np.nan), where=area_grid > 0)

    maps = [("legacy baseline (linear)", grid(baseline), "viridis"),
            ("fitted residual (linear)", grid(residual), "coolwarm"),
            ("applied total (clipped)", grid(clipped), "viridis")]
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True, sharey=True)
    for axis, (title, values, cmap) in zip(axes, maps):
        finite = values[np.isfinite(values)]
        if cmap == "coolwarm":
            limit = np.percentile(np.abs(finite), 99.0) if len(finite) else 1.0
            limit = max(limit, np.finfo(float).tiny)
            vmin, vmax = -limit, limit
        else:
            vmin = np.percentile(finite, 1.0) if len(finite) else 0.0
            vmax = np.percentile(finite, 99.0) if len(finite) else 1.0
            if vmax <= vmin:
                vmax = vmin + 1.0
        image = axis.pcolormesh(lon_edges, lat_edges, values, shading="auto",
                                cmap=cmap, vmin=vmin, vmax=vmax)
        fig.colorbar(image, ax=axis, label="dimensionless activity")
        axis.set_title(title, fontsize=9)
        axis.set_ylabel("body-frame latitude [deg]")
        axis.grid(True, alpha=0.15)
    axes[-1].set_xlabel("body-frame longitude [deg]")
    fig.suptitle("Area-weighted facet activity maps\n%s" % identity_label,
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    map_name = "activity_maps.png"
    fig.savefig(directory / map_name, dpi=130,
                metadata=_metadata(identity, identity_label))
    plt.close(fig)
    return [facet_path.name, degree_path.name, template_path.name,
            region_path.name, map_name], {
        "regions": {row[0]: {
            "area_fraction": float(row[1]),
            "residual_mean_square": float(row[7]),
            "clipped_area_fraction": float(row[9]),
        } for row in rows},
        "degree_count": len(degree_components),
        "template_count": len(template_components),
    }


def identity_label(payload):
    """Compact visible label shared by every report plot."""
    profile = payload.get("profile") or "legacy CLI selection"
    model = ("enhancement L=%s" % payload.get("enhancement_degree")
             if payload.get("enhancement_enabled") else
             "legacy activity degree=%s" % payload.get("legacy_degree"))
    return ("profile=%s; %s; lambda=%.4g; selection=%s; report=%s" %
            (profile, model, payload["lambda"],
             payload["selection_sha256"][:12],
             report_identity(payload)[:12]))


def write_report_descriptor(directory, identity, payload, products):
    """Write a human/machine-readable index before output hashing."""
    path = Path(directory) / "fit_report.json"
    path.write_text(json.dumps({
        "schema": REPORT_SCHEMA,
        "schema_version": REPORT_SCHEMA_VERSION,
        "report_identity": identity,
        "identity_payload": payload,
        "products": list(products),
    }, indent=2, sort_keys=True) + "\n")
    return path.name
