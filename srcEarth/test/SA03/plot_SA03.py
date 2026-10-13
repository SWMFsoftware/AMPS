#!/usr/bin/env python3
"""Generate independent PNG/EPS SA03 figures from saved comparison artifacts.

Plotting is deliberately downstream of comparison and has no authority to rescore a
gate. Every figure is derived from persisted comparison/context tables, and the plot
manifest hashes those inputs, this implementation, and every requested output. This
makes post-run plotting reproducible while preventing a figure-only success from being
mistaken for physical execution or scientific qualification.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sa03_common import (CHANNELS, SPACECRAFT, atomic_json, read_csv, resolve_data_root,
                         sha256_file, utc_now, validate_spacecraft)


COLORS = dict(zip((name for name, _, _ in CHANNELS), plt.cm.viridis(np.linspace(0.05, 0.95, 7))))


def _spacecraft(data: pd.DataFrame) -> tuple[str, ...]:
    """Respect explicit single-spacecraft runs in every figure family."""
    return tuple(item for item in SPACECRAFT if item in set(data.spacecraft))


def _series(axis: object, rows: pd.DataFrame, column: str, **style: object) -> None:
    """Break lines at gaps exceeding 1.5 nominal sample intervals.

    Preserve all supplied values. Logarithmic callers separately annotate physical
    zeros; this utility also serves linear driver/orbit plots without changing data.
    """
    times = _time_values(rows.time)
    values = rows[column].to_numpy(dtype=float).copy()
    if len(times) > 1:
        delta = np.diff(times).astype("timedelta64[s]").astype(float)
        cadence = float(np.min(delta[delta > 0])) if np.any(delta > 0) else 0
        breaks = np.flatnonzero(delta > 1.5 * cadence) + 1 if cadence else []
        times = np.insert(times, breaks, times[breaks])
        values = np.insert(values, breaks, np.nan)
    style.setdefault("marker", ".")  # retain isolated observations across gaps
    style.setdefault("markersize", 2)
    axis.plot(times, values, **style)


def _time_values(values: object) -> np.ndarray:
    """Return timezone-naive UTC values accepted by this host's Matplotlib 3.1."""

    converted = pd.to_datetime(values, utc=True)
    if isinstance(converted, pd.Series):
        return converted.dt.tz_convert(None).to_numpy()
    return converted.tz_convert(None).to_numpy()


def _save(fig: plt.Figure, directory: Path, stem: str, formats: tuple[str, ...]) -> list[str]:
    """Save each requested format and restore artists after EPS opacity conversion."""
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for extension in formats:
        path = directory / f"{stem}.{extension}"
        # EPS has no transparency. Make its existing opaque fallback explicit,
        # then restore artists so PNG appearance never depends on format order.
        opaque = [(artist, artist.get_alpha()) for artist in fig.findobj()
                  if extension == "eps" and hasattr(artist, "get_alpha")
                  and artist.get_alpha() is not None and 0 < artist.get_alpha() < 1]
        try:
            for artist, _ in opaque:
                artist.set_alpha(1)
            fig.savefig(path, dpi=180 if extension == "png" else None, bbox_inches="tight")
        finally:
            for artist, alpha in opaque:
                artist.set_alpha(alpha)
        paths.append(str(path))
    plt.close(fig)
    return paths


def _load(output: Path, spacecraft: tuple[str, ...], product: str) -> pd.DataFrame:
    """Load saved comparison rows strictly; plotting never regenerates model values."""
    frames = []
    for item in spacecraft:
        path = output / item / "comparison" / f"{product}_comparison.csv"
        frame = pd.DataFrame(read_csv(path))
        for column in frame.columns:
            if column not in ("spacecraft", "utc", "channel", "reference_interpretation"):
                frame[column] = pd.to_numeric(frame[column], errors="raise")
        frame["time"] = pd.to_datetime(frame["utc"], utc=True)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _time_series_channels(data: pd.DataFrame, figures: Path, formats: tuple[str, ...]) -> list[str]:
    """Show all seven channels, with gaps/zeros visible and no favorable selection."""
    ships = _spacecraft(data)
    fig, axes = plt.subplots(len(ships), 1, figsize=(12, 4 * len(ships)), sharex=True)
    axes = np.atleast_1d(axes)
    for axis, spacecraft in zip(axes, ships):
        subset = data[data.spacecraft == spacecraft]
        for channel, _, _ in CHANNELS:
            rows = subset[subset.channel == channel].sort_values("time")
            _series(axis, rows, "reference", color=COLORS[channel], lw=1.1, label=f"{channel} REPT")
            _series(axis, rows, "model", color=COLORS[channel], lw=0.8, ls="--", label=f"{channel} AMPS")
        axis.text(0.01, 0.02, f"Physical model zeros: {int((subset.model == 0).sum())}",
                  transform=axis.transAxes, bbox=dict(facecolor="white", edgecolor="none"))
        axis.set_yscale("log"); axis.set_ylabel(f"{spacecraft}\nJ [cm$^{{-2}}$ s$^{{-1}}$ sr$^{{-1}}$ MeV$^{{-1}}$]")
        axis.grid(True, which="both", alpha=0.25)
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    handles, labels = axes[0].get_legend_handles_labels()
    # Anchor to the axes: old Matplotlib repositions figure-level legends during
    # bbox_inches='tight', which can move them onto an otherwise separated title.
    axes[0].legend(handles[::2] + handles[1::2], labels[::2] + labels[1::2], ncol=7,
                   fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.015))
    axes[0].set_title("SA03 seven-channel differential intensity: solid REPT, dashed AMPS", pad=48)
    return _save(fig, figures, "sa03_channel_time_series", formats)


def _representative_spectra(data: pd.DataFrame, figures: Path, formats: tuple[str, ...]) -> list[str]:
    """Choose first/middle/last matched epochs without examining their agreement."""
    ships = _spacecraft(data)
    fig, axes = plt.subplots(len(ships), 3, figsize=(12, 3.5 * len(ships)), sharex=True, sharey=True, squeeze=False)
    for row_axes, spacecraft in zip(axes, ships):
        subset = data[data.spacecraft == spacecraft]
        times = sorted(subset.utc.unique())
        selected = [times[0], times[len(times) // 2], times[-1]]
        for axis, utc in zip(row_axes, selected):
            rows = subset[subset.utc == utc].sort_values("energy_min_mev")
            centers = np.sqrt(rows.energy_min_mev * rows.energy_max_mev)
            axis.plot(centers.to_numpy(), rows.reference.to_numpy(), "o-", label="REPT")
            axis.plot(centers.to_numpy(), rows.model.to_numpy(), "s--", label="AMPS folded")
            axis.set(xscale="log", yscale="log", title=f"{spacecraft} {utc[5:16]} UTC")
            axis.set_xlim(18.5, 116)
            axis.set_xticks([20, 40, 80, 116]); axis.set_xticklabels(["20", "40", "80", "116"])
            axis.tick_params(axis="x", which="minor", labelbottom=False)
            zeros = int((rows.model == 0).sum())
            if zeros:
                axis.text(0.02, 0.04, f"Physical zeros: {zeros} (off log scale)",
                          transform=axis.transAxes, fontsize=8)
            axis.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend()
    # Figure.supxlabel/supylabel postdate the supported Matplotlib 3.1 runtime.
    fig.text(0.5, 0.01, "Proton energy [MeV]", ha="center")
    fig.text(0.01, 0.5, "Differential intensity [cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$]", va="center", rotation=90)
    fig.suptitle("Deterministic first, middle, and last matched spectra")
    return _save(fig, figures, "sa03_representative_spectra", formats)


def _scalar_series(data: pd.DataFrame, figures: Path, formats: tuple[str, ...],
                   stem: str, title: str, ylabel: str) -> list[str]:
    """Plot native integrated quantities against their correlated REPT references."""
    ships = _spacecraft(data)
    fig, axes = plt.subplots(len(ships), 1, figsize=(12, 3 * len(ships)), sharex=True)
    axes = np.atleast_1d(axes)
    for axis, spacecraft in zip(axes, ships):
        rows = data[data.spacecraft == spacecraft].sort_values("time")
        _series(axis, rows, "reference", lw=1.2, label="REPT-derived reference")
        _series(axis, rows, "model", lw=1.0, ls="--", label="AMPS")
        axis.text(0.01, 0.02, f"Physical model zeros: {int((rows.model == 0).sum())}", transform=axis.transAxes)
        axis.set_yscale("log"); axis.set_ylabel(f"{spacecraft}\n{ylabel}"); axis.grid(True, which="both", alpha=0.25)
    axes[0].legend(); axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    fig.suptitle(title)
    return _save(fig, figures, stem, formats)


def _scatter(datasets: dict[str, pd.DataFrame], figures: Path, formats: tuple[str, ...]) -> list[str]:
    """Display zeros explicitly below the log range without changing saved values."""
    fig, axes = plt.subplots(2, 2, figsize=(9, 9))
    for axis, (product, frame) in zip(axes.flat, datasets.items()):
        positive_values = np.r_[frame.reference.to_numpy(), frame.model[frame.model > 0].to_numpy()]
        low, high = np.nanmin(positive_values), np.nanmax(positive_values)
        for spacecraft in _spacecraft(frame):
            # RBSP-B is the only source-quality-eligible SA03 branch. Keep its
            # marker explicit instead of retaining an unreachable RBSP-A branch:
            # plotting code must not imply that archived RBSP-A audit rows can be
            # promoted into the scientific comparison population.
            marker = "^"
            rows = frame[frame.spacecraft == spacecraft]
            positive = rows.model > 0
            axis.scatter(rows.loc[positive, "reference"].to_numpy(), rows.loc[positive, "model"].to_numpy(), s=9, alpha=0.5,
                         marker=marker, label=spacecraft)
            if (~positive).any():
                axis.scatter(rows.loc[~positive, "reference"].to_numpy(), np.full((~positive).sum(), low / 3), marker="x",
                             label=f"{spacecraft} zero (shown below positive range)")
        guide = np.geomspace(low, high, 100)
        axis.plot(guide, guide, "k-", lw=1); axis.plot(guide, 2 * guide, "k:", lw=0.8); axis.plot(guide, guide / 2, "k:", lw=0.8)
        unit = {"spectrum": "cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$",
                "channel_flux": "cm$^{-2}$ s$^{-1}$", "total_flux": "cm$^{-2}$ s$^{-1}$",
                "density": "cm$^{-3}$"}[product]
        axis.set(xscale="log", yscale="log", title=product.replace("_", " ").title(),
                 xlabel=f"Reference [{unit}]", ylabel=f"Model [{unit}]")
        axis.grid(True, which="both", alpha=0.2)
    axes[0, 0].legend(); fig.suptitle("SA03 model/reference scatter (unity and factor-two guides)")
    return _save(fig, figures, "sa03_model_reference_scatter", formats)


def _time_energy_maps(data: pd.DataFrame, figures: Path, formats: tuple[str, ...]) -> list[str]:
    """Use a labelled chronological index so missing times do not imply exposure."""
    outputs = []
    for spacecraft in _spacecraft(data):
        subset = data[data.spacecraft == spacecraft].copy()
        table_ref = subset.pivot(index="channel", columns="time", values="reference").reindex([x[0] for x in CHANNELS])
        table_mod = subset.pivot(index="channel", columns="time", values="model").reindex([x[0] for x in CHANNELS])
        residual = np.log10(np.where(table_mod > 0, table_mod / table_ref, np.nan))
        fig, axes = plt.subplots(3, 1, figsize=(12, 7), sharex=True)
        for axis, values, label, cmap in zip(
            axes, (np.log10(table_ref), np.log10(table_mod.where(table_mod > 0)), residual),
            ("log10 REPT intensity", "log10 AMPS intensity", "log10(AMPS/REPT)"),
            ("viridis", "viridis", "coolwarm"),
        ):
            image = axis.imshow(values, aspect="auto", origin="lower", interpolation="nearest", cmap=cmap,
                                extent=(0, values.shape[1], 0.5, 7.5))
            axis.set_yticks(range(1, 8)); axis.set_yticklabels([x[0] for x in CHANNELS]); axis.set_ylabel(label)
            fig.colorbar(image, ax=axis, pad=0.01)
        # A physical zero and an absent reference are scientifically different.
        # Both are off the logarithmic color scale, so mark only actual zero-model
        # cells with crosses. Blank unmarked cells remain missing, not exclusions
        # manufactured to improve numerical or observational agreement scores.
        zero_row, zero_column = np.where(table_mod.to_numpy() == 0)
        for axis in axes[1:]:
            axis.scatter(zero_column + 0.5, zero_row + 1, marker="x", c="black", s=12, linewidths=0.6)
        axes[-1].set_xlabel("Matched epoch index (chronological)")
        fig.suptitle(f"{spacecraft} observed/model/residual maps; crosses: {len(zero_row)} physical zeros; blanks: missing")
        outputs.extend(_save(fig, figures, f"sa03_{spacecraft.lower()}_time_energy_maps", formats))
    return outputs


def _residual_context(data: pd.DataFrame, figures: Path, formats: tuple[str, ...]) -> list[str]:
    """Show orbit dependence; count negative-infinite zero-model residuals explicitly."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    residual = np.log10(np.where(data.model > 0, data.model / data.reference, np.nan))
    scatter = axes[0].scatter(data.l_shell.to_numpy(), residual, c=data.energy_min_mev.to_numpy(), s=10, cmap="viridis")
    axes[1].scatter(data.mlt_hour.to_numpy(), residual, c=data.energy_min_mev.to_numpy(), s=10, cmap="viridis")
    for axis, label in zip(axes, ("L shell", "MLT [hour]")):
        axis.axhline(0, color="k", lw=0.8); axis.axhline(math.log10(2), color="k", ls=":", lw=0.7)
        axis.axhline(-math.log10(2), color="k", ls=":", lw=0.7); axis.set(xlabel=label, ylabel="log10(AMPS/REPT)")
        axis.grid(True, alpha=0.2)
    fig.colorbar(scatter, ax=axes, label="Channel lower energy [MeV]")
    zeros = int((data.model == 0).sum())
    fig.suptitle("SA03 spectral residuals versus orbit coordinates" +
                 (f"\n{zeros} physical-zero predictions have residual −∞ (not plotted)" if zeros else ""))
    return _save(fig, figures, "sa03_residuals_l_mlt", formats)


def _upstream_context(validation_root: Path, spectrum: pd.DataFrame, figures: Path,
                      formats: tuple[str, ...]) -> list[str]:
    """Read frozen upstream/driver context, with a deterministic near-57-MeV channel."""
    processed = validation_root / "processed" / "SA03"
    erne = pd.DataFrame(read_csv(processed / "soho_erne_5min.csv"))
    qd = pd.DataFrame(read_csv(processed / "qin_denton_5min.csv"))
    erne["time"] = pd.to_datetime(erne["utc"], utc=True)
    qd["time"] = pd.to_datetime(qd["utc"], utc=True)
    flux_col = "flux_cm2_s_sr_mev"
    erne[flux_col] = pd.to_numeric(erne[flux_col])
    center_col = "energy_center_mev"
    erne[center_col] = pd.to_numeric(erne[center_col])
    channel = erne.iloc[(erne[center_col] - 57.0).abs().argsort()[:1]][center_col].iloc[0]
    upstream = erne[erne[center_col] == channel].sort_values("time")
    dst_col = next(name for name in qd.columns if name.lower() in ("dst_nt", "dst"))
    qd[dst_col] = pd.to_numeric(qd[dst_col])
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    _series(axes[0], upstream, flux_col); axes[0].set(yscale="log", ylabel=f"ERNE {channel:g} MeV\nJ [cm$^{{-2}}$ s$^{{-1}}$ sr$^{{-1}}$ MeV$^{{-1}}$]")
    _series(axes[1], qd.sort_values("time"), dst_col); axes[1].set_ylabel("Qin-Denton Dst [nT]")
    for spacecraft in _spacecraft(spectrum):
        rows = spectrum[spectrum.spacecraft == spacecraft].drop_duplicates("utc").sort_values("time")
        _series(axes[2], rows, "l_shell", label=spacecraft)
    axes[2].set_ylabel("L shell"); axes[2].legend(); axes[2].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    for axis in axes: axis.grid(True, alpha=0.25)
    fig.suptitle("SA03 upstream boundary, storm driver, and sampled orbit context")
    return _save(fig, figures, "sa03_orbit_upstream_storm_context", formats)


def _diagnostics(output: Path, spacecraft: tuple[str, ...], spectrum: pd.DataFrame,
                 figures: Path, formats: tuple[str, ...]) -> list[str]:
    """Keep native reconstruction residuals separate from resolution convergence."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for item in spacecraft:
        rows = spectrum[spectrum.spacecraft == item].drop_duplicates("utc").sort_values("time")
        axes[0].plot(_time_values(rows.time), rows.max_unresolved_fraction_at_epoch.to_numpy(), label=item)
        reconstruction = pd.DataFrame(read_csv(output / item / "comparison" / "native_reconstruction.csv"))
        axes[1].plot(_time_values(reconstruction.utc), pd.to_numeric(reconstruction.total_flux_fractional_residual).to_numpy(), label=f"{item} flux")
        axes[1].plot(_time_values(reconstruction.utc), pd.to_numeric(reconstruction.density_fractional_residual).to_numpy(), ls="--", label=f"{item} density")
    axes[0].axhline(0.01, color="r", ls=":", label="unresolved gate"); axes[0].set(ylabel="Maximum unresolved fraction", title="Trajectory support")
    axes[1].set(yscale="log", ylabel="Native/reconstructed fractional residual", title="Independent reconstruction")
    for axis in axes: axis.grid(True, alpha=0.25); axis.legend(fontsize=8)
    fig.suptitle("SA03 unresolved-support and numerical diagnostics")
    return _save(fig, figures, "sa03_unresolved_refinement_diagnostics", formats)


def plot_all(data_path: Path, output: Path, spacecraft: tuple[str, ...],
             formats: tuple[str, ...]) -> list[str]:
    """Render saved evidence only; add convergence figures only when evidence exists."""
    _, validation_root = resolve_data_root(data_path)
    output = output.expanduser().resolve(); figures = output / "figures"
    datasets = {name: _load(output, spacecraft, name) for name in ("spectrum", "channel_flux", "total_flux", "density")}
    paths: list[str] = []
    paths += _time_series_channels(datasets["spectrum"], figures, formats)
    paths += _representative_spectra(datasets["spectrum"], figures, formats)
    paths += _scalar_series(datasets["total_flux"], figures, formats, "sa03_total_flux_time_series",
                            "SA03 18.5–116 MeV omnidirectional integral intensity", "F [cm$^{-2}$ s$^{-1}$]")
    paths += _scalar_series(datasets["density"], figures, formats, "sa03_density_time_series",
                            "SA03 energetic-proton density (REPT reference is isotropic-equivalent)", "n [cm$^{-3}$]")
    paths += _scatter(datasets, figures, formats)
    paths += _time_energy_maps(datasets["spectrum"], figures, formats)
    paths += _residual_context(datasets["spectrum"], figures, formats)
    paths += _upstream_context(validation_root, datasets["spectrum"], figures, formats)
    paths += _diagnostics(output, spacecraft, datasets["spectrum"], figures, formats)
    convergence_path = output / "convergence_summary.csv"
    if convergence_path.is_file():
        convergence = pd.read_csv(convergence_path)
        convergence["time"] = pd.to_datetime(convergence.utc, utc=True)
        fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
        for axis, product in zip(axes, ("channel_flux", "total_flux", "density")):
            for item in spacecraft:
                rows = convergence[(convergence.spacecraft == item) & (convergence.product == product)]
                # Show the worst channel at each common epoch, retaining zeros.
                rows = rows.groupby("time", as_index=False).fractional_change.max()
                _series(axis, rows, "fractional_change", label=item)
            axis.axhline(0.02, color="r", ls=":", label="fixed 2% gate")
            axis.set(ylabel="Fractional change", title=product)
            axis.legend(); axis.grid(True, alpha=0.25)
        fig.suptitle("SA03 numerical refinement on common epochs")
        paths += _save(fig, figures, "sa03_convergence", formats)
    summary = json.loads((output / "comparison_summary.json").read_text(encoding="utf-8"))
    fig, axis = plt.subplots(figsize=(10, 5)); axis.axis("off")
    lines = [f"SA03 observational status: {summary['status']}", "", "Pooled product metrics:"]
    for product, metric in summary["pooled"].items():
        lines.append(f"{product}: N={metric['count']}, median |log10 ratio|={metric['median_abs_log10_ratio']}, factor-two={metric['factor_two_fraction']:.3f}")
    lines += ["", "REPT density is isotropic-equivalent and derived from the same seven spectra.",
              "Static T05 snapshots do not demonstrate dynamic electric-field transport."]
    axis.text(0.02, 0.95, "\n".join(lines), va="top", family="monospace")
    paths += _save(fig, figures, "sa03_validation_overview", formats)
    # Post-processing can occur long after execution, using a newer plotting
    # implementation. Bind these figures to the actual saved comparisons, frozen
    # context and plotting sources instead of relying on the physical run's earlier
    # repository identity. Nothing in this manifest constitutes scientific PASS.
    inputs = [output / "comparison_summary.json"]
    inputs += [path for item in spacecraft for path in (output / item / "comparison").glob("*.csv")]
    inputs += [validation_root / "processed/SA03" / name
               for name in ("soho_erne_5min.csv", "qin_denton_5min.csv")]
    if convergence_path.is_file():
        inputs.append(convergence_path)
    atomic_json(output / "plot_manifest.json", {
        "created_utc": utc_now(), "interpretation": "saved-evidence rendering, not qualification",
        "source_sha256": {name: sha256_file(Path(__file__).with_name(name))
                          for name in ("plot_SA03.py", "sa03_common.py")},
        "input_sha256": {str(path.resolve()): sha256_file(path) for path in inputs},
        "figure_sha256": {str(Path(path).resolve()): sha256_file(Path(path)) for path in paths},
    })
    return paths


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    source = result.add_mutually_exclusive_group(required=True)
    source.add_argument("--data-root", type=Path)
    source.add_argument(
        "--data-archive", type=Path,
        help=("Alias for an extracted package, its model_inputs directory, or its "
              "validation_data directory; not a tar archive."),
    )
    result.add_argument("--output-dir", required=True, type=Path)
    result.add_argument("--spacecraft", nargs="+", default=list(SPACECRAFT))
    result.add_argument("--formats", nargs="+", choices=("png", "eps"), default=["png", "eps"])
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    paths = plot_all(args.data_root or args.data_archive, args.output_dir,
                     validate_spacecraft(args.spacecraft), tuple(args.formats))
    print(f"[PASS] wrote {len(paths)} figure files under {args.output_dir.resolve() / 'figures'}")
    return 0


if __name__ == "__main__":
    try: raise SystemExit(main())
    except (FileNotFoundError, KeyError, ValueError) as error:
        print(f"ERROR: {error}"); raise SystemExit(2)
