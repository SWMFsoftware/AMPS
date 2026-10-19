#!/usr/bin/env python3
"""C00: freeze the CPU scientific and numerical reference for the CUDA port.

Run this tool with the CPU reference interpreter (the one that produced the
campaign, currently /usr/bin/python3 3.8 with NumPy 1.24 / SciPy 1.10), never
with the CUDA environment: the fixtures must be what the reference itself
computes.  Each stage writes one new subdirectory of a capture directory and
refuses to overwrite anything.

Stages
------
new          create ``<COPS_OUTPUT_ROOT>/cuda_reference/c00_<UTC>`` and print it
identities   reference source tree, interpreter, packages, base cache, every
             kernel file, SPICE dependencies and omega, campaign inventory
masks        exact NG-only visibility masks (reference function) plus the
             per-point predicate sets the GPU union must reproduce
fit          pressures, residuals, regularization and analytic Jacobians of
             the reference evaluator for a deterministic production subset,
             representative parameter vectors, branch coverage, exact-sum
             references and branch-aware central differences
joint        C00b (ROADMAP rev. 2): the same fixtures for the joint families
             B29/N30/C30/NC31 and the labelled R30 against the repaired CPU
             tree, with the real NG-only mask, a CPU replica of the RG-only
             mask, a synthetic cliff stand-in and the archived B29 state
gravity      polyhedral_gravity 3.3 reference values on synthetic shapes and
             on the exact 788k gravity body (DISABLE after an independent
             topology audit; the package HEAL is defective at this scale and
             is only replicated on a sample as evidence), summation-order
             repeatability and the offset-halving slope audit

The stages read only; they never modify the reference tree, the base cache, or
running campaign products.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))
from cops_cuda import provenance as pv  # noqa: E402
from cops_cuda import reference  # noqa: E402


CAPTURE_SCHEMA = "liouville_cops_cuda.c00_capture"
CAPTURE_VERSION = 1
DEFAULT_BASE = "output_fit_shap5_788k_base"
DEFAULT_PROFILE = "all_available"
# Campaign late-NG interval (ROADMAP 3.2); overridable on the command line.
LATE_INTERVAL = ("2016-09-30T10:33:10", "2016-09-30T10:39:26")
# A second, earlier interval exercises the same predicate on a geometry with
# both gauges active.  It is a test case, not a scientific template.
CONTROL_INTERVAL = ("2016-09-30T09:00:00", "2016-09-30T09:10:00")
DEFAULT_DENSITY = 533.0
DEFAULT_OFFSET = 0.1
EPS = np.finfo(float).eps


# ----------------------------------------------------------------------------- common
def log(message):
    print("[c00 %s] %s" % (time.strftime("%H:%M:%S"), message), flush=True)


def stage_directory(capture, name):
    capture = Path(capture).resolve()
    if not (capture / "capture.json").is_file():
        raise SystemExit("Error: %s is not a C00 capture directory (run 'new' first)" % capture)
    return pv.make_new_directory(capture / name)


def finish_stage(out, stage, summary, started):
    """Write the stage manifest with the identity of every product file."""
    root = reference.bind()
    products = [pv.file_identity(path) for path in sorted(Path(out).iterdir())
                if path.is_file()]
    manifest = {
        "schema": CAPTURE_SCHEMA, "schema_version": CAPTURE_VERSION,
        "stage": stage, "started_utc": started, "finished_utc": pv.utc_now(),
        "process": pv.process_context(),
        "reference_tree_sha256": pv.tree_identity(
            root, reference.SOURCE_EXCLUDES)["tree_sha256"],
        "capture_tool": pv.file_identity(__file__),
        "products": products, "summary": summary,
    }
    pv.write_json_exclusive(Path(out) / "stage_manifest.json", manifest)
    log("stage %s complete: %s" % (stage, out))


def save_npz(path, **arrays):
    """Compressed, pickle-free arrays; refuses an existing path."""
    with Path(path).open("xb") as stream:
        np.savez_compressed(stream, **arrays)


def load_reference_modules():
    return {name: reference.load(name) for name in (
        "cops_paths", "cops_io", "cops_selection", "fit_cops", "gravity_cache",
        "build_physical_activity_templates", "activity_basis")}


def base_directory(modules, base):
    return Path(modules["cops_paths"].product_input(base)).resolve()


def read_base_table(modules, run):
    headers, names, rows, columns = modules["cops_io"].read_model_table(
        run / "COPS_model.dat")
    return headers, names, rows, columns


def kernel_index(capture):
    """Per-file kernel identities recorded by the identities stage, if present."""
    path = Path(capture) / "identities" / "kernel_files.json"
    if not path.is_file():
        return None
    return {entry["point"]: entry for entry in json.loads(path.read_text())["kernels"]}


def verify_used_kernels(capture, run, points):
    """Check every kernel a stage reads against the frozen per-file digests."""
    index = kernel_index(capture)
    if index is None:
        raise SystemExit("Error: run the identities stage first; kernel digests are needed")
    for point in points:
        path = run / "kernels" / ("point=%04d.bin" % point)
        if pv.sha256_file(path) != index[point]["sha256"]:
            raise SystemExit("Error: kernel %s changed after the identities stage" % path)
    return len(points)


# ----------------------------------------------------------------------------- new
def stage_new(args):
    modules = load_reference_modules()
    root = modules["cops_paths"].output_root()
    capture = pv.make_new_directory(
        Path(args.capture) if args.capture else root / "cuda_reference" / ("c00_" + pv.utc_stamp()))
    for forbidden in (reference.bind(), APP_ROOT):
        if forbidden == capture or forbidden in capture.parents:
            raise SystemExit("Error: capture directory must be outside the source trees")
    pv.write_json_exclusive(capture / "capture.json", {
        "schema": CAPTURE_SCHEMA, "schema_version": CAPTURE_VERSION,
        "created_utc": pv.utc_now(), "reference_root": str(reference.bind()),
        "application_root": str(APP_ROOT), "output_root": str(root)})
    print(capture)


# ----------------------------------------------------------------------------- identities
def running_campaign_processes(paths_module):
    """Map campaign process command lines to their resolved --out products."""
    found = []
    ticks = os.sysconf("SC_CLK_TCK")
    boot_time = next(int(line.split()[1]) for line in Path("/proc/stat").read_text().splitlines()
                     if line.startswith("btime "))
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            argv = (entry / "cmdline").read_bytes().split(b"\0")
            # Field 22 of /proc/<pid>/stat is the start time in clock ticks
            # after boot; the command name (field 2) may contain spaces.
            fields = (entry / "stat").read_text().rsplit(")", 1)[1].split()
            started = boot_time + int(fields[19]) / ticks
        except (OSError, IndexError, ValueError):
            continue
        argv = [item.decode("utf-8", "replace") for item in argv if item]
        tool = next((Path(item).name for item in argv if Path(item).name in (
            "fit_cops.py", "build_physical_activity_templates.py",
            "run_fit_profile.py", "liouville_cops")), None)
        if tool is None:
            continue
        out = None
        if "--out" in argv and argv.index("--out") + 1 < len(argv):
            out = str(paths_module.output_path(argv[argv.index("--out") + 1]))
        found.append({"pid": int(entry.name), "tool": tool, "out": out,
                      "argv": argv, "start_unix": started})
    return found


def campaign_inventory(modules, processes, reference_root):
    """Classify every campaign product; running products are never references."""
    root = modules["cops_paths"].output_root()
    running = {}
    for process in processes:
        if process["out"]:
            running.setdefault(process["out"], []).append(process["pid"])
    items = []
    for path in sorted(root.iterdir()):
        entry = {"path": str(path), "running_pids": sorted(running.get(str(path), []))}
        if path.is_dir():
            for name in ("fit_manifest.json", "run_manifest.json"):
                manifest_path = path / name
                if manifest_path.is_file():
                    manifest = json.loads(manifest_path.read_text())
                    entry.update({"manifest": pv.file_identity(manifest_path),
                                  "kind": manifest.get("fit_kind") or manifest.get("run_kind"),
                                  "states": manifest.get("states"),
                                  "command": manifest.get("command"),
                                  "parameters": manifest.get("parameters"),
                                  "optimizer_settings": manifest.get("optimizer_settings")})
                    break
            states = entry.get("states") or {}
            if entry["running_pids"]:
                status = "running_not_a_reference"
            elif "manifest" not in entry:
                status = "not_a_campaign_product"
            elif states.get("process") != "complete":
                status = "incomplete_not_a_reference"
            elif entry.get("kind") in ("forward",):
                status = "complete_forward_base"
            elif states.get("optimizer") in ("converged",):
                status = "complete_fit_converged"
            else:
                status = "complete_fit_needs_inspection"
            entry["status"] = status
            parameters = entry.get("parameters")
            if isinstance(parameters, dict):
                entry["free_parameter_count"] = len(parameters.get("names", []) or [])
        else:
            entry["file"] = pv.file_identity(path)
            entry["status"] = ("rejected_product" if "rejected" in path.name
                               else "file_product")
        items.append(entry)
    # Reference sources edited after a running job started mean that job used
    # different code than the current tree; flag it rather than assume.
    newest_source = max((p.stat().st_mtime for p in (reference_root / "tools").glob("*.py")))
    for process in processes:
        process["reference_tools_modified_after_start"] = newest_source > process["start_unix"]
    return {"output_root": str(root), "products": items, "processes": processes}


def stage_identities(args):
    started = pv.utc_now()
    out = stage_directory(args.capture, "identities")
    modules = load_reference_modules()
    root = reference.bind()
    cops_io = modules["cops_io"]
    builder = modules["build_physical_activity_templates"]
    run = base_directory(modules, args.base)

    log("hashing reference and application trees")
    reference_tree = pv.tree_identity(root, reference.SOURCE_EXCLUDES)
    application_tree = pv.tree_identity(APP_ROOT)

    log("validating the base cache (all kernels; this reads ~10 GB)")
    t0 = time.monotonic()
    audit = cops_io.validate_fit_base(run)
    validation_seconds = time.monotonic() - t0
    log("hashing every kernel file individually")
    kernels, digest = [], hashlib.sha256()
    for point in audit["point_ids"]:
        path = run / "kernels" / ("point=%04d.bin" % point)
        sha = pv.sha256_file(path)
        digest.update(path.name.encode("ascii"))
        digest.update(sha.encode("ascii"))
        kernels.append({"point": point, "file": path.name,
                        "size_bytes": path.stat().st_size, "sha256": sha})
    if digest.hexdigest() != audit["kernels"]["set_sha256"]:
        raise SystemExit("Error: per-file kernel digests disagree with the validated set digest")
    pv.write_json_exclusive(out / "kernel_files.json", {
        "base": str(run), "set_sha256": audit["kernels"]["set_sha256"],
        "set_digest_definition": "sha256 over (filename, hex sha256) in table point order",
        "kernels": kernels})

    manifest = audit["manifest"] or {}
    dependencies = []
    for entry in manifest.get("dependencies", []):
        recorded = entry.get("file", {})
        try:
            actual = pv.file_identity(recorded["path"])
            status = ("verified" if actual["sha256"] == recorded.get("sha256") and
                      actual["size_bytes"] == recorded.get("size_bytes") else "MISMATCH")
        except OSError as exc:
            actual, status = None, "unreadable: %s" % exc
        dependencies.append({"role": entry.get("role"), "recorded": recorded,
                             "status": status})
    bad = [d for d in dependencies if d["status"] != "verified"]

    log("evaluating SPICE omega with the reference helper")
    epochs = [args.late_start, builder.utc_midpoint(args.late_start, args.late_stop),
              args.late_stop]
    omega_executable = root / "spice_angular_velocity"
    omega, variation, command = builder.spice_omega(
        omega_executable, run / "input.snapshot.in", epochs)

    import importlib.util
    import numpy
    import scipy
    gravity_spec = importlib.util.find_spec("polyhedral_gravity")
    gravity_binary = (pv.file_identity(gravity_spec.origin)
                      if gravity_spec and gravity_spec.origin else None)
    wheel = sorted((APP_ROOT / "third_party" / "download_polyhedral_gravity_v3.3").glob(
        "polyhedral_gravity-3.3-cp38-*.whl"))
    wheel_so_sha = None
    if wheel:
        import zipfile
        with zipfile.ZipFile(str(wheel[0])) as archive:
            for name in archive.namelist():
                if name.endswith(".so"):
                    wheel_so_sha = hashlib.sha256(archive.read(name)).hexdigest()

    processes = running_campaign_processes(modules["cops_paths"])
    inventory = campaign_inventory(modules, processes, root)
    pv.write_json_exclusive(out / "campaign_inventory.json", inventory)

    identities = {
        "reference_tree": reference_tree,
        "reference_git": pv.git_state(root),
        "application_tree": application_tree,
        "application_git": pv.git_state(APP_ROOT),
        "reference_executables": {
            name: pv.file_identity(root / name) for name in (
                "liouville_cops", "spice_angular_velocity") if (root / name).is_file()},
        "interpreter": {"python": sys.version, "executable": sys.executable,
                        "numpy": numpy.__version__, "scipy": scipy.__version__},
        "polyhedral_gravity": {
            "binary": gravity_binary,
            "pypi_cp38_wheel": pv.file_identity(wheel[0]) if wheel else None,
            "binary_equals_pypi_wheel": (gravity_binary is not None and
                                         gravity_binary["sha256"] == wheel_so_sha),
            "pinned_requirement": (root / "requirements-gravity.txt").read_text()},
        "base": {
            "directory": str(run), "legacy_manifest": audit["legacy_manifest"],
            "validation_seconds": validation_seconds,
            "model_table": audit["model_table"], "elements": audit["elements"],
            "kernels": audit["kernels"], "point_count": len(audit["point_ids"]),
            "run_manifest": pv.file_identity(run / "run_manifest.json"),
            "input_snapshot": pv.file_identity(run / "input.snapshot.in"),
            "rebuild_record": (pv.file_identity(run / "rebuild.json")
                               if (run / "rebuild.json").is_file() else None),
            "manifest_states": manifest.get("states"),
            "effective_settings": manifest.get("effective_settings"),
            "dependencies": dependencies},
        "spice_omega": {"epochs_utc": epochs, "median_rad_s": omega.tolist(),
                        "max_sample_deviation_rad_s": variation,
                        "command": [str(item) for item in command]},
        "campaign_settings": {"late_interval_utc": [args.late_start, args.late_stop],
                              "density_kg_m3": DEFAULT_DENSITY,
                              "gravity_offset_m": DEFAULT_OFFSET,
                              "cliff_thresholds_deg": [45.0, 60.0]},
    }
    pv.write_json_exclusive(out / "identities.json", identities)
    summary = {"dependency_failures": len(bad), "kernels": len(kernels),
               "kernel_records": audit["kernels"]["records"],
               "reference_tree_sha256": reference_tree["tree_sha256"],
               "omega_rad_s": omega.tolist(),
               "running_processes": [(p["pid"], p["tool"], p["out"]) for p in processes],
               "product_status": {item["path"]: item["status"] for item in inventory["products"]}}
    finish_stage(out, "identities", summary, started)
    if bad:
        raise SystemExit("Error: %d manifest dependencies failed verification" % len(bad))


# ----------------------------------------------------------------------------- masks
def stage_masks(args):
    started = pv.utc_now()
    out = stage_directory(args.capture, "masks")
    modules = load_reference_modules()
    fit_cops = modules["fit_cops"]
    builder = modules["build_physical_activity_templates"]
    parse_utc = modules["cops_selection"].parse_utc
    run = base_directory(modules, args.base)
    header, names, rows, columns = read_base_table(modules, run)
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    nfaces = len(elements["area"])
    summary = {}
    for label, (start, stop) in (("late", (args.late_start, args.late_stop)),
                                 ("control", CONTROL_INTERVAL)):
        log("mask %s %s..%s" % (label, start, stop))
        mask, ng_response, rg_response, meta = builder.late_gauge_template(
            str(run), header, names, rows, nfaces, start, stop)
        points = meta["selected_point_ids"]
        verify_used_kernels(args.capture, run, points)
        # Per-point predicate sets, in CSR form.  The GPU must reproduce each
        # set and then their union; the union is checked against the
        # reference function's own result before anything is saved.
        offsets, faces, ng_counts, rg_counts = [0], [], [], []
        union = np.zeros(nfaces, dtype=bool)
        for point in points:
            saved, _, records = fit_cops.read_kernel(
                str(run / "kernels" / ("point=%04d.bin" % point)))
            assert saved == point
            ng = records["Kn"] > 0.0
            rg = (records["KPhi"] != 0.0) | (records["KV"] != 0.0)
            selected = np.unique(records["face"][ng & ~rg]).astype(np.int32)
            union[selected] = True
            faces.append(selected)
            offsets.append(offsets[-1] + len(selected))
            ng_counts.append(int(np.count_nonzero(ng)))
            rg_counts.append(int(np.count_nonzero(rg)))
        if not np.array_equal(union, mask.astype(bool)):
            raise SystemExit("Error: per-point predicate union differs from the reference mask")
        save_npz(out / ("mask_%s.npz" % label),
                 mask=mask.astype(np.uint8), ng_response=ng_response,
                 rg_response=rg_response, point_ids=np.asarray(points, dtype=np.int32),
                 point_offsets=np.asarray(offsets, dtype=np.int64),
                 point_faces=(np.concatenate(faces) if faces else np.empty(0, np.int32)),
                 point_ng_records=np.asarray(ng_counts, dtype=np.int64),
                 point_rg_records=np.asarray(rg_counts, dtype=np.int64))
        meta = dict(meta, interval_seconds=parse_utc(stop) - parse_utc(start),
                    mask_sha256=pv.array_sha256(mask.astype(np.uint8)),
                    area_fraction=float(np.sum(elements["area"][mask > 0]) /
                                        np.sum(elements["area"])))
        pv.write_json_exclusive(out / ("mask_%s.json" % label), meta)
        summary[label] = {key: meta[key] for key in (
            "selected_point_count", "ng_not_rg_union_face_count",
            "never_rg_over_interval_face_count", "area_fraction", "mask_sha256")}
    finish_stage(out, "masks", summary, started)


# ----------------------------------------------------------------------------- fit
def fit_selection(modules, run, profile_name):
    """Apply a named profile to the FULL table (profiles assert full row counts)."""
    fit_cops, selection = modules["fit_cops"], modules["cops_selection"]
    root = reference.bind()
    header, names, rows, columns = read_base_table(modules, run)
    _, intervals = fit_cops.run_intervals(header, str(root), str(run))
    hk = modules["cops_paths"].product_input("data/COPS_HK_20160930.txt")
    saturated = fit_cops.ng_saturated(hk, [row[columns["UTC"]] for row in rows])
    profile, profile_identity = selection.load_profile(
        profile_name, selection.profile_directory(str(root)))
    rows, report = selection.apply_selection_profile(rows, columns, profile,
                                                     intervals, saturated)
    return header, names, rows, columns, profile_identity, report, saturated


def choose_subset(rows, columns, count, late, parse_utc):
    """Deterministic production subset: evenly spaced fit rows plus late rows.

    Evenly spaced rows span altitude, illumination and both gauges; up to
    eight late-interval rows add the NG-only geometry.  The returned list is
    the production-subset index list, saved verbatim.
    """
    use = [k for k, row in enumerate(rows)
           if int(row[columns["useNG"]]) or int(row[columns["useRG"]])]
    spaced = [use[int(i)] for i in np.unique(np.linspace(0, len(use) - 1, count).astype(int))]
    t0, t1 = parse_utc(late[0]), parse_utc(late[1])
    late_rows = [k for k in use if t0 <= parse_utc(rows[k][columns["UTC"]]) <= t1]
    late_pick = [late_rows[int(i)] for i in np.unique(
        np.linspace(0, len(late_rows) - 1, min(8, len(late_rows))).astype(int))] if late_rows else []
    return sorted(set(spaced) | set(late_pick))


def synthetic_templates(centroid):
    """Two deterministic nonnegative binary fields used only to test the evaluator."""
    radius = np.linalg.norm(centroid, axis=1)
    longitude = np.degrees(np.arctan2(centroid[:, 1], centroid[:, 0]))
    cap = (centroid[:, 2] > 0.5 * radius).astype(float)
    wedge = ((longitude >= 0.0) & (longitude < 60.0)).astype(float)
    definitions = {"syn_north_cap": "z > 0.5*|r| (body frame)",
                   "syn_lon_wedge": "0 <= atan2(y,x) < 60 deg"}
    return ["syn_north_cap", "syn_lon_wedge"], np.column_stack((cap, wedge)), definitions


def make_problem(fit_cops, M, elements, coef, config, templates):
    args = argparse.Namespace(
        degree=config["degree"], enhancement_degree=config["enhancement_degree"],
        fix_legacy_activity=config["fix_legacy_activity"],
        fix_temperature=config["fix_temperature"], fix_kappa=config["fix_kappa"],
        ordered_temperature=config["ordered_temperature"], tbounds=(30.0, 400.0),
        enhancement_scale=None, max_nfev=1000, ftol=1e-10, xtol=1e-10, gtol=1e-10,
        jobs=1)
    if config["templates"]:
        names, values, _ = templates
        args.activity_template_names = list(names)
        args.activity_template_sources = [{"name": n, "synthetic": True} for n in names]
        args.activity_template_basis = np.ascontiguousarray(values)
    return fit_cops.Problem(M, elements, coef, args)


FIT_CONFIGS = (
    {"name": "legacy_d4_free_T_kappa", "degree": 4, "enhancement_degree": 4,
     "fix_legacy_activity": False, "fix_temperature": False, "fix_kappa": False,
     "ordered_temperature": False, "templates": False},
    {"name": "legacy_d3_ordered_T", "degree": 3, "enhancement_degree": 4,
     "fix_legacy_activity": False, "fix_temperature": False, "fix_kappa": False,
     "ordered_temperature": True, "templates": False},
    {"name": "legacy_d2_fixed_T", "degree": 2, "enhancement_degree": 4,
     "fix_legacy_activity": False, "fix_temperature": True, "fix_kappa": False,
     "ordered_temperature": False, "templates": False},
    {"name": "enhancement_l4_control", "degree": 4, "enhancement_degree": 4,
     "fix_legacy_activity": True, "fix_temperature": True, "fix_kappa": True,
     "ordered_temperature": False, "templates": False},
    {"name": "enhancement_l6_harmonic", "degree": 4, "enhancement_degree": 6,
     "fix_legacy_activity": True, "fix_temperature": True, "fix_kappa": True,
     "ordered_temperature": False, "templates": False},
    {"name": "enhancement_l5_templates", "degree": 4, "enhancement_degree": 5,
     "fix_legacy_activity": True, "fix_temperature": True, "fix_kappa": True,
     "ordered_temperature": False, "templates": True},
)


def parameter_vectors(P, config, seed):
    """Prior plus seeded vectors that exercise clipping, bounds and T branches."""
    rng = np.random.default_rng(seed)
    x0 = P.x0.copy()
    vectors = [("prior", x0.copy())]
    if not len(x0):
        return vectors
    names = list(P.names)

    def with_temperature(x, values):
        x = x.copy()
        if "Tmax" in names:
            for key, value in zip(("Tmax", "T75", "Tmin"), values):
                x[names.index(key)] = value
        elif "Tmin" in names:          # ordered: (Tmin, T75-Tmin, Tmax-T75)
            tmax, t75, tmin = values
            x[names.index("Tmin")] = tmin
            x[names.index("dT75")] = max(t75 - tmin, 0.0)
            x[names.index("dTmax")] = max(tmax - t75, 0.0)
        return x

    if P.na:
        mild = x0.copy()
        mild[:P.na] += rng.normal(0.0, 0.05, P.na)
        strong = x0.copy()
        strong[:P.na] += rng.normal(0.0, 0.4, P.na)
        vectors += [("legacy_mild", mild), ("legacy_strong_clip", strong)]
        if not P.args.fix_temperature:
            # The lit law falls from Tmax at SZA=0 to exactly T75 at SZA=75 deg,
            # so the Tmin floor can only bind for lit facets when T75 < Tmin.
            # Tmax=150, T75=140, Tmin=145 splits lit contributions between both
            # branches.  The ordered parameterization enforces T75 >= Tmin and
            # therefore cannot reach that split; it gets a T75 == Tmin edge.
            if P.ordered:
                vectors.append(("temperature_ordered_edge",
                                with_temperature(strong, (150.0, 140.0, 140.0))))
            else:
                vectors.append(("temperature_floor_mix",
                                with_temperature(strong, (150.0, 140.0, 145.0))))
            vectors.append(("temperature_wide", with_temperature(mild, (220.0, 130.0, 60.0))))
        if "kappa" in names:
            last = vectors[-1][1].copy()
            last[names.index("kappa")] = 0.37
            vectors.append(("kappa_scaled", last))
    if P.nb:
        nh = P.nh
        for label, sigma in (("enhancement_mild", 0.02), ("enhancement_strong_clip", 0.4)):
            x = x0.copy()
            x[P.bslice.start:P.bslice.start + nh] = rng.normal(0.0, sigma, nh)
            if P.nt:
                x[P.bslice.start + nh:P.bslice.stop] = np.abs(rng.normal(0.0, sigma, P.nt))
            vectors.append((label, x))
        if P.nt:
            x = vectors[-1][1].copy()
            x[P.bslice.start + nh] = 0.0     # one template amplitude exactly on its bound
            vectors.append(("template_on_bound", x))
    for label, x in vectors:
        if np.any(x < P.lo) or np.any(x > P.hi):
            raise SystemExit("Error: generated vector %s violates the bounds" % label)
    return vectors


def contribution_terms(M, act_lin, T):
    """The five per-contribution products of Model.evaluate, in its exact order."""
    fit_cops = sys.modules["fit_cops"]
    Tmax, T75, Tmin = T
    af = np.maximum(act_lin, 0.0)[M.face]
    Tl = Tmax + fit_cops.KT * (Tmax - T75) * M.xlit
    use_lit = M.lit & (Tl > Tmin)
    u = 1.0 / np.sqrt(np.where(use_lit, Tl, Tmin))
    return {"ng_h2o": M.wn["H2O"] * af * u, "rg_phi_h2o": M.wp["H2O"] * af,
            "rg_v_h2o": M.wv["H2O"] * af * u, "ng_co2": M.wn["CO2"] * u,
            "rg_v_co2": M.wv["CO2"] * u}, use_lit, Tl


def exact_segment_sums(values, indptr):
    """Correctly rounded per-point sums (math.fsum) and sums of magnitudes."""
    exact = np.empty(len(indptr) - 1)
    magnitude = np.empty(len(indptr) - 1)
    for k in range(len(indptr) - 1):
        segment = values[indptr[k]:indptr[k + 1]]
        exact[k] = math.fsum(segment.tolist())
        magnitude[k] = math.fsum(np.abs(segment).tolist())
    return exact, magnitude


def exact_pressures(M, terms, kappa):
    """Pressures assembled from correctly rounded sums of the reference terms."""
    sums = {name: exact_segment_sums(values, M.indptr) for name, values in terms.items()}
    phi_co2 = exact_segment_sums(M.Q["CO2"] * M.ACO2[M.face] * M.KPhi, M.indptr)
    cNG, cRG, bm = M.cNG, M.cRG, M.bm
    pNG = (cNG["H2O"] * bm["H2O"] * sums["ng_h2o"][0] +
           kappa * cNG["CO2"] * bm["CO2"] * sums["ng_co2"][0])
    pRG = (cRG["H2O"] * (sums["rg_phi_h2o"][0] + bm["H2O"] * sums["rg_v_h2o"][0]) +
           kappa * cRG["CO2"] * (phi_co2[0] + bm["CO2"] * sums["rg_v_co2"][0]))
    # Sum of magnitudes of every contribution entering each pressure: the
    # scale for a summation-error bound that is independent of cancellation.
    mNG = (abs(cNG["H2O"] * bm["H2O"]) * sums["ng_h2o"][1] +
           abs(kappa * cNG["CO2"] * bm["CO2"]) * sums["ng_co2"][1])
    mRG = (np.abs(cRG["H2O"]) * (sums["rg_phi_h2o"][1] + bm["H2O"] * sums["rg_v_h2o"][1]) +
           np.abs(kappa * cRG["CO2"]) * (phi_co2[1] + bm["CO2"] * sums["rg_v_co2"][1]))
    return pNG, pRG, mNG, mRG


def branch_signature(P, x):
    """Activity active set on contributing faces and the lit-temperature branch."""
    fit_cops = sys.modules["fit_cops"]
    a, b, T, _ = P.unpack(x)
    act = P.M.Y @ a
    if P.nb:
        act = act + P.enhancement_scale * (P.enhancement_basis @ b)
    Tmax, T75, Tmin = T
    Tl = Tmax + fit_cops.KT * (Tmax - T75) * P.M.xlit
    return (act[P.M.face] > 0.0), (P.M.lit & (Tl > Tmin))


def central_differences(P, x, step=1e-6):
    """Central differences of the data residuals with branch-crossing flags."""
    base_active, base_lit = branch_signature(P, x)
    columns, crossed, steps = [], [], []
    for k in range(len(x)):
        h = step * max(1.0, abs(x[k]))
        plus, minus = x.copy(), x.copy()
        plus[k] += h
        minus[k] -= h
        # Respect bounds: one-sided at an active bound, as the optimizer would.
        if minus[k] < P.lo[k]:
            minus[k] = x[k]
        if plus[k] > P.hi[k]:
            plus[k] = x[k]
        width = plus[k] - minus[k]
        columns.append((P.data_residuals(plus) - P.data_residuals(minus)) / width)
        flags = []
        for probe in (plus, minus):
            active, lit = branch_signature(P, probe)
            flags.append(bool(np.any(active != base_active) or np.any(lit != base_lit)))
        crossed.append(any(flags))
        steps.append(width)
    return np.column_stack(columns), np.asarray(crossed), np.asarray(steps)


def stage_fit(args):
    started = pv.utc_now()
    out = stage_directory(args.capture, "fit")
    modules = load_reference_modules()
    fit_cops = modules["fit_cops"]
    parse_utc = modules["cops_selection"].parse_utc
    root = reference.bind()
    run = base_directory(modules, args.base)
    header, names, rows, columns, profile_identity, report, saturated = fit_selection(
        modules, run, args.profile)
    subset = choose_subset(rows, columns, args.points, (args.late_start, args.late_stop),
                           parse_utc)
    subset_rows = [rows[k] for k in subset]
    point_ids = [int(row[columns["iPoint"]]) for row in subset_rows]
    verify_used_kernels(args.capture, run, point_ids)
    log("subset: %d rows of %d (profile %s)" % (len(subset), len(rows), args.profile))

    coef = fit_cops.read_activity_coefficients(str(root), header)
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    M = fit_cops.Model(str(run), header, columns, subset_rows, elements, coef)
    M.load_kernels(M.sel)
    # C++ replay agreement on the subset (the reference fit performs the same
    # check on all rows before fitting).
    pNG0, pRG0 = M.evaluate(coef[0], M.T0, 1.0)
    ok = (M.pNG_run[M.rows] > 0) & (M.pRG_run[M.rows] > 0)
    replay = {"NG": float(np.max(np.abs(pNG0[ok] / M.pNG_run[M.rows][ok] - 1.0))),
              "RG": float(np.max(np.abs(pRG0[ok] / M.pRG_run[M.rows][ok] - 1.0)))}
    log("C++ replay agreement on subset: %s" % replay)

    templates = synthetic_templates(elements["centroid"])
    save_npz(out / "subset.npz", subset_row_index=np.asarray(subset, dtype=np.int64),
             point_ids=np.asarray(point_ids, dtype=np.int32),
             use_ng=np.asarray([int(r[columns["useNG"]]) for r in subset_rows], np.int8),
             use_rg=np.asarray([int(r[columns["useRG"]]) for r in subset_rows], np.int8),
             ref_ng=np.asarray([float(r[columns["refNG"]]) for r in subset_rows]),
             ref_rg=np.asarray([float(r[columns["refRG"]]) for r in subset_rows]),
             saturated_ng=np.asarray([bool(saturated[k]) for k in subset]),
             contribution_count=np.diff(M.indptr),
             synthetic_templates=templates[1].astype(np.uint8),
             replay_prior_png=pNG0, replay_prior_prg=pRG0)

    summary = {"profile": profile_identity, "selection_counts": report["counts"],
               "subset_rows": len(subset), "subset_point_ids": point_ids,
               "contributions": int(M.indptr[-1]), "cpp_replay_max_rel": replay,
               "synthetic_templates": templates[2], "configs": {}}
    fd_targets = {"legacy_d4_free_T_kappa": ("legacy_strong_clip", "temperature_floor_mix"),
                  "legacy_d3_ordered_T": ("temperature_ordered_edge",),
                  "enhancement_l6_harmonic": ("enhancement_strong_clip",),
                  "enhancement_l5_templates": ("template_on_bound",)}
    for number, config in enumerate(FIT_CONFIGS):
        log("config %s" % config["name"])
        P = make_problem(fit_cops, M, elements, coef, config, templates)
        vectors = parameter_vectors(P, config, seed=1000 + number)
        record = {"names": list(P.names), "x0": P.x0.tolist(), "lo": P.lo.tolist(),
                  "hi": P.hi.tolist(), "na": P.na, "nh": P.nh, "nt": P.nt,
                  "enhancement_scale": P.enhancement_scale, "N0": P.N0,
                  "ascale": float(P.ascale), "NNG": int(P.NNG), "NRG": int(P.NRG),
                  "lambda": args.lam, "vectors": {}}
        arrays = {"Lt": P.Lt, "use_ng_mask": P.useNG, "use_rg_mask": P.useRG}
        for label, x in vectors:
            a, b, T, kappa = P.unpack(x)
            pNG, pRG = M.evaluate(a, T, kappa, b=b if P.nb else None,
                                  enhancement_basis=P.enhancement_basis,
                                  enhancement_scale=P.enhancement_scale)
            P.lam = 0.0
            data_res = P.data_residuals(x)
            data_jac = P.data_jacobian(x)
            P.lam = args.lam
            full_res = P.residuals(x)
            full_jac = P.jacobian(x)
            act_lin, _ = P.activity_linear_and_derivative(x, len(x))
            terms, use_lit, Tl = contribution_terms(M, act_lin, T)
            xpNG, xpRG, mNG, mRG = exact_pressures(M, terms, kappa)
            seen = np.unique(M.face)
            coverage = {
                "contributing_faces": int(len(seen)),
                "contributing_faces_clipped": int(np.count_nonzero(act_lin[seen] <= 0.0)),
                "contributing_faces_exact_zero": int(np.count_nonzero(act_lin[seen] == 0.0)),
                "contributions": int(len(M.face)),
                "night_contributions": int(np.count_nonzero(M.night)),
                "lit_le75_contributions": int(np.count_nonzero(M.lit)),
                "lit_branch_contributions": int(np.count_nonzero(use_lit)),
                "temperature_floor_lit_contributions": int(np.count_nonzero(M.lit & ~use_lit)),
                "rg_velocity_negative_terms": int(np.count_nonzero(terms["rg_v_h2o"] < 0.0)),
                "pressure_floor_hits": int(np.count_nonzero(pNG <= fit_cops.FLOOR) +
                                           np.count_nonzero(pRG <= fit_cops.FLOOR)),
                "parameters_on_bound": [P.names[k] for k in range(len(x))
                                        if x[k] == P.lo[k] or x[k] == P.hi[k]],
                "reference_summation_error_rel_max": {
                    "NG": float(np.max(np.abs(pNG - xpNG) / np.maximum(mNG, 1e-300))),
                    "RG": float(np.max(np.abs(pRG - xpRG) / np.maximum(mRG, 1e-300)))},
                "rg_cancellation_ratio_max": float(np.max(mRG / np.maximum(np.abs(xpRG), 1e-300))),
            }
            key = "%s__%s" % (config["name"], label)
            arrays.update({key + "__x": x, key + "__png": pNG, key + "__prg": pRG,
                           key + "__png_exact": xpNG, key + "__prg_exact": xpRG,
                           key + "__png_magnitude": mNG, key + "__prg_magnitude": mRG,
                           key + "__data_residuals": data_res, key + "__data_jacobian": data_jac,
                           key + "__residuals": full_res, key + "__jacobian": full_jac})
            if label in fd_targets.get(config["name"], ()):
                log("  central differences for %s" % label)
                P.lam = 0.0
                fd, crossed, steps = central_differences(P, x)
                denominator = np.maximum(np.max(np.abs(data_jac), axis=0), 1e-300)
                error = np.max(np.abs(fd - data_jac), axis=0) / denominator
                arrays.update({key + "__fd_jacobian": fd, key + "__fd_branch_crossed": crossed,
                               key + "__fd_steps": steps})
                coverage["fd_max_rel_error_smooth_columns"] = float(
                    np.max(error[~crossed])) if np.any(~crossed) else None
                coverage["fd_branch_crossed_columns"] = [P.names[k] for k in np.nonzero(crossed)[0]]
            record["vectors"][label] = coverage
        save_npz(out / ("fit_%s.npz" % config["name"]), **arrays)
        summary["configs"][config["name"]] = record
    pv.write_json_exclusive(out / "fit_fixtures.json", summary)
    finish_stage(out, "fit", {k: summary[k] for k in (
        "subset_rows", "contributions", "cpp_replay_max_rel")}, started)


# ----------------------------------------------------------------------------- joint (C00b)
# ROADMAP revision 2, C00b: re-freeze the evaluator reference against the
# repaired CPU tree for the primary joint families.  This stage writes a new
# subdirectory of an existing capture; its stage manifest records the CURRENT
# reference-tree digest, so the historical 'fit' stage (older digest) and this
# versioned set coexist.  Kernel files are checked against the capture's frozen
# per-file digests (they are unchanged inputs).
ARCHIVE_FIT = ("/home/vtenishe/ROS-END-OF-MISSION/"
               "RES-COPS-101826-6pm-7459fd3c350733bb38c87799df62a3cdcbdf9b34/"
               "output_fit_full_current16")
JOINT_CONFIGS = (
    # label, templates (names resolved in stage_joint), role
    ("B29", (), "primary"),
    ("N30", ("late_ng_not_rg",), "primary"),
    ("C30", ("syn_cliff_standin",), "primary (synthetic stand-in until C05)"),
    ("NC31", ("late_ng_not_rg", "syn_cliff_standin"), "primary (synthetic stand-in until C05)"),
    ("R30", ("late_rg_not_ng",), "labelled extension (user decision 2026-10-19)"),
)
JOINT_FD = {"B29": ("legacy_strong_clip", "archived_b29_state"),
            "NC31": ("enhancement_strong_clip", "template_on_bound"),
            "R30": ("enhancement_mild",)}


def archived_b29(path):
    """Fitted values of the archived 29-parameter fit (fit_result.txt table)."""
    values = {}
    pattern = re.compile(r"^(a\d\d|Tmax|T75|Tmin|kappa)\s+(\S+)\s+(\S+)\s+\S+\s*$")
    for line in Path(path, "fit_result.txt").read_text().splitlines():
        match = pattern.match(line)
        if match:
            values[match.group(1)] = (float(match.group(2)), float(match.group(3)))
    expected = ["a%02d" % j for j in range(25)] + ["Tmax", "T75", "Tmin", "kappa"]
    if sorted(values) != sorted(expected):
        raise SystemExit("Error: archived fit_result.txt lacks parameters: %s" %
                         sorted(set(expected) - set(values)))
    return values


def rg_only_mask(fit_cops, run, point_ids, nfaces):
    """CPU replica of late_rg_not_ng: union over points of (KPhi!=0 | KV!=0) & !(Kn>0).

    Uses the builder's own per-point predicates (point_rg, point_ng) in the
    opposite direction; this is the C05 test oracle for the GPU mask.
    """
    mask = np.zeros(nfaces, dtype=bool)
    per_point = []
    for point in point_ids:
        saved, _, records = fit_cops.read_kernel(str(run / "kernels" / ("point=%04d.bin" % point)))
        if saved != point:
            raise SystemExit("Error: kernel point mismatch %d" % point)
        point_ng = records["Kn"] > 0.0
        point_rg = (records["KPhi"] != 0.0) | (records["KV"] != 0.0)
        faces = records["face"][point_rg & ~point_ng]
        mask[faces] = True
        per_point.append(np.sort(faces))
    return mask, per_point


def stage_joint(args):
    started = pv.utc_now()
    out = stage_directory(args.capture, "joint")
    modules = load_reference_modules()
    fit_cops = modules["fit_cops"]
    parse_utc = modules["cops_selection"].parse_utc
    root = reference.bind()
    run = base_directory(modules, args.base)
    header, names, rows, columns, profile_identity, report, saturated = fit_selection(
        modules, run, args.profile)
    subset = choose_subset(rows, columns, args.points, (args.late_start, args.late_stop),
                           parse_utc)
    subset_rows = [rows[k] for k in subset]
    point_ids = [int(row[columns["iPoint"]]) for row in subset_rows]
    verify_used_kernels(args.capture, run, point_ids)
    previous = Path(args.capture) / "fit" / "fit_fixtures.json"
    same_subset = (json.loads(previous.read_text())["subset_point_ids"] == point_ids
                   if previous.is_file() else None)
    log("subset: %d rows (identical to C00 fit subset: %s)" % (len(subset), same_subset))

    coef = fit_cops.read_activity_coefficients(str(root), header)
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    nfaces = len(elements["area"])
    M = fit_cops.Model(str(run), header, columns, subset_rows, elements, coef)
    M.load_kernels(M.sel)

    # Templates: the real NG-only mask (C00 masks stage, equal to the reference
    # builder's), the RG-only replica, and a synthetic cliff stand-in.
    late = np.load(str(Path(args.capture) / "masks" / "mask_late.npz"))
    ng_only = late["mask"] > 0.0
    late_points = [int(p) for p in late["point_ids"]]
    verify_used_kernels(args.capture, run, late_points)
    rg_only, rg_per_point = rg_only_mask(fit_cops, run, late_points, nfaces)
    syn_names, syn_values, syn_defs = synthetic_templates(elements["centroid"])
    library = {"late_ng_not_rg": ng_only.astype(float),
               "late_rg_not_ng": rg_only.astype(float),
               "syn_cliff_standin": syn_values[:, 0]}
    definitions = {
        "late_ng_not_rg": "C00 masks stage (equals reference builder): union of Kn>0 & KPhi==0 & KV==0",
        "late_rg_not_ng": "union over the same %d late points of (KPhi!=0 | KV!=0) & !(Kn>0)" % len(late_points),
        "syn_cliff_standin": "synthetic %s; stands in for cliff_ge_60deg until the C05 bundle" %
                             syn_defs["syn_north_cap"]}
    area = elements["area"]
    template_stats = {name: {"faces": int(np.count_nonzero(v)),
                             "area_fraction": float(np.sum(area[v > 0]) / np.sum(area))}
                      for name, v in library.items()}
    template_stats["overlap_ng_rg_faces"] = int(np.count_nonzero(ng_only & rg_only))
    log("templates: %s" % template_stats)
    archived = archived_b29(args.archive_fit)

    summary = {"stage": "C00b joint", "profile": profile_identity,
               "subset_point_ids": point_ids, "identical_to_c00_fit_subset": same_subset,
               "contributions": int(M.indptr[-1]), "templates": template_stats,
               "template_definitions": definitions,
               "archive_fit": str(args.archive_fit), "configs": {}}
    save_npz(out / "templates.npz", late_ng_not_rg=ng_only, late_rg_not_ng=rg_only,
             syn_cliff_standin=library["syn_cliff_standin"] > 0.0,
             late_point_ids=np.asarray(late_points, dtype=np.int32),
             rg_only_point_offsets=np.concatenate([[0], np.cumsum([len(f) for f in rg_per_point])]).astype(np.int64),
             rg_only_point_faces=(np.concatenate(rg_per_point) if rg_per_point
                                  else np.zeros(0)).astype(np.int32),
             subset_point_ids=np.asarray(point_ids, dtype=np.int32))

    for number, (label, template_names, role) in enumerate(JOINT_CONFIGS):
        log("config %s %s" % (label, template_names))
        config = {"name": label, "degree": 4, "enhancement_degree": 4,
                  "fix_legacy_activity": False, "fix_temperature": False, "fix_kappa": False,
                  "ordered_temperature": False, "templates": bool(template_names)}
        templates = (list(template_names),
                     np.column_stack([library[n] for n in template_names]) if template_names
                     else np.empty((nfaces, 0)), definitions)
        P = make_problem(fit_cops, M, elements, coef, config, templates)
        vectors = parameter_vectors(P, config, seed=2000 + number)
        # The archived 29-parameter optimum (Tmax=30 K, T75=400 K at bounds;
        # heavy clipping) as a production-like state; amplitudes stay zero.
        x = P.x0.copy()
        for k, name in enumerate(P.names):
            if name in archived:
                fitted = archived[name][1]
                x[k] = fitted / P.ascale if name.startswith("a") else fitted
        vectors.append(("archived_b29_state", np.clip(x, P.lo, P.hi)))
        record = {"role": role, "templates": list(template_names), "names": list(P.names),
                  "x0": P.x0.tolist(), "lo": P.lo.tolist(), "hi": P.hi.tolist(),
                  "na": P.na, "nh": P.nh, "nt": P.nt, "reg_slice": [P.reg_slice.start, P.reg_slice.stop],
                  "enhancement_scale": P.enhancement_scale, "N0": P.N0,
                  "ascale": float(P.ascale), "NNG": int(P.NNG), "NRG": int(P.NRG),
                  "lambda": args.lam, "vectors": {}}
        arrays = {"Lt": P.Lt, "use_ng_mask": P.useNG, "use_rg_mask": P.useRG,
                  "enhancement_basis_columns": P.enhancement_basis if P.nb else np.empty((0,))}
        for vlabel, x in vectors:
            a, b, T, kappa = P.unpack(x)
            pNG, pRG = M.evaluate(a, T, kappa, b=b if P.nb else None,
                                  enhancement_basis=P.enhancement_basis,
                                  enhancement_scale=P.enhancement_scale)
            P.lam = 0.0
            data_res, data_jac = P.data_residuals(x), P.data_jacobian(x)
            P.lam = args.lam
            full_res, full_jac = P.residuals(x), P.jacobian(x)
            act_lin, _ = P.activity_linear_and_derivative(x, len(x))
            terms, use_lit, _ = contribution_terms(M, act_lin, T)
            xpNG, xpRG, mNG, mRG = exact_pressures(M, terms, kappa)
            seen = np.unique(M.face)
            coverage = {
                "contributing_faces_clipped": int(np.count_nonzero(act_lin[seen] <= 0.0)),
                "contributing_faces": int(len(seen)),
                "lit_branch_contributions": int(np.count_nonzero(use_lit)),
                "temperature_floor_lit_contributions": int(np.count_nonzero(M.lit & ~use_lit)),
                "parameters_on_bound": [P.names[k] for k in range(len(x))
                                        if x[k] == P.lo[k] or x[k] == P.hi[k]],
                "reference_summation_error_rel_max": {
                    "NG": float(np.max(np.abs(pNG - xpNG) / np.maximum(mNG, 1e-300))),
                    "RG": float(np.max(np.abs(pRG - xpRG) / np.maximum(mRG, 1e-300)))},
                "cost_data": float(0.5 * np.sum(data_res ** 2))}
            key = "%s__%s" % (label, vlabel)
            arrays.update({key + "__x": x, key + "__png": pNG, key + "__prg": pRG,
                           key + "__png_exact": xpNG, key + "__prg_exact": xpRG,
                           key + "__png_magnitude": mNG, key + "__prg_magnitude": mRG,
                           key + "__data_residuals": data_res, key + "__data_jacobian": data_jac,
                           key + "__residuals": full_res, key + "__jacobian": full_jac})
            if vlabel in JOINT_FD.get(label, ()):
                log("  central differences for %s" % vlabel)
                P.lam = 0.0
                fd, crossed, steps = central_differences(P, x)
                denominator = np.maximum(np.max(np.abs(data_jac), axis=0), 1e-300)
                error = np.max(np.abs(fd - data_jac), axis=0) / denominator
                arrays.update({key + "__fd_jacobian": fd, key + "__fd_branch_crossed": crossed,
                               key + "__fd_steps": steps})
                coverage["fd_max_rel_error_smooth_columns"] = float(
                    np.max(error[~crossed])) if np.any(~crossed) else None
                coverage["fd_branch_crossed_columns"] = [P.names[k] for k in np.nonzero(crossed)[0]]
            record["vectors"][vlabel] = coverage
        save_npz(out / ("joint_%s.npz" % label), **arrays)
        summary["configs"][label] = record
    pv.write_json_exclusive(out / "joint_fixtures.json", summary)
    finish_stage(out, "joint", {k: summary[k] for k in (
        "contributions", "templates", "identical_to_c00_fit_subset")}, started)


# ----------------------------------------------------------------------------- gravity
def cube(side):
    h = 0.5 * side
    v = np.array([[x, y, z] for x in (-h, h) for y in (-h, h) for z in (-h, h)], dtype=float)
    # Outward (counter-clockwise seen from outside) triangulation.
    f = np.array([[0, 1, 3], [0, 3, 2], [4, 6, 7], [4, 7, 5], [0, 4, 5], [0, 5, 1],
                  [2, 3, 7], [2, 7, 6], [0, 2, 6], [0, 6, 4], [1, 5, 7], [1, 7, 3]],
                 dtype=np.int64)
    return orient_outward(v, f)


def tetrahedron(edge):
    s = edge / (2.0 * math.sqrt(2.0))
    v = s * np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
    f = np.array([[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]], dtype=np.int64)
    return orient_outward(v, f)


def icosphere(radius, level):
    t = (1.0 + math.sqrt(5.0)) / 2.0
    v = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t],
         [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]]
    f = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4],
         [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8],
         [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]]
    v = [np.array(p, dtype=float) / np.linalg.norm(p) for p in v]
    for _ in range(level):
        cache, faces = {}, []

        def midpoint(i, j):
            key = (min(i, j), max(i, j))
            if key not in cache:
                p = v[i] + v[j]
                v.append(p / np.linalg.norm(p))
                cache[key] = len(v) - 1
            return cache[key]
        for a, b, c in f:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            faces += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        f = faces
    return orient_outward(radius * np.asarray(v), np.asarray(f, dtype=np.int64))


def orient_outward(v, f):
    """Flip any face whose normal points toward the (convex) body's centroid."""
    f = f.copy()
    center = v.mean(axis=0)
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    outward = np.einsum("ij,ij->i", np.cross(b - a, c - a), (a + b + c) / 3.0 - center) > 0.0
    f[~outward] = f[~outward][:, [1, 0, 2]]
    return v, f


def rotation_matrix(alpha, beta, gamma):
    ca, sa, cb, sb, cg, sg = (math.cos(alpha), math.sin(alpha), math.cos(beta),
                              math.sin(beta), math.cos(gamma), math.sin(gamma))
    rz = np.array([[ca, -sa, 0], [sa, ca, 0], [0, 0, 1]])
    ry = np.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]])
    rx = np.array([[1, 0, 0], [0, cg, -sg], [0, sg, cg]])
    return rz @ ry @ rx


def probe_points(v, f):
    """Far, near-face/edge/vertex, interior and exactly-singular evaluation points."""
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    n = np.cross(b - a, c - a)
    n /= np.linalg.norm(n, axis=1)[:, None]
    centroid = (a + b + c) / 3.0
    scale = np.max(np.linalg.norm(v, axis=1))
    labels, points = [], []

    def add(label, point):
        labels.append(label)
        points.append(np.asarray(point, dtype=float))
    for d in (5.0, 50.0):
        for direction in ([1, 0, 0], [0.3, -0.5, 0.8], [-0.2, -0.9, -0.4]):
            u = np.asarray(direction, dtype=float)
            add("far_%gR" % d, d * scale * u / np.linalg.norm(u))
    for k in (0, len(f) // 2):
        for h in (1e-3, 0.1, 1.0):
            add("near_face_%d_%gm" % (k, h), centroid[k] + h * n[k])
        add("on_face_centroid_%d" % k, centroid[k])
        midpoint = 0.5 * (a[k] + b[k])
        add("on_edge_midpoint_%d" % k, midpoint)
        add("near_edge_%d_0.1m" % k, midpoint + 0.1 * n[k])
        add("on_vertex_%d" % k, a[k])
        add("near_vertex_%d_0.1m" % k, a[k] + 0.1 * a[k] / np.linalg.norm(a[k]))
    add("interior_origin", [0.0, 0.0, 0.0])
    add("interior_offset", 0.2 * scale * np.array([0.5, 1.0, -1.5]) / math.sqrt(3.5))
    return labels, np.asarray(points)


def evaluate_package(pg, v, f, density, points, integrity, parallel):
    polyhedron = pg.Polyhedron(polyhedral_source=(v, f), density=density,
                               normal_orientation=pg.NormalOrientation.OUTWARDS,
                               integrity_check=integrity, metric_unit=pg.MetricUnit.METER)
    evaluator = pg.GravityEvaluable(polyhedron=polyhedron)
    result = evaluator(computation_points=np.asarray(points).tolist(), parallel=parallel)
    potential = np.asarray([r[0] for r in result])
    acceleration = np.asarray([r[1] for r in result])
    tensor = np.asarray([r[2] for r in result])
    return potential, acceleration, tensor, polyhedron


def mesh_topology(v, f):
    """Closedness/orientation facts independent of the package's ray casting."""
    directed = np.concatenate((f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]))
    undirected = np.sort(directed, axis=1)
    _, und_counts = np.unique(undirected, axis=0, return_counts=True)
    _, dir_counts = np.unique(directed, axis=0, return_counts=True)
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    signed_volume = float(np.sum(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0)
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    return {"vertices": int(len(v)), "faces": int(len(f)),
            "edges": int(len(und_counts)),
            "edges_not_shared_by_exactly_two_faces": int(np.count_nonzero(und_counts != 2)),
            "directed_edges_repeated": int(np.count_nonzero(dir_counts != 1)),
            "consistently_oriented_closed_manifold": bool(
                np.all(und_counts == 2) and np.all(dir_counts == 1)),
            "euler_characteristic": int(len(v) - len(und_counts) + len(f)),
            "signed_volume_m3": signed_volume, "area_min_m2": float(area.min()),
            "zero_area_faces": int(np.count_nonzero(area == 0.0))}


def package_ray_flags(V, F, sample, eps=1e-14):
    """Exact NumPy replica of Polyhedron::countRayPolyhedronIntersections (v3.3).

    Same operation order as the C++ (no FMA): ray from centroid + eps*N along
    N, Moller-Trumbore against every face, odd number of *distinct*
    intersection points means "inward".  Returns the flags and how many of
    them come from the face hitting itself, which happens when eps is below
    one rounding unit of the coordinates (GRAVITY_REFERENCE.md section 5).
    """
    A, B, C = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    e1, e2 = B - A, C - A

    def cross(a, b):
        return np.stack((a[..., 1] * b[..., 2] - a[..., 2] * b[..., 1],
                         a[..., 2] * b[..., 0] - a[..., 0] * b[..., 2],
                         a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]), axis=-1)

    def dot(a, b):
        return ((0.0 + a[..., 0] * b[..., 0]) + a[..., 1] * b[..., 1]) + a[..., 2] * b[..., 2]
    flags, own = [], []
    for i in sample:
        a, b, c = A[i], B[i], C[i]
        n = cross(b - a, c - b)
        n = n / np.sqrt(dot(n, n))
        origin = ((a + b) + c) / 3.0 + n * eps
        h = cross(np.broadcast_to(n, e2.shape), e2)
        det = dot(e1, h)
        with np.errstate(divide="ignore", invalid="ignore"):
            f = 1.0 / det
            s = origin - A
            u = f * dot(s, h)
            q = cross(s, e1)
            v = f * dot(np.broadcast_to(n, q.shape), q)
            t = f * dot(e2, q)
        hit = (~((det > -eps) & (det < eps)) & ~((u < 0.0) | (u > 1.0)) &
               ~((v < 0.0) | (u + v > 1.0)) & (t > eps))
        distinct = {tuple(point) for point in origin + n * t[hit][:, None]}
        flags.append(len(distinct) % 2 != 0)
        own.append(bool(hit[i]))
    return np.asarray(flags), np.asarray(own)


def stage_gravity(args):
    started = pv.utc_now()
    if args.cpus:
        cpus = set()
        for part in args.cpus.split(","):
            lo, _, hi = part.partition("-")
            cpus.update(range(int(lo), int(hi or lo) + 1))
        # TBB sizes its pool from the affinity mask at first use, so this must
        # happen before polyhedral_gravity is imported.
        os.sched_setaffinity(0, cpus)
    out = stage_directory(args.capture, "gravity")
    modules = load_reference_modules()
    import polyhedral_gravity as pg
    fit_cops, builder = modules["fit_cops"], modules["build_physical_activity_templates"]
    gravity_cache = modules["gravity_cache"]
    package = {"version": pg.__version__, "parallelization": getattr(pg, "__parallelization__", None),
               "commit": getattr(pg, "__commit__", None), "binary": pv.file_identity(pg.__file__)}
    summary = {"package": package, "synthetic": {}, "cpu_affinity": sorted(os.sched_getaffinity(0))}

    # Synthetic shapes: exact-position singular cases, symmetry/equivariance,
    # density scaling and refinement convergence.  Orientation is proven by the
    # independent topology audit; the package VERIFY cannot be used because it
    # misjudges km-scale meshes (see package_orientation_check_vs_scale).
    shapes = [("cube_1km", cube(1000.0)), ("tetrahedron_1km", tetrahedron(1000.0))]
    shapes += [("icosphere_r1km_l%d" % level, icosphere(1000.0, level)) for level in range(5)]
    arrays = {}
    shift = np.array([1234.5, -987.25, 456.125])
    rotation = rotation_matrix(0.3, -0.7, 1.1)
    for name, (v, f) in shapes:
        topology = mesh_topology(v, f)
        if not (topology["consistently_oriented_closed_manifold"] and
                topology["signed_volume_m3"] > 0.0):
            raise SystemExit("Error: synthetic shape %s is not a closed outward manifold" % name)
        labels, points = probe_points(v, f)
        pot, acc, ten, _ = evaluate_package(pg, v, f, DEFAULT_DENSITY, points,
                                            pg.PolyhedronIntegrity.DISABLE, False)
        pot_t, acc_t, ten_t, _ = evaluate_package(pg, v + shift, f, DEFAULT_DENSITY,
                                                  points + shift,
                                                  pg.PolyhedronIntegrity.DISABLE, False)
        pot_r, acc_r, ten_r, _ = evaluate_package(pg, v @ rotation.T, f, DEFAULT_DENSITY,
                                                  points @ rotation.T,
                                                  pg.PolyhedronIntegrity.DISABLE, False)
        _, acc_2, _, _ = evaluate_package(pg, v, f, 2.0 * DEFAULT_DENSITY, points,
                                          pg.PolyhedronIntegrity.DISABLE, False)
        for suffix, value in (("vertices", v), ("faces", f), ("points", points),
                              ("potential", pot), ("acceleration", acc), ("tensor", ten),
                              ("acceleration_translated", acc_t),
                              ("acceleration_rotated", acc_r),
                              ("acceleration_double_density", acc_2)):
            arrays["%s__%s" % (name, suffix)] = value
        arrays["%s__labels" % name] = np.asarray(labels, dtype="U48")
        a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
        volume = float(np.sum(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0)
        far = [i for i, label in enumerate(labels) if label.startswith("far_50")]
        mass_gm = 6.67430e-11 * DEFAULT_DENSITY * volume
        far_rel = [abs(np.linalg.norm(acc[i]) * np.dot(points[i], points[i]) / mass_gm - 1.0)
                   for i in far]
        scale = np.max(np.linalg.norm(acc, axis=1))
        summary["synthetic"][name] = {
            "faces": int(len(f)), "points": len(labels), "volume_m3": volume,
            "translation_max_abs_diff_rel": float(np.max(np.abs(acc_t - acc)) / scale),
            "rotation_max_abs_diff_rel": float(np.max(np.abs(acc_r - acc @ rotation.T)) / scale),
            "density_scaling_max_rel": float(np.max(np.abs(acc_2 - 2.0 * acc)) / scale),
            "far_field_point_mass_rel_at_50R": far_rel,
            "nonfinite_values": int(np.count_nonzero(~np.isfinite(acc)))}
        if name.startswith("icosphere"):
            # Analytic sphere of the same radius: a convergence check only.
            inside = [i for i, label in enumerate(labels) if label == "interior_offset"][0]
            r = np.linalg.norm(points[inside])
            analytic = 4.0 / 3.0 * math.pi * 6.67430e-11 * DEFAULT_DENSITY * r
            summary["synthetic"][name]["interior_vs_sphere_rel"] = float(
                abs(np.linalg.norm(acc[inside]) / analytic - 1.0))
    # HEAL repair on a deliberately flipped face of the cube.
    v, f = cube(1000.0)
    broken = f.copy()
    broken[3] = broken[3][[1, 0, 2]]
    healed = pg.Polyhedron(polyhedral_source=(v, broken), density=DEFAULT_DENSITY,
                           normal_orientation=pg.NormalOrientation.OUTWARDS,
                           integrity_check=pg.PolyhedronIntegrity.HEAL,
                           metric_unit=pg.MetricUnit.METER)
    healed_faces = np.asarray(healed.faces, dtype=np.int64)
    summary["synthetic"]["heal_flipped_cube_face"] = {
        "flipped_index": 3, "healed_face": healed_faces[3].tolist(),
        "original_face": f[3].tolist(), "orientation": str(healed.normal_orientation),
        "healed_matches_outward_winding_cyclically": bool(
            any(np.array_equal(np.roll(healed_faces[3], s), f[3]) for s in range(3)))}
    scale_test = {}
    for radius in (1.0, 100.0, 1000.0, 2000.0):
        v, f = icosphere(radius, 3)
        probe = pg.Polyhedron(polyhedral_source=(v, f), density=DEFAULT_DENSITY,
                              normal_orientation=pg.NormalOrientation.OUTWARDS,
                              integrity_check=pg.PolyhedronIntegrity.DISABLE,
                              metric_unit=pg.MetricUnit.METER)
        orientation, violating = probe.check_normal_orientation()
        scale_test["icosphere_l3_r%gm" % radius] = {
            "independent_topology_outward": bool(
                mesh_topology(v, f)["consistently_oriented_closed_manifold"]),
            "package_majority_orientation": str(orientation),
            "package_flagged_faces": len(violating), "faces": int(len(f))}
    summary["synthetic"]["package_orientation_check_vs_scale"] = scale_test
    save_npz(out / "gravity_synthetic.npz", **arrays)
    log("synthetic gravity fixtures written")

    # Exact 788k gravity body: the reference builder's reader and options.
    run = base_directory(modules, args.base)
    identities = json.loads((Path(args.capture) / "identities" / "identities.json").read_text())
    shape_dep = [d for d in identities["base"]["dependencies"] if d["role"] == "shape_model"]
    shape_path = Path(args.gravity_shape or shape_dep[0]["recorded"]["path"])
    shape = pv.file_identity(shape_path)
    log("reading gravity body %s" % shape_path)
    V, F = builder.read_nastran_surface(shape_path)
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    centroid, normal, area = elements["centroid"], elements["normal"], elements["area"]
    topology = mesh_topology(V, F)
    # The transport facets (targets) and the gravity body (sources) are both
    # read from the same BDF in the current campaign; record whether they are
    # the same triangles in the same order rather than assume it.
    same_order = None
    if len(F) == len(area):
        same_order = float(np.max(np.abs(V[F].mean(axis=1) - centroid)))
    log("topology: %s" % topology)

    # The package's HEAL/VERIFY ray test is unreliable at this coordinate
    # scale (GRAVITY_REFERENCE.md section 5): it would flip about half of the
    # correctly oriented faces.  The reference therefore evaluates the input
    # winding with DISABLE, admitted only by the independent topology audit.
    # An exact replica of the package test on a sample documents the defect.
    if not (topology["consistently_oriented_closed_manifold"] and
            topology["signed_volume_m3"] > 0.0 and topology["zero_area_faces"] == 0 and
            topology["euler_characteristic"] == 2):
        raise SystemExit("Error: gravity body fails the topology audit; refusing DISABLE")
    log("replicating the package orientation test on %d sampled faces" % args.heal_sample)
    flag_sample = np.sort(np.random.default_rng(5).choice(len(F), args.heal_sample, replace=False))
    t0 = time.monotonic()
    flags, own = package_ray_flags(V, F, flag_sample)
    healed = F                       # the evaluated face table (input winding)
    heal = {"policy": "DISABLE after independent topology audit; package HEAL not used",
            "replica_sample_faces": int(len(flag_sample)),
            "replica_flagged_inward": int(np.count_nonzero(flags)),
            "replica_flagged_fraction": float(np.mean(flags)),
            "replica_flags_from_own_face_hit": int(np.count_nonzero(flags & own)),
            "replica_seconds": time.monotonic() - t0,
            "evaluated_faces_sha256": pv.array_sha256(F)}
    log("package orientation test replica: %s" % heal)
    polyhedron = pg.Polyhedron(polyhedral_source=(V, F), density=DEFAULT_DENSITY,
                               normal_orientation=pg.NormalOrientation.OUTWARDS,
                               integrity_check=pg.PolyhedronIntegrity.DISABLE,
                               metric_unit=pg.MetricUnit.METER)
    evaluator = pg.GravityEvaluable(polyhedron=polyhedron)

    n = len(area)
    sample = np.unique(np.linspace(0, n - 1, min(512, n), dtype=int))
    rng = np.random.default_rng(788)
    radius = np.linalg.norm(centroid, axis=1)
    order_area, order_radius = np.argsort(area, kind="stable"), np.argsort(radius, kind="stable")
    difficult = {"smallest_area": order_area[:64], "largest_area": order_area[-64:],
                 "smallest_radius": order_radius[:64], "largest_radius": order_radius[-64:],
                 "random": np.sort(rng.choice(n, 64, replace=False))}
    targets = np.unique(np.concatenate([sample] + list(difficult.values())))
    on_surface = targets[:: max(1, len(targets) // 32)][:32]

    def run_points(points):
        result = evaluator(computation_points=np.asarray(points).tolist(), parallel=True)
        return (np.asarray([r[0] for r in result]), np.asarray([r[1] for r in result]),
                np.asarray([r[2] for r in result]))
    log("evaluating %d targets at two offsets plus %d on-surface points" %
        (len(targets), len(on_surface)))
    t0 = time.monotonic()
    p_full, g_full, t_full = run_points(centroid[targets] + args.offset * normal[targets])
    eval_seconds = time.monotonic() - t0
    p_half, g_half, t_half = run_points(centroid[targets] + 0.5 * args.offset * normal[targets])
    p_on, g_on, t_on = run_points(centroid[on_surface])
    omega = np.asarray(identities["spice_omega"]["median_rad_s"])
    slope_full = gravity_cache.derive_effective_fields(
        centroid[targets], normal[targets], g_full, omega)[2]
    slope_half = gravity_cache.derive_effective_fields(
        centroid[targets], normal[targets], g_half, omega)[2]
    in_sample = np.isin(targets, sample)
    halving = np.abs(slope_full - slope_half)

    # Summation-order repeatability: identical faces in reverse order through
    # the same package path (DISABLE skips a second O(n^2) audit; the healed
    # orientation is already established).
    log("repeatability: reversed face order on 64 targets")
    repeat = targets[:: max(1, len(targets) // 64)][:64]
    reversed_poly = pg.Polyhedron(polyhedral_source=(V, healed[::-1].copy()),
                                  density=DEFAULT_DENSITY,
                                  normal_orientation=pg.NormalOrientation.OUTWARDS,
                                  integrity_check=pg.PolyhedronIntegrity.DISABLE,
                                  metric_unit=pg.MetricUnit.METER)
    reversed_eval = pg.GravityEvaluable(polyhedron=reversed_poly)
    rev = reversed_eval(computation_points=(centroid[repeat] + args.offset * normal[repeat]).tolist(),
                        parallel=True)
    g_rev = np.asarray([r[1] for r in rev])
    g_ref = g_full[np.searchsorted(targets, repeat)]
    rel = np.linalg.norm(g_rev - g_ref, axis=1) / np.linalg.norm(g_ref, axis=1)
    slope_rev = gravity_cache.derive_effective_fields(
        centroid[repeat], normal[repeat], g_rev, omega)[2]
    slope_ref = slope_full[np.searchsorted(targets, repeat)]

    save_npz(out / "gravity_788k.npz", targets=targets.astype(np.int64),
             in_offset_halving_sample=in_sample,
             points_offset=centroid[targets] + args.offset * normal[targets],
             points_half_offset=centroid[targets] + 0.5 * args.offset * normal[targets],
             potential_offset=p_full, acceleration_offset=g_full, tensor_offset=t_full,
             potential_half_offset=p_half, acceleration_half_offset=g_half,
             tensor_half_offset=t_half, on_surface_targets=on_surface.astype(np.int64),
             potential_on_surface=p_on, acceleration_on_surface=g_on, tensor_on_surface=t_on,
             omega_rad_s=omega, slope_offset_deg=slope_full, slope_half_offset_deg=slope_half,
             repeat_targets=repeat.astype(np.int64), acceleration_reversed_order=g_rev,
             evaluated_faces=healed, flag_sample=flag_sample, replica_flags=flags,
             **{"difficult_" + k: v.astype(np.int64)
                                     for k, v in difficult.items()})
    gravity = {
        "shape": shape, "topology": topology, "orientation": heal,
        "elements_match_bdf_triangle_order_max_centroid_diff_m": same_order,
        "elements_faces": int(n), "density_kg_m3": DEFAULT_DENSITY, "offset_m": args.offset,
        "targets": int(len(targets)), "evaluation_seconds_per_target_offset": eval_seconds / len(targets),
        "offset_halving_sample_count": int(np.count_nonzero(in_sample)),
        "offset_halving_slope_deg_max_sample": float(np.max(halving[in_sample])),
        "offset_halving_slope_deg_max_all_targets": float(np.max(halving)),
        "repeatability_reversed_order": {
            "targets": int(len(repeat)),
            "acceleration_rel_max": float(np.max(rel)),
            "acceleration_rel_median": float(np.median(rel)),
            "slope_deg_max": float(np.max(np.abs(slope_rev - slope_ref)))},
        "nonfinite_values": int(np.count_nonzero(~np.isfinite(g_full)) +
                                np.count_nonzero(~np.isfinite(g_on)))}
    summary["exact_788k"] = gravity
    pv.write_json_exclusive(out / "gravity_fixtures.json", summary)
    finish_stage(out, "gravity", gravity, started)


# ----------------------------------------------------------------------------- main
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stage", choices=("new", "identities", "masks", "fit", "joint", "gravity"))
    parser.add_argument("--capture", help="capture directory (from 'new'; for 'new': optional explicit path)")
    parser.add_argument("--reference", help="CPU LiouvilleCOPS tree (default: $%s or ../LiouvilleCOPS)" %
                        reference.ENVIRONMENT_VARIABLE)
    parser.add_argument("--base", default=DEFAULT_BASE, help="manifested fitting base")
    parser.add_argument("--late-start", default=LATE_INTERVAL[0])
    parser.add_argument("--late-stop", default=LATE_INTERVAL[1])
    parser.add_argument("--profile", default=DEFAULT_PROFILE, help="selection profile for 'fit'")
    parser.add_argument("--points", type=int, default=40, help="evenly spaced fit rows for 'fit'")
    parser.add_argument("--lambda", dest="lam", type=float, default=1e-4,
                        help="regularization weight for the augmented residual/Jacobian fixtures")
    parser.add_argument("--archive-fit", default=ARCHIVE_FIT,
                        help="archived 29-parameter fit directory for 'joint' (fit_result.txt)")
    parser.add_argument("--gravity-shape", help="override the base's shape_model as gravity body")
    parser.add_argument("--offset", type=float, default=DEFAULT_OFFSET)
    parser.add_argument("--cpus", help="CPU list for 'gravity', e.g. 30-35 (limits TBB threads)")
    parser.add_argument("--heal-sample", type=int, default=300,
                        help="faces for the exact replica of the package orientation test")
    args = parser.parse_args()
    reference.bind(args.reference)
    if args.stage != "new" and not args.capture:
        parser.error("--capture is required for stage %s" % args.stage)
    {"new": stage_new, "identities": stage_identities, "masks": stage_masks,
     "fit": stage_fit, "joint": stage_joint, "gravity": stage_gravity}[args.stage](args)


if __name__ == "__main__":
    main()
