#!/usr/bin/env python3
"""Read-only preparation of the archived-physics/788k comparison campaign.

Print JSON containing exact commands and lightweight provenance checks. Never
launch a fit, read the multi-gigabyte transport kernels, or write a product.
The fitter remains responsible for its full cache and template validation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import sys

import numpy as np

from cops_io import read_model_table
from cops_selection import apply_selection_profile, load_profile, profile_directory
from gravity_cache import ALGORITHM
from gravity_topology import ORIENTATION_POLICY

APP = Path(__file__).resolve().parents[1]


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def settings_audit(archived, current):
    """Compare actual effective settings and dependency content identities.

    A changed pathname is not itself a changed HK table. Manifest hashes are
    compared by role, including multiplicity, so moved SPICE inputs are also
    handled without inventing filename equivalence.
    """
    left, right = archived["effective_settings"], current["effective_settings"]
    changes = {key: {"archived": left.get(key), "current": right.get(key)}
               for key in sorted(set(left) | set(right)) if left.get(key) != right.get(key)}
    roles = sorted({d["role"] for m in (archived, current) for d in m["dependencies"]})
    dependencies = {}
    for role in roles:
        a, b = [sorted(d["file"]["sha256"] for d in m["dependencies"] if d["role"] == role)
                for m in (archived, current)]
        dependencies[role] = {"same_recorded_content": bool(a and b and a == b),
                              "archived_count": len(a), "current_count": len(b)}
    return {"effective_setting_changes": changes, "dependency_roles": dependencies,
            "unexpected_setting_changes": sorted(set(changes) - {
                "OutputDirectory", "ShapeModel", "TemperatureRGTable"})}


def compare_selections(archived_table, current_table):
    """Compare selected observations by point ID *and* canonical UTC.

    The archived fit table gives the actual used pressures, not its base mask
    (which excludes late NG). The current all_available policy is evaluated by
    the shared selector. Saturation is included, so an all-false flag vector
    changes no selection decision; it is not a saturation count estimate.
    """
    from cops_selection import parse_utc
    _, _, old, oc = read_model_table(archived_table, ("UTC", "pNG_used_Pa", "pRG_used_Pa"))
    _, _, new, nc = read_model_table(current_table, ("UTC", "inside", "refNG", "refRG"))
    profile, _ = load_profile("all_available", profile_directory(APP))
    new, _ = apply_selection_profile(new, nc, profile, {"FitInterval": []}, [False] * len(new))
    def key(row, columns):
        return (int(row[columns["iPoint"]]), parse_utc(row[columns["UTC"]]))
    result = {}
    for gauge in ("NG", "RG"):
        a = {key(r, oc): float(r[oc["p" + gauge + "_used_Pa"]]) for r in old
             if float(r[oc["p" + gauge + "_used_Pa"]]) > 0}
        b = {key(r, nc): float(r[nc["ref" + gauge]]) for r in new if int(r[nc["use" + gauge]])}
        common = sorted(set(a) & set(b))
        result[gauge] = {"archived_selected": len(a), "current_selected": len(b),
                         "common_selected": len(common),
                         "common_ids": [k[0] for k in common],
                         "current_only_ids": [k[0] for k in sorted(set(b) - set(a))],
                         "archived_only_ids": [k[0] for k in sorted(set(a) - set(b))],
                         "changed_common_pressure_ids": [k[0] for k in common if a[k] != b[k]]}
    return result


def template_precheck(path, names):
    """Reject missing/old cliff provenance before emitting an eligible plan.

    This inspects only small NPZ metadata members. It deliberately does not
    assert that the per-face arrays match the active mesh; load_bundle and the
    shared selection gate will do that when the real fitter is launched.
    """
    if path is None or not Path(path).is_file():
        return ["template bundle missing; regenerate/re-export a validated bundle"]
    with np.load(str(path), allow_pickle=False) as archive:
        available = archive["names"].tolist()
        metadata = json.loads(str(archive["metadata_json"].item()))
    issues = ["missing template " + name for name in names if name not in available]
    # NG visibility was computed from transport kernels, independently of
    # gravity. Selecting only that column is permitted from the quarantined
    # bundle; cliff columns require corrected gravity provenance.
    if any(name.startswith("cliff_") for name in names):
        gravity = metadata.get("gravity", {})
        definition = gravity.get("cache", {}).get("definition", {})
        if (definition.get("algorithm") != ALGORITHM or
                definition.get("orientation_policy") != ORIENTATION_POLICY or
                gravity.get("integrity_check") != "DISABLE after topology audit"):
            issues.append("cliff provenance predates audited DISABLE gravity; HEAL products are ineligible")
    return issues


def fit_argv(run, output, templates=(), bundle=None):
    """Use the archived 29-parameter scheme plus zero, one, or two masks.

    Absence of --fix-* and --ordered-temperature is intentional: all 25
    degree-0..4 coefficients, three independent temperatures and kappa vary.
    One lambda gives one optimizer task; 36 workers would not speed this fit.
    """
    command = [sys.executable, str(APP / "tools/run_fit_profile.py"),
               "--run", str(run), "--profile", "all_available", "--degree", "4",
               "--enhancement-degree", "4", "--lambda", "0", "--jobs", "1",
               "--threads-per-process", "1", "--tbounds", "30", "400",
               "--max-nfev", "1000", "--ftol", "1e-10", "--xtol", "1e-10",
               "--gtol", "1e-10", "--check-tolerance", "1e-6", "--out", str(output)]
    for name in templates:
        command.extend(["--activity-template", str(bundle), name])
    return command


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-fit", type=Path, required=True)
    parser.add_argument("--archive-base", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--template-bundle", type=Path)
    parser.add_argument("--ng-bundle", type=Path,
                        help="optional independent NG-only bundle (may be the quarantined mixed bundle; no cliff column is selected)")
    args = parser.parse_args(argv)
    try:
        old = read_json(args.archive_base / "run_manifest.json")
        new = read_json(args.run / "run_manifest.json")
        audit = settings_audit(old, new)
        # FitInterval is absent in the archived full-mission scheme. Refuse to
        # label a restricted run equivalent rather than silently ignoring it.
        if new["effective_settings"].get("FitInterval"):
            raise ValueError("campaign requires the archived full FitInterval")
        selection = compare_selections(args.archive_fit / "fit_data_used.txt", args.run / "COPS_model.dat")
        dependency_checks = []
        for dep in new["dependencies"]:
            identity = dep["file"]
            path = Path(identity["path"])
            check = {"role": dep["role"], "path": str(path), "exists": path.is_file()}
            # Bound the preparatory audit's I/O; the authoritative base audit
            # hashes all dependencies and transport records on actual launch.
            if path.is_file() and path.stat().st_size <= 1024 * 1024:
                check["current_small_file_hash_matches"] = hashlib.sha256(path.read_bytes()).hexdigest() == identity["sha256"]
            dependency_checks.append(check)
        shared_blockers = []
        if audit["unexpected_setting_changes"]:
            shared_blockers.append("review unexpected effective-setting differences")
        if any(not value["same_recorded_content"] for role, value in audit["dependency_roles"].items() if role != "shape_model"):
            shared_blockers.append("non-shape dependency identities differ")
        if any(not d["exists"] or d.get("current_small_file_hash_matches") is False for d in dependency_checks):
            shared_blockers.append("a current dependency is missing or has changed")
        if any(value["changed_common_pressure_ids"] for value in selection.values()):
            shared_blockers.append("common pressure observations differ")
        plans = []
        for label, names in (("B", []), ("N", ["late_ng_not_rg"]),
                             ("C", ["cliff_ge_60deg"]),
                             ("NC", ["late_ng_not_rg", "cliff_ge_60deg"])):
            output = args.output_root.resolve() / ("fit_shap5_788k_restored_" + label)
            blockers = list(shared_blockers)
            if output.exists():
                blockers.append("output already exists")
            bundle = args.ng_bundle if label == "N" and args.ng_bundle else args.template_bundle
            if names:
                blockers.extend(template_precheck(bundle, names))
            command = fit_argv(args.run.resolve(), output, names, bundle)
            plans.append({"label": label, "free_parameters": 29 + len(names),
                          "precheck_blockers": blockers, "argv": command,
                          "shell_command": "env COPS_OUTPUT_ROOT=" + shlex.quote(str(args.output_root.resolve())) + " " + shlex.join(command)})
        print(json.dumps({"schema": "liouville_cops.restored_campaign_preparation_v1",
                          "scope": "read-only precheck, not full cache validation or execution approval",
                          "audit": audit, "dependency_checks": dependency_checks,
                          "selection_comparison": selection, "plans": plans,
                          "cpp_replay_threads": 36,
                          "launch_note": "Execute only after shared fitter tests pass; use distinct outputs for additional starts. Full cache/template validation remains mandatory."}, indent=2))
    except (OSError, ValueError, KeyError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
