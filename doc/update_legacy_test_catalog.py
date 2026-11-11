#!/usr/bin/env python3
"""Document legacy application tests without running or regenerating them.

The Table is a maintained declaration; Makefile.test is the host-filtered
effective recipe snapshot. Read both as text, reconcile every generated
compile/run/check owner with a real Table block, and preserve omitted records.
This tool never runs Config.pl, the table parser, make recipes, Git, or MPI.
Normal PDF compilation consumes the supplied local TeX/JSON snapshot and does
not require the surrounding source tree.
"""

import argparse
import collections
import hashlib
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def escape(value):
    replacements = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%",
                    "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}"}
    return "".join(replacements.get(c, c) for c in str(value))


def purpose(name):
    """Editorial descriptions of the named fixture, never a physical PASS.

    The group rules describe the configuration variant explicitly named by
    the maintained recipes; numerical observations are extracted separately.
    An unknown future generated owner requires an editorial description.
    """
    exact = {
        "Saturn-MOP": "Saturn MOP water-population regression; reduced H2O mesh output.",
        "Moon": "Lunar exosphere particle-moment baseline and restart producer.",
        "Moon-tetra-cut-cells": "Lunar tetrahedral cut-cell mesh and particle-output regression.",
        "Moon-restart": "Lunar continuation using particle restart files from test_Moon.",
        "Enceladus": "Legacy Enceladus application moment regression; distinct from the newer E3 forward campaign.",
        "Mercury-test": "Mercury application particle-output baseline.",
        "unit_test": "Legacy GoogleTest application; inspect assertions/exit and stage log, not an empty numeric diff alone.",
        "unit_test_field-lines": "Legacy GoogleTest field-line suite and summary-report comparison.",
        "unit_test_particle_distribution": "Legacy GoogleTest particle-distribution suite and summary-report comparison.",
        "unit_test_collisions": "Legacy GoogleTest collision suite and summary-report comparison.",
        "linear-system-solver": "Standalone linear-system solver diagnostic and reference-output regression.",
        "linear-system-solver-openmp": "OpenMP linear-system solver variant with its reference diagnostic.",
        "periodic-bc-gmres": "GMRES/periodic-boundary solver regression.",
        "Periodic-BC": "Uniform periodic-boundary wave/particle outputs and tracked-trajectory regression.",
        "Periodic-BC-NonUniform": "Nonuniform periodic-boundary wave/particle outputs and tracked trajectories.",
        "SEP3D-ARMS-coupling": "SEP application with ARMS background data; compares the configured SEP3D diagnostic.",
        "CG": "Comet CG gas-coma application baseline.",
        "CG-dust-test": "Comet gas/dust regression including water and dust sample outputs.",
        "CG-dust-test-tetra-cut-mesh": "Comet gas/dust regression on the tetrahedral cut-cell mesh.",
        "model-dust": "Legacy dust charging/size model application regression.",
        "CCMC-Individual_Trajectories": "CCMC individual charged-particle trajectory regression.",
        "Mover-gca-relativistic": "Relativistic guiding-centre mover application baseline.",
        "Europa-test1": "First maintained Europa source/particle regression configuration.",
        "Europa-test2": "Second maintained Europa source/particle regression configuration.",
        "Collision-Test-NTC": "No-time-counter collision fixture with multiple distribution/reference subchecks.",
        "Collision-Test-MF": "Majorant-frequency collision fixture with multiple distribution/reference subchecks.",
        "Chemistry-Test": "Chemical-solver application comparing two maintained kinetic reference outputs.",
        "Venus-Hot-Oxygen": "Venus hot-oxygen application particle-moment regression.",
        "InterpolateAMR-test": "AMR interpolation reference values through the maintained test adapter.",
        "InterpolateAMR-AMPS-test": "Native AMPS centre- and corner-based AMR interpolation against one maintained reference.",
        "cdap_test": "CDAP diagnostic application and its configured numeric baseline.",
        "const-field-bc": "Constant-field boundary application; compares its declared numerical boundary products.",
        "Geospace-Flux-Gridless": "Gridless geospace density/flux regression for directional, spatial and dipole fixtures.",
        "Geospace-Flux-Gridless-hybrid": "MPI/OpenMP variant of the gridless density/flux fixture set.",
        "Geospace-SWMF-AMPS": "Live SWMF geospace/PT setup followed by proton trajectory and volume-output comparisons.",
        "Earth-cutoff-rigidity-gridless-steadystate_t96": "Steady-state gridless cutoff-rigidity regression with T96 field.",
        "Earth-vertical-cutoff-rigidity-gridless": "Gridless vertical cutoff-rigidity baseline.",
        "Earth-vertical-cutoff-rigidity-gridless-hybrid": "MPI/OpenMP gridless vertical cutoff-rigidity baseline.",
        "Earth-cutoff-rigidity-gridless-trajectory": "Gridless cutoff-rigidity products along the declared trajectory.",
        "Earth-vertical-cutoff-rigidity-gridless-trajectory": "Gridless vertical cutoff along a spacecraft trajectory.",
        "Earth-vertical-cutoff-rigidity-gridless-trajectory_dipole_only": "Dipole-only gridless vertical cutoff trajectory fixture.",
        "Earth-cutoff-rigidity-3d": "Mesh-backed three-dimensional cutoff-rigidity application regression.",
        "Earth-cutoff-rigidity-3d-analytic": "Three-dimensional cutoff-rigidity analytic-fixture configuration.",
        "Earth-3d_forward-dipole_analytic": "Forward three-dimensional analytic-dipole fixture; trajectory and volume products.",
        "Earth-3d_forward-t96_analytic": "Forward three-dimensional T96/analytic-fixture configuration.",
        "Earth-3d_forward-t96": "Forward three-dimensional T96 trajectory/volume regression.",
    }
    if name in exact:
        return exact[name]
    if name.startswith("fast-wave"):
        variants = {
            "fast-wave": "baseline",
            "fast-wave--legacy-ecsim": "legacy ECSIM",
            "fast-wave__correction": "correction-enabled",
            "fast-wave_rebalancing": "load-rebalancing",
            "fast-wave_AVX256": "AVX256",
            "fast-wave_AVX512": "AVX512",
            "fast-wave_multithread_UpdateJMassMatrix": "threaded current/mass-matrix update",
            "fast-wave_multithread_mover_UpdateJMassMatrix": "threaded mover and current/mass-matrix update",
            "fast-wave_multithread_mover_matmul_UpdateJMassMatrix": "threaded mover/matrix multiply/current update",
            "fast-wave-openmp": "OpenMP",
            "fast-wave-openmp__shared_temp_list": "OpenMP shared temporary particle list",
        }
        return ("Electromagnetic fast-wave " + variants[name] +
                " regression, with proton/electron diagnostics.")
    if name.startswith("lightwave-current-"):
        return "Electromagnetic light-wave regression with configuration " + name[len("lightwave-"):] + "."
    if name.startswith("magnetic-reconnection"):
        return "Magnetic-reconnection diagnostic regression" + (
            " using AVX256." if name.endswith("avx256") else ".")
    for prefix in ("Individual_Trajectories--", "Mover-"):
        if name.startswith(prefix):
            return "Charged-particle " + name[len(prefix):] + " trajectory/mover fixture."
    if name.startswith("SEP--Parker_spiral--"):
        return "Prescribed Parker-spiral SEP transport; " + name.split("--")[-1] + " formulation."
    if name == "SEP--Parker3D--ParkerIMF":
        return "Three-dimensional Parker-background SEP transport regression."
    if name == "SEP--Parker3D--SWMF_IMF":
        return "Three-dimensional SEP transport with an SWMF-supplied IMF."
    if name.startswith("SEP--focused_transport-IH-"):
        return "SWMF inner-heliosphere focused-transport configuration " + name.split("transport-IH-", 1)[1] + "."
    if name.startswith("SEP-SC-IH--"):
        return "Solar-corona/inner-heliosphere SEP attachment configuration " + name[len("SEP-SC-IH--"):] + "."
    if name in ("SEP_SWMF_test19", "SEP_SWMF_test19_fl"):
        return "SWMF test19 SEP coupling regression" + (
            " with field-line configuration." if name.endswith("_fl") else ".")
    if name.startswith("AMPS-BATSRUS"):
        return "Live BATS-R-US/AMPS coupling configuration " + name[len("AMPS-BATSRUS"):].lstrip("-_") + "."
    raise ValueError("Add an editorial purpose for generated owner: " + name)


def family(name):
    if name.startswith(("SEP", "AMPS-BATSRUS")):
        return "Heliospheric transport and coupled plasma"
    if name.startswith(("Earth", "Geospace")):
        return "Geospace cutoff, access and flux"
    if name.startswith(("fast-wave", "lightwave", "magnetic-reconnection",
                        "linear-system", "periodic-bc", "Periodic-BC", "const-field")):
        return "Fields, waves, solvers and periodic boundaries"
    if name.startswith(("Individual_", "Mover-", "CCMC-Individual")):
        return "Individual trajectories and mover applications"
    if name.startswith(("unit_test", "Collision-", "Chemistry-", "Interpolate", "cdap")):
        return "Numerical and process fixtures"
    return "Planetary, comet and dust applications"


def read_inputs():
    table = (ROOT / "MakefileTest/Table").read_text()
    effective = (ROOT / "Makefile.test").read_text()
    records = []
    for match in re.finditer(r"^#>\s*\n(.*?)^<#\s*$", table, re.M | re.S):
        body = match.group(1)
        values = re.findall(r"^Name=(.+)$", body, re.M)
        if len(values) != 1:
            raise ValueError("A Table block must contain exactly one Name")
        raw = values[0].strip()
        application, _, selector = raw.partition("@@")
        application = application.strip()
        # Match the maintained generator's removal of the test/ prefix.
        target = "test_" + application.replace("test/", "")
        records.append({
            "application": application, "target": target,
            "selector": selector or "all machines",
            "table_line": table.count("\n", 0, match.start()) + 1,
            "body": body,
        })
    phase_pattern = re.compile(
        r"^test_([^\s:]+)_(compile|rundir|run|check):[^\n]*\n(.*?)"
        r"(?=^test_[^\s:]+:|\Z)", re.M | re.S)
    owners = collections.OrderedDict()
    for match in phase_pattern.finditer(effective):
        target = "test_" + match.group(1)
        owner = owners.setdefault(target, {"target": target, "phases": {}})
        owner["phases"][match.group(2)] = {
            "line": effective.count("\n", 0, match.start()) + 1,
            "recipe": match.group(3),
        }
    for target, owner in owners.items():
        if set(owner["phases"]) != {"compile", "rundir", "run", "check"}:
            raise ValueError("Incomplete generated phase registry: " + target)
        declarations = [x for x in records if x["target"] == target]
        if not declarations:
            raise ValueError("Generated owner absent from Table: " + target)
        name = target[len("test_"):]
        check = owner["phases"]["check"]["recipe"]
        combined = "\n".join(x["recipe"] for x in owner["phases"].values())
        key_match = re.search(r"^TEST" + re.escape(name) + r"KEYS=(.*)$",
                              effective, re.M)
        keys = key_match.group(1).strip() if key_match else ""
        requirements = []
        if re.search(r"\bSWMF[/; ]", combined):
            requirements.append("SWMF component tree")
        if any(x in (combined + keys).lower() for x in ("kameleon",)):
            requirements.append("Kameleon and its data")
        spice_values = re.findall(r"-spice-(?:path|kernels)=([^\s]+)", combined + " " + keys)
        if any(x != "nospice" for x in spice_values):
            requirements.append("configured SPICE toolkit/kernels")
        if "gtest" in combined:
            requirements.append("GoogleTest")
        if re.search(r"openmp|OMP_NUM_THREADS|OMPIRUN", combined + keys, re.I):
            requirements.append("declared OpenMP setup")
        if re.search(r"-avx[= ]", combined + keys):
            requirements.append("matching CPU/compiler vector target")
        if "-cplr-data-path=" in combined + keys:
            requirements.append("external background data")
        owner.update({
            "family": family(name), "purpose": purpose(name),
            "table": [{k: v for k, v in x.items() if k != "body"} for x in declarations],
            "diffnum_calls": check.count("DiffNum.pl"),
            "plain_diff_calls": len(re.findall(r"(?:^|\n)\s*[-@]*diff\s", check)),
            "explicit_tolerances": sorted(set(re.findall(r"-[ar]=[0-9.eE+-]+", check))),
            "git_instructions_present": bool(re.search(r"\bgit\s+(?:clone|checkout|pull|fetch|switch)\b", combined)),
            "native_preflight_helper": "prepare_legacy_native.py" in combined,
            "applications": list(dict.fromkeys(re.findall(r"-application=([^\s]+)", combined))),
            "phase_line_numbers": {k: v["line"] for k, v in owner["phases"].items()},
            "requirements_detected": requirements,
            "hardcoded_mpi_ranks": sorted(set(int(x) for x in re.findall(
                r"\b(?:mpirun|mpiexec)\s+-(?:np|n)\s+(\d+)",
                owner["phases"]["run"]["recipe"]))),
        })
        # Save only compact metadata in the document snapshot, not executable
        # shell recipes. The actual generated file remains the recipe authority.
        owner.pop("phases")
    grouped = collections.OrderedDict()
    for record in records:
        if record["target"] in owners:
            continue
        grouped.setdefault(record["target"], []).append(
            {k: v for k, v in record.items() if k != "body"})
    host = re.search(r"^HOSTNAME=(.*)$", effective, re.M)
    snapshot = {
        "schema_version": 1,
        "scope": "Registration/recipe documentation only; no tests executed",
        "generated_hostname": host.group(1).strip() if host else "unspecified",
        "table_records": len(records),
        "table_distinct_targets": len({x["target"] for x in records}),
        "generated_targets": len(owners), "omitted_targets": len(grouped),
        "source_sha256": {
            "MakefileTest/Table": hashlib.sha256(table.encode()).hexdigest(),
            "Makefile.test": hashlib.sha256(effective.encode()).hexdigest(),
        },
        "generated": list(owners.values()), "not_generated": list(grouped.values()),
    }
    return snapshot


def tex(snapshot):
    lines = [
        "% Generated by update_legacy_test_catalog.py; source-only editorial snapshot.",
        r"\chapter{Legacy application-regression catalog}\label{ch:legacy-catalog}",
        "{} Table records describe {} distinct targets. The current generated "
        "Makefile contains {} owners for host {}; {} targets are absent from "
        "this host snapshot. Presence is not execution, support or PASS.".format(
            snapshot["table_records"], snapshot["table_distinct_targets"],
            snapshot["generated_targets"], escape(snapshot["generated_hostname"]),
            snapshot["omitted_targets"]),
        r"Read Chapter~\ref{ch:application-regression} before executing a recipe.",
        "The local companion " + r"\path{legacy_test_catalog.json}" +
        " records all phase line numbers, selectors and source hashes.",
    ]
    groups = collections.OrderedDict()
    for owner in snapshot["generated"]:
        groups.setdefault(owner["family"], []).append(owner)
    for name, rows in groups.items():
        lines += [
            r"\section{" + escape(name) + "}",
            r"\begin{longtable}{L{0.40\textwidth}L{0.49\textwidth}}",
            r"\caption{Generated targets: " + escape(name) + r".}\\",
            r"\toprule Target and recipe location & Purpose and actual check contract\\\midrule\endfirsthead",
            r"\toprule Target and recipe location & Purpose and actual check contract\\\midrule\endhead",
        ]
        for owner in rows:
            left = r"\path{" + owner["target"] + r"}\newline "
            left += "Table lines " + ", ".join(str(x["table_line"]) for x in owner["table"])
            left += r".\newline "
            left += "Makefile.test compile/check lines {}/{}.".format(
                owner["phase_line_numbers"]["compile"], owner["phase_line_numbers"]["check"])
            right = escape(owner["purpose"]) + r"\newline "
            if owner["diffnum_calls"]:
                right += "{} DiffNum call(s); ".format(owner["diffnum_calls"])
                right += escape(", ".join(owner["explicit_tolerances"])) if owner["explicit_tolerances"] else "default comparator tolerances"
                right += ". "
            elif owner["plain_diff_calls"]:
                right += "{} plain diff call(s). ".format(owner["plain_diff_calls"])
            else:
                right += "No DiffNum/plain-diff invocation detected; inspect custom check and process result. "
            if owner["git_instructions_present"]:
                right += r"\textbf{Git instructions in recipe.} "
            if owner["native_preflight_helper"]:
                right += "Native preflight helper present. "
            if owner["requirements_detected"]:
                right += r"\newline " + escape("; ".join(owner["requirements_detected"])) + "."
            lines.append(left + " & " + right + r"\\")
        lines += [r"\bottomrule", r"\end{longtable}"]
    lines += [
        r"\section{Declared but absent from the generated host snapshot}",
        "These are not executed SKIPs: they are omitted targets. Filters are "
        "shown as declared, without pretending this editorial tool implements "
        "the complete production generator.",
        r"\begin{longtable}{L{0.48\textwidth}L{0.41\textwidth}}",
        r"\caption{Additional Table declarations and machine selectors.}\\",
        r"\toprule Target & Declaration and selector\\\midrule\endfirsthead",
        r"\toprule Target & Declaration and selector\\\midrule\endhead",
    ]
    for entries in snapshot["not_generated"]:
        left = r"\path{" + entries[0]["target"] + "}"
        right = "; ".join("Table line {}: {}".format(x["table_line"], x["selector"]) for x in entries)
        lines.append(left + " & " + escape(right) + r"\\")
    lines += [r"\bottomrule", r"\end{longtable}", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    snapshot = read_inputs()
    outputs = {
        HERE / "LegacyApplicationCatalog.tex": tex(snapshot),
        HERE / "legacy_test_catalog.json": json.dumps(snapshot, indent=2) + "\n",
    }
    for path, content in outputs.items():
        if args.check:
            if not path.is_file() or path.read_text() != content:
                sys.exit("ERROR: legacy application catalog is stale; run doc/update_legacy_test_catalog.py")
        else:
            path.write_text(content)
    print("PASS legacy catalog{}: {} generated / {} declared targets".format(
        " check" if args.check else " refresh",
        snapshot["generated_targets"], snapshot["table_distinct_targets"]))
