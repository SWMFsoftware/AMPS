#!/usr/bin/env python3
"""Verify source-explicit reduced kinematics using independent mathematics.

No AMPS configuration, PIC execution, observation, physical source activation
or qualification flag is produced. The adapter supplies analytical fixtures
and random normals; every product solution is the maintained production code.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[2]
DEPENDENCIES = ("srcMoon/species/ReducedKinematics.h", "srcMoon/species/ReducedKinematics.cpp",
    "srcMoon/species/ReactionExecution.h", "srcMoon/test/reduced_kinematics/probe.cpp",
    "srcMoon/test/reduced_kinematics/CONTRACT.md", "srcMoon/test/reduced_kinematics/README.md",
    "srcMoon/test/test_m3_reduced_kinematics.py")
TRIALS = 40000
SEED = 731993
STATISTICS = 17
FAMILY_ALPHA = .01
ARITHMETIC_MULTIPLIER = 128
DECIMAL_PRECISION = 400


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hex_decimal(text):
    """Exact native binary input, without first converting through float64."""
    mantissa, exponent = text.lstrip("+-").lower().split("p")
    integer, _, fraction = mantissa[2:].partition(".")
    value = Decimal(int(integer+fraction, 16))*Decimal(2)**(int(exponent)-4*len(fraction))
    return -value if text.startswith("-") else value


def norm(v): return sum((x*x for x in v), Decimal(0)).sqrt()


def two_body(masses, total_momentum, available, direction):
    """Independent reduced-mass analytical oracle, not a production call."""
    total = sum(masses)
    reduced = masses[0]*masses[1]/total
    relative = (2*reduced*available).sqrt()
    p = [[m/total*total_momentum[d]+sign*relative*direction[d] for d in range(3)]
         for m, sign in zip(masses, (1, -1))]
    return p, [[x/m for x in row] for row, m in zip(p, masses)]


def jacobi_three(masses, normals, energy):
    """Independent orthogonal Jacobi basis: pair products1/2, spectator0.

    The production implementation instead projects sqrt(m)*g momenta and
    selects a dependent heavy coordinate. This basis diagonalizes the two
    relative modes, providing a separate algebraic construction of the same
    selected microcanonical measure.
    """
    m0, m1, m2 = masses; pair = m1+m2; total = m0+pair
    z1 = [(m2/pair).sqrt()*normals[1][d]-(m1/pair).sqrt()*normals[2][d] for d in range(3)]
    z2 = [(m1*m0/(pair*total)).sqrt()*normals[1][d]+
          (m2*m0/(pair*total)).sqrt()*normals[2][d]-(pair/total).sqrt()*normals[0][d] for d in range(3)]
    scale = (2*energy/(norm(z1)**2+norm(z2)**2)).sqrt()
    k = [(m1*m2/pair).sqrt()*x*scale for x in z1]
    q = [(pair*m0/total).sqrt()*x*scale for x in z2]
    p = [[-x for x in q], [m1/pair*q[d]+k[d] for d in range(3)],
         [m2/pair*q[d]-k[d] for d in range(3)]]
    return p, [[x/m for x in row] for row, m in zip(p, masses)]


def beta_cdf(x):
    """Closed-form Beta(3/2,3/2) CDF; no fitted histogram/reference output."""
    return 2/math.pi*(math.asin(math.sqrt(x))-math.sqrt(x*(1-x))*(1-2*x))


def run(output, compiler="g++"):
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("output must be a new or empty evidence directory: "+str(output))
    output.mkdir(parents=True, exist_ok=True)
    checks = []; comparisons = []; ensemble_rows = []
    result = dict(schema="srcMoon-reduced-kinematics-controls-v1", status="ERROR", results=[],
        counts={s: 0 for s in ("PASS", "FAIL", "ERROR", "SKIPPED")},
        scope="ANALYTICAL SOFTWARE CONTROLS; approved reduced models remain EXPERIMENTAL, not physical M3 qualification",
        physical_provider_activation=False, formal_I21_I24_status="NOT EVALUATED")
    manifest = dict(schema="srcMoon-reduced-kinematics-manifest-v1", source_revision=subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(), source_hashes={p: sha(REPO/p) for p in DEPENDENCIES},
        external_data=[], MPI="none", OMP="none", random_seed=SEED,
        fixture_notice="Convenient analytical kg/J controls, never isotope data or observations",
        controls=dict(trials=TRIALS, simultaneous_statistics=STATISTICS, family_alpha=FAMILY_ALPHA,
            arithmetic_multiplier=ARITHMETIC_MULTIPLIER, decimal_precision=DECIMAL_PRECISION), commands=[])
    def check(name, passed, **details):
        checks.append(dict(name=name, status="PASS" if passed else "FAIL", **details))
    try:
        cxx = shutil.which(compiler)
        if not cxx:
            checks.append(dict(name="compiler-availability", status="SKIPPED", reason="C++17 compiler unavailable"))
            result.update(status="SKIPPED", reason="compiler unavailable: "+compiler)
            return result
        version = subprocess.run([cxx, "--version"], capture_output=True, text=True, timeout=15)
        (output/"compiler-version.log").write_text(version.stdout+version.stderr)
        if version.returncode: raise RuntimeError("compiler identity query failed")
        manifest["compiler"] = dict(path=cxx, version=version.stdout)
        executable = output/"reduced-kinematics-probe"
        command = [cxx, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror", "-pedantic", "-I", str(REPO),
            str(REPO/"srcMoon/species/ReducedKinematics.cpp"), str(REPO/"srcMoon/test/reduced_kinematics/probe.cpp"), "-o", str(executable)]
        with (output/"compile.log").open("w") as log:
            built = subprocess.run(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT, timeout=120)
        manifest["commands"].append(dict(argv=command, returncode=built.returncode))
        if built.returncode: raise RuntimeError("production kernel compilation failed; see compile.log")
        manifest["executable"] = dict(path=str(executable), sha256=sha(executable))
        completed = subprocess.run([str(executable)], capture_output=True, text=True, timeout=180)
        (output/"stdout.log").write_text(completed.stdout); (output/"stderr.log").write_text(completed.stderr)
        manifest["commands"].append(dict(argv=[str(executable)], returncode=completed.returncode))
        rows = [(status, name) for status, name in re.findall(r"^(PASS|FAIL) (.+)$", completed.stdout, re.M)]
        summary = re.findall(r"^SUMMARY checks=(\d+) failed=(\d+) seed=(\d+) trials=(\d+)$", completed.stdout, re.M)
        if len(summary) != 1 or not rows or len({name for _, name in rows}) != len(rows):
            raise RuntimeError("missing/duplicate native control records")
        total, failed, seed, trials = map(int, summary[0])
        if total != len(rows) or failed != sum(status == "FAIL" for status, _ in rows) or seed != SEED or trials != TRIALS or completed.returncode != int(failed > 0):
            raise RuntimeError("native counts/exit/seed/trials contradict actual controls")
        checks.extend(dict(name="native/"+name, status=status) for status, name in rows)
        events = {}; products = {}; platform = None; model_inputs = None; model_errors = None; ensemble = None
        for line in completed.stdout.splitlines():
            fields = line.split("\t")
            if fields[0] == "PLATFORM": platform = (Decimal(fields[1]), int(fields[2]), int(fields[3]))
            elif fields[0] == "EVENT":
                if fields[1] in events: raise RuntimeError("duplicate event receipt")
                events[fields[1]] = list(map(Decimal, fields[2:]))
            elif fields[0] == "PRODUCT":
                name, index = fields[1], int(fields[2])
                if index in products.setdefault(name, {}): raise RuntimeError("duplicate product occurrence")
                products[name][index] = list(map(Decimal, fields[3:]))
            elif fields[0] == "MODELINPUT": model_inputs = fields[1:]
            elif fields[0] == "ERRORS": model_errors = list(map(Decimal, fields[1:]))
            elif fields[0] == "ENSEMBLE": ensemble = fields[1:]
        if platform is None or model_inputs is None or model_errors is None or ensemble is None:
            raise RuntimeError("missing platform/model/ensemble receipt")
        manifest["long_double"] = dict(epsilon=str(platform[0]), digits=platform[1], maximum_exponent=platform[2])
        with localcontext() as context:
            context.prec = DECIMAL_PRECISION; context.Emin = -999999; context.Emax = 999999
            eps = platform[0]; relative = ARITHMETIC_MULTIPLIER*eps
            def compare(name, actual, reference, scale=None):
                actual, reference = Decimal(actual), Decimal(reference)
                denominator = abs(reference) if scale is None else abs(Decimal(scale))
                limit = relative*denominator
                passed = actual.is_finite() and abs(actual-reference) <= limit
                check(name, passed, actual=str(actual), reference=str(reference), absolute_bound=str(limit))
                comparisons.append(dict(metric=name, actual=str(actual), reference=str(reference), absolute_bound=str(limit),
                    relative_error=str(abs(actual-reference)/denominator) if denominator else "0" if passed else "Infinity"))
            physical, speed = map(hex_decimal, model_inputs)
            c = Decimal(299792458); beta = speed/c; gamma = 1/(1-beta*beta).sqrt()
            expected_errors = [(physical-1)/physical, physical-1, ((gamma-1)-beta*beta/2)/(gamma-1)]
            for name, actual, reference in zip(("mass-inertia", "mass-acceleration", "nonrelativistic-kinetic"), model_errors, expected_errors):
                compare("model-error/"+name, actual, reference)
            direction = [(Decimal(3)/8).sqrt(), (Decimal(3)/8).sqrt(), -Decimal(1)/2]
            references = {}
            masses = [Decimal(1), Decimal(2)]; momentum = list(map(Decimal, (1, -2, 3)))
            references["general-two-body"] = (masses, *two_body(masses, momentum, Decimal(4), list(map(Decimal, (0, 0, 1)))))
            for name, masses, threshold, mean in (("neutral-photon", [Decimal(1), Decimal(1)], Decimal(3), Decimal(2)),
                    ("atomic-ion-electron-photon", [Decimal(".999"), Decimal(".001")], Decimal(5), Decimal(2)),
                    ("tiny-open-polar", [Decimal(1), Decimal(1)], Decimal(3), Decimal(2)),
                    ("separate-small-excess-preserved", [Decimal(1), Decimal(1)], Decimal("1e-100"), Decimal("1e-200"))):
                photon = threshold+mean; momentum = [photon/c, Decimal(0), Decimal(0)]
                available = mean-(photon/c)**2/(2*sum(masses))
                selected_direction = direction
                if name == "tiny-open-polar":
                    u = Decimal("1e-25"); transverse = (2*u*(1-u)).sqrt()
                    selected_direction = [transverse, transverse, 2*u-1]
                references[name] = (masses, *two_body(masses, momentum, available, selected_direction))
                if name not in events: check(name+"/event-receipt", False); continue
                compare(name+"/available-COM-energy", events[name][0], available)
                compare(name+"/photon-energy", events[name][3], photon)
                compare(name+"/state-Q", events[name][4], -threshold)
                compare(name+"/source-threshold", events[name][5], threshold)
                compare(name+"/source-mean", events[name][6], mean)
                compare(name+"/source-available", events[name][7], mean)
            normals = [[Decimal(x) for x in values] for values in [(1, 2, -1), (2, -1, ".5"), ("-.2", ".7", -2)]]
            # Solve for *velocities* from a prescribed relative speed in each
            # mode, independently of the native momentum construction. The
            # electron and the ion each use their own direction in the same
            # axes. A nonzero boost tests that neither local axis is applied
            # to a velocity in a different frame.
            charged_masses = [Decimal(".001"), Decimal(".999"), Decimal(2)]
            electron_mass, ion_mass, neutral_mass = charged_masses
            pair_mass = ion_mass + neutral_mass
            total_mass = electron_mass + pair_mass
            for name, fraction, kinetic, p in (
                    ("charged-rest", ".25", 4, (0, 0, 0)),
                    ("charged-drift", ".25", 4, (3, -2, 1)),
                    ("charged-heavy-only", "0", 4, (3, -2, 1)),
                    ("charged-electron-only", "1", 4, (3, -2, 1)),
                    ("charged-at-threshold", "1", 0, (3, -2, 1)),
                    ("charged-photon", ".25", 4-(12/c)**2/(2*total_mass), (12/c, 0, 0))):
                electron_k = Decimal(fraction) * kinetic
                heavy_k = kinetic - electron_k
                e_speed = (2 * electron_k * (1/electron_mass + 1/pair_mass)).sqrt()
                heavy_speed = (2 * heavy_k * (1/ion_mass + 1/neutral_mass)).sqrt()
                e_rel = [Decimal(0), e_speed, Decimal(0)]
                ion_rel = [Decimal(".6") * heavy_speed, Decimal(0), Decimal(".8") * heavy_speed]
                pair_v = [(p[d] - electron_mass * e_rel[d])/total_mass for d in range(3)]
                velocities = [[pair_v[d] + e_rel[d] for d in range(3)],
                              [pair_v[d] + neutral_mass/pair_mass * ion_rel[d] for d in range(3)],
                              [pair_v[d] - ion_mass/pair_mass * ion_rel[d] for d in range(3)]]
                momenta = [[m * v for v in vector] for m, vector in zip(charged_masses, velocities)]
                references[name] = (charged_masses, momenta, velocities)
                if name not in events or set(products.get(name, {})) != {0, 1, 2}:
                    check(name+"/three-product-receipt", False)
                    continue
                if name == "charged-photon":
                    compare(name+"/source-mean-photon-energy",events[name][3],12)
                    compare(name+"/ordered-formation-Q",events[name][4],-8)
                    compare(name+"/outgoing-lab-kinetic-energy",events[name][2],4)
                compare(name+"/COM-energy", events[name][1], Decimal(kinetic),
                        Decimal(kinetic) if kinetic else sum(Decimal(x)**2 for x in p) / (2*total_mass))
                ve, vi, vn = (products[name][i][4:] for i in range(3))
                pair_actual = [(ion_mass*a+neutral_mass*b)/pair_mass for a, b in zip(vi, vn)]
                relative_energy = norm([a-b for a,b in zip(ve, pair_actual)])**2/(2*(1/electron_mass+1/pair_mass))
                heavy_energy = norm([a-b for a,b in zip(vi, vn)])**2/(2*(1/ion_mass+1/neutral_mass))
                # For a zero mode the residual scales with the finite boosted
                # total energy, while a positive mode is checked on its own scale.
                lab_energy = kinetic + sum(Decimal(x)**2 for x in p)/(2*total_mass)
                compare(name+"/electron-pair-relative-energy",relative_energy,electron_k,electron_k or lab_energy)
                compare(name+"/ion-neutral-relative-energy",heavy_energy,heavy_k,heavy_k or lab_energy)
            for name, masses, energy in (("microcanonical-deterministic", [Decimal(16), Decimal(1), Decimal(1)], Decimal(1)),
                    ("extreme-mass-energy-scales", [Decimal("1e100"), Decimal("1e-100"), Decimal("1e-100")], Decimal("1e-200"))):
                if name == "extreme-mass-energy-scales" and platform[2] < 1400:
                    checks.append(dict(name="extreme-mass-energy-scales/positive-domain", status="SKIPPED", reason="individual kinetic energy requires extended exponent range; actual refusal separately checked"))
                    continue
                references[name] = (masses, *jacobi_three(masses, normals, energy))
            references["zero-energy-deterministic"] = ([Decimal(16), Decimal(1), Decimal(1)],
                [[Decimal(16), Decimal(0), Decimal(0)], [Decimal(1), Decimal(0), Decimal(0)], [Decimal(1), Decimal(0), Decimal(0)]],
                [[Decimal(1), Decimal(0), Decimal(0)]]*3)
            check("complete-deterministic-event-set", set(events) == set(references) == set(products))
            for name, (masses, momenta, velocities) in references.items():
                if name not in products: check(name+"/products-present", False); continue
                check(name+"/product-count", set(products[name]) == set(range(len(masses))))
                pscale = max(map(norm, momenta)); vscale = max(map(norm, velocities))
                for index, mass in enumerate(masses):
                    if index not in products[name]: continue
                    row = products[name][index]
                    if len(row) != 7: raise RuntimeError("truncated product phase-space receipt")
                    compare(f"{name}/{index}/mass", row[0], mass)
                    for d in range(3):
                        compare(f"{name}/{index}/momentum-{d}", row[d+1], momenta[index][d], pscale)
                        compare(f"{name}/{index}/velocity-{d}", row[d+4], velocities[index][d], vscale)
                    if name == "tiny-open-polar":
                        compare(f"{name}/{index}/tiny-transverse-preserved", row[2], momenta[index][1])
        if len(ensemble) != 18: raise RuntimeError("truncated ensemble metrics")
        completed_trials = int(ensemble[0]); check("ensemble/full-predeclared-count", completed_trials == TRIALS)
        observed = list(map(float, ensemble[1:9])); counts = list(map(int, ensemble[9:]))
        expected = [.5, 5/16, 1/18, 17/36, 17/36, 1/3, 1/3, 1/3]
        bound = math.sqrt(math.log(2*STATISTICS/FAMILY_ALPHA)/(2*TRIALS))
        names = ["Jacobi-mean", "Jacobi-second-moment", "O-energy-share", "H1-energy-share", "H2-energy-share", "isotropy-x", "isotropy-y", "isotropy-z"]
        for name, value, target in zip(names, observed, expected):
            check("ensemble/"+name, math.isfinite(value) and abs(value-target) <= bound, actual=value, reference=target, simultaneous_99pct_bound=bound)
        for index, count in enumerate(counts, 1):
            x = index/10; target = beta_cdf(x); actual = count/TRIALS
            check("ensemble/Beta-CDF-"+str(index), 0 <= count <= TRIALS and abs(actual-target) <= bound,
                actual=actual, reference=target, simultaneous_99pct_bound=bound)
            ensemble_rows.append(dict(jacobi_energy_fraction=x, count=count, trials=TRIALS,
                empirical_cdf=actual, beta_3_2_3_2_cdf=target, simultaneous_99pct_bound=bound))
        if len(observed)+len(counts) != STATISTICS: raise RuntimeError("frozen statistical family size changed")
        if {p: sha(REPO/p) for p in DEPENDENCIES} != manifest["source_hashes"]:
            raise RuntimeError("source changed during compilation/execution/scoring")
        with (output/"comparison.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=("metric", "actual", "reference", "absolute_bound", "relative_error")); writer.writeheader(); writer.writerows(comparisons)
        with (output/"microcanonical_cdf.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(ensemble_rows[0])); writer.writeheader(); writer.writerows(ensemble_rows)
        result["visualization"] = dict(status="SKIPPED", reason="Matplotlib unavailable")
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
            x = [r["jacobi_energy_fraction"] for r in ensemble_rows]
            axes[0].plot(x, [r["beta_3_2_3_2_cdf"] for r in ensemble_rows], label="independent Beta(3/2,3/2)")
            axes[0].plot(x, [r["empirical_cdf"] for r in ensemble_rows], ".", label="production sampler")
            axes[0].set_xlabel("Jacobi relative-mode energy fraction"); axes[0].set_ylabel("CDF"); axes[0].legend()
            axes[1].errorbar(range(3), observed[2:5], yerr=bound, fmt=".", label="production mean / simultaneous bound")
            axes[1].plot(range(3), expected[2:5], "x", label="independent mass-weighted identity")
            axes[1].set_xticks(range(3)); axes[1].set_xticklabels(["O fixture", "H1 fixture", "H2 fixture"])
            axes[1].set_ylabel("Mean kinetic-energy fraction"); axes[1].legend()
            fig.suptitle("ANALYTICAL REDUCED-MODEL CONTROLS — NO OBSERVATIONS")
            for extension in ("png", "eps"): fig.savefig(output/("reduced_kinematics_comparison."+extension))
            plt.close(fig); result["visualization"] = dict(status="PASS", scope="mathematical diagnostic only")
        except ImportError: pass
        failed = any(c["status"] == "FAIL" for c in checks)
        result.update(status="FAIL" if failed else "PASS", metrics=dict(native_checks=total,
            independent_numeric_comparisons=len(comparisons), ensemble_trials=TRIALS,
            statistical_family=STATISTICS, simultaneous_99pct_bound=bound), artifacts=dict(
                numerical_csv=str(output/"comparison.csv"), statistical_csv=str(output/"microcanonical_cdf.csv")))
    except (OSError, ValueError, ArithmeticError, RuntimeError, subprocess.SubprocessError) as error:
        checks.append(dict(name="infrastructure", status="ERROR", reason=str(error)))
        result.update(status="ERROR", reason=str(error))
    finally:
        result["results"] = checks
        result["counts"] = {s: sum(row["status"] == s for row in checks) for s in ("PASS", "FAIL", "ERROR", "SKIPPED")}
        result["run_manifest"] = str(output/"run_manifest.json")
        (output/"run_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
        (output/"result.json").write_text(json.dumps(result, indent=2)+"\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path); parser.add_argument("--cxx", default="g++")
    args = parser.parse_args()
    output = args.output_dir or REPO/"test_output/srcMoon"/("reduced-kinematics-"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ"))
    try: result = run(output, args.cxx)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(str(error), file=sys.stderr); return 2
    print(json.dumps(result, indent=2))
    return {"PASS":0,"FAIL":1,"ERROR":2,"SKIPPED":77}[result["status"]]


if __name__ == "__main__": raise SystemExit(main())
