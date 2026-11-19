#!/usr/bin/env python3
"""Independent mathematical controls for production piecewise-affine hazards.

No observation, physical provider, PIC build or random selection is present.
Decimal implements the analytical reference independently; the native adapter
only invokes the maintained production source and serializes its results.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
DEPENDENCIES = ("srcMoon/species/PiecewiseHazard.h", "srcMoon/species/PiecewiseHazard.cpp",
    "srcMoon/test/reaction_execution/piecewise_hazard/probe.cpp",
    "srcMoon/test/reaction_execution/piecewise_hazard/CONTRACT.md",
    "srcMoon/test/reaction_execution/piecewise_hazard/README.md",
    "srcMoon/test/reaction_execution/piecewise_hazard/run.py")
DECIMAL_PRECISION = 160
MAXIMUM_ROOT_REFINEMENTS = 65536  # explicit arithmetic work bound, not physics
ROUNDING_MULTIPLIER = 128  # predeclared relative arithmetic comparison bound


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def segment(begin, end, *channels):
    return {"begin": str(begin), "end": str(end),
            "channels": [[int(identity), str(first), str(last)] for identity, first, last in channels]}


def fixtures():
    """Frozen analytical definitions: changing them after FAIL needs review."""
    rows = []
    def add(name, history, end, depth="1", begin="0", operation="locate",
            expected=None, exact_binary=True, limit=MAXIMUM_ROOT_REFINEMENTS):
        rows.append(dict(name=name, history=history, begin=str(begin), end=str(end),
            depth=str(depth), operation=operation, expected=expected,
            exact_binary=exact_binary, maximum_root_refinements=limit))
    constant = [segment(0, 2, (11, 2, 2))]
    add("constant", constant, 2)
    add("constant_exposure", constant, 2, operation="integrate")
    add("zero_hazard_no_event", [segment(0, 2, (11, 0, 0))], 2)
    add("rising", [segment(0, 2, (11, 0, 4))], 2)
    add("falling", [segment(0, 2, (11, 2, 0))], 2)
    competing = [segment(0, 2, (11, 0, 4), (29, 1, 1))]
    add("instantaneous_not_integrated_branches", competing, 2, 2)
    add("competing_exposure", competing, 2, operation="integrate")
    jump = [segment(0, 2, (11, 1, 1)), segment(2, 4, (11, 10, 10))]
    add("left_limit_jump", jump, 4, 2)
    add("strictly_future_into_right_segment", jump, 4, "up:2")
    add("strictly_future_beyond_window", jump, 2, "up:2")
    add("exact_requested_endpoint", jump, 2, 2)
    left_zero = [segment(0, 1, (11, 2, 0)), segment(1, 2, (11, 7, 7))]
    add("left_zero_right_positive", left_zero, 2)
    add("left_zero_strictly_future", left_zero, 2, "up:1")
    add("flat_derivative_near_endpoint", left_zero, 1, "down:1")
    add("rounded_prefix_is_not_consumed_depth", [segment(0, 1, (11, 1, "down:1")),
        segment(1, 2, (11, 1, 1))], 2)
    add("underflowed_prefix_preserved_exactly", [segment(0, "1e-4000", (11, "1e-4000", "1e-4000")),
        segment("1e-4000", 1, (11, 1, 1))], 1, ".5", exact_binary=False)
    add("zero_width_covered", constant, 1, begin=1)
    add("partial_window", competing, 2, 1, begin=1)
    unpartitioned = [segment(0, 3, (11, 0, 6), (29, 4, 1))]
    partitioned = [segment(i, i+1, (11, 2*i, 2*(i+1)), (29, 4-i, 3-i)) for i in range(3)]
    for name, history in (("whole", unpartitioned), ("partitioned", partitioned)):
        add(name+"_event", history, 3, 10)
        add(name+"_exposure", history, 3, operation="integrate")
    for scale in ("1e-100", "1e100", "1e-4000", "1e4000"):
        exponent = int(scale.split("e")[1]); duration = "1e"+str(-exponent)
        add("scale_"+str(exponent).replace("-", "minus"),
            [segment(0, duration, (11, scale, scale))], duration, ".5", exact_binary=False)
    add("overflowing_horizon_finite_event", [segment(0, "1e4000", (11, "1e4000", "1e4000"))],
        "1e4000", "1e4000", exact_binary=False)
    for index in range(1, 33):
        # Powers-of-two depths and integer affine inputs support exact binary
        # endpoint/predecessor admission without a fitted numerical window.
        add(f"quadratic_grid_{index:02d}", competing, 2, str(Decimal(index)/Decimal(8)))
    gap = [segment(0, 1, (11, 1, 1)), segment(2, 3, (11, 1, 1))]
    add("gap_even_after_possible_event", gap, 3, ".5", expected="IncompleteSupport")
    add("zero_width_gap", gap, "1.5", begin="1.5", expected="IncompleteSupport")
    add("beyond_support", constant, 3, ".5", expected="IncompleteSupport")
    add("before_support", [segment(1, 2, (11, 1, 1))], 2, expected="IncompleteSupport")
    add("empty_history", [], 1, expected="IncompleteSupport")
    add("overlap", [segment(0, 2, (11, 1, 1)), segment(1, 3, (11, 1, 1))], 3, expected="InvalidInput")
    add("reversed_segment", [segment(2, 1, (11, 1, 1))], 1, expected="InvalidInput")
    add("zero_duration_segment", [segment(0, 0, (11, 1, 1))], 0, expected="InvalidInput")
    add("negative_rate", [segment(0, 1, (11, -1, 1))], 1, expected="InvalidInput")
    add("nonfinite_rate", [segment(0, 1, (11, "nan", 1))], 1, expected="InvalidInput")
    add("infinite_rate", [segment(0, 1, (11, 1, "inf"))], 1, expected="InvalidInput")
    add("duplicate_channel", [segment(0, 1, (11, 1, 1), (11, 2, 2))], 1, expected="InvalidInput")
    add("changed_channel", [segment(0, 1, (11, 1, 1)), segment(1, 2, (29, 1, 1))], 2, expected="InvalidInput")
    add("missing_channel", [segment(0, 1, (11, 1, 1), (29, 1, 1)), segment(1, 2, (11, 1, 1))], 2, expected="InvalidInput")
    add("empty_channels", [segment(0, 1)], 1, expected="InvalidInput")
    add("reversed_query", constant, 0, begin=1, expected="InvalidInput")
    add("negative_query", constant, 1, begin=-1, expected="InvalidInput")
    add("nonfinite_query", constant, "inf", expected="InvalidInput")
    add("zero_depth", constant, 1, 0, expected="InvalidInput")
    add("negative_depth", constant, 1, -1, expected="InvalidInput")
    add("nonfinite_depth", constant, 1, "nan", expected="InvalidInput")
    add("zero_work_limit", constant, 1, limit=0, expected="InvalidInput")
    add("root_limit_exhaustion", [segment(0, 2, (11, 2, 0))], 2,
        limit=1, expected="NumericalError")
    add("integral_overflow", [segment(0, "1e4000", (11, "1e4000", "1e4000"))],
        "1e4000", operation="integrate", expected="NumericalError", exact_binary=False)
    add("integral_underflow", [segment(0, "1e-4000", (11, "1e-4000", "1e-4000"))],
        "1e-4000", operation="integrate", expected="NumericalError", exact_binary=False)
    return rows


def hex_decimal(text):
    """Decode the adapter's exact native binary value, without float64 loss."""
    if text in ("nan", "-nan", "inf", "-inf"):
        return Decimal(text)
    negative = text.startswith("-")
    mantissa, exponent = text.lstrip("+-").lower().split("p")
    integer, dot, fraction = (mantissa[2:] if mantissa.startswith("0x") else mantissa).partition(".")
    value = Decimal(int(integer+fraction, 16)) * (Decimal(2)**(int(exponent)-4*len(fraction)))
    return -value if negative else value


def analytical(case, begin, end, depth):
    """Independent affine integral and closed-form quadratic first crossing."""
    channel_areas = {}; cumulative = Decimal(0)
    for index, s in enumerate(case["history"]):
        first, last = Decimal(s["begin"]), Decimal(s["end"])
        a, b = max(first, begin), min(last, end)
        if b <= a:
            continue
        slopes = {k: (Decimal(r1)-Decimal(r0))/(last-first) for k, r0, r1 in s["channels"]}
        starts = {k: Decimal(r0)+slopes[k]*(a-first) for k, r0, r1 in s["channels"]}
        length = b-a
        areas = {k: starts[k]*length+slopes[k]*length*length/2 for k in starts}
        area = sum(areas.values(), Decimal(0))
        if case["operation"] == "locate" and cumulative+area >= depth:
            remaining = depth-cumulative
            rate, slope = sum(starts.values(), Decimal(0)), sum(slopes.values(), Decimal(0))
            offset = length if remaining == area else remaining/rate if slope == 0 else (
                2*remaining/(rate+(rate*rate+2*slope*remaining).sqrt()))
            return dict(status="Event", time=a+offset, segment=index, exposure=depth)
        cumulative += area
        for k, value in areas.items():
            channel_areas[k] = channel_areas.get(k, Decimal(0))+value
    return dict(status="Complete", exposure=cumulative, channels=channel_areas)


def cumulative_at(case, begin, end):
    copy = dict(case, operation="integrate")
    return analytical(copy, begin, end, Decimal(1))["exposure"]


def resolved_case(case, digits):
    """Resolve fixture-adjacent binary rates without duplicating any kernel."""
    def number(text):
        if text.startswith(("up:", "down:")):
            up = text.startswith("up:"); base = Decimal(text.split(":", 1)[1])
            if base not in (1, 2):
                raise RuntimeError("adjacent-value fixture oracle only declares bases1 and2")
            power = (1 if base == 1 else 2)-digits
            return str(base+(Decimal(2)**power if up else -(Decimal(2)**(power-1))))
        return text
    return dict(case, history=[dict(s, channels=[[k, number(a), number(b)] for k, a, b in s["channels"]]) for s in case["history"]])


def run(output, compiler="g++"):
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("output directory must be new or empty: "+str(output))
    output.mkdir(parents=True, exist_ok=True)
    result = dict(schema="srcMoon-piecewise-hazard-controls-v1", status="ERROR", results=[],
        counts={s: 0 for s in ("PASS", "FAIL", "ERROR", "SKIPPED")},
        scope="PURE ANALYTICAL SOFTWARE CONTROLS; no observations, physical provider, PIC wiring or formal qualification",
        formal_I21_I24_status="NOT EVALUATED")
    manifest = dict(schema="srcMoon-piecewise-hazard-manifest-v1", source_revision=subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
        source_hashes={p: sha(REPO/p) for p in DEPENDENCIES}, external_data=[], MPI="none", random_seed=None,
        roadmap_hashes={str(p): sha(p) for p in [REPO/"srcMoon/development_and_validation.md",
            Path("/data/vtenishe/moon_validation_data/development_and_validation.md")] if p.exists()},
        numeric_controls=dict(maximum_root_refinements=MAXIMUM_ROOT_REFINEMENTS,
            comparison_rounding_multiplier=ROUNDING_MULTIPLIER, decimal_precision=DECIMAL_PRECISION), commands=[])
    comparisons = []; checks = []
    def check(name, passed, **details):
        checks.append(dict(name=name, status="PASS" if passed else "FAIL", **details))
    try:
        cxx = shutil.which(compiler)
        if not cxx:
            reason = "C++17 compiler unavailable: "+compiler
            checks.append(dict(name="compiler-availability", status="SKIPPED", reason=reason))
            result.update(status="SKIPPED", reason=reason, results=checks,
                counts={s: int(s == "SKIPPED") for s in ("PASS", "FAIL", "ERROR", "SKIPPED")})
            return result
        version = subprocess.run([cxx, "--version"], capture_output=True, text=True, timeout=15)
        (output/"compiler-version.log").write_text(version.stdout+version.stderr)
        if version.returncode:
            raise RuntimeError("compiler identity query failed")
        manifest["compiler"] = dict(path=cxx, version=version.stdout)
        executable = output/"piecewise-hazard-probe"
        command = [cxx, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror", "-pedantic", "-I", str(REPO),
            str(REPO/"srcMoon/species/PiecewiseHazard.cpp"), str(HERE/"probe.cpp"), "-o", str(executable)]
        with (output/"compile.log").open("w") as log:
            completed = subprocess.run(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT, timeout=120)
        manifest["commands"].append(dict(argv=command, returncode=completed.returncode))
        if completed.returncode:
            raise RuntimeError("production compilation failed; see compile.log")
        manifest["executable"] = dict(path=str(executable), sha256=sha(executable))
        cases = fixtures(); (output/"fixtures.json").write_text(json.dumps(cases, indent=2)+"\n")
        lines = [str(len(cases))]
        for case in cases:
            lines.append(" ".join(map(str, [case["name"], case["operation"], case["begin"], case["end"],
                case["depth"], case["maximum_root_refinements"], len(case["history"])])))
            for s in case["history"]:
                lines.append(f'{s["begin"]} {s["end"]} {len(s["channels"])}')
                lines.extend(" ".join(map(str, channel)) for channel in s["channels"])
        text = "\n".join(lines)+"\n"; (output/"fixtures.txt").write_text(text)
        completed = subprocess.run([str(executable)], input=text, capture_output=True, text=True, timeout=300)
        (output/"stdout.log").write_text(completed.stdout); (output/"stderr.log").write_text(completed.stderr)
        manifest["commands"].append(dict(argv=[str(executable)], returncode=completed.returncode))
        if completed.returncode:
            raise RuntimeError("native adapter failed; see stdout/stderr")
        native = {}; platform = None
        for line in completed.stdout.splitlines():
            fields = line.split("\t"); tag = fields[0]
            if tag == "PLATFORM":
                if platform is not None or len(fields) != 5:
                    raise RuntimeError("duplicate/malformed platform receipt")
                platform = dict(digits=int(fields[1]), epsilon=fields[2], maximum_exponent=int(fields[3]), minimum_exponent=int(fields[4]))
                continue
            if len(fields) < 2:
                raise RuntimeError("malformed native row")
            row = native.setdefault(fields[1], {"rates": {}, "channels": {}})
            if tag == "RATE": row["rates"][int(fields[2])] = fields[3]
            elif tag == "EXPOSURE": row["channels"][int(fields[2])] = fields[3]
            elif tag in ("RESULT", "INPUT", "INPUTHEX", "EVENTHEX", "DIAGNOSTIC"):
                if tag in row: raise RuntimeError("duplicate native receipt")
                row[tag] = fields[2:]
            else: raise RuntimeError("unknown native protocol record: "+tag)
        if not platform or set(native) != {c["name"] for c in cases}:
            raise RuntimeError("missing/extra native cases")
        manifest["long_double_platform"] = platform
        with localcontext() as context:
            context.prec = DECIMAL_PRECISION; context.Emax = 999999; context.Emin = -999999
            epsilon = Decimal(platform["epsilon"]); relative_bound = ROUNDING_MULTIPLIER*epsilon
            for case in cases:
                name = case["name"]; row = native[name]
                if any(key not in row for key in ("RESULT", "INPUTHEX", "EVENTHEX", "DIAGNOSTIC")) or len(row["RESULT"]) != 7:
                    raise RuntimeError("incomplete native control receipt: "+name)
                status, unchanged, history_same, actual_exposure, _, segment_index, total_rate = row["RESULT"]
                check(name+"/history-unchanged", history_same == "1")
                if any("4000" in token for s in case["history"] for token in [s["begin"], s["end"]]+[x for ch in s["channels"] for x in ch[1:]]) and platform["maximum_exponent"] < 14000:
                    checks.append(dict(name=name+"/platform-domain", status="SKIPPED", reason="extended exponent fixture exceeds native long-double range"))
                    continue
                if case["expected"]:
                    check(name+"/status", status == case["expected"], expected=case["expected"], actual=status)
                    check(name+"/output-preserved", unchanged == "1")
                    check(name+"/diagnostic", bool(row["DIAGNOSTIC"][0]))
                    continue
                begin, end, depth = map(hex_decimal, row["INPUTHEX"])
                mathematical_case = resolved_case(case, platform["digits"])
                reference = analytical(mathematical_case, begin, end, depth)
                check(name+"/status", status == reference["status"], expected=reference["status"], actual=status)
                if status != reference["status"]:
                    continue
                def number(metric, actual, expected):
                    actual = Decimal(actual); expected = Decimal(expected)
                    passed = actual.is_finite() and (actual == 0 if expected == 0 else abs(actual-expected) <= relative_bound*abs(expected))
                    check(name+"/"+metric, passed, actual=str(actual), expected=str(expected))
                    comparisons.append(dict(case=name, metric=metric, actual=str(actual), reference=str(expected),
                        relative_error=str(abs((actual-expected)/expected)) if expected else "0" if actual == 0 else "Infinity"))
                number("exposure", actual_exposure, reference["exposure"])
                if status == "Complete":
                    if case["operation"] == "integrate":
                        expected_channels = reference["channels"] or {int(k): Decimal(0) for k, _, _ in case["history"][0]["channels"]}
                        check(name+"/channel-identities", set(row["channels"]) == set(expected_channels))
                        for k, expected in expected_channels.items(): number("channel-exposure-"+str(k), row["channels"].get(k, "NaN"), expected)
                    else: check(name+"/no-event-no-rates", not row["rates"])
                    continue
                event, predecessor = map(hex_decimal, row["EVENTHEX"])
                number("event-time", event, reference["time"])
                check(name+"/left-limit-segment", int(segment_index) == reference["segment"])
                s = mathematical_case["history"][reference["segment"]]; first, last = Decimal(s["begin"]), Decimal(s["end"])
                rates = {int(k): Decimal(a)+(Decimal(b)-Decimal(a))*(event-first)/(last-first) for k, a, b in s["channels"]}
                check(name+"/instantaneous-channel-identities", set(row["rates"]) == set(rates))
                for k, expected in rates.items(): number("instantaneous-rate-"+str(k), row["rates"].get(k, "NaN"), expected)
                number("instantaneous-total", total_rate, sum(rates.values(), Decimal(0)))
                if case["exact_binary"]:
                    at = cumulative_at(mathematical_case, begin, event); before = cumulative_at(mathematical_case, begin, predecessor)
                    check(name+"/first-representable-crossing", before < depth <= at,
                        predecessor_exposure=str(before), depth=str(depth), event_exposure=str(at))
                if name == "underflowed_prefix_preserved_exactly":
                    check(name+"/strict-future-representable-time", event == Decimal(".5")+Decimal(2)**(-platform["digits"]))
        if {p: sha(REPO/p) for p in DEPENDENCIES} != manifest["source_hashes"]:
            raise RuntimeError("source changed during mathematical control execution")
        with (output/"comparison.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=("case", "metric", "actual", "reference", "relative_error"))
            writer.writeheader(); writer.writerows(comparisons)
        result["visualization"] = dict(status="SKIPPED", reason="Matplotlib unavailable")
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            points = [r for r in comparisons if r["metric"] == "event-time" and r["case"].startswith("quadratic_grid_")]
            fig, axes = plt.subplots(2, 1, figsize=(8, 6), constrained_layout=True)
            axes[0].plot(range(len(points)), [float(r["reference"]) for r in points], "-", label="independent Decimal analytical root")
            axes[0].plot(range(len(points)), [float(r["actual"]) for r in points], ".", label="production")
            axes[0].set_ylabel("Event relative time [s]"); axes[0].legend()
            axes[1].plot(range(len(points)), [float(r["relative_error"]) for r in points], ".")
            axes[1].set_ylabel("Relative root difference"); axes[1].set_xlabel("Predeclared analytical fixture index")
            fig.suptitle("MATHEMATICAL SOFTWARE CONTROLS — NO OBSERVATIONS")
            for extension in ("png", "eps"): fig.savefig(output/("affine_hazard_comparison."+extension))
            plt.close(fig); result["visualization"] = dict(status="PASS", scope="diagnostic mathematical comparison only")
        except ImportError:
            pass
        counts = {s: sum(row["status"] == s for row in checks) for s in ("PASS", "FAIL", "ERROR", "SKIPPED")}
        result.update(status="FAIL" if counts["FAIL"] else "PASS", results=checks, counts=counts,
            metrics=dict(cases=len(cases), arithmetic_relative_bound=str(relative_bound)),
            artifacts=dict(comparison_csv=str(output/"comparison.csv"), stdout=str(output/"stdout.log"), stderr=str(output/"stderr.log")))
    except (OSError, ValueError, ArithmeticError, RuntimeError, subprocess.SubprocessError) as error:
        checks.append(dict(name="infrastructure", status="ERROR", reason=str(error)))
        result.update(status="ERROR", reason=str(error), results=checks,
            counts={s: sum(row["status"] == s for row in checks) for s in ("PASS", "FAIL", "ERROR", "SKIPPED")})
    finally:
        result["run_manifest"] = str(output/"run_manifest.json")
        (output/"run_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
        (output/"result.json").write_text(json.dumps(result, indent=2)+"\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--cxx", default="g++")
    args = parser.parse_args()
    output = args.output_dir or REPO/"test_output/srcMoon"/("piecewise-hazard-"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ"))
    try: result = run(output, args.cxx)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(str(error), file=sys.stderr); return 2
    print(json.dumps(result, indent=2))
    return {"PASS": 0, "FAIL": 1, "ERROR": 2, "SKIPPED": 77}[result["status"]]


if __name__ == "__main__": raise SystemExit(main())
