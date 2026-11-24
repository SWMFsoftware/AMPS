#!/usr/bin/env python3
"""Exercise the I24 single-point production provider and independently audit actual native plans."""
from __future__ import annotations
import argparse
import csv
from decimal import Decimal as D
from hashlib import sha256
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
TEST = HERE.parents[1]
ROOT = TEST.parents[1]
CASE = HERE / "reference/sodium_case.json"
DEFAULT_INVENTORY = ROOT / "test_output/srcMoon/m3-approved-products-20261118T163142Z/source-acquisition-v2/reference_mass_inventory.json"
sys.path.insert(0, str(TEST / "sodium_cohort"))
from verify_numerical import source_uniform
EXIT = {"PASS": 0, "FAIL": 1, "ERROR": 2, "SKIPPED": 77}


class MissingPrerequisite(FileNotFoundError):
    """Only a missing external input/tool can skip the case; lost run output is an error."""


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def prepare_inputs(lock, inventory_path, output):
    """Obtain CLI inputs for either MPI adapter from the same frozen sources.

    Authentication of the immutable release is also performed by both the
    producer export and the native reader. The independent rate oracle uses
    the original printed area (in Å²) recorded before this test was executed.
    """
    lock, inventory_path, output = Path(lock).resolve(), Path(inventory_path).resolve(), Path(output).resolve()
    controls = json.loads(CASE.read_text())
    for path, field in ((lock, "lock_sha256"), (inventory_path, "mass_inventory_sha256")):
        if not path.is_file():
            raise MissingPrerequisite("required external input is unavailable: "+str(path))
        if sha(path) != controls[field]:
            raise ValueError("changed source identity: " + str(path))
    release = json.loads(lock.read_text())
    inventory = json.loads(inventory_path.read_text())
    extractor = TEST / "physical_state_sources/prepare.py"
    spec = importlib.util.spec_from_file_location("i24_mass_evidence", extractor)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    raw = inventory_path.parent / "raw"
    for filename in module.SOURCES:
        if not (raw/filename).is_file():
            raise MissingPrerequisite("required NIST source packet member is unavailable: "+str(raw/filename))
    if module.extract(raw) != inventory:
        raise ValueError("NIST isotope/energy/constant native-source replay differs from the retained inventory")
    hashes = {str(p): sha(p) for p in (lock, inventory_path, CASE, extractor, *(raw / name for name in module.SOURCES))}
    # Only the declared r14 selected source can generate the compiled projection.
    exporter = ROOT / "srcMoon/tools/export_ground_capture.py"
    point = output / "point.csv"
    command = [sys.executable, "-B", str(exporter), "--lock", str(lock), "--output", str(point)]
    result = subprocess.run(command, capture_output=True, text=True, timeout=30, check=False)
    (output / "export.stdout.log").write_text(result.stdout)
    (output / "export.stderr.log").write_text(result.stderr)
    dump(output / "export-command.json", command)
    if result.returncode:
        raise ValueError("could not export r14 model: " + result.stderr.strip())
    hashes.update({str(p): sha(p) for p in (point, exporter)})
    neutral = {x["identity"]: D(x["reference_neutral_mass_kg"]) for x in inventory["isotopes"]}
    constants = {key: D(value["value"]) for key, value in inventory["constants"].items()}
    thresholds = {x["element"]: D(x["ionization_energy_eV"]) for x in inventory["atomic_ground_thresholds"]}
    inputs = [neutral["1H"], neutral["23Na"], constants["electron mass"], constants["atomic mass constant"],
              constants["elementary charge"], thresholds["H"], thresholds["Na"]]
    argv = [str(lock.parent), release["export_path"], release["export_manifest_sha256"], release["release_id"],
            release["release_manifest_sha256"], release["profile_id"], release["profile_sha256"],
            str(point), *map(str, inputs)]
    h, na, electron, mu, ev, h_ie, na_ie = inputs
    speed = (2 * D(controls["energy_ev_per_u"]) * ev / mu).sqrt()
    sigma = D(controls["printed_area_angstrom_squared"]) * D("1e-20")
    rate = D(controls["proton_density_m3"]) * sigma * speed
    oracle = dict(h=h, na=na, proton=h-electron, ion=na-electron, speed=speed, sigma=sigma,
                  q=(h_ie-na_ie)*ev, rate=rate, horizon=D(controls["optical_depth"])/rate)
    return argv, oracle, hashes


def audit(directory, ranks, case, reference):
    """Reconstruct each first hit and product momentum independently of the C++ solver.

    Computing the energies in the source frame from both recorded velocities
    detects mass swaps and wrong-Q signs. The incident projection then tests
    the declared angular convention in the pair COM; the transverse azimuth
    test detects collapse of the trajectory onto one plane.
    """
    count = case["parents"]
    paths = [directory / ("parents-rank-%d.csv" % rank) for rank in range(ranks)]
    rows = {}
    for rank, path in enumerate(paths):
        with path.open() as stream:
            for row in csv.DictReader(stream):
                if None in row or any(v is None for v in row.values()):
                    raise ValueError("malformed native parent CSV record")
                i = int(row["parent"])
                if not 0 <= i < count or i % ranks != rank or i in rows:
                    raise ValueError("duplicate or misplaced native parent ID")
                rows[i] = row
    if len(rows) != count:
        raise ValueError("native execution omitted a parent history")
    with (directory / "rate.csv").open() as stream:
        rate, = csv.DictReader(stream)
    if any(not math.isfinite(float(rate[k])) for k in ("rate_s", "sigma_m2", "speed_m_s")):
        raise ValueError("nonfinite native rate receipt")
    frequency, speed, q = map(float, (reference["rate"], reference["speed"], reference["q"]))
    mp, mi, mh = map(float, (reference["proton"], reference["ion"], reference["h"]))
    # NIST model uses a common reference mass scheme with explicit formation
    # energies. Total incoming/outgoing masses match under this convention.
    px, py = .6*mp*speed, .8*mp*speed
    kinetic = .5*mp*speed**2
    vx_com, vy_com = px/(mh+mi), py/(mh+mi)
    horizon = float(reference["horizon"])
    mu_samples, az_samples, times = [], [], []
    worst_time = worst_momentum = worst_energy = worst_q = 0.
    wrong_ids = wrong_fates = 0
    for i in range(count):
        row = rows[i]
        t_first = -math.log(source_uniform(case["seed"], case["cohort"], i))/frequency
        hit = row["fate"] == "capture"
        if row["fate"] not in ("capture", "survived"):
            raise ValueError("unknown native fate")
        wrong_fates += (t_first <= horizon) != hit
        time_s = float(row["time_s"])
        if not math.isfinite(time_s):
            raise ValueError("nonfinite native event time")
        worst_time = max(worst_time, abs(time_s-min(t_first, horizon))/horizon)
        if not hit:
            continue
        times.append(time_s/horizon)
        ion = [float(row["na"+j]) for j in "xyz"]
        hydrogen = [float(row["h"+j]) for j in "xyz"]
        q_out = float(row["Q_J"])
        if not all(math.isfinite(x) for x in ion+hydrogen+[q_out]):
            raise ValueError("nonfinite product kinematics")
        wrong_ids += (row["identity_ok"] != "1" or
                      any(D(row[k]) != D(case["parent_weight"]) for k in ("na_weight","h_weight","bath_weight")))
        residual = [mi*a + mh*b - c for a,b,c in zip(ion,hydrogen,(px,py,0.))]
        worst_momentum = max(worst_momentum, math.sqrt(sum(x*x for x in residual))/(mp*speed))
        k_after = .5*(mi*sum(x*x for x in ion) + mh*sum(x*x for x in hydrogen))
        worst_energy = max(worst_energy, abs(k_after-kinetic-q)/kinetic)
        worst_q = max(worst_q, abs(q_out-q)/abs(q))
        hx,hy,hz = hydrogen[0]-vx_com, hydrogen[1]-vy_com, hydrogen[2]
        norm = math.sqrt(hx*hx+hy*hy+hz*hz)
        mu_samples.append((.6*hx+.8*hy)/norm)
        az_samples.append((math.atan2(hz,-.8*hx+.6*hy)/(2*math.pi)) % 1)
    captures = len(times)
    f_survive = (count-captures)/count
    tolerance = case["arithmetic_relative_tolerance"]
    checks = []
    def check(name, passed, actual, reference_value):
        checks.append(dict(name=name, status="PASS" if passed else "FAIL",
                           actual=actual, reference=reference_value))
    check("direct n sigma g [1/s]", abs(float(rate["rate_s"])/frequency-1)<=.01, float(rate["rate_s"]), frequency)
    for field, ref in (("sigma_m2", float(reference["sigma"])),("speed_m_s", speed)):
        check(field, abs(float(rate[field])/ref-1)<=tolerance, float(rate[field]), ref)
    for guard in ("off_energy_refusal","negative_density_refusal","zero_density_zero_rate"):
        check(guard, rate[guard]=="1", rate[guard], "1")
    # These two controls use the first native collision. A run with zero
    # collisions fails its lifetime and weight criteria and has no guard row.
    if captures:
        with (directory / "event-controls.csv").open() as stream:
            guard_row, = csv.DictReader(stream)
        for name in ("closed_endpoint_refused_prior_plan_preserved","rescaled_daughter_weights_invariant"):
            check(name, guard_row[name]=="1", guard_row[name], "1")
    check("all parent first-hit classifications", wrong_fates==0, wrong_fates, 0)
    check("all event-time errors / horizon", worst_time<=tolerance, worst_time, tolerance)
    lifetime = -horizon/math.log(f_survive) if 0 < f_survive < 1 else None
    check("end-point lifetime [s]", lifetime is not None and abs(lifetime*frequency-1)<=.03, lifetime, 1/frequency)
    check("product identity, charge and represented weights", captures>0 and wrong_ids==0, wrong_ids, 0)
    for label, value in (("momentum error / incoming momentum",worst_momentum),
                         ("energy error / incoming kinetic energy",worst_energy),
                         ("reaction Q relative error",worst_q)):
        check(label, captures>0 and value<=tolerance, value, tolerance)
    # Dvoretzky–Kiefer–Wolfowitz simultaneous CDF checks. Bonferroni divides
    # the predeclared 1% diagnostic false-positive allowance across two CDFs.
    limit = math.sqrt(math.log(4/case["stochastic_family_alpha"])/(2*captures)) if captures else 0.
    for name, values in (("pair-COM cos(theta)",[(v+1)/2 for v in mu_samples]),("pair-COM azimuth",az_samples)):
        distance = max((max((j+1)/captures-x, x-j/captures) for j,x in enumerate(sorted(values))),default=1.)
        check(name+" CDF", distance<=limit, distance, limit)
    return rows, checks, sorted(mu_samples), sorted(times)


def plot(output, times, mu, parents, exposure):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,2,figsize=(9.5,4), constrained_layout=True)
    grid=[j/100 for j in range(101)]
    axes[0].step([0]+times+[1], [1]+[1-(j+1)/parents for j in range(len(times))]+[1-len(times)/parents],
                 where="post", label="Native parents")
    axes[0].plot(grid,[math.exp(-v*exposure) for v in grid],"k--",label="Exponential")
    axes[0].set(xlabel="Elapsed time / horizon", ylabel="Surviving Na fraction")
    axes[1].step(mu,[(j+1)/len(mu) for j in range(len(mu))],where="post",label="Native H recoil")
    axes[1].plot([-1,1],[0,1],"k--",label="Isotropic COM")
    axes[1].set(xlabel="cos θ from incident beam",ylabel="Cumulative probability")
    for ax in axes: ax.legend(framealpha=1)
    for suffix in ("png","eps"): fig.savefig(output/("comparison."+suffix))
    plt.close(fig)


def run(lock, inventory, output, ranks=8, wall_seconds=3600):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    case = json.loads(CASE.read_text())
    manifest = dict(schema="srcMoon-I24-ground-capture-manifest-v1", scope="compiled MPI provider and event planner",
                    source_revision=subprocess.check_output(["git","rev-parse","HEAD"], cwd=ROOT, text=True).strip(),
                    seed=case["seed"], cohorts=case["parents"], mpi_layouts=list(dict.fromkeys((1,ranks))),
                    openmp_threads=1, commands=[], model=case["model_id"], formal_status="SKIPPED / NOT QUALIFIED")
    checks = []
    result = dict(status="ERROR", scope=manifest["scope"], formal_status=manifest["formal_status"], checks=checks)
    stdout, stderr = output/"stdout.log", output/"stderr.log"
    stdout.touch(); stderr.touch()
    try:
        if not 0 < wall_seconds <= 3600 or not math.isfinite(wall_seconds) or ranks < 1:
            raise ValueError("positive MPI ranks and a total wall allowance of at most 3600 s are required")
        def invoke(argv, directory):
            remaining = wall_seconds - (time.monotonic()-start)
            if remaining <= 0:
                raise TimeoutError("I24 total time allowance exhausted")
            argv = list(map(str, argv))
            manifest["commands"].append({"cwd": str(directory), "argv": argv})
            with stdout.open("a") as out, stderr.open("a") as err:
                out.write(json.dumps(argv)+"\n"); out.flush()
                result = subprocess.run(argv, cwd=directory, stdout=out, stderr=err, timeout=remaining,
                                        env=dict(os.environ, OMP_NUM_THREADS="1", OMPI_MCA_hwloc_base_binding_policy="none"))
            manifest["commands"][-1]["returncode"] = result.returncode
            if result.returncode:
                raise RuntimeError("native command exited %d (see %s)" % (result.returncode, stderr))
        bundle_args, ref, hashes = prepare_inputs(lock, inventory, output)
        manifest["input_sha256"] = hashes
        shutil.copyfile(CASE, output/"case.json")
        dump(output/"source_reference.json", {k:str(v) for k,v in ref.items()})
        mpi_compiler, launcher = shutil.which("mpicxx"), shutil.which("mpirun")
        if not mpi_compiler or not launcher:
            raise MissingPrerequisite("mpicxx and mpirun required for the MPI provider test")
        manifest["compiler"] = subprocess.check_output([mpi_compiler, "--version"], text=True).splitlines()[0]
        manifest["mpi_launcher"] = subprocess.check_output([launcher, "--version"], text=True).splitlines()[0]
        species = ROOT/"srcMoon/species"
        shared = ROOT/"src/models/chemical_database"
        sources = [HERE/"ground_capture_probe.cpp",
                   *(species/(name+".cpp") for name in ("SodiumGroundCapture","ReducedKinematics","QualificationRandom","ReactionExecution","PiecewiseHazard")),
                   shared/"Runtime.cpp", shared/"DataIO.cpp"]
        binary = output/"ground-capture-probe"
        # -MMD below is intentionally not passed on the multi-unit compiler
        # invocation. Record the actual local include DAG with the same flags.
        flags = ["-std=c++17","-O2","-g","-Wall","-Wextra","-pedantic","-DOMPI_SKIP_MPICXX","-I",str(ROOT)]
        dep = subprocess.run([mpi_compiler, *flags, "-MM", *map(str,sources)], cwd=ROOT,
                             capture_output=True, text=True, check=True, timeout=45)
        (output/"compiler-dependencies.txt").write_text(dep.stdout)
        dependencies = set()
        for line in dep.stdout.replace("\\\n"," ").splitlines():
            _, member_paths = line.split(":",1)
            dependencies.update((ROOT/path).resolve() for path in member_paths.split())
        dependencies.update((HERE/"ground_capture_setup.h", Path(__file__), TEST/"sodium_cohort/verify_numerical.py"))
        manifest["source_sha256"] = {str(path):sha(path) for path in sorted(dependencies)}
        invoke([mpi_compiler,*flags,*sources,"-o",binary], ROOT)
        manifest["binary_sha256"] = sha(binary)
        first = None
        for size in manifest["mpi_layouts"]:
            directory = output/("mpi-%d" % size)
            directory.mkdir()
            invoke([launcher,"-np",size,binary,*bundle_args,case["proton_density_m3"],str(ref["speed"]),
                    str(ref["horizon"]),case["seed"],case["cohort"],case["parents"],directory,
                    case["parent_weight"],*case["daughter_weights"]], directory)
            rows, metrics, angles, times = audit(directory, size, case, ref)
            checks.extend(dict(check, mpi_ranks=size) for check in metrics)
            if first is None:
                first = rows
                plot(output, times, angles, case["parents"], float(case["optical_depth"]))
            else:
                differences = sum(row != first[i] for i,row in rows.items())
                checks.append(dict(name="per-parent rank invariance", status="PASS" if not differences else "FAIL",
                                   actual=differences, reference=0, mpi_ranks=size))
        with (output/"comparison.csv").open("w",newline="") as stream:
            writer = csv.DictWriter(stream,fieldnames=("mpi_ranks","name","status","actual","reference"))
            writer.writeheader();writer.writerows(checks)
        if any(sha(Path(path)) != identity for path,identity in manifest["source_sha256"].items()):
            raise RuntimeError("production provider or adapter changed during execution")
        # A native evidence index is written for evaluated runs (including FAIL).
        # This indexes the actual MPI histories, separately from the formal PIC
        # gate which has no qualifying particle trajectory receipt.
        dump(output/"native_receipt.json", {"histories": {str(p.relative_to(output)):sha(p)
             for p in sorted(output.glob("mpi-*/*.csv"))}, "scope":manifest["scope"]})
        failures = [c["name"] for c in checks if c["status"]=="FAIL"]
        result.update(status="FAIL" if failures else "PASS",
                      reason=", ".join(failures) if failures else "all source, clock, MPI and product-plan checks met")
    except MissingPrerequisite as error:
        result.update(status="SKIPPED", reason=str(error))
    except (OSError, ValueError, ArithmeticError, KeyError, TypeError, RuntimeError, ImportError, subprocess.SubprocessError) as error:
        result.update(status="ERROR", reason=str(error))
        with stderr.open("a") as stream: stream.write(str(error)+"\n")
    result["duration_s"] = time.monotonic()-start
    result["counts"] = {status:sum(c["status"]==status for c in checks) for status in EXIT}
    dump(output/"run_manifest.json",manifest);dump(output/"result.json",result)
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True, type=Path)
    parser.add_argument("--mass-inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--test-runtime-seconds", type=float, default=3600)
    args=parser.parse_args()
    value=run(args.lock,args.mass_inventory,args.output_dir,args.mpi_ranks,args.test_runtime_seconds)
    print(json.dumps(value,indent=2));sys.exit(EXIT[value["status"]])
