#!/usr/bin/env python3
"""Exercise the I24 total-capture production provider and audit its native plans independently.

Approved 2026-11-25: H+ + Na(3s) -> Na+ + H(n), n = 1, 2, 3, at 1-5 keV/u;
Morgan 1985 total rate, Shingal 1986 shell split (log-log partials), NIST ASD
shell centroids, zero-angle COM exchange and radiated energy booked at the
event. See TOTAL_CAPTURE_CONTRACT.md. The oracle below is written in Python
Decimal from the printed source values; it shares no code with the C++ path.
"""
from __future__ import annotations
import argparse
import csv
from decimal import Context, Decimal as D, localcontext
import functools
from hashlib import sha256
import hmac
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import time



def precise(function):
    """Run the oracle at 50 Decimal digits without changing the caller's context."""
    @functools.wraps(function)
    def wrapper(*args, **kwargs):
        with localcontext(Context(prec=50)):
            return function(*args, **kwargs)
    return wrapper

HERE = Path(__file__).resolve().parent
TEST = HERE.parents[1]
ROOT = TEST.parents[1]
CASE = HERE / "reference/total_capture_case.json"
sys.path.insert(0, str(TEST))  # srcMoon/test: shared common/ helpers
from common.python_environment import describe as python_environment, require as require_python
EXIT = {"PASS": 0, "FAIL": 1, "ERROR": 2, "SKIPPED": 77}
LIGHT_SPEED = D(299792458)


class MissingPrerequisite(FileNotFoundError):
    """Only a missing external input/tool can skip the case; lost run output is an error."""


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def keyed_uniform(seed, cohort, parent, purpose, index):
    """RFC 2104 HMAC-SHA256 reproduction of the production keyed open uniform.

    Purposes: 1 source clock, 2 channel, 3 product clock, 4 kinematics. A
    source parent has an empty genealogy. The 52 leading digest bits map to
    the binary64 midpoint grid exactly.
    """
    u64 = lambda n: struct.pack(">Q", n)
    label = cohort.encode("utf-8")
    message = (b"AMPS-MOON-M3-HMAC-SHA256-COUNTER-v1" + u64(len(label)) + label +
               u64(parent) + u64(0) + u64(purpose) + u64(index))
    bits = int.from_bytes(hmac.new(u64(seed), message, "sha256").digest()[:7], "big") >> 4
    return (bits + .5) * 2. ** -52


@precise
def loglog(x, x0, x1, y0, y1):
    t = (x / x0).ln() / (x1 / x0).ln()
    return (y0.ln() + t * (y1.ln() - y0.ln())).exp()


class Oracle:
    """Independent Decimal reference for rates, shells and two-body products."""

    @precise
    def __init__(self, case, inventory, table_csv, levels_csv):
        neutral = {x["identity"]: D(x["reference_neutral_mass_kg"]) for x in inventory["isotopes"]}
        constants = {k: D(v["value"]) for k, v in inventory["constants"].items()}
        ie = {x["element"]: D(x["ionization_energy_eV"]) for x in inventory["atomic_ground_thresholds"]}
        self.u, self.e = constants["atomic mass constant"], constants["elementary charge"]
        electron = constants["electron mass"]
        self.h, self.na = neutral["1H"], neutral["23Na"]
        self.electron = electron
        self.proton, self.ion = self.h - electron, self.na - electron
        self.q1 = (ie["H"] - ie["Na"]) * self.e
        self.ie_h = ie["H"] * self.e
        self.inputs = [neutral["1H"], neutral["23Na"], electron, self.u, self.e, ie["H"], ie["Na"]]
        knots = case["morgan_oracle"]
        self.morgan = [(D(E) * 1000, D(s) * D("1e-20")) for E, s in zip(knots["energy_kev_u"], knots["sigma_1e-16_cm2"])]
        with table_csv.open() as stream:
            rows = list(csv.DictReader(stream))
        self.shingal = [(D(r["energy_kev_u"]) * 1000,
                         [D(r["sigma_1s_A2"]), D(r["sigma_2s_A2"]) + D(r["sigma_2p_A2"]),
                          D(r["sigma_3s_A2"]) + D(r["sigma_3p_A2"]) + D(r["sigma_3d_A2"])]) for r in rows]
        weighted, weights = [D(0)] * 3, [D(0)] * 3
        with levels_csv.open() as stream:
            levels = list(csv.DictReader(stream))
        for r in levels:
            n, g = int(r["n"]), D(r["g"])
            weighted[n - 1] += g * D(r["level_eV"]); weights[n - 1] += g
        self.excitation = [w / g * self.e for w, g in zip(weighted, weights)]
        self.theta = D(case["theta_max_com_deg"])

    @precise
    def speed(self, energy_ev):
        return (2 * D(energy_ev) * self.e / self.u).sqrt()

    @precise
    def sigma(self, energy):
        for (e0, s0), (e1, s1) in zip(self.morgan, self.morgan[1:]):
            if e0 <= energy <= e1:
                return s0 if energy == e0 else s1 if energy == e1 else loglog(energy, e0, e1, s0, s1)
        raise ValueError("oracle energy outside the Morgan knots")

    @precise
    def fractions(self, energy):
        for (e0, p0), (e1, p1) in zip(self.shingal, self.shingal[1:]):
            if e0 <= energy <= e1:
                part = p0 if energy == e0 else p1 if energy == e1 else [loglog(energy, e0, e1, a, b) for a, b in zip(p0, p1)]
                total = sum(part)
                return [x / total for x in part]
        raise ValueError("oracle energy outside the Shingal rows")

    @precise
    def products(self, g, shell):
        """Zero-angle exchange for a stationary Na target and beam along (0.6, 0.8, 0)."""
        mass = self.proton + self.na
        mu_i, mu_f = self.proton * self.na / mass, self.ion * self.h / mass
        q = self.q1 - self.excitation[shell - 1]
        available = mu_i * g * g / 2 + q
        pf = (2 * mu_f * available).sqrt()
        unit = (D("0.6"), D("0.8"), D(0))
        vcm = [self.proton * g * c / mass for c in unit]
        hydrogen = [v + pf / self.h * c for v, c in zip(vcm, unit)]
        ion = [v - pf / self.ion * c for v, c in zip(vcm, unit)]
        impulse = 2 * pf * D(math.sin(math.radians(float(self.theta)) / 2))
        return dict(q=q, pf=pf, hydrogen=hydrogen, ion=ion, impulse=impulse,
                    recoil=impulse * impulse / (2 * self.ion),
                    photon=self.excitation[shell - 1] / (self.h * LIGHT_SPEED),
                    bath=-(self.proton * g * g / 2 + self.ie_h) + self.excitation[shell - 1])


def binomial_two_sided(k, n, p):
    """Exact two-sided binomial p-value (sum of outcomes no more likely than k)."""
    if n == 0:
        return 1.
    logp = lambda j: (math.lgamma(n + 1) - math.lgamma(j + 1) - math.lgamma(n - j + 1) +
                      (j * math.log(p) if j else 0.) + ((n - j) * math.log1p(-p) if n - j else 0.))
    reference = logp(k) + 1e-12
    return min(1., sum(math.exp(logp(j)) for j in range(n + 1) if logp(j) <= reference))


@precise
def prepare_inputs(case, output):
    """Authenticate sources, regenerate the pinned CSVs and build CLI inputs."""
    for name, entry in case["sources"].items():
        path = ROOT / entry["path"]
        if name == "morgan1985_directory":
            path = path / "nist_morgan_1985_retry2.pdf"
        if not path.exists():
            raise MissingPrerequisite("required external input is unavailable: " + str(path))
        if sha(path) != entry["sha256"]:
            raise ValueError("changed source identity: " + str(path))
    inventory_path = ROOT / case["sources"]["mass_inventory"]["path"]
    inventory = json.loads(inventory_path.read_text())
    extractor = TEST / "physical_state_sources/prepare.py"
    spec = importlib.util.spec_from_file_location("i24_mass_evidence", extractor)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    raw = inventory_path.parent / "raw"
    for filename in module.SOURCES:
        if not (raw / filename).is_file():
            raise MissingPrerequisite("required NIST source packet member is unavailable: " + str(raw / filename))
    if module.extract(raw) != inventory:
        raise ValueError("NIST isotope/energy/constant native-source replay differs from the retained inventory")
    generator = HERE / "prepare_total_capture.py"
    spec = importlib.util.spec_from_file_location("i24_total_generator", generator)
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    try:
        outputs = gen.generate(case)
    except gen.NotApproved as error:
        raise MissingPrerequisite(str(error))
    derived = output / "derived"
    derived.mkdir()
    for name, data in outputs.items():
        (derived / name).write_bytes(data)
        if data != (HERE / "reference" / name).read_bytes() or sha256(data).hexdigest() != case["derived_files"][name]:
            raise ValueError("regenerated %s differs from the frozen reviewed bytes" % name)
    hashes = {str(p): sha(p) for p in (CASE, generator, extractor, inventory_path, HERE / "reference/total_capture_approval_20261125.json",
                                       *(raw / name for name in module.SOURCES), *derived.iterdir())}
    oracle = Oracle(case, inventory, derived / "shingal_table4.csv", derived / "hydrogen_levels.csv")
    sources = [str(ROOT / case["sources"]["morgan1985_directory"]["path"]), str(derived / "shingal_table4.csv"),
               str(derived / "hydrogen_levels.csv"), *map(str, oracle.inputs)]
    return sources, oracle, hashes


@precise
def audit(directory, ranks, case, oracle, energy_ev, cohort):
    """Score one energy's native histories against the independent oracle."""
    ensemble = case["ensemble"]
    count = ensemble["parents"]
    rows = {}
    for rank in range(ranks):
        with (directory / ("parents-rank-%d.csv" % rank)).open() as stream:
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
        native, = csv.DictReader(stream)
    g = oracle.speed(energy_ev)
    sigma = oracle.sigma(D(energy_ev))
    fractions = oracle.fractions(D(energy_ev))
    rate = D(ensemble["proton_density_m3"]) * sigma * g
    horizon = D(ensemble["optical_depth"]) / rate
    tol = case["arithmetic_relative_tolerance"]
    checks = []

    def check(name, passed, actual, reference):
        checks.append(dict(name=name, status="PASS" if passed else "FAIL", actual=actual, reference=reference))

    rel = lambda a, b: abs(D(a) / D(b) - 1)
    check("direct n sigma g [1/s] (contract <= 1%)", rel(native["rate_s"], rate) <= D("0.01"), native["rate_s"], str(rate))
    for field, ref in (("rate_s", rate), ("sigma_m2", sigma), ("speed_m_s", g), ("energy_ev_per_u", D(energy_ev))):
        check("arithmetic " + field, rel(native[field], ref) <= D(tol), native[field], str(ref))
    for n in range(3):
        check("shell fraction f%d" % (n + 1), abs(D(native["f%d" % (n + 1)]) - fractions[n]) <= D(tol),
              native["f%d" % (n + 1)], str(fractions[n]))
    for n in (2, 3):
        check("NIST shell centroid E_exc(%d) [J]" % n, rel(native["excitation%d_J" % n], oracle.excitation[n - 1]) <= D(tol),
              native["excitation%d_J" % n], str(oracle.excitation[n - 1]))
    for guard in ("low_energy_refused", "high_energy_refused", "negative_density_refused",
                  "nonfinite_density_refused", "zero_density_zero_rate"):
        check(guard, native[guard] == "1", native[guard], "1")
    frequency = float(rate)
    hor = float(horizon)
    wrong_fates = wrong_shells = wrong_ids = 0
    worst = dict(time=0., hydrogen=0., ion_recoil=0., q=0., radiated=0., bath=0., impulse=0., recoil=0., photon=0.)
    shells = [0, 0, 0]
    times = []
    reference_products = {n: oracle.products(g, n) for n in (1, 2, 3)}
    momentum_scale = float(oracle.proton * g)
    electron_scale = float(oracle.electron * g)
    for i in range(count):
        row = rows[i]
        t_first = -math.log(keyed_uniform(ensemble["seed"], cohort, i, 1, 0)) / frequency
        hit = row["fate"] == "capture"
        if row["fate"] not in ("capture", "survived"):
            raise ValueError("unknown native fate")
        wrong_fates += (t_first <= hor) != hit
        worst["time"] = max(worst["time"], abs(float(row["time_s"]) - min(t_first, hor)) / hor)
        if not hit:
            continue
        times.append(float(row["time_s"]) / hor)
        u = keyed_uniform(ensemble["seed"], cohort, i, 4, 0)
        expected = 1 if u < fractions[0] else 2 if u < fractions[0] + fractions[1] else 3
        shell = int(row["shell"])
        wrong_shells += shell != expected
        shells[shell - 1] += 1
        ref = reference_products[shell]
        wrong_ids += (row["identity_ok"] != "1" or
                      any(D(row[k]) != D(ensemble["parent_weight"]) for k in ("na_weight", "h_weight", "bath_weight")))
        h = [D(row["h" + j]) for j in "xyz"]
        ion = [D(row["na" + j]) for j in "xyz"]
        diff = lambda a, b: math.sqrt(sum(float(x - y) ** 2 for x, y in zip(a, b)))
        # H on the incoming momentum scale. The small longitudinal Na+ recoil
        # (terms -Q/v and -v/2 nearly cancel) on the electron-transfer scale m_e g.
        worst["hydrogen"] = max(worst["hydrogen"], float(oracle.h) * diff(h, ref["hydrogen"]) / momentum_scale)
        worst["ion_recoil"] = max(worst["ion_recoil"], float(oracle.ion) * diff(ion, ref["ion"]) / electron_scale)
        represented = D(row["bath_weight"])
        for key, native_value, reference in (("q", row["Q_J"], ref["q"]),
                                             ("radiated", row["radiated_J"], oracle.excitation[shell - 1]),
                                             ("bath", D(row["bath_energy_change_J"]) / represented, ref["bath"]),
                                             ("impulse", row["impulse_bound_kg_m_s"], ref["impulse"]),
                                             ("recoil", row["na_recoil_bound_J"], ref["recoil"]),
                                             ("photon", row["photon_recoil_m_s"], ref["photon"])):
            if reference == 0:
                worst[key] = max(worst[key], float(abs(D(native_value))))
            else:
                worst[key] = max(worst[key], float(rel(native_value, reference)))
    captures = sum(shells)
    check("all parent first-hit classifications", wrong_fates == 0, wrong_fates, 0)
    check("all event-time errors / horizon", worst["time"] <= tol, worst["time"], tol)
    survive = (count - captures) / count
    lifetime = -hor / math.log(survive) if 0 < survive < 1 else None
    check("end-point lifetime [s] (contract <= 3%)", lifetime is not None and abs(lifetime * frequency - 1) <= .03,
          lifetime, 1 / frequency)
    check("per-event shell equals oracle selection", captures > 0 and wrong_shells == 0, wrong_shells, 0)
    check("product identity, charge and represented weights", captures > 0 and wrong_ids == 0, wrong_ids, 0)
    for key, label in (("hydrogen", "H velocity vs zero-angle oracle (momentum scale)"),
                       ("q", "Q(n) relative error"), ("radiated", "radiated energy E_exc(n) relative error"),
                       ("bath", "bath energy change per event relative error"),
                       ("impulse", "omitted impulse bound relative error"),
                       ("recoil", "Na+ recoil-energy bound relative error"),
                       ("photon", "photon-recoil speed bound relative error")):
        check(label, captures > 0 and worst[key] <= tol, worst[key], tol)
    check("Na+ recoil momentum vs oracle / (m_e g)", captures > 0 and worst["ion_recoil"] <= case["recoil_relative_tolerance"],
          worst["ion_recoil"], case["recoil_relative_tolerance"])
    if captures:
        with (directory / "event-controls.csv").open() as stream:
            guard, = csv.DictReader(stream)
        for name in ("closed_endpoints_refused_prior_plan_preserved", "rescaled_daughter_weights_invariant"):
            check(name, guard[name] == "1", guard[name], "1")
    stats = dict(captures=captures, shells=shells, fractions=[float(f) for f in fractions], times=sorted(times),
                 sigma=float(sigma), rate=frequency, native_sigma=float(native["sigma_m2"]))
    return rows, checks, stats


def plot(output, results, oracle):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    grid = [j / 100 for j in range(101)]
    for label, stats in results.items():
        n = stats["parents"]
        t = stats["times"]
        axes[0].step([0] + t + [1], [1] + [1 - (j + 1) / n for j in range(len(t))] + [1 - len(t) / n],
                     where="post", label=label + " keV/u")
    axes[0].plot(grid, [math.exp(-v) for v in grid], "k--", label="exp(−t/τ)")
    axes[0].set(xlabel="Elapsed time / horizon", ylabel="Surviving Na fraction")
    energies = [1000 * (1 + 4 * j / 200) for j in range(201)]
    for n, style in ((1, ":"), (2, "-"), (3, "--")):
        axes[1].plot([e / 1000 for e in energies], [float(oracle.fractions(D(e))[n - 1]) for e in energies],
                     "k" + style, label="Shingal log-log, n=%d" % n)
        xs, ys, err = [], [], []
        for label, stats in results.items():
            k, total = stats["shells"][n - 1], stats["captures"]
            p = k / total
            xs.append(float(label)); ys.append(p); err.append(1.96 * math.sqrt(max(p * (1 - p), 1 / total) / total))
        axes[1].errorbar(xs, ys, yerr=err, fmt="o", capsize=3, label="native n=%d" % n)
    axes[1].set(yscale="log", xlabel="Projectile energy [keV/u]", ylabel="Shell fraction")
    axes[2].plot([e / 1000 for e in energies], [float(oracle.sigma(D(e))) * 1e4 for e in energies], "k-",
                 label="Morgan log-log")
    axes[2].plot([float(e) / 1000 for e, _ in oracle.morgan], [float(s) * 1e4 for _, s in oracle.morgan], "ks",
                 label="Morgan knots")
    axes[2].plot([float(k) for k in results], [r["native_sigma"] * 1e4 for r in results.values()], "o",
                 label="native provider")
    axes[2].set(xlim=(0.9, 5.1), xlabel="Projectile energy [keV/u]", ylabel="σ_tot [cm²]")
    for ax in axes:
        ax.legend(framealpha=1, fontsize=7)
    for suffix in ("png", "eps"):
        fig.savefig(output / ("comparison." + suffix))
    plt.close(fig)


@precise
def run(output, ranks=8, wall_seconds=3600):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    case = json.loads(CASE.read_text())
    ensemble = case["ensemble"]
    manifest = dict(schema="srcMoon-I24-total-capture-manifest-v1", scope="compiled MPI provider and event planner",
                    source_revision=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                    seed=ensemble["seed"], parents_per_energy=ensemble["parents"],
                    mpi_layouts=list(dict.fromkeys((1, ranks))), openmp_threads=1, commands=[],
                    model=case["model_id"], formal_status="SKIPPED / NOT QUALIFIED",
                    python=python_environment())
    checks = []
    result = dict(status="ERROR", scope=manifest["scope"], formal_status=manifest["formal_status"], checks=checks)
    stdout, stderr = output / "stdout.log", output / "stderr.log"
    stdout.touch(); stderr.touch()
    try:
        if not 0 < wall_seconds <= 3600 or not math.isfinite(wall_seconds) or ranks < 1:
            raise ValueError("positive MPI ranks and a total wall allowance of at most 3600 s are required")

        def invoke(argv, directory):
            remaining = wall_seconds - (time.monotonic() - start)
            if remaining <= 0:
                raise TimeoutError("I24 total time allowance exhausted")
            argv = list(map(str, argv))
            manifest["commands"].append({"cwd": str(directory), "argv": argv})
            with stdout.open("a") as out, stderr.open("a") as err:
                out.write(json.dumps(argv) + "\n"); out.flush()
                completed = subprocess.run(argv, cwd=directory, stdout=out, stderr=err, timeout=remaining,
                                           env=dict(os.environ, OMP_NUM_THREADS="1",
                                                    OMPI_MCA_hwloc_base_binding_policy="none"))
            manifest["commands"][-1]["returncode"] = completed.returncode
            if completed.returncode:
                raise RuntimeError("native command exited %d (see %s)" % (completed.returncode, stderr))
        sources, oracle, hashes = prepare_inputs(case, output)
        # After input authentication (a corrupt input stays ERROR) and before any
        # compilation or MPI launch: a missing package is SKIPPED, not a late ERROR.
        require_python(("matplotlib",), MissingPrerequisite)  # PNG/EPS figures
        manifest["input_sha256"] = hashes
        shutil.copyfile(CASE, output / "case.json")
        mpi_compiler, launcher = shutil.which("mpicxx"), shutil.which("mpirun")
        if not mpi_compiler or not launcher:
            raise MissingPrerequisite("mpicxx and mpirun required for the MPI provider test")
        manifest["compiler"] = subprocess.check_output([mpi_compiler, "--version"], text=True).splitlines()[0]
        manifest["mpi_launcher"] = subprocess.check_output([launcher, "--version"], text=True).splitlines()[0]
        species, shared = ROOT / "srcMoon/species", ROOT / "src/models/chemical_database"
        units = [HERE / "total_capture_probe.cpp",
                 *(species / (n + ".cpp") for n in ("SodiumTotalCapture", "ReducedKinematics", "QualificationRandom",
                                                    "ReactionExecution", "PiecewiseHazard")),
                 *(shared / (n + ".cpp") for n in ("Runtime", "DataIO", "ChargeExchange"))]
        binary = output / "total-capture-probe"
        flags = ["-std=c++17", "-O2", "-g", "-Wall", "-Wextra", "-pedantic", "-DOMPI_SKIP_MPICXX", "-I", str(ROOT)]
        dep = subprocess.run([mpi_compiler, *flags, "-MM", *map(str, units)], cwd=ROOT,
                             capture_output=True, text=True, check=True, timeout=60)
        (output / "compiler-dependencies.txt").write_text(dep.stdout)
        dependencies = set()
        for line in dep.stdout.replace("\\\n", " ").splitlines():
            dependencies.update((ROOT / p).resolve() for p in line.split(":", 1)[1].split())
        dependencies.update((HERE / "total_capture_setup.h", Path(__file__), HERE / "prepare_total_capture.py"))
        manifest["source_sha256"] = {str(p): sha(p) for p in sorted(dependencies)}
        invoke([mpi_compiler, *flags, *units, "-o", binary], ROOT)
        manifest["binary_sha256"] = sha(binary)
        low, high = (oracle.speed(e) for e in ensemble["refusal_energies_ev_per_u"])
        results = {}
        for energy in ensemble["energies_kev_u"]:
            energy_ev = D(energy) * 1000
            cohort = "%s-%skeVu" % (ensemble["cohort_prefix"], energy)
            g = oracle.speed(energy_ev)
            horizon = D(ensemble["optical_depth"]) / (D(ensemble["proton_density_m3"]) * oracle.sigma(energy_ev) * g)
            first = None
            for size in manifest["mpi_layouts"]:
                directory = output / ("E%s-mpi-%d" % (energy, size))
                directory.mkdir()
                invoke([launcher, "-np", size, binary, *sources, ensemble["proton_density_m3"], g, horizon,
                        ensemble["seed"], cohort, ensemble["parents"], directory, ensemble["parent_weight"],
                        *ensemble["daughter_weights"], low, high], directory)
                rows, metrics, stats = audit(directory, size, case, oracle, energy_ev, cohort)
                checks.extend(dict(c, energy_kev_u=energy, mpi_ranks=size) for c in metrics)
                if first is None:
                    first = rows
                    results[energy] = dict(stats, parents=ensemble["parents"])
                else:
                    differences = sum(row != first[i] for i, row in rows.items())
                    checks.append(dict(name="per-parent rank invariance", status="PASS" if not differences else "FAIL",
                                       actual=differences, reference=0, energy_kev_u=energy, mpi_ranks=size))
        # Branch frequencies: exact binomial tests, Bonferroni over all
        # shells and energies for the declared 1% family false-alarm rate.
        tests = 3 * len(results)
        alpha = case["stochastic_family_alpha"] / tests
        for energy, stats in results.items():
            for n in range(3):
                p_value = binomial_two_sided(stats["shells"][n], stats["captures"], stats["fractions"][n])
                checks.append(dict(name="shell n=%d frequency binomial p-value" % (n + 1),
                                   status="PASS" if p_value > alpha else "FAIL", actual=p_value, reference=alpha,
                                   energy_kev_u=energy, mpi_ranks=manifest["mpi_layouts"][0],
                                   counts="%d/%d" % (stats["shells"][n], stats["captures"])))
        plot(output, results, oracle)
        with (output / "comparison.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=("energy_kev_u", "mpi_ranks", "name", "status", "actual", "reference", "counts"))
            writer.writeheader()
            writer.writerows({k: c.get(k, "") for k in writer.fieldnames} for c in checks)
        if any(sha(Path(path)) != identity for path, identity in manifest["source_sha256"].items()):
            raise RuntimeError("production provider or adapter changed during execution")
        dump(output / "native_receipt.json", {"histories": {str(p.relative_to(output)): sha(p)
             for p in sorted(output.glob("E*-mpi-*/*.csv"))}, "scope": manifest["scope"]})
        failures = sorted({c["name"] for c in checks if c["status"] == "FAIL"})
        result.update(status="FAIL" if failures else "PASS",
                      reason=", ".join(failures) if failures else "all rate, shell, kinematic, ledger, statistical and MPI checks met")
    except MissingPrerequisite as error:
        result.update(status="SKIPPED", reason=str(error))
    except (OSError, ValueError, ArithmeticError, KeyError, TypeError, RuntimeError, ImportError,
            subprocess.SubprocessError) as error:
        result.update(status="ERROR", reason=str(error))
        with stderr.open("a") as stream:
            stream.write(str(error) + "\n")
    result["duration_s"] = time.monotonic() - start
    result["counts"] = {status: sum(c["status"] == status for c in checks) for status in EXIT}
    dump(output / "run_manifest.json", manifest)
    dump(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--test-runtime-seconds", type=float, default=3600)
    args = parser.parse_args()
    value = run(args.output_dir, args.mpi_ranks, args.test_runtime_seconds)
    print(json.dumps({k: v for k, v in value.items() if k != "checks"}, indent=2))
    sys.exit(EXIT[value["status"]])
