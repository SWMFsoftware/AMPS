#!/usr/bin/env python3
"""Independent reference checks of the pitch-angle library (PAR01-PAR11).

The library is driven through build/pitch_angle_driver (one command per line)
and compared with values that do not come from the library:

  PAR01  I(q,H) by 40-digit mpmath quadrature vs the closed form and the
         library quadrature (Section 5.6, Eqs. 27, 37, 38);
  PAR02  phi(eps) by mpmath vs Eq. (39) and the library quadrature;
  PAR03  D and dD/dmu by mpmath (definition and numerical derivative) for all
         four Release A shapes (Eqs. 24, 26, 29, 34, 163, 164);
  PAR04  lambda_parallel by mpmath quadrature of Eq. (23) of the mpmath D
         vs the library's normalized coefficient, including PA-EPREM's
         factor 2 (MFP v1.2 D-2);
  PAR05  finite-volume transmissibility Eq. (178) by mpmath quadrature;
  PAR06  mean-free-path companion bundle F-PA-01/02/04
         (mean_free_path/MEAN_FREE_PATH_MODEL_DATA/benchmark_points.json);
  PAR07  parallel-diffusion companion bundle: particle kinematics, the
         pitch-angle constants and kappa = v lambda/3 for the isotropic law
         (parallel_diffusion/parallel_diffusion_model_data);
  PAR08  mean-free-path bundle F-KIN-01 rigidities;
  PAR09  the v1.2 supplement's independent mathematical checks
         (pitch_angle_diffusion_model_data/check_physics.py, needs mpmath);
  PAR10  the compiled registry equals the specification inventory
         (generate_registry.py --check, specification SHA-256);
  PAR11  the supplement's C++17 interface contract compiles
         (compile_contract.cpp, -Wall -Wextra -Wpedantic -Werror).

Tolerances: closed forms 1e-12 relative; quadrature-class values 1e-9;
library quadrature vs mpmath 1e-12 (it is converged to 1e-14).  Output lines
'[ID] PASS|FAIL|SKIP message'; exit status 1 on any failure.  A missing
dependency (mpmath, a bundle file) is SKIP with the reason, never PASS.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEP_COMMON = HERE.parent
DATA = HERE / "pitch_angle_diffusion_model_data"
MFP_BUNDLE = SEP_COMMON / "mean_free_path" / "MEAN_FREE_PATH_MODEL_DATA" / "benchmark_points.json"
PARALLEL_BUNDLE = (SEP_COMMON / "parallel_diffusion" / "parallel_diffusion_model_data" /
                   "benchmark_points.json")

C = 299792458.0
E = 1.602176634e-19
MP = 1.67262192369e-27
AU = 149597870700.0

try:
    import mpmath as mp
    mp.mp.dps = 40
except ImportError:  # reported as SKIP by the checks that need it
    mp = None

RESULTS: list[tuple[str, str, str]] = []


def report(test_id: str, ok: bool | None, message: str) -> None:
    status = "SKIP" if ok is None else ("PASS" if ok else "FAIL")
    RESULTS.append((test_id, status, message))
    print(f"[{test_id}] {status} {message}", flush=True)


class Driver:
    def __init__(self, path: Path):
        self.process = subprocess.Popen([str(path)], stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, text=True)

    def __call__(self, command: str) -> list[float]:
        self.process.stdin.write(command + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline().strip()
        if line.startswith("error"):
            raise RuntimeError(f"{command!r}: {line}")
        return [float(x) for x in line.split()]

    def close(self) -> None:
        self.process.stdin.close()
        self.process.wait()


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(b), 1e-300)


# ----------------------------------------------------------------------------
def mp_I(q, h):
    """I(q,H) by mpmath's own 2F1 (Eq. 37, second form; exact at H = 0),
    confirmed by an mpmath quadrature after mu = t^(1/(2-q)) (removes the
    |mu|^(1-q) endpoint singularity that a plain mp.quad does not resolve at
    q -> 2).  The two must agree to 1e-25 or the reference is rejected."""
    q, h = mp.mpf(q), mp.mpf(h)
    if h == 0:
        closed = 2 / ((2 - q) * (4 - q))
    else:
        z = 1 / (1 + h)
        closed = (mp.hyp2f1(1, 1, q / (q - 1), z) - mp.hyp2f1(1, 1, (q + 2) / (q - 1), z) / 3) / (1 + h)
    m = 1 / (2 - q)
    quad = mp.quad(lambda t: m * t**(m - 1) * (1 - t**(2 * m)) / (t**(m * (q - 1)) + h), [0, 1])
    if abs(closed - quad) > mp.mpf("1e-25") * abs(closed):
        raise RuntimeError(f"mpmath reference for I({q},{h}) inconsistent: {closed} vs {quad}")
    return closed


def mp_phi(eps):
    eps = mp.mpf(eps)
    a = 1 + eps
    return mp.mpf(3) / 2 * mp.quad(lambda m: (1 - m) * (1 + m)**2 / (a * m + eps),
                                   [0, eps, 10 * eps, 1] if eps < 0.1 else [0, 1])


def mp_h(model, m, q, h, eps):
    a = abs(m)
    if model in ("PA-ISO", "PA-EPREM"):
        return mp.mpf(1)
    if model == "PA-QFORM":
        return a**(q - 1) + h
    return a / (1 + a) + eps


def mp_D(model, d0, m, q, h, eps):
    m = mp.mpf(m)
    return d0 * (1 - m**2) * mp_h(model, m, mp.mpf(q), mp.mpf(h), mp.mpf(eps))


def check_par01(driver):
    if mp is None:
        return report("PAR01", None, "mpmath not installed")
    worst_closed = worst_quad = 0.0
    cases = [(q, h) for q in ("5/3", "1.67", "1.41", "1.2", "1.9") for h in
             ("0", "0.0001", "0.01", "0.05", "0.5", "3", "40")]
    for q, h in cases:
        qf = float(mp.mpf(q)) if "/" not in q else 5.0 / 3.0
        qm = mp.mpf(5) / 3 if q == "5/3" else mp.mpf(q)
        ref = float(mp_I(qm, h))
        closed, quad = driver(f"qform {qf!r} {float(h)!r}")
        worst_closed = max(worst_closed, rel(closed, ref))
        worst_quad = max(worst_quad, rel(quad, ref))
    ok = worst_closed <= 1e-12 and worst_quad <= 1e-12
    report("PAR01", ok, f"I(q,H) {len(cases)} points vs mpmath: closed {worst_closed:.2e}, "
                        f"quadrature {worst_quad:.2e} (budget 1e-12)")


def check_par02(driver):
    if mp is None:
        return report("PAR02", None, "mpmath not installed")
    worst = 0.0
    eps_values = ("1e-12", "1e-6", "0.001", "0.01", "0.045", "0.048", "0.1", "1", "10", "1000")
    for e in eps_values:
        ref = float(mp_phi(e))
        closed, quad = driver(f"phi {float(e)!r}")
        worst = max(worst, rel(closed, ref), rel(quad, ref))
    report("PAR02", worst <= 1e-12, f"phi(eps) {len(eps_values)} points vs mpmath: worst "
                                    f"{worst:.2e} (budget 1e-12)")


def check_par03(driver):
    if mp is None:
        return report("PAR03", None, "mpmath not installed")
    worst_d = worst_dp = 0.0
    shapes = [("PA-ISO", 0, 0, 0), ("PA-EPREM", 0, 0, 0), ("PA-QFORM", 5 / 3, 0.05, 0),
              ("PA-QFORM", 1.3, 0.0, 0), ("PA-EPS", 0, 0, 0.048), ("PA-EPS", 0, 0, 0.0)]
    mus = (-0.999, -0.7, -0.31, -0.02, 0.004, 0.25, 0.6, 0.97, 1.0, -1.0)
    for model, q, h, eps in shapes:
        for m in mus:
            d, dp = driver(f"shape {model} 0.7 {m!r} {q!r} {h!r} {eps!r}")
            ref = mp_D(model, mp.mpf("0.7"), m, q, h, eps)
            worst_d = max(worst_d, abs(d - float(ref)) / 0.7)
            if abs(m) < 1.0:
                refp = mp.diff(lambda x: mp_D(model, mp.mpf("0.7"), x, q, h, eps), mp.mpf(m))
                worst_dp = max(worst_dp, abs(dp - float(refp)) / max(1.0, abs(float(refp))))
    ok = worst_d <= 1e-15 and worst_dp <= 1e-12
    report("PAR03", ok, f"D and dD/dmu for 6 shapes x {len(mus)} mu vs mpmath: "
                        f"{worst_d:.1e} / {worst_dp:.1e}")


def check_par04(driver):
    if mp is None:
        return report("PAR04", None, "mpmath not installed")
    worst = 0.0
    lam = 0.2 * AU
    t = mp.mpf(10) * 10**6 * mp.mpf(E)
    mc2 = mp.mpf(MP) * mp.mpf(C)**2
    p = float(mp.sqrt(t * (t + 2 * mc2)) / mp.mpf(C))
    gamma = mp.sqrt(1 + (mp.mpf(p) / (mp.mpf(MP) * C))**2)
    v = mp.mpf(p) * C**2 / (gamma * mp.mpf(MP) * C**2)
    for model, q, h, eps in [("PA-ISO", 0, 0, 0), ("PA-EPREM", 0, 0, 0),
                             ("PA-QFORM", 5 / 3, 0.05, 0), ("PA-QFORM", 1.5, 0.0, 0),
                             ("PA-EPS", 0, 0, 0.01), ("PA-EPS", 0, 0, 0.7)]:
        lam23, d0, _ = driver(f"lambda {model} {q!r} {h!r} {eps!r} {lam!r} {MP!r} {E!r} {p!r}")
        integral = mp.quad(lambda m: (1 - m**2)**2 / mp_D(model, mp.mpf(d0), m, q, h, eps),
                           [-1, -0.01, 0, 0.01, 1])
        ref = float(3 * v / 8 * integral)
        expected = 2 * lam if model == "PA-EPREM" else lam
        worst = max(worst, rel(lam23, ref), rel(lam23, expected))
    report("PAR04", worst <= 1e-10, f"Eq. (23) of the normalized D (6 shapes, 10 MeV proton) "
                                    f"vs mpmath and the law (EPREM x2): worst {worst:.2e}")


def check_par05(driver):
    if mp is None:
        return report("PAR05", None, "mpmath not installed")
    worst = 0.0
    for model, q, h, eps, edges in [("PA-QFORM", 5 / 3, 0.0, 0, (-0.1, 0.0, 0.1)),
                                    ("PA-QFORM", 1.4, 0.2, 0, (-0.3, 0.1, 0.2)),
                                    ("PA-EPS", 0, 0, 0.048, (-0.05, 0.0, 0.05)),
                                    ("PA-ISO", 0, 0, 0, (0.5, 0.8, 0.9))]:
        (t,) = driver(f"transmissibility {model} {q!r} {h!r} {eps!r} 0.9 {edges[0]!r} "
                      f"{edges[1]!r} {edges[2]!r}")
        a = (edges[0] + edges[1]) / 2
        b = (edges[1] + edges[2]) / 2
        points = [a, 0, b] if a < 0 < b else [a, b]
        ref = 1 / mp.quad(lambda m: 1 / mp_D(model, mp.mpf("0.9"), m, q, h, eps), points)
        worst = max(worst, rel(t, float(ref)))
    report("PAR05", worst <= 1e-9, f"Eq. (178) transmissibility vs mpmath: worst {worst:.2e}")


def check_par06(driver):
    if not MFP_BUNDLE.exists():
        return report("PAR06", None, f"bundle not found: {MFP_BUNDLE}")
    fixtures = {f["id"]: f for f in json.loads(MFP_BUNDLE.read_text())["fixtures"]}
    worst = 0.0
    for h, values in fixtures["F-PA-01"]["outputs"]["values"].items():
        closed, _ = driver(f"qform {5.0 / 3.0!r} {float(h)!r}")
        worst = max(worst, rel(closed, values["I"]), rel(0.75 * closed, values["D0_lambda_over_v"]))
    for e, values in fixtures["F-PA-02"]["outputs"]["values"].items():
        closed, _ = driver(f"phi {float(e)!r}")
        worst = max(worst, rel(closed, values["phi_exact"]))
    lam23, _, _ = driver(f"lambda PA-EPREM 0 0 0 1e10 {MP!r} {E!r} 1e-19")
    worst = max(worst, rel(lam23 / 1e10, fixtures["F-PA-04"]["outputs"]["lambda_effective_over_lambda"]))
    report("PAR06", worst <= 1e-12, f"MFP bundle F-PA-01/02/04: worst {worst:.2e} (budget 1e-12)")


def check_par07(driver):
    if not PARALLEL_BUNDLE.exists():
        return report("PAR07", None, f"bundle not found: {PARALLEL_BUNDLE}")
    data = json.loads(PARALLEL_BUNDLE.read_text())
    worst = 0.0
    for name, particle in data["particles"].items():
        mass = float(particle["mass_kg"])
        charge = particle["charge_number"] * E
        p = particle["momentum_kg_m_s"]
        gamma, speed, rigidity, _ = driver(f"kinematics {mass!r} {charge!r} {p!r}")
        worst = max(worst, rel(gamma, particle["gamma"]), rel(speed, particle["speed_m_s"]),
                    rel(rigidity, particle["rigidity_V"]))
        lam23, _, _ = driver(f"lambda PA-ISO 0 0 0 {particle['lambda_parallel_m']!r} "
                             f"{mass!r} {charge!r} {p!r}")
        worst = max(worst, rel(speed * lam23 / 3.0, particle["kappa_parallel_m2_s"]))
    constants = data["spectrum_and_pitch_angle"]
    closed, _ = driver(f"qform {5.0 / 3.0!r} 0.01")
    # The PARALLEL bundle integrates over [-1,1] (twice Section 5's I).
    worst = max(worst, rel(2 * closed, constants["I_q5_3_h0_01"]),
                rel(0.75 * closed, constants["D0_lambda_over_v_q5_3_h0_01"]))
    closed, _ = driver(f"qform {5.0 / 3.0!r} 0")
    worst = max(worst, rel(2 * closed, constants["I_q5_3_h0_exact_Eq24"]))
    report("PAR07", worst <= 1e-12, f"PARALLEL bundle kinematics, kappa = v lambda/3, "
                                    f"pitch-angle constants: worst {worst:.2e}")


def check_par08(driver):
    if not MFP_BUNDLE.exists() or mp is None:
        return report("PAR08", None, "bundle or mpmath not available")
    out = {f["id"]: f for f in json.loads(MFP_BUNDLE.read_text())["fixtures"]}["F-KIN-01"]["outputs"]
    worst = 0.0
    for key, mass, charge, t_mev in [("proton_1MeV_R_MV", MP, E, "1"),
                                     ("electron_0.94MeV_R_MV", 9.1093837015e-31, -E, "0.94"),
                                     ("electron_2.0MeV_R_MV", 9.1093837015e-31, -E, "2.0")]:
        t = mp.mpf(t_mev) * 10**6 * mp.mpf(E)
        mc2 = mp.mpf(mass) * mp.mpf(C)**2
        p = float(mp.sqrt(t * (t + 2 * mc2)) / C)
        _, _, rigidity, _ = driver(f"kinematics {mass!r} {charge!r} {p!r}")
        worst = max(worst, rel(rigidity / 1e6, out[key]))
    report("PAR08", worst <= 1e-12, f"MFP bundle F-KIN-01 rigidities: worst {worst:.2e}")


def check_par09():
    if mp is None:
        return report("PAR09", None, "mpmath not installed (python3 -m pip install --user mpmath)")
    output = HERE / "build" / "check_physics_results.json"
    run = subprocess.run([sys.executable, str(DATA / "check_physics.py"), "--output", str(output)],
                         capture_output=True, text=True)
    if run.returncode != 0 or not output.exists():
        return report("PAR09", False, f"check_physics.py exit {run.returncode}: {run.stdout[-300:]}")
    result = json.loads(output.read_text())
    report("PAR09", result["failed"] == 0,
           f"supplement check_physics.py: {result['passed']}/{result['checks']} passed")


def check_par10():
    run = subprocess.run([sys.executable, str(HERE / "generate_registry.py"), "--check"],
                         capture_output=True, text=True)
    report("PAR10", run.returncode == 0, run.stdout.strip() or run.stderr.strip())


def check_par11():
    run = subprocess.run(["g++", "-std=c++17", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
                          "-fsyntax-only", str(DATA / "compile_contract.cpp")],
                         capture_output=True, text=True, cwd=DATA)
    report("PAR11", run.returncode == 0,
           "compile_contract.cpp compiles (C++17, warnings as errors)" if run.returncode == 0
           else run.stderr[-300:])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--driver", default=str(HERE / "build" / "pitch_angle_driver"))
    args = parser.parse_args()
    driver = Driver(Path(args.driver))
    try:
        for check in (check_par01, check_par02, check_par03, check_par04, check_par05,
                      check_par06, check_par07, check_par08):
            try:
                check(driver)
            except Exception as error:  # a driver error is a failure, not a skip
                report(check.__name__.replace("check_", "").upper(), False, str(error))
    finally:
        driver.close()
    check_par09()
    check_par10()
    check_par11()
    failed = sum(status == "FAIL" for _, status, _ in RESULTS)
    passed = sum(status == "PASS" for _, status, _ in RESULTS)
    skipped = sum(status == "SKIP" for _, status, _ in RESULTS)
    print(f"pitch-angle reference summary: {passed} passed, {failed} failed, {skipped} skipped")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
