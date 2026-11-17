#!/usr/bin/env python3
"""Validation-only full parallel Vlasov--Maxwell dispersion for isotropic ep plasma.

No production domain limits or physical parameters are supplied implicitly.
The input is an explicit controlled-study manifest, including both temperatures,
densities, rest masses, charge, c, mu0 and kB. The nonrelativistic Maxwellian
model is declared in that manifest; finite electron inertia, displacement current
and the Landau continuation of both cyclotron resonances are retained.

For exp[i(kz-omega*t)] and e_x+i*h*e_y, h=+1 is the spatial Right convention
used by sep_resonant_collision. At k>0, chi_s,h=omega_ps^2 Z(zeta_s)/
(omega*k*v_ts), zeta_s=(omega+h*Omega_s)/(k*v_ts), signed Omega_s.
This is the k_perp=0, isotropic/no-drift reduction of the full Stix tensor;
see the primary NHDS paper https://verscharen.com/papers/NHDS.pdf and author
implementation https://github.com/danielver02/NHDS/blob/master/src/suscep.f90.

Arbitrary-precision entire Z, precision refinement and numerical winding checks
are independent of the production nondispersive MHD/scattering implementation.
Unresolved roots/damping are reported, never replaced with Hall or MHD values.
"""
import argparse
import dataclasses
import hashlib
import json
import os
from pathlib import Path
import tempfile
import platform
from typing import Any, Dict, List, Tuple

import mpmath as mp
from mpmath.libmp.libhyper import NoConvergence

NUMERICAL_EXCEPTIONS = (ZeroDivisionError, OverflowError, ValueError, NoConvergence)


class Unresolved(RuntimeError):
    """The requested numerical branch could not be resolved without a fallback."""


def scalar(value: Any, name: str) -> mp.mpf:
    if not isinstance(value, str):
        raise ValueError(name + " must be an explicit decimal string")
    result = mp.mpf(value)
    if not mp.isfinite(result):
        raise ValueError(name + " must be finite")
    return result


@dataclasses.dataclass
class Plasma:
    values: Dict[str, Any]

    def numbers(self) -> Dict[str, mp.mpf]:
        required = ("proton_density_m3", "electron_density_m3", "magnetic_field_t",
                    "proton_temperature_k", "electron_temperature_k", "proton_mass_kg",
                    "electron_mass_kg", "charge_magnitude_c", "light_speed_mps",
                    "vacuum_permeability_hpm", "boltzmann_jpk")
        allowed = set(required) | {"distribution", "propagation", "provenance",
                                   "proton_drift_mps", "electron_drift_mps", "context"}
        if set(self.values)-allowed:
            raise ValueError("unrecognized equilibrium input: " + ", ".join(sorted(set(self.values)-allowed)))
        values = {name: scalar(self.values[name], name) for name in required}
        if any(value <= 0 for value in values.values()):
            raise ValueError("all physical plasma inputs must be strictly positive")
        if values["proton_density_m3"] != values["electron_density_m3"]:
            raise ValueError("this uniform proton/electron equilibrium requires explicit quasineutral densities")
        if self.values.get("distribution") != "isotropic-nonrelativistic-Maxwellian":
            raise ValueError("distribution must explicitly select isotropic nonrelativistic Maxwellians")
        if self.values.get("propagation") != "parallel" or not isinstance(self.values.get("provenance"), str) or not self.values["provenance"].strip():
            raise ValueError("parallel propagation and controlled-study provenance must be explicit")
        for name in ("proton_drift_mps", "electron_drift_mps"):
            if scalar(self.values[name], name) != 0:
                raise ValueError("the declared independent solver supports zero-drift controlled cases")
        return values

    def scales(self) -> Dict[str, mp.mpf]:
        p = self.numbers()
        n, b = p["proton_density_m3"], p["magnetic_field_t"]
        mi, me = p["proton_mass_kg"], p["electron_mass_kg"]
        mu, kb, q = p["vacuum_permeability_hpm"], p["boltzmann_jpk"], p["charge_magnitude_c"]
        va = b / mp.sqrt(mu*n*mi)
        omega = q*b/mi
        return {"mass_ratio": mi/me, "alfven_ion_mps": va,
                "alfven_total_mass_mps": b/mp.sqrt(mu*n*(mi+me)),
                "alfven_ion_over_c": va/p["light_speed_mps"],
                "ion_inertial_length_m": va/omega, "ion_gyrofrequency_radps": omega,
                "beta_proton": 2*mu*n*kb*p["proton_temperature_k"]/b**2,
                "beta_electron": 2*mu*n*kb*p["electron_temperature_k"]/b**2,
                "proton_thermal_gyroradius_m": mp.sqrt(2*kb*p["proton_temperature_k"]/mi)/omega,
                "electron_thermal_over_c": mp.sqrt(2*kb*p["electron_temperature_k"]/me)/p["light_speed_mps"]}


def plasma_z(zeta: mp.mpc) -> mp.mpc:
    """Entire Landau function, including the lower-half-plane residue."""
    return mp.j*mp.sqrt(mp.pi)*mp.exp(-zeta*zeta)*mp.erfc(-mp.j*zeta)


def z_contour(zeta: mp.mpc, offset: mp.mpf) -> mp.mpc:
    """Independent defining integral along a contour strictly BELOW the pole.

    This verification path uses complex quadrature, not erfc/Faddeeva. No
    finite integration cutoff or principal-value-only damping shortcut is used.
    Two different offsets and precision refinement should agree with plasma_z.
    """
    if offset <= max(mp.mpf(0), -mp.im(zeta)):
        raise ValueError("Landau verification contour must lie strictly below the pole")
    def integrand(x: mp.mpf) -> mp.mpc:
        t = x-mp.j*offset
        return mp.exp(-t*t)/(t-zeta)/mp.sqrt(mp.pi)
    return mp.quad(integrand, [-mp.inf, -8, -4, 0, 4, 8, mp.inf])


class Dispersion:
    def __init__(self, plasma: Plasma):
        self.plasma = plasma

    def evaluate(self, frequency: mp.mpc, kdi: mp.mpf, helicity: int
                 ) -> Tuple[mp.mpc, mp.mpc, mp.mpc]:
        """Return F,F_z,F_K, with F=(omega/omega_pi)^2 D.

        z=omega/Omega_p and K=k*d_p. Multiplication removes the apparent
        omega=0 poles without adding a root at finite k: F(0,K)=-K^2.
        Both thermal species and displacement current are retained exactly.
        """
        if kdi <= 0 or helicity not in (-1, 1):
            raise ValueError("positive k*d_p and an explicit circular helicity are required")
        s = self.plasma.scales()
        mass = s["mass_ratio"]
        a2 = s["alfven_ion_over_c"]**2
        value = frequency*frequency*a2-kdi*kdi
        dz = 2*frequency*a2
        dk = -2*kdi
        for weight, gyro, thermal in ((mp.mpf(1), mp.mpf(1), mp.sqrt(s["beta_proton"])),
                                      (mass, -mass, mp.sqrt(mass*s["beta_electron"]))):
            kv = kdi*thermal
            argument = (frequency+helicity*gyro)/kv
            zz = plasma_z(argument)
            derivative = -2*(1+argument*zz)
            value += weight*frequency*zz/kv
            dz += weight*(zz+frequency*derivative/kv)/kv
            dk -= weight*frequency*(zz+argument*derivative)/(kdi*kv)
        return value, dz, dk

    def weak_damping_energy_ratio(self, real_frequency: mp.mpf, kdi: mp.mpf,
                                 helicity: int) -> mp.mpf:
        """Hermitian dielectric/Brillouin response W_total/E_B at real omega.

        This is a weak-damping wave-energy expression, not an exact conserved
        energy for heavily damped modes. No damping threshold is selected here.
        For a circular phasor it is 1+(omega_r/(kc))^2 Re[d(omega*epsilon_h)/domega].
        The same Hermitian-dielectric definition is used by the primary NHDS
        waveenergy.f90; the derivative below is analytic, not a finite stencil.
        """
        s = self.plasma.scales()
        mass, a2 = s["mass_ratio"], s["alfven_ion_over_c"]**2
        derivative = mp.mpc(1)
        for weight, gyro, thermal in ((mp.mpf(1), mp.mpf(1), mp.sqrt(s["beta_proton"])),
                                      (mass, -mass, mp.sqrt(mass*s["beta_electron"]))):
            kv = kdi*thermal
            argument = (real_frequency+helicity*gyro)/kv
            derivative += weight*(-2*(1+argument*plasma_z(argument)))/(a2*kv*kv)
        return 1+a2*(real_frequency/kdi)**2*mp.re(derivative)


def newton(operator: Dispersion, kdi: mp.mpf, helicity: int, seed: mp.mpc,
           digits: int) -> Tuple[mp.mpc, int]:
    tolerance = mp.power(10, -(digits-12))
    root = mp.mpc(seed)
    for iteration in range(80):
        try:
            value, derivative, _ = operator.evaluate(root, kdi, helicity)
        except NUMERICAL_EXCEPTIONS as error:
            raise Unresolved("complex dispersion evaluation failed: " + str(error)) from error
        if derivative == 0 or not mp.isfinite(value) or not mp.isfinite(derivative):
            raise Unresolved("singular/nonfinite complex dispersion derivative")
        step = value/derivative
        if abs(step) <= tolerance*max(abs(root), kdi):
            return root, iteration+1
        accepted = False
        for _ in range(24):
            trial = root-step
            try:
                smaller = abs(operator.evaluate(trial, kdi, helicity)[0]) < abs(value)
            except NUMERICAL_EXCEPTIONS:
                smaller = False
            if smaller:
                root = trial
                accepted = True
                break
            step /= 2
        if not accepted:
            raise Unresolved("complex Newton residual did not decrease; no model substitution performed")
    raise Unresolved("complex Newton iteration limit reached")


def winding(operator: Dispersion, center: mp.mpc, kdi: mp.mpf, helicity: int,
            radius: mp.mpf, samples: int) -> Tuple[int, mp.mpf, mp.mpf, mp.mpc]:
    offsets = [radius*mp.exp(2*mp.pi*mp.j*i/samples) for i in range(samples)]
    try:
        evaluations = [operator.evaluate(center+offset, kdi, helicity) for offset in offsets]
    except NUMERICAL_EXCEPTIONS as error:
        raise Unresolved("root-isolation evaluation failed: " + str(error)) from error
    values = [value[0] for value in evaluations]
    if any(value == 0 or not mp.isfinite(value) for value in values):
        raise Unresolved("root-isolation contour passes a zero or unresolved value")
    steps = [mp.arg(values[(i+1) % samples]/values[i]) for i in range(samples)]
    turns = mp.fsum(steps)/(2*mp.pi)
    count = int(mp.nint(turns))
    if abs(turns-count) > mp.power(10, -12):
        raise Unresolved("root-isolation winding is not an integer")
    # Independent argument-principle integral, with no modulo-2*pi phase
    # unwrapping. This detects e.g. z^65 aliasing on a64-point phase grid.
    integral = mp.fsum(offset*evaluation[1]/evaluation[0]
                      for offset, evaluation in zip(offsets, evaluations))/samples
    return count, max(abs(step) for step in steps), min(abs(value) for value in values), integral


def isolate(operator: Dispersion, root: mp.mpc, kdi: mp.mpf, helicity: int
            ) -> Dict[str, Any]:
    # These are numerical contour/refinement controls, never physical limits.
    radius = abs(root)/8
    for _ in range(12):
        for samples, refined in ((43, 67), (67, 131), (131, 257), (257, 521)):
            try:
                low = winding(operator, root, kdi, helicity, radius, samples)
                high = winding(operator, root, kdi, helicity, radius, refined)
            except Unresolved:
                break  # Shrink away from a contour zero, without moving the root.
            integral_tolerance = mp.power(10, -min(25, mp.mp.dps//2))
            integral_agrees = abs(low[3]-low[0]) < integral_tolerance and abs(high[3]-high[0]) < integral_tolerance
            resolved = low[0] == high[0] and integral_agrees and high[1] < mp.pi/4
            if resolved and high[0] == 1:
                return {"winding_count": 1, "contour_samples": refined,
                        "independent_log_derivative_count_real": mp.nstr(mp.re(high[3]), mp.mp.dps),
                        "independent_log_derivative_count_imaginary": mp.nstr(mp.im(high[3]), mp.mp.dps),
                        "contour_radius": mp.nstr(radius, mp.mp.dps),
                        "maximum_contour_phase_step": mp.nstr(high[1], mp.mp.dps),
                        "minimum_contour_residual": mp.nstr(high[2], mp.mp.dps)}
            if resolved:
                break  # More than one actual root: seek a smaller isolating disk.
        radius /= 2
    raise Unresolved("could not numerically isolate a single branch root under contour refinement")


def numerical_anchor(plasma: Plasma, kdi: mp.mpf) -> mp.mpf:
    s = plasma.scales()
    # Seed inside the long-wavelength expansion even in high-beta controlled
    # cases. This selects the branch; it never truncates the requested domain.
    thermal = max(mp.mpf(1), s["beta_proton"], s["beta_electron"]/s["mass_ratio"])
    return min(kdi/8, mp.mpf("0.0001")/thermal)


def arithmetic_guard(plasma: Plasma, kdi: mp.mpf) -> int:
    s = plasma.scales()
    anchor = numerical_anchor(plasma, kdi)
    # Direct exp(-zeta^2)*erfc(-i*zeta) may lose ~2 log10|zeta|
    # digits when its enormous opposite exponents cancel. Z'=-2(1+zeta*Z)
    # loses another ~2 log10|zeta|, and neutral species sums ~log10(1/K).
    # Increase working precision instead of declaring a small Newton step
    # converged when its derivative or neutral cancellation is unresolved.
    largest = max(1/(anchor*mp.sqrt(s["beta_proton"])),
                  s["mass_ratio"]/(anchor*mp.sqrt(s["mass_ratio"]*s["beta_electron"])), mp.mpf(1))
    loss = 4*int(mp.ceil(mp.log10(largest)))+int(mp.ceil(max(mp.mpf(0), -mp.log10(anchor))))
    guard = max(0, loss-15)
    if guard > 10000:
        raise Unresolved("required cancellation guard exceeds the explicit numerical working-precision resource bound")
    return guard


def trace_branch(plasma: Plasma, kdi: mp.mpf, helicity: int, digits: int,
                 maximum_log_step: mp.mpf) -> Tuple[mp.mpc, int, int]:
    op = Dispersion(plasma)
    anchor = numerical_anchor(plasma, kdi)
    scales = plasma.scales()
    seed = anchor/mp.sqrt(1+1/scales["mass_ratio"]+scales["alfven_ion_over_c"]**2)
    root, iterations = newton(op, anchor, helicity, seed, digits)
    if mp.re(root) <= 0:
        raise Unresolved("MHD anchor did not select the positive-frequency branch")
    current, step, continuation_steps = anchor, maximum_log_step, 0
    while current < kdi:
        target = min(kdi, current*mp.exp(step))
        try:
            _, dz, dk = op.evaluate(root, current, helicity)
            if dz == 0 or not mp.isfinite(dz) or not mp.isfinite(dk):
                raise Unresolved("continuation predictor derivative is zero or nonfinite")
            predictor = root-dk/dz*(target-current)
            candidate, used = newton(op, target, helicity, predictor, digits)
            prediction_error = abs(candidate-predictor)/abs(root)
            if mp.re(candidate) <= 0 or prediction_error > mp.mpf(1)/3:
                raise Unresolved("continuation left the local positive-frequency branch")
        except (Unresolved,) + NUMERICAL_EXCEPTIONS:
            step /= 2
            if step < mp.mpf("0.000001"):
                raise Unresolved("branch continuation exhausted numerical refinement")
            continue
        root, current = candidate, target
        iterations += used
        continuation_steps += 1
        # Recover efficiency after a difficult interval without relaxing the
        # physical equation, Newton criterion or independent path comparison.
        if used <= 8 and prediction_error < mp.mpf("0.05"):
            step = min(maximum_log_step, step*mp.mpf("1.25"))
        if continuation_steps > 200000:
            raise Unresolved("branch continuation work limit reached")
    return root, iterations, continuation_steps


def solve_point(plasma: Plasma, kdi_text: str, helicity: int, digits: int = 40,
                verification_digits: int = 70) -> Dict[str, Any]:
    if digits < 25 or verification_digits < digits+15:
        raise ValueError("precision study needs >=25 base digits and >=15 additional verification digits")
    with mp.workdps(verification_digits+15):
        kdi = scalar(kdi_text, "kdi")
        if kdi <= 0:
            raise ValueError("positive k*d_p required")
        plasma.numbers()
        guard = arithmetic_guard(plasma, kdi)
        anchor_text = mp.nstr(numerical_anchor(plasma, kdi), 25)
    # Track the selected mode continuously from its analytic MHD limit.
    # A continuation failure remains unresolved rather than jumping to a Hall,
    # cold-plasma, nearby heavily damped, or arbitrary alternative root.
    with mp.workdps(digits+15+guard):
        coarse, iterations, coarse_steps = trace_branch(plasma, kdi, helicity, digits, mp.mpf("0.1"))
        base, refined_iterations, refined_steps = trace_branch(plasma, kdi, helicity, digits, mp.mpf("0.05"))
        path_error = abs(base-coarse)
        if path_error > mp.power(10, -(digits-15))*abs(base):
            raise Unresolved("mode identity did not converge under independent ln(k) continuation refinement")
        iterations += refined_iterations
    with mp.workdps(verification_digits+15+guard):
        op = Dispersion(plasma)
        root, used = newton(op, kdi, helicity, base, verification_digits)
        residual, dz, dk = op.evaluate(root, kdi, helicity)
        correction = abs(residual/dz)
        numerical_error = path_error+abs(root-base)+correction+mp.power(10, -(digits-8))*abs(root)
        if abs(root-base) > mp.power(10, -(digits-15))*abs(root):
            raise Unresolved("complex frequency did not converge under precision refinement")
        isolation = isolate(op, root, kdi, helicity)
        scales = plasma.scales()
        omega = root*scales["ion_gyrofrequency_radps"]
        output_digits = verification_digits+5
        def number(value: mp.mpf) -> str:
            return mp.nstr(value, output_digits)
        real_text, imaginary_text = number(mp.re(omega)), number(mp.im(omega))
        serialized = mp.mpc(real_text, imaginary_text)
        rounding_error = abs(omega-serialized)
        total_error = numerical_error*scales["ion_gyrofrequency_radps"]+rounding_error
        damping_resolved = abs(mp.im(serialized)) > total_error
        if damping_resolved and mp.im(serialized) > 0:
            raise Unresolved("resolved growth contradicts this isotropic zero-drift equilibrium")
        record = {"status": "RESOLVED", "spatial_helicity": "Right" if helicity == 1 else "Left",
                  "positive_phase_speed_branch": True,
                  "numerical_anchor_kdi": anchor_text, "additional_cancellation_guard_digits": guard,
                  "requested_k_di": kdi_text, "serialized_significant_digits": output_digits,
                  "k_di": number(kdi), "k_per_m": number(kdi/scales["ion_inertial_length_m"]),
                  "omega_real_radps": real_text,
                  "omega_imaginary_raw_radps": imaginary_text,
                  "frequency_numerical_error_radps": number(total_error),
                  "frequency_serialization_roundoff_radps": number(rounding_error),
                  "numerical_error_is_convergence_estimate": True,
                  "damping_resolved": damping_resolved,
                  "phase_speed_over_model_va_ion": number(mp.re(root)/kdi),
                  "group_speed_over_model_va_ion": number(mp.re(-dk/dz)),
                  "phase_speed_relative_mhd_error": number(abs(mp.re(root)/kdi-1)),
                  "group_speed_relative_mhd_error": number(abs(mp.re(-dk/dz)-1)),
                  "weak_damping_hermitian_energy_over_magnetic_energy":
                      number(op.weak_damping_energy_ratio(mp.re(root), kdi, helicity)),
                  "damping_over_real_frequency_raw": number(abs(mp.im(root)/mp.re(root))),
                  "newton_correction_relative": number(correction/abs(root)),
                  "precision_refinement_relative": number(abs(root-base)/abs(root)),
                  "continuation_refinement_relative": number(path_error/abs(root)),
                  "continuation_coarse_steps": coarse_steps, "continuation_refined_steps": refined_steps,
                  "newton_iterations": iterations+used,
                  "numerical_isolation": isolation}
        return record


def main() -> int:
    implementation_source = Path(__file__).read_bytes()
    implementation_sha256 = hashlib.sha256(implementation_source).hexdigest()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--kdi", nargs="+", required=True, help="explicit positive k*d_p study points")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--digits", type=int, default=40)
    parser.add_argument("--verification-digits", type=int, default=70)
    parser.add_argument("--helicities", nargs="+", choices=("Right", "Left"), default=("Right", "Left"))
    args = parser.parse_args()
    input_bytes = args.case.read_bytes()
    document = json.loads(input_bytes.decode("utf-8"))
    plasma = Plasma(document)
    with mp.workdps(args.verification_digits+15):
        plasma.numbers()
        for point in args.kdi:
            if scalar(point, "kdi") <= 0:
                raise ValueError("all study k*d_p values must be positive")
        scales = {name: mp.nstr(value, args.verification_digits+5) for name, value in plasma.scales().items()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    archive = args.output.parent / ("warm_parallel_implementation-"+implementation_sha256+".py")
    if archive.exists():
        if archive.read_bytes() != implementation_source:
            raise OSError("content-addressed implementation archive differs from its declared digest")
    else:
        descriptor, temporary = tempfile.mkstemp(prefix=archive.name+".", dir=str(archive.parent))
        try:
            with os.fdopen(descriptor, "wb") as stream:
                stream.write(implementation_source)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, archive)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    rows: List[Dict[str, Any]] = []
    for label in args.helicities:
        helicity = 1 if label == "Right" else -1
        for kdi in args.kdi:
            try:
                row = solve_point(plasma, kdi, helicity, args.digits, args.verification_digits)
            except (Unresolved,) + NUMERICAL_EXCEPTIONS as error:
                row = {"status": "UNRESOLVED", "k_di": kdi,
                       "spatial_helicity": "Right" if helicity == 1 else "Left", "reason": str(error)}
            rows.append(row)
            print("warm-parallel {} kdi={} {}".format(label, kdi, row["status"]), flush=True)
    output = {"schema": "srcsep-controlled-warm-parallel-study-v1", "input": document,
              "mpmath_version": mp.__version__, "base_digits": args.digits,
              "python_version": platform.python_version(),
              "implementation_sha256": implementation_sha256,
              "implementation_archive": archive.name,
              "input_sha256": hashlib.sha256(input_bytes).hexdigest(),
              "requested_k_di": list(args.kdi), "requested_helicities": list(args.helicities),
              "algorithm": "complex Landau-erfc; analytic-derivative Newton; two-continuation-meshes; precision-refinement; phase-plus-log-derivative-winding-v3",
              "verification_digits": args.verification_digits, "scales": scales, "rows": rows,
              "physical_production_qualified": False, "domain_limits_selected": False,
              "physical_error_acceptance_budget": None,
              "root_isolation": "numerical winding/refinement check; no interval-arithmetic proof",
              "counterpropagating_symmetry": "omega(branch=-,h,k)=-conjugate(omega(branch=+,-h,k)); k>0; damping sign retained",
              "limitations": ["linear homogeneous isotropic nonrelativistic Maxwellian proton/electron plasma",
                              "strictly parallel propagation; no amplitude, obliquity or inhomogeneous WKB qualification",
                              "imaginary frequencies below numerical uncertainty are unresolved, not exact zero",
                              "Hermitian dielectric wave-energy ratio is a weak-damping diagnostic; no physical damping acceptance threshold selected",
                              "no observational validation or automatic spectral truncation"]}
    descriptor, temporary = tempfile.mkstemp(prefix=args.output.name+".", dir=str(args.output.parent))
    try:
        with os.fdopen(descriptor, "w") as stream:
            json.dump(output, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, args.output)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return 1 if any(row["status"] != "RESOLVED" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
