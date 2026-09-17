#!/usr/bin/env python3
"""Generate the frozen SHK05 v1 oblique-shock reference header.

This validation-only program is intentionally independent of swcme_shock.hpp.
It solves the complete eight-equation planar ideal-MHD Rankine-Hugoniot system
with Python's standard-library Decimal arithmetic at 80-digit precision.  The
unknowns are rho2, p2, three shock-frame velocity components, and three magnetic
components; no production compression reduction, candidate reconstruction, or
root-bracketing code is imported or translated here.

Run this script from any directory and redirect its stdout to a candidate file.
Review the diff and validation report before replacing the versioned header:

  python3 test/reference/generate_shk05_oblique_v1.py > /tmp/shk05_v1.hpp

The checked-in header is the release fixture.  It is never regenerated during
normal builds, so a production algebra change cannot silently update both sides
of the benchmark.
"""

from decimal import Decimal as D
from decimal import getcontext


# Eighty decimal digits leave more than the requested 50-digit guard margin
# after numerical differentiation and Gaussian elimination.  The output is
# rounded to 17 significant decimal digits because the production API uses
# IEEE binary64; the unrounded residual is retained as fixture metadata.
getcontext().prec = 80
PI = D("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899")
MU0 = D("1.25663706127e-6")
PROTON_MASS = D("1.67262192595e-27")


def sin_decimal(angle):
    """Return sin(angle) without depending on a binary64 math library."""
    angle %= 2 * PI
    if angle > PI:
        angle -= 2 * PI
    term = angle
    result = angle
    order = 1
    while True:
        term *= -angle * angle / D((2 * order) * (2 * order + 1))
        previous = result
        result += term
        if result == previous:
            return result
        order += 1


def cos_decimal(angle):
    """Return cos(angle) through the same high-precision sine series."""
    return sin_decimal(angle + PI / 2)


def vector_square(vector):
    return sum(value * value for value in vector)


def jump_fluxes(rho, pressure, velocity, magnetic, gamma):
    """Evaluate all eight conserved quantities for an x-normal discontinuity."""
    ux, uy, uz = velocity
    bx, by, bz = magnetic
    magnetic_square = vector_square(magnetic)
    velocity_square = vector_square(velocity)
    mass = rho * ux
    momentum = (
        rho * ux * ux + pressure + magnetic_square / 2 - bx * bx,
        rho * ux * uy - bx * by,
        rho * ux * uz - bx * bz,
    )
    # These are the two tangential components of u x B.  Using them directly
    # avoids choosing the same 2x2 reduced system as the production solver.
    electric_tangent = (uz * bx - ux * bz, ux * by - uy * bx)
    energy = ux * (
        D("0.5") * rho * velocity_square
        + gamma / (gamma - 1) * pressure
        + magnetic_square
    ) - bx * sum(velocity[i] * magnetic[i] for i in range(3))
    return (mass, *momentum, *electric_tangent, energy, bx)


def solve_linear(matrix, rhs):
    """Solve an 8x8 Decimal system with independent partial pivoting."""
    size = len(rhs)
    augmented = [matrix[row][:] + [rhs[row]] for row in range(size)]
    minimum_pivot = None
    for column in range(size):
        pivot_row = max(range(column, size),
                        key=lambda row: abs(augmented[row][column]))
        pivot = abs(augmented[pivot_row][column])
        if pivot < D("1e-65"):
            raise ArithmeticError("singular high-precision Newton Jacobian")
        minimum_pivot = pivot if minimum_pivot is None else min(minimum_pivot, pivot)
        augmented[column], augmented[pivot_row] = (
            augmented[pivot_row], augmented[column]
        )
        divisor = augmented[column][column]
        for entry in range(column, size + 1):
            augmented[column][entry] /= divisor
        for row in range(size):
            if row == column:
                continue
            multiplier = augmented[row][column]
            for entry in range(column, size + 1):
                augmented[row][entry] -= multiplier * augmented[column][entry]
    return [augmented[row][size] for row in range(size)], minimum_pivot


def fast_speed(rho, pressure, magnetic, gamma):
    """Compute the x-normal fast characteristic speed in normalized units."""
    magnetic_square = vector_square(magnetic)
    alfven_square = magnetic_square / rho
    sound_square = gamma * pressure / rho
    cosine_square = (magnetic[0] * magnetic[0] / magnetic_square
                     if magnetic_square > 0 else D(0))
    discriminant = ((alfven_square + sound_square) ** 2
                    - 4 * alfven_square * sound_square * cosine_square)
    return ((alfven_square + sound_square + discriminant.sqrt()) / 2).sqrt()


def solve_case(case):
    """Find and cross-check the fast root with seeds plus Mach continuation."""
    gamma = D(case["gamma"])
    beta = D(case["beta"])
    theta = D(case["theta_deg"]) * PI / 180
    phi = D(case["phi_deg"]) * PI / 180
    polarity = D(case["polarity"])
    magnetic1 = (
        polarity * cos_decimal(theta),
        polarity * sin_decimal(theta) * cos_decimal(phi),
        polarity * sin_decimal(theta) * sin_decimal(phi),
    )
    pressure1 = beta / 2
    fast1 = fast_speed(D(1), pressure1, magnetic1, gamma)
    def residual_norm(values):
        return max(abs(value) for value in values)

    def system_at_mach(fast_mach):
        """Build a full-system residual without reusing production reductions."""
        velocity = (
            -fast_mach * fast1,
            D(case["vt_y_alfven"]),
            D(case["vt_z_alfven"]),
        )
        upstream_flux = jump_fluxes(
            D(1), pressure1, velocity, magnetic1, gamma
        )

        def residual(state):
            downstream_flux = jump_fluxes(
                state[0], state[1], state[2:5], state[5:8], gamma
            )
            return [
                downstream_flux[index] - upstream_flux[index]
                for index in range(8)
            ]

        return velocity, residual

    def initial_state(compression, velocity):
        """Construct a deliberately approximate primitive-state Newton seed."""
        bx, by, bz = magnetic1
        magnetic2 = (bx, compression * by, compression * bz)
        pressure2 = (
            pressure1 + velocity[0] ** 2 * (1 - 1 / compression)
            + (by * by + bz * bz
               - magnetic2[1] ** 2 - magnetic2[2] ** 2) / 2
        )
        if pressure2 <= 0:
            pressure2 = D("0.1")
        return [
            compression,
            pressure2,
            velocity[0] / compression,
            velocity[1],
            velocity[2],
            *magnetic2,
        ]

    def newton_solve(state, residual):
        """Solve one complete eight-variable system and report conditioning."""
        solve_minimum_pivot = None
        for _iteration in range(100):
            values = residual(state)
            initial_norm = residual_norm(values)
            if initial_norm < D("1e-55"):
                return state, True, solve_minimum_pivot
            jacobian = [[D(0) for _ in range(8)] for _ in range(8)]
            for column in range(8):
                delta = D("1e-28") * max(D(1), abs(state[column]))
                plus = state[:]
                minus = state[:]
                plus[column] += delta
                minus[column] -= delta
                plus_values = residual(plus)
                minus_values = residual(minus)
                for row in range(8):
                    jacobian[row][column] = (
                        (plus_values[row] - minus_values[row]) / (2 * delta)
                    )
            try:
                correction, pivot = solve_linear(
                    jacobian, [-value for value in values]
                )
            except ArithmeticError:
                return state, False, solve_minimum_pivot
            solve_minimum_pivot = (
                pivot if solve_minimum_pivot is None
                else min(solve_minimum_pivot, pivot)
            )

            # Backtracking is a numerical globalization device only.  It does
            # not encode shock physics or select a branch; admissibility is
            # evaluated independently after the full system converges.
            damping = D(1)
            accepted = False
            while damping > D("1e-12"):
                candidate = [
                    state[index] + damping * correction[index]
                    for index in range(8)
                ]
                if (candidate[0] > 0 and candidate[1] > 0
                        and residual_norm(residual(candidate)) < initial_norm):
                    state = candidate
                    accepted = True
                    break
                damping /= 2
            if not accepted:
                return state, False, solve_minimum_pivot
        return state, False, solve_minimum_pivot

    target_mach = D(case["fast_mach"])
    velocity1, residual = system_at_mach(target_mach)

    roots = []
    converged_seed_count = 0
    physical_seed_count = 0
    minimum_pivot = None
    # Multiple compression guesses expose the trivial root and any additional
    # branches.  Only independently classified evolutionary fast roots enter
    # the frozen fixture; the number of unique admissible roots is recorded.
    for seed_text in ("1.2", "1.5", "2.0", "2.8", "3.5", "3.9"):
        seed = D(seed_text)
        state, converged, pivot = newton_solve(
            initial_state(seed, velocity1), residual
        )
        if pivot is not None:
            minimum_pivot = (
                pivot if minimum_pivot is None else min(minimum_pivot, pivot)
            )

        if not converged:
            continue
        converged_seed_count += 1
        compression = state[0]
        fast2 = fast_speed(state[0], state[1], state[5:8], gamma)
        normal_alfven2 = abs(state[5]) / state[0].sqrt()
        downstream_fast_mach = abs(state[2]) / fast2
        entropy_ratio = ((state[1] / pressure1)
                         / (state[0] ** gamma))
        evolutionary_fast = (
            compression > 1 + D("1e-20")
            and compression < (gamma + 1) / (gamma - 1) + D("1e-20")
            and target_mach > 1
            and downstream_fast_mach < 1
            and abs(state[2]) > normal_alfven2
            and entropy_ratio > 1
        )
        if not evolutionary_fast:
            continue
        physical_seed_count += 1
        if not any(abs(compression - old[0]) < D("1e-30") for old in roots):
            roots.append(state)

    if len(roots) != 1 or physical_seed_count < 2:
        raise RuntimeError(
            f"{case['id']}: expected one well-conditioned fast root; "
            f"found {len(roots)} from {physical_seed_count} seeds"
        )

    root = roots[0]

    # Trace the selected solution away from the target Mach number and back in
    # eight ten-percent steps.  This is a genuine continuation check: every
    # new solve starts from the preceding root, not from a compression formula.
    # Agreement after the round trip demonstrates that the fixture did not
    # jump to a different nonlinear branch despite the successful direct seeds.
    continuation_state = root[:]
    continuation_mach = target_mach
    continuation_step_count = 0
    for scale_index in (*range(1, 5), *range(3, -1, -1)):
        next_mach = target_mach * (D(1) + D(scale_index) / 10)
        next_velocity, next_residual = system_at_mach(next_mach)
        # Mass conservation makes u2n scale approximately with u1n at fixed
        # compression.  Applying only that predictor keeps continuation
        # independent of the production solver's pressure/field reductions.
        continuation_state[2] *= next_velocity[0] / (
            -continuation_mach * fast1
        )
        continuation_state, converged, pivot = newton_solve(
            continuation_state, next_residual
        )
        if not converged:
            raise RuntimeError(
                f"{case['id']}: Mach continuation failed at {next_mach}"
            )
        # Direct-seed and continuation conditioning are intentionally kept
        # separate: the fixture's stored minimum pivot describes the six
        # branch-search solves, while continuation is a pass/fail cross-check.
        continuation_mach = next_mach
        continuation_step_count += 1

    continuation_root_difference = max(
        abs(continuation_state[index] - root[index])
        / max(D(1), abs(root[index]))
        for index in range(8)
    )
    if continuation_root_difference > D("1e-45"):
        raise RuntimeError(
            f"{case['id']}: continuation returned a different root "
            f"({continuation_root_difference})"
        )

    residual_max = residual_norm(residual(root))
    gamma = D(case["gamma"])
    density1 = D(case["number_density_cm3"]) * D("1e6") * PROTON_MASS
    magnetic_scale = D(case["magnetic_nT"]) * D("1e-9")
    alfven_scale = magnetic_scale / (MU0 * density1).sqrt()
    pressure_scale = magnetic_scale * magnetic_scale / MU0
    normal_velocity1 = D(case["upstream_normal_kms"]) * D("1e3")
    shock_speed = normal_velocity1 + D(case["fast_mach"]) * fast1 * alfven_scale

    upstream_velocity_si = (
        normal_velocity1,
        velocity1[1] * alfven_scale,
        velocity1[2] * alfven_scale,
    )
    downstream_velocity_si = (
        shock_speed + root[2] * alfven_scale,
        root[3] * alfven_scale,
        root[4] * alfven_scale,
    )
    upstream_magnetic_si = tuple(value * magnetic_scale for value in magnetic1)
    downstream_magnetic_si = tuple(value * magnetic_scale for value in root[5:8])
    fast2 = fast_speed(root[0], root[1], root[5:8], gamma)
    entropy_ratio = (root[1] / pressure1) / (root[0] ** gamma)

    return {
        "id": case["id"],
        "gamma": gamma,
        "upstream_rho": density1,
        "upstream_pressure": pressure1 * pressure_scale,
        "upstream_velocity": upstream_velocity_si,
        "upstream_magnetic": upstream_magnetic_si,
        "shock_speed": shock_speed,
        "fast_mach": D(case["fast_mach"]),
        "compression": root[0],
        "downstream_rho": root[0] * density1,
        "downstream_pressure": root[1] * pressure_scale,
        "downstream_velocity": downstream_velocity_si,
        "downstream_magnetic": downstream_magnetic_si,
        "entropy_ratio": entropy_ratio,
        "downstream_fast_mach": abs(root[2]) / fast2,
        "reference_residual": residual_max,
        "minimum_newton_pivot": minimum_pivot,
        "converged_seed_count": converged_seed_count,
        "physical_seed_count": physical_seed_count,
        "physical_root_count": len(roots),
    }


# The matrix varies every quantity named by the validation plan while staying
# away from the separately tested near-Mach-one and exactly parallel/
# perpendicular limits.  phi and two tangential-flow components ensure the
# benchmark exercises the full vector system rather than a coplanar shortcut.
CASES = (
    dict(id="low_beta_15deg", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667", beta="0.1", theta_deg="15", phi_deg="0", polarity="1", fast_mach="1.6", vt_y_alfven="0", vt_z_alfven="0", number_density_cm3="3", magnetic_nT="3", upstream_normal_kms="350"),
    dict(id="moderate_30deg", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667", beta="0.5", theta_deg="30", phi_deg="20", polarity="1", fast_mach="2.0", vt_y_alfven="0.15", vt_z_alfven="-0.05", number_density_cm3="5", magnetic_nT="5", upstream_normal_kms="400"),
    dict(id="unit_beta_45deg", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667", beta="1.0", theta_deg="45", phi_deg="45", polarity="1", fast_mach="3.0", vt_y_alfven="0.20", vt_z_alfven="0.10", number_density_cm3="7", magnetic_nT="7", upstream_normal_kms="425"),
    dict(id="high_beta_60deg", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667", beta="2.0", theta_deg="60", phi_deg="70", polarity="1", fast_mach="4.0", vt_y_alfven="-0.10", vt_z_alfven="0.25", number_density_cm3="10", magnetic_nT="8", upstream_normal_kms="300"),
    dict(id="hot_75deg", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667", beta="5.0", theta_deg="75", phi_deg="110", polarity="1", fast_mach="2.5", vt_y_alfven="0.30", vt_z_alfven="-0.20", number_density_cm3="2", magnetic_nT="4", upstream_normal_kms="500"),
    dict(id="gamma_1p4_25deg", gamma="1.4", beta="0.2", theta_deg="25", phi_deg="35", polarity="1", fast_mach="3.5", vt_y_alfven="0.10", vt_z_alfven="0.15", number_density_cm3="4", magnetic_nT="10", upstream_normal_kms="275"),
    dict(id="gamma_1p5_50deg", gamma="1.5", beta="1.5", theta_deg="50", phi_deg="80", polarity="1", fast_mach="1.5", vt_y_alfven="-0.25", vt_z_alfven="0.05", number_density_cm3="12", magnetic_nT="6", upstream_normal_kms="450"),
    dict(id="strong_70deg", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667", beta="0.8", theta_deg="70", phi_deg="125", polarity="1", fast_mach="6.0", vt_y_alfven="0.05", vt_z_alfven="-0.30", number_density_cm3="1.5", magnetic_nT="12", upstream_normal_kms="325"),
    dict(id="negative_polarity", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667", beta="1.0", theta_deg="40", phi_deg="15", polarity="-1", fast_mach="3.0", vt_y_alfven="0.20", vt_z_alfven="-0.10", number_density_cm3="5", magnetic_nT="5", upstream_normal_kms="400"),
    dict(id="dense_35deg", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667", beta="0.7", theta_deg="35", phi_deg="150", polarity="1", fast_mach="2.8", vt_y_alfven="-0.15", vt_z_alfven="-0.20", number_density_cm3="25", magnetic_nT="15", upstream_normal_kms="550"),
    dict(id="weak_field_55deg", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667", beta="3.0", theta_deg="55", phi_deg="210", polarity="1", fast_mach="2.2", vt_y_alfven="0.35", vt_z_alfven="0.20", number_density_cm3="8", magnetic_nT="2", upstream_normal_kms="375"),
    dict(id="gamma_1p6_65deg", gamma="1.6", beta="0.35", theta_deg="65", phi_deg="300", polarity="-1", fast_mach="4.5", vt_y_alfven="-0.05", vt_z_alfven="0.35", number_density_cm3="6", magnetic_nT="9", upstream_normal_kms="225"),
)


def literal(value):
    """Format one Decimal as an auditable binary64 C++ literal."""
    # Values below 1e-300 are numerical zero at the scale of every fixture and
    # may underflow a compiler's unsuffixed decimal-literal conversion.  Emit
    # an exact zero so strict warning builds remain clean; the high-precision
    # residual metadata still records the solve accuracy separately.
    if value == 0 or abs(value) < D("1e-300"):
        return "0.0"
    return f"{value:.17e}"


def vector_literal(vector):
    # Case stores C-style double[3] members, so one brace pair initializes the
    # array.  A second pair would be appropriate for std::array but is rejected
    # by strict aggregate initialization of these deliberately dependency-free
    # fixture records.
    return "{" + ", ".join(literal(value) for value in vector) + "}"


def main():
    results = [solve_case(case) for case in CASES]
    print("#pragma once")
    print()
    print("// Generated by generate_shk05_oblique_v1.py at Decimal precision=80.")
    print("// Do not edit numerical values by hand; create a new fixture version.")
    print("namespace swcme_test { namespace shk05_reference_v1 {")
    print("struct Case {")
    print("  const char* id; double gamma; double upstream_rho_kg_m3;")
    print("  double upstream_pressure_Pa; double upstream_velocity_m_s[3];")
    print("  double upstream_magnetic_T[3]; double shock_speed_m_s;")
    print("  double fast_mach; double compression; double downstream_rho_kg_m3;")
    print("  double downstream_pressure_Pa; double downstream_velocity_m_s[3];")
    print("  double downstream_magnetic_T[3]; double entropy_ratio;")
    print("  double downstream_fast_mach; double reference_max_residual;")
    print("  double minimum_newton_pivot; int converged_seed_count;")
    print("  int physical_seed_count; int physical_root_count; const char* branch;")
    print("};")
    print("inline constexpr Case CASES[] = {")
    for result in results:
        print("  {")
        print(f"    \"{result['id']}\", {literal(result['gamma'])},")
        print(f"    {literal(result['upstream_rho'])}, {literal(result['upstream_pressure'])},")
        print(f"    {vector_literal(result['upstream_velocity'])},")
        print(f"    {vector_literal(result['upstream_magnetic'])}, {literal(result['shock_speed'])},")
        print(f"    {literal(result['fast_mach'])}, {literal(result['compression'])},")
        print(f"    {literal(result['downstream_rho'])}, {literal(result['downstream_pressure'])},")
        print(f"    {vector_literal(result['downstream_velocity'])},")
        print(f"    {vector_literal(result['downstream_magnetic'])},")
        print(f"    {literal(result['entropy_ratio'])}, {literal(result['downstream_fast_mach'])},")
        print(f"    {literal(result['reference_residual'])}, {literal(result['minimum_newton_pivot'])},")
        print(f"    {result['converged_seed_count']}, {result['physical_seed_count']},")
        print(f"    {result['physical_root_count']}, \"EVOLUTIONARY_FAST\"")
        print("  },")
    print("};")
    print("inline constexpr int DECIMAL_PRECISION_DIGITS = 80;")
    print("inline constexpr int FIXTURE_VERSION = 1;")
    print("} }  // namespace swcme_test::shk05_reference_v1")


if __name__ == "__main__":
    main()
