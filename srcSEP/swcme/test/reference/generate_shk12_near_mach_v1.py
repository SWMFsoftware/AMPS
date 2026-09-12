#!/usr/bin/env python3
"""Generate the frozen SHK12 v1 near-Mach-one reference fixture.

This validation-only generator uses the complete eight-variable ideal-MHD
Rankine-Hugoniot residual from the independent SHK05 Decimal implementation.
It traces the weak evolutionary-fast branch from Mfast-1=0.5 down to 1e-5
with 80-digit Newton continuation.  No production compression scan, candidate
reconstruction, branch-selection rule, or binary64 result is used.

Normal builds consume the reviewed header and never execute this program.
Regenerate a candidate explicitly and review its diff before versioning it:

  python3 test/reference/generate_shk12_near_mach_v1.py > /tmp/shk12_v1.hpp
"""

from generate_shk05_oblique_v1 import (
    D,
    MU0,
    PI,
    PROTON_MASS,
    cos_decimal,
    fast_speed,
    jump_fluxes,
    literal,
    sin_decimal,
    solve_linear,
    vector_literal,
)


# These families cover four decades of plasma beta, four genuinely oblique
# geometries, two gamma values, non-coplanar B, and nonzero tangential flow.
# The separate runtime sweep also includes the more singular one-degree cases;
# they are deliberately not frozen as well-conditioned numerical references.
FAMILIES = (
    dict(id="b0p01_t30_g1p4", beta="0.01", theta_deg="30", gamma="1.4"),
    dict(id="b0p01_t60_g5p3", beta="0.01", theta_deg="60", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667"),
    dict(id="b0p1_t89_g1p4", beta="0.1", theta_deg="89", gamma="1.4"),
    dict(id="b0p1_t30_g5p3", beta="0.1", theta_deg="30", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667"),
    dict(id="b1_t15_g1p4", beta="1", theta_deg="15", gamma="1.4"),
    dict(id="b1_t60_g5p3", beta="1", theta_deg="60", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667"),
    dict(id="b10_t30_g1p4", beta="10", theta_deg="30", gamma="1.4"),
    dict(id="b10_t89_g5p3", beta="10", theta_deg="89", gamma="1.6666666666666666666666666666666666666666666666666666666666666666666666666666667"),
)

# Small logarithmic ratios make Newton remain on the same weak branch as the
# compression approaches one.  The emitted subset keeps the header compact;
# all intermediate continuation points are still mandatory generation checks.
CONTINUATION_EXCESSES = tuple(
    D(value) for value in (
        "0.5", "0.2", "0.1", "0.05", "0.02", "0.01", "0.005",
        "0.002", "0.001", "0.0005", "0.0002", "0.0001",
        "0.00005", "0.00002", "0.00001",
    )
)
EMITTED_EXCESSES = {
    D("0.5"), D("0.1"), D("0.01"), D("0.001"),
    D("0.0001"), D("0.00001"),
}


def residual_norm(values):
    """Return the infinity norm without converting Decimal values to float."""
    return max(abs(value) for value in values)


def newton_solve(state, residual):
    """Solve all eight jump equations with a numerical full Decimal Jacobian."""
    minimum_pivot = None
    for iteration in range(100):
        values = residual(state)
        old_norm = residual_norm(values)
        if old_norm < D("1e-55"):
            return state, True, iteration, old_norm, minimum_pivot

        # Central differences are evaluated at far higher precision than the
        # binary64 production calculation.  The full Jacobian is intentionally
        # unrelated to the production solver's one-dimensional reduction.
        jacobian = [[D(0) for _ in range(8)] for _ in range(8)]
        for column in range(8):
            step = D("1e-28") * max(D(1), abs(state[column]))
            plus = state[:]
            minus = state[:]
            plus[column] += step
            minus[column] -= step
            plus_values = residual(plus)
            minus_values = residual(minus)
            for row in range(8):
                jacobian[row][column] = (
                    (plus_values[row] - minus_values[row]) / (2 * step)
                )

        try:
            correction, pivot = solve_linear(
                jacobian, [-value for value in values]
            )
        except ArithmeticError:
            return state, False, iteration, old_norm, minimum_pivot
        minimum_pivot = (
            pivot if minimum_pivot is None else min(minimum_pivot, pivot)
        )

        # Positivity-preserving backtracking globalizes Newton without encoding
        # any production bracket or branch preference.
        damping = D(1)
        while damping > D("1e-14"):
            candidate = [
                state[index] + damping * correction[index]
                for index in range(8)
            ]
            if (candidate[0] > 0 and candidate[1] > 0
                    and residual_norm(residual(candidate)) < old_norm):
                state = candidate
                break
            damping /= 2
        else:
            return state, False, iteration, old_norm, minimum_pivot

    return state, False, 100, residual_norm(residual(state)), minimum_pivot


def solve_family(family):
    """Trace one independently classified weak fast-shock family."""
    gamma = D(family["gamma"])
    beta = D(family["beta"])
    theta = D(family["theta_deg"]) * PI / 180
    azimuth = D("25") * PI / 180
    magnetic1 = (
        cos_decimal(theta),
        sin_decimal(theta) * cos_decimal(azimuth),
        sin_decimal(theta) * sin_decimal(azimuth),
    )
    pressure1 = beta / 2
    fast1 = fast_speed(D(1), pressure1, magnetic1, gamma)
    tangential_velocity = (D("0.2"), D("-0.1"))

    def system(excess):
        velocity1 = (
            -(D(1) + excess) * fast1,
            tangential_velocity[0],
            tangential_velocity[1],
        )
        flux1 = jump_fluxes(D(1), pressure1, velocity1, magnetic1, gamma)

        def residual(state):
            flux2 = jump_fluxes(
                state[0], state[1], state[2:5], state[5:8], gamma
            )
            return [flux2[index] - flux1[index] for index in range(8)]

        return velocity1, residual

    def seed(compression, velocity1):
        # This intentionally rough MHD seed is used only at the strong end of
        # the family.  Every weaker point starts from the preceding full root.
        bx, by, bz = magnetic1
        magnetic2 = (bx, compression * by, compression * bz)
        pressure2 = (
            pressure1 + velocity1[0] ** 2 * (1 - 1 / compression)
            + (by * by + bz * bz
               - magnetic2[1] ** 2 - magnetic2[2] ** 2) / 2
        )
        if pressure2 <= 0:
            pressure2 = D("0.1")
        return [
            compression, pressure2, velocity1[0] / compression,
            velocity1[1], velocity1[2], *magnetic2,
        ]

    first_excess = CONTINUATION_EXCESSES[0]
    velocity1, residual = system(first_excess)
    roots = []
    minimum_pivot = None
    # Multiple starting compressions prove that the initial continuation root
    # is not an artifact of one Newton guess.  The evolutionary-fast classifier
    # is based only on independent characteristic and entropy conditions.
    for guess in (D("1.1"), D("1.3"), D("1.7"), D("2.3"), D("3.0"), D("4.0")):
        root, converged, _, _, pivot = newton_solve(seed(guess, velocity1), residual)
        if not converged:
            continue
        if pivot is not None:
            minimum_pivot = pivot if minimum_pivot is None else min(minimum_pivot, pivot)
        fast2 = fast_speed(root[0], root[1], root[5:8], gamma)
        normal_alfven2 = abs(root[5]) / root[0].sqrt()
        entropy_ratio = (root[1] / pressure1) / (root[0] ** gamma)
        evolutionary_fast = (
            root[0] > 1 + D("1e-30")
            and abs(root[2]) / fast2 < 1
            and abs(root[2]) > normal_alfven2
            and entropy_ratio > 1
        )
        if evolutionary_fast and not any(
            abs(root[0] - old[0]) < D("1e-28") for old in roots
        ):
            roots.append(root)

    if len(roots) != 1:
        raise RuntimeError(
            f"{family['id']}: expected one evolutionary root, found {len(roots)}"
        )

    state = roots[0]
    previous_excess = first_excess
    results = []
    for excess in CONTINUATION_EXCESSES:
        velocity1, residual = system(excess)
        if excess != first_excess:
            # Mass continuity supplies only a normal-velocity predictor.  All
            # other primitive components are obtained by the independent full
            # solve, not by production reconstruction formulas.
            state[2] *= (D(1) + excess) / (D(1) + previous_excess)
            state, converged, iterations, residual_max, pivot = newton_solve(
                state, residual
            )
            if not converged:
                raise RuntimeError(
                    f"{family['id']}: continuation failed at excess={excess}"
                )
            if pivot is not None:
                minimum_pivot = pivot if minimum_pivot is None else min(minimum_pivot, pivot)
        else:
            iterations = 0
            residual_max = residual_norm(residual(state))

        if not (state[0] > 1 and residual_max < D("1e-50")):
            raise RuntimeError(
                f"{family['id']}: invalid weak root at excess={excess}"
            )

        if excess in EMITTED_EXCESSES:
            density1 = D("5e6") * PROTON_MASS
            magnetic_scale = D("5e-9")
            alfven_scale = magnetic_scale / (MU0 * density1).sqrt()
            pressure_scale = magnetic_scale * magnetic_scale / MU0
            normal_velocity1 = D("400e3")
            shock_speed = (
                normal_velocity1 + (D(1) + excess) * fast1 * alfven_scale
            )
            results.append(dict(
                id=f"{family['id']}_dm_{str(excess).replace('.', 'p')}",
                beta=beta,
                theta_deg=D(family["theta_deg"]),
                gamma=gamma,
                mach_excess=excess,
                upstream_rho=density1,
                upstream_pressure=pressure1 * pressure_scale,
                upstream_velocity=(
                    normal_velocity1,
                    tangential_velocity[0] * alfven_scale,
                    tangential_velocity[1] * alfven_scale,
                ),
                upstream_magnetic=tuple(value * magnetic_scale for value in magnetic1),
                shock_speed=shock_speed,
                compression=state[0],
                downstream_rho=state[0] * density1,
                downstream_pressure=state[1] * pressure_scale,
                downstream_velocity=(
                    shock_speed + state[2] * alfven_scale,
                    state[3] * alfven_scale,
                    state[4] * alfven_scale,
                ),
                downstream_magnetic=tuple(value * magnetic_scale for value in state[5:8]),
                reference_residual=residual_max,
                continuation_iterations=iterations,
                minimum_pivot=minimum_pivot,
            ))
        previous_excess = excess

    return results


def main():
    """Emit a dependency-free C++17 aggregate fixture to standard output."""
    results = [result for family in FAMILIES for result in solve_family(family)]
    print("#pragma once")
    print()
    print("// Generated by generate_shk12_near_mach_v1.py at Decimal precision=80.")
    print("// Do not edit numerical values by hand; create a new fixture version.")
    print("namespace swcme_test { namespace shk12_reference_v1 {")
    print("struct Case {")
    print("  const char* id; double beta; double theta_deg; double gamma;")
    print("  double mach_excess; double upstream_rho_kg_m3;")
    print("  double upstream_pressure_Pa; double upstream_velocity_m_s[3];")
    print("  double upstream_magnetic_T[3]; double shock_speed_m_s;")
    print("  double compression; double downstream_rho_kg_m3;")
    print("  double downstream_pressure_Pa; double downstream_velocity_m_s[3];")
    print("  double downstream_magnetic_T[3]; double reference_max_residual;")
    print("  int continuation_iterations; double minimum_newton_pivot;")
    print("};")
    print("inline constexpr Case CASES[] = {")
    for result in results:
        print("  {")
        print(f"    \"{result['id']}\", {literal(result['beta'])}, {literal(result['theta_deg'])}, {literal(result['gamma'])},")
        print(f"    {literal(result['mach_excess'])}, {literal(result['upstream_rho'])}, {literal(result['upstream_pressure'])},")
        print(f"    {vector_literal(result['upstream_velocity'])},")
        print(f"    {vector_literal(result['upstream_magnetic'])}, {literal(result['shock_speed'])},")
        print(f"    {literal(result['compression'])}, {literal(result['downstream_rho'])}, {literal(result['downstream_pressure'])},")
        print(f"    {vector_literal(result['downstream_velocity'])},")
        print(f"    {vector_literal(result['downstream_magnetic'])}, {literal(result['reference_residual'])},")
        print(f"    {result['continuation_iterations']}, {literal(result['minimum_pivot'])}")
        print("  },")
    print("};")
    print("inline constexpr int DECIMAL_PRECISION_DIGITS = 80;")
    print("inline constexpr int FIXTURE_VERSION = 1;")
    print("} }  // namespace swcme_test::shk12_reference_v1")


if __name__ == "__main__":
    main()
