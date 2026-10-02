#!/usr/bin/env python3
"""Audit the frozen SSE3D08 weak-flank compression with 80-digit fluxes.

Run from any directory with Python's standard library. This script neither
imports SWCME nor evaluates its multiplied/deflated compression polynomial.
It reconstructs downstream states from direct mass, tangential momentum,
induction and normal momentum conditions, then bisects TOTAL ENERGY FLUX.
It does not modify a fixture or any application source.

Inputs below are the binary64 upstream state, normal and normal surface speed
at the reported negative-Z ray, t=180 s, from the complete SSE example.
The positive-Z case has identical compression by reflection symmetry.
"""
from decimal import Decimal, localcontext


def binary64(value):
    """Preserve the supplied double exactly before the high-precision solve."""
    return Decimal.from_float(float(value))


def verify():
    with localcontext() as context:
        context.prec = 80
        dot = lambda a, b: sum(x * y for x, y in zip(a, b))
        scale = lambda a, factor: [x * factor for x in a]
        add = lambda a, b: [x + y for x, y in zip(a, b)]
        subtract = lambda a, b: [x - y for x, y in zip(a, b)]
        rho = binary64("1.2544509358396501e-18")
        pressure = binary64("1.0354739485627494e-09")
        surface_speed = binary64("720637.90072562301")
        gamma = binary64("1.6666666666666667")
        mu0 = binary64("1.2566370612700001e-06")
        normal = list(map(binary64, ("0.54183421148481936",
            "-0.83664771383675862", "-0.08022649310763437")))
        velocity = list(map(binary64, ("371301.35830900434",
            "-148100.04152198954", "-14201.373844574337")))
        magnetic = list(map(binary64, ("4.42138828602069e-07",
            "-2.2104459636040566e-07", "-1.7498930221389042e-08")))
        # Restore an exactly unit vector in Decimal; the logged normalized
        # binary64 vector has an O(epsilon) norm error. The C++ comparison uses
        # a 5e-13 absolute compression tolerance, far above this tiny effect.
        normal = scale(normal, 1 / dot(normal, normal).sqrt())
        u1 = subtract(velocity, scale(normal, surface_speed))
        un = dot(u1, normal)
        ut = subtract(u1, scale(normal, un))
        bn = dot(magnetic, normal)
        bt = subtract(magnetic, scale(normal, bn))

        def energy_flux(density, thermal_pressure, speed, field):
            # Direct ideal-MHD total-energy flux normal to the surface in
            # the normal-moving shock frame. No divided energy residual.
            return dot(speed, normal) * (
                density * dot(speed, speed) / 2
                + gamma * thermal_pressure / (gamma - 1)
                + dot(field, field) / mu0
            ) - dot(field, normal) * dot(speed, field) / mu0

        upstream_flux = energy_flux(rho, pressure, u1, magnetic)

        def residual(compression):
            # This closed tangential solution follows the two linear RH
            # conditions. Normal velocity/density follow mass conservation;
            # pressure follows normal momentum. Energy is tested separately.
            b2t = scale(bt, compression * (rho * un * un - bn * bn / mu0)
                / (rho * un * un - compression * bn * bn / mu0))
            u2t = add(ut, scale(subtract(b2t, bt), bn / (mu0 * rho * un)))
            p2 = pressure + rho * un * un * (1 - 1 / compression)
            p2 += (dot(bt, bt) - dot(b2t, b2t)) / (2 * mu0)
            u2 = add(scale(normal, un / compression), u2t)
            b2 = add(scale(normal, bn), b2t)
            return energy_flux(compression * rho, p2, u2, b2) - upstream_flux

        # Isolate the nontrivial weak root above unity. Eighty-digit arithmetic
        # separates its energy residual from that of the trivial r=1 state.
        lower, upper = Decimal(1) + Decimal("1e-12"), Decimal("1.00001")
        lower_residual = residual(lower)
        assert lower_residual * residual(upper) < 0, "reference bracket lost"
        for _ in range(240):
            middle = (lower + upper) / 2
            middle_residual = residual(middle)
            if middle_residual * lower_residual > 0:
                lower, lower_residual = middle, middle_residual
            else:
                upper = middle
        compression = (lower + upper) / 2
        expected = Decimal("1.0000011572267355924125831759717077013695476763206588150089224851523829659872940")
        assert abs(compression - expected) < Decimal("1e-70"), "frozen root changed"
        assert abs(residual(compression)) < Decimal("1e-65"), "energy flux does not close"
        print("PASS: independent 80-digit conserved-flux weak-flank reference")
        print("compression=" + str(compression))
        print("energy_flux_residual=" + str(residual(compression)))


if __name__ == "__main__":
    verify()
