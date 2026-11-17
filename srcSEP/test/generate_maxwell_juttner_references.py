#!/usr/bin/env python3
"""Print independent high-precision rate references; never updates fixtures.

Requires mpmath 1.3.0, uses 90 decimal digits, and compares its special-function
normalization to an original positive rest-volume integral. The printed values
are manually reviewable inputs for MJ01; invoking the production C++ kernel is
deliberately impossible here. This is a reference generator, not a native
source/observational validation program.
"""
import mpmath as mp


def main():
    mp.mp.dps = 90
    if mp.__version__ != "1.3.0":
        raise RuntimeError("Reference provenance pins mpmath 1.3.0; review a version change explicitly")
    temperatures = ("1e-12", "1e-6", "0.01", "0.03125", "0.032", "0.1", "1", "10", "1e6")
    for text in temperatures:
        theta = mp.mpf(text)
        bessel = 2 * theta * (1 + theta) / (mp.exp(1 / theta) * mp.besselk(2, 1 / theta))
        # x=theta*y, scaled original rest radial normalization. Infinite
        # positive integration has no prescribed momentum/temperature cutoff.
        integral = mp.quad(
            lambda y: (1 + theta*y) * mp.sqrt(y*(theta*y+2)) * mp.exp(-y),
            [0, 1, 4, 16, mp.inf],
        )
        from_integral = 2 * mp.sqrt(theta) * (1 + theta) / integral
        relative = abs(from_integral/bessel-1)
        if relative > mp.mpf("1e-70"):
            raise RuntimeError(f"Independent normalization mismatch at theta={text}: {relative}")
        print(f"theta={text} mean_v_over_c={mp.nstr(bessel, 42)} relative_cross_check={mp.nstr(relative, 4)}")


if __name__ == "__main__":
    main()
