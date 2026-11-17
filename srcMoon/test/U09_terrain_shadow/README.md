# U09 — Topographic illumination and horizon shadowing

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U09_terrain_shadow/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

The linked production probe also covers the four analytical branches of the
I17 uniform-disc eclipse kernel: no overlap, totality, a concentric smaller
occulter, and partial overlap of equal radii. This establishes the local circle
geometry only. It does not verify SPICE wiring, terrestrial-atmosphere light,
terrain/eclipsing-body composition, thermal response, or Diviner agreement;
those remain linked I17 responsibilities.

The November 12 extension adds closed-surface radial rays through shared
vertices of a translated/scaled icosphere, in addition to planar edge/vertex
and distant-coordinate cases. The independent oracle is geometric: a ray
aimed through a closed convex surface vertex must hit, at its prescribed
radius. A pre-fix run missed 49 vertices. Production BVH slab and barycentric
bounds now include floating-point arithmetic enclosures scaled by coordinate
magnitude and triangle conditioning. These are numerical roundoff bounds,
not an adjustable physical clearance or relaxed acceptance threshold.

The M2 extension also runs `psd-illumination` through the compiled shared-core
photon-stimulated desorption kernel. Its geometry-only callback invokes the
production finite-disc and terrain BVH service; the actual PSD kernel retains
the physical population, compiled effective cross section and 1-AU photon
flux. It excludes thermal Bond albedo. Analytical fixtures verify full and
half projected incidence, a back-facing facet, actual-point inverse-square
dilution, opaque-Earth totality and ridge umbra. Half projected incidence uses
a complete disc above a facet with central cosine 0.5; it does not confuse
visible-disc area with photon irradiance. This is local source/geometry
wiring coverage, not an observational determination of Na source strength.

The same probe verifies that a NULL callback retains the original expression
and cached-total policy, zero physical stock makes no geometry call, and an
enabled callback refreshes total and per-area maximum from current physical
stocks. Refresh must preserve pending adsorption/desorption counters and
the reservoir itself. Geometry/formula comparisons require `2e-13` relative
agreement; exact zero, legacy arithmetic and counter checks use zero tolerance.
The fixture never edits compiled source coefficients. A build without positive
Na PSD coefficients reports `SKIPPED` for that capability.

`psd-illumination-invalid` is a separate expected-error probe: with positive
stock it supplies a nonfinite geometry factor and must terminate through
the AMPS diagnostic. It is deliberately absent from the normal passing probe
list. For an already built production-kernel probe executable, run:

```bash
/absolute/path/to/production_kernel_probe psd-illumination
/absolute/path/to/production_kernel_probe psd-illumination-invalid
```

The second command must fail with the message that the direct photon-flux
factor must be finite and nonnegative. No real observation or production
campaign is executed by either command. The full source entry-point, source
interval timing and cache contract is documented in
[PSD_ILLUMINATION_CONTRACT.md](../PSD_ILLUMINATION_CONTRACT.md).
