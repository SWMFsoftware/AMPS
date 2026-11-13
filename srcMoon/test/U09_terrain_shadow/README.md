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
