# Stage 10: field-line export and `srcSEP` integration

Stage 10 establishes a non-lossy, versioned exchange between a prepared 3-D
coronal/heliospheric state and one or more field-aligned `srcSEP` simulations.
The exchange is divided into three ownership layers:

1. `sep_common` owns neutral records, validation, bounded interpolation, and
   transactional bundle I/O.
2. `sep_coronal_cme` owns both-sign tracing, Sun-branch selection, outward arc
   length, finite-tube measures, and front intersections.
3. `srcSEP3D` and `srcSEP` contain only thin export/import adapters. They do not
   copy shared headers or implement a second tracing/source model.

## Neutral state and topology

`NodeState` stores position, outward tangent, magnetic and plasma vectors,
density, pressure, temperature, focusing length, directional wave energies,
tube area, region, both topology classifications, magnetic sector, source
label, longitude-map Jacobians, mapping validity, and interface identity in SI
units.  Arc length is strictly increasing outward and is independent of the
sign of **B**.

Interpolation is bounded and one-sided. Positive primitives use log-linear
interpolation; geometry and continuous signed vectors use linear interpolation
with tangent renormalization. Region, topology, sector, source label, mapping
validity, and interface identity must match at both ends and are copied rather
than interpolated. A line crossing an unresolved `R_i`, `R_scs`, separatrix,
HCS, or transition clearance is rejected instead of being truncated into an
apparently complete Sun-to-observer line.

History coverage, intersection presence, and source-measure status are
orthogonal enums. `FrontIntersection` likewise keeps geometric connection,
source eligibility, rejection-bit set, front generation, and time. The
non-lossy `LineConnectionSummary` reports all root counts, first/last
geometric and source-active times, the union of rejection causes, and the
evaluated horizon.

## Tracing and orientation

`TraceSunConnectedLine` integrates both signs of

`dx/dell = +/- B/|B|`

from an explicit Cartesian/observer seed using RK4. Sphere crossings are
event-located on the final segment. Exactly one branch must reach the solar
surface and the other the configured outer sphere. The solar branch is
reversed, joined to the outer branch without duplicating the seed, and
reparameterized by outward geometric arc length. Thus inward magnetic
polarity never reverses the 1-D coordinate.

The exported line has one stable string ID. Request ordering, vector position,
rank ownership, and adaptive point count are not physical identity.

## Absolute tube and quadrature measures

For a finite traced footprint carrying unsigned magnetic flux `Phi`, every
node uses

`A(s)|B(s)| = Phi`.

The quadrature weight derives from the same positive flux, not an unrelated
user scalar. Multi-line sets validate unique IDs and positive measures and are
canonically sorted before identity calculation. Characteristic-only lines are
typed separately and cannot be used for count-based observer sampling.

`FindFrontIntersections` evaluates a prepared implicit front on every line
segment, brackets sign-changing roots, retains exactly resolved tangent-root
creation, and deduplicates a root shared by adjacent segments. Stable root IDs
bind the line, front generation, and segment lineage. Source eligibility is a
separate field and never overwrites geometric existence.

## Observer exchange

Each repeated observer mapping retains its stable observer ID, exact energy
edges, detector-frame declaration, and cadence-resolved overlap components.
Components store stable IDs, topology-constant time intervals, closest arc
coordinate, separation, represented volume, and integrated exposure. A
component ends at tangency/creation/annihilation; interpolation through such a
topology event is prohibited. `ObserverMappingCoverage` must span the run
horizon.

The same finite-tube volume/exposure therefore normalizes 3-D and 1-D
residence estimators. Empty-bin validity remains explicit through the Stage-8
finite-output contract.

## Transactional bundle and provenance

`WriteBundleTransactional` validates the complete `FieldLineSet`, serializes
lines in stable-ID order, writes one tabular member per line, computes SHA-256
for every member, and writes a canonical JSON manifest. Publication occurs by
renaming the complete sibling temporary directory. An existing target or any
member/manifest failure leaves no partial new bundle.

The canonical bundle identity covers schema, frame, epoch, history interval,
resolved rotation provenance, and all line/member content. Import verifies
the schema, safe serialized IDs, expected member names, every checksum, every
record, and the recomputed identity before particle allocation. Changes to
geometry, physical state, front history, observer mapping, or provenance
therefore change restart identity.

## Application adapters

`srcSEP3D/adapters/field_line_export_adapter` gathers parsed requests and an
owning immutable-state evaluator, calls the shared tracing/reduction API, and
publishes the returned bundle. It contains no tracing or source physics.

`srcSEP/adapters/field_line_bundle_adapter` reads and validates the immutable
bundle, calls an application-provided publisher for the existing background
snapshot/coefficient registry, and only then replaces its owning active
handle. Failed publication preserves the prior bundle. Evaluation selects an
exact stable line ID and performs bounded shared interpolation; there is no
nearest-line fallback. The baseline `srcSEP` adapter links `sep_common` only.

Repeated line meshes call `ResampleLine(begin,end,count)`. Requested endpoints
and point count are exact; positive primitives remain positive, and any range
outside the bundle or interval crossing an unresolved discontinuity fails.
`ParseLineMeshRequests` consumes the normalized `srcSEP` assignments and
requires this complete repeated record for every selected line:

```ini
[line_mesh.earth]
line_id = line-A
start_s_m = 6.957e8
end_s_m = 1.495978707e11
point_count = 4097
resampling = conservative-positive-one-sided
```

The parser rejects unknown or missing keys, unit-bearing/non-finite numbers,
reversed intervals, fewer than two points, and more than one mesh for the same
physical line. `BuildLineMeshes` then resolves the exact serialized line ID
and rejects absent lines and any requested extrapolation before particle
allocation.

## Verification

`make test-stage10` runs 196 cumulative tests. Stage 10 adds `HCS3D07`,
`FLX3D01--09`, `FLX1D01--03`, `OBS1D01--02`, `NAT1D01`, `TIM1D01`,
`POP1D01--02`, `RUN1D01`, `RST1D01`, and `XM3D01--02`. They verify radial
both-sign tracing, outward orientation, flux conservation, bounded one-sided
interpolation, checksummed transactional round trip and corruption rejection,
request-order determinism, ordinary/tangent front roots, source/geometric
history separation, observer coverage, 1-D resampling, finite values,
population conservation, restart identity mutation, and Parker/focused
3-D-to-1-D state parity.
