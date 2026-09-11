# Fix 12 — shock-surface mesh topology and quality

## Purpose

This update removes the topological defects in the original shock-surface
triangulation and establishes one production area-weighting path for future SEP
source injection.

The previous `build_shock_mesh()` treated the polar surface as a rectangular
`(theta,phi)` array with `(nTheta+1)*(nPhi+1)` vertices.  That representation is
not a valid unique-vertex mesh on a polar surface:

1. every `phi` value at `theta=0` represents the same physical apex;
2. `phi=0` and `phi=2*pi` represent the same physical seam point on every ring;
3. for a closed full surface, all `phi` values at `theta=pi` represent the same
   rear pole;
4. triangulating those repeated vertices creates zero-area polar cells and
   duplicate seam topology.

Those defects cannot be repaired correctly by simply discarding zero-area
triangles after meshing.  Filtering leaves duplicated vertices, ambiguous edge
incidence, and incorrect node/source accounting.  The topology therefore has to
be correct when the mesh is constructed.

## New topology

`ShockMesh` is now an explicit triangular manifold.

### Finite SSE cap

For `nTheta` polar intervals and `nPhi` azimuthal samples, the cap contains:

- one unique apex;
- `nTheta` non-polar rings;
- exactly `nPhi` unique vertices on each ring;
- the last ring at the configured half width `lambda`;
- no `phi=2*pi` duplicate.

The node and cell counts are

```text
Nv = 1 + nTheta*nPhi
Ne = nPhi*(2*nTheta - 1)
```

The first `nPhi` cells form the apex fan.  Each adjacent ring pair contributes
`2*nPhi` cells.  The final ring is the one physical open boundary of the cap.

### Sphere and ellipsoid

Closed Sphere/Ellipsoid surfaces contain:

- one unique apex/north pole;
- `nTheta-1` non-polar rings;
- one unique rear/south pole;
- no duplicate periodic-seam vertices.

Their counts are

```text
Nv = 2 + (nTheta - 1)*nPhi
Ne = 2*nPhi*(nTheta - 1)
```

The first and last rings are connected to the two unique poles by triangle fans.
All other rings use two triangles per angular quadrilateral.

## Periodic seam

A ring contains only

```text
phi_k = 2*pi*k/nPhi,   k = 0,...,nPhi-1.
```

The seam is represented only by connectivity:

```text
next = (k + 1) % nPhi
```

There is therefore no second node at `phi=2*pi`.  This is important for both
surface integration and source-patch selection: a seam node can no longer be
double counted simply because it had two angular-coordinate representations.

## Triangle winding and quality

The apex-aligned basis is right handed.  Triangles are generated in the same
orientation as the parameterization `dX/dtheta x dX/dphi`, which points outward.
The rear-pole fan uses the complementary ordering required to preserve outward
winding.

`compute_triangle_metrics()` now treats malformed geometry as an error rather
than a benign zero-area cell.  It rejects:

- repeated indices within a triangle;
- out-of-range connectivity;
- non-finite vertex coordinates;
- zero or numerically unresolved edge scales;
- area below a scale-aware threshold
  `128*epsilon*L_max^2`;
- a cell normal whose dot product with the averaged analytical nodal normal is
  not positive;
- non-finite derived area/centroid/source quantities.

No normal is manufactured for a degenerate triangle and no bad cell is removed
silently.

## Mesh metadata

`ShockMesh` now records:

```cpp
std::size_t n_theta_intervals;
std::size_t n_phi;
bool closed_surface;
```

These values let source/output code audit the construction without inferring
resolution from vertex counts.  Existing 1-based `tri_i/tri_j/tri_k`
connectivity is retained for Tecplot compatibility.

## Area-weighted source-patch sampling

A uniform source intensity per unit shock area must not select mesh cells
uniformly by index because the angular grid has nonuniform physical triangle
areas.

The production API now provides:

```cpp
AreaSamplingTable Model::build_area_sampling_table(const TriMetrics&) const;
std::size_t Model::sample_triangle_by_area(
    const AreaSamplingTable&, double unit_uniform) const;
```

`build_area_sampling_table()`:

1. requires every triangle area to be finite and strictly positive;
2. accumulates total area in `long double`;
3. creates a strictly increasing normalized CDF;
4. sets the mathematical endpoint exactly to one after normalization.

`sample_triangle_by_area()` accepts a caller-provided `U[0,1)` variate and uses
`upper_bound` on that CDF.  SWCME deliberately does **not** own an RNG.  AMPS or
another caller remains responsible for seed selection, MPI/thread RNG policy,
and reproducibility while the SWCME mesh layer owns the physical area
probabilities.

## Analytical area references

`MSH03` uses independent area references.

For a sphere of radius `R`,

```text
A = 4*pi*R^2.
```

A true SSE front is part of a translated sphere.  With half width `lambda`,

```text
c = R/(1 + sin(lambda))
a = c*sin(lambda)
```

and the ray-tangent boundary corresponds to polar angle
`beta_max = pi/2 + lambda` about the translated sphere center.  Therefore

```text
A_SSE = 2*pi*a^2*(1 - cos(beta_max))
      = 2*pi*R^2*sin(lambda)^2/(1 + sin(lambda)).
```

This reference is independent of the triangle-area implementation.

## Validation tests

Five deterministic mesh tests were added.

### MSH01 — nondegenerate cells

Builds multiple Sphere and SSE meshes at minimum, typical, and fine angular
resolution and at narrow/broad SSE widths.  Every cell must have finite positive
area.  The test reports minimum, median, and maximum cell area.

### MSH02 — outward orientation

Builds rotated Sphere, anisotropic Ellipsoid, and SSE meshes.  Every geometric
triangle normal is compared with the analytical shock normal at the centroid
direction.  The dot product must be positive for every cell.

### MSH03 — area convergence

Builds `1x`, `2x`, `4x`, and `8x`-style refinement levels for Sphere and SSE.
Integrated planar-triangle area must converge monotonically to the independent
analytical area.  The final observed order must exceed 1.5 and the fine-grid
relative error must be below 0.5%.

### MSH04 — unique apex and periodic seam

Inspects topology directly for a finite SSE cap.  It verifies:

- exact vertex and triangle counts;
- no duplicate coordinate pairs;
- every edge has incidence one or two;
- exactly `nPhi` boundary edges;
- apex valence `nPhi`;
- Euler characteristic one for the disk topology;
- explicit wrapped last-to-first seam adjacency.

### MSH05 — area-weighted stochastic patch sampling

Builds the production area CDF and verifies it is strictly increasing and ends
at one.  For each of three fixed RNG seeds, 200,000 production cell selections
are aggregated into 16 bins.  The observed counts are compared with exact
physical-area probabilities using a chi-square statistic and an independent
regularized-incomplete-gamma p-value implementation.  The deterministic CI
criterion is `p > 1e-3`.

## Compatibility

The public `build_shock_mesh(S,nTheta,nPhi)` signature is unchanged, as are the
nodal fields and 1-based triangle arrays.  The number/order of generated nodes
and cells changes because duplicate topological vertices and degenerate cells no
longer exist.  Consumers that assumed a rectangular `(nTheta+1)*(nPhi+1)` node
layout must stop indexing the surface as a structured grid and instead use the
explicit triangle connectivity.

This change is intentional: a shock surface used for area integration and SEP
source sampling is an unstructured triangular manifold, not a logically
rectangular array with duplicate polar/seam points.
