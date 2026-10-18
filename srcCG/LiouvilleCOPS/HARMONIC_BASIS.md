# Localized-activity real spherical-harmonic basis, version 1

This document is the normative mathematical specification for the optional
higher-degree activity basis in LiouvilleCOPS. The implementation
identifier is:

```text
liouville_cops.real_orthonormal_cs_v1
```

Stage A2 implemented and tested this basis independently. Stage A3 now lets
the C++ source model apply a strictly versioned, optional H2O residual using
this evaluator. The default remains disabled. The Python fitter does not yet
estimate these coefficients; that integration belongs to Stage A4.

## 1. Coordinate system

For a nonzero facet-centroid vector

```text
r = (x, y, z),       R = sqrt(x^2 + y^2 + z^2),
mu = cos(theta) = z/R,
phi = atan2(y, x),   -pi <= phi <= pi,
theta = acos(mu),     0 <= theta <= pi.
```

The implementation uses `hypot(hypot(x,y),z)` and clamps `mu` to `[-1,1]`
to protect against floating-point roundoff. At either pole (`x = y = 0`), it
sets `phi = 0`; all modes with nonzero order are zero there, so the convention
does not affect the field.

The coordinates are the axes of the mesh supplied to LiouvilleCOPS. The
application interprets those axes as its configured `BodyFrame`: SPICE
spacecraft and Sun vectors are transformed into that frame before they are
combined with the mesh. The supplied manifested SHAP5 fitting base records
`BodyFrame = 67P/C-G_CK`. A NASTRAN file contains no frame metadata, so this
association is an input/provenance requirement, not a fact inferred from the
mesh file itself.

Multiplying `(x,y,z)` by any positive scalar must leave every basis value
unchanged. `(0,0,0)`, a nonfinite coordinate, or a nonfinite angle is invalid.

## 2. Associated Legendre functions and phase

`P_l^m(mu)` includes the Condon--Shortley phase `(-1)^m`. For nonnegative
`m`, both implementations use the following unnormalized recurrence:

```text
P_0^0(mu) = 1

P_m^m(mu) = -(2m-1) sqrt(max(0, 1-mu^2)) P_(m-1)^(m-1)(mu),  m >= 1

P_(m+1)^m(mu) = (2m+1) mu P_m^m(mu)

P_l^m(mu) = [(2l-1) mu P_(l-1)^m(mu)
             - (l+m-1) P_(l-2)^m(mu)] / (l-m),             l >= m+2.
```

The square-root guard handles a `mu` value rounded to a pole. This recurrence
is used only through degree 8; no claim is made here about its conditioning at
arbitrarily high degree.

## 3. Fully normalized real harmonics

Define

```text
N_lm = sqrt[(2l+1)/(4 pi) * (l-m)!/(l+m)!],    m >= 0.
```

The real functions are

```text
Y_(l,0)(theta,phi) = N_l0 P_l^0(cos(theta)),

Y_(l,+m)(theta,phi) = sqrt(2) N_lm P_l^m(cos(theta)) cos(m phi),  m > 0,

Y_(l,-m)(theta,phi) = sqrt(2) N_lm P_l^m(cos(theta)) sin(m phi),  m > 0.
```

Thus positive order is the cosine branch and negative order is the sine
branch. Because the Legendre function includes the Condon--Shortley phase,
for example:

```text
Y_(0,0) = 1/sqrt(4 pi)
Y_(1,0) at the north pole = sqrt(3/(4 pi))
Y_(1,+1) on the +x axis = -sqrt(3/(4 pi))
Y_(1,-1) on the +y axis = -sqrt(3/(4 pi)).
```

The normalization is orthonormal over solid angle:

```text
integral_0^(2pi) integral_0^pi
    Y_(l,m) Y_(l',m') sin(theta) dtheta dphi
  = delta_(l,l') delta_(m,m').
```

Normalization is evaluated from log-gamma factorials in both languages:

```text
(l-m)!/(l+m)! = exp[lgamma(l-m+1) - lgamma(l+m+1)].
```

## 4. Degree range and coefficient order

The evaluator supports degrees 0 through 8 so low-degree identities and
normalization can be tested. The optional activity enhancement uses only
degrees 5 through `L`, with `5 <= L <= 8`.

For an inclusive range `l_min..l_max`, coefficients are ordered by increasing
degree and, within a degree, increasing signed order:

```text
(l_min,-l_min), ..., (l_min,+l_min),
(l_min+1,-l_min-1), ..., (l_max,+l_max).
```

The zero-based array index is

```text
index(l,m;l_min) = l^2 + l + m - l_min^2.
```

The inclusive coefficient count is

```text
(l_max+1)^2 - l_min^2.
```

Consequently, degrees 5 only contain 11 coefficients; degrees 5--6 contain
24; degrees 5--7 contain 39; and degrees 5--8 contain 56.

## 5. Relationship to the legacy srcCG basis

The 25 legacy activity functions in `src/Source.cpp` are retained unchanged.
They are explicit, unnormalized expressions indexed 0--24 through degree 4.
The basis defined here is independently normalized and versioned. It must not
be described as the next entries of the legacy table and must not be written
into the version-1 `elements.bin` 25-function field.

Stage A3 evaluates this basis from mesh centroids in C++ when a validated
enhancement file is enabled. Stage A4 may evaluate it from the centroids
already stored in `elements.bin`. The two fields are combined only through the
separately specified enhancement model, with the legacy degree-0--4 field
remaining distinguishable. `ACTIVITY_ENHANCEMENT_FILE.md` defines that model's
file and compatibility contract.

## 6. Numerical contract

- Inputs and outputs are IEEE-754 double precision.
- Supported degree bounds satisfy `0 <= l_min <= l_max <= 8`.
- Output contains exactly the count and ordering specified above.
- Cartesian inputs must be finite and have positive radius.
- Python and C++ must agree at deterministic poles, axes, near-pole points,
  and generic points to a maximum absolute difference of `2e-14` for degrees
  0--8. This is much tighter than the downstream pressure-agreement target of
  `1e-6`.
- Numerical quadrature must reproduce unit norm and cross-mode orthogonality
  through degree 8 to the tolerance recorded in the tests.
- No transport cache format changes in Stages A2 or A3.

The authoritative implementations are `tools/activity_basis.py` and
`src/ActivityEnhancementBasis.{h,cpp}`. Tests compile a small independent C++
probe and compare every returned mode with Python.
