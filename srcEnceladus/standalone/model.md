# Enceladus Stage-0 neutral plume model

## Purpose and scope

This document describes the physics, numerical implementation, inference
method, assumptions, inputs, outputs, and execution of the standalone
Enceladus plume application in this directory. The application implements
Stage 0 of the Enceladus AMPS development plan: infer an H2O neutral-source
description from selected Cassini INMS and UVIS observations before that
source is coupled to the full AMPS and BATS-R-US numerical machinery.

The application is intentionally independent of AMPS. It is a steady,
collisionless, gravity-free analytic/semianalytic neutral model with linear
source amplitudes. Its main uses are:

- reproducing E3 and E5 sticking-corrected INMS density profiles;
- predicting H2O columns along UVIS occultation rays;
- comparing three published source geometries;
- determining which combinations of jets, tiger stripes, and empirical
  background are observable;
- separating spatial source shape from epoch-dependent total production; and
- generating a source prescription and validation diagnostics for later AMPS
  stages.

It does **not** simulate neutral-neutral collisions, gravity, plasma, pickup
ions, dust, chemistry, or time-dependent vent flow. A good fit is therefore a
conditional result under the assumptions documented below, not a validation
of the omitted physics.

## Software and data flow

```text
model JSON + stripe geometry + optional external jet catalog
                              |
                              v
                  validated physical sources
                              |
                 +------------+-------------+
                 |                          |
                 v                          v
       local-density design matrix   ray-integrated column matrix
              (INMS)                         (UVIS)
                 |                          |
                 +------------+-------------+
                              |
                 scaled nonnegative fitting
                              |
                              v
       rates, backgrounds, profiles, ranks, residuals, hashes, figures
```

The primary modules are:

| File | Responsibility |
|---|---|
| `plume_model.py` | Coordinate transformations, shifted-Maxwellian moments, flux normalization, point and line sources, local-density matrices, and ray-integrated column matrices. |
| `fitting.py` | Configuration loading, source-case materialization, weighted nonnegative least squares, optional bounded shape optimization, and identifiability diagnostics. |
| `plume_fit.py` | Strict INMS density-data interface and one-observation command-line application. |
| `joint_fitting_study.py` | E3/E5 and UVIS independent fits, instrument-combined fits, all-data fits, fixed-shape rate fits, ablations, cross-predictions, plots, and reports. |
| `compare_source_case_studies.py` | Checked synthesis of completed studies for the three source cases. |
| `config/tenishev2014_stage0.json` | Reference physics, source kinematics, source-case options, and background selection. |
| `data/published_stripes_2014.json` | Seven tiger-stripe polylines containing 102 IAU_ENCELADUS surface nodes. |
| `tools/spice_ground_tracks.cpp` | CSPICE calculation of flyby states and ellipsoidal ground tracks. |
| `tools/spice_observer_ranges.cpp` | CSPICE reconstruction of finite Cassini-to-Sun UVIS ray geometry. |

Observation data and derived catalogs remain outside the code tree under
`/data/vtenishe/cassini_validation_data`. No observation directory is inferred
by the code; paths are supplied explicitly by CLI options or the documented
environment variables.

## Units, frame, and geometry conventions

All physical calculations use SI units.

- Positions and distances: metres.
- Local number density: m^-3.
- Column density: m^-2.
- Temperature: kelvin.
- Speed: m/s.
- Production rate `Q`: molecules/s.
- Model frame: `IAU_ENCELADUS`.
- Reference spherical radius: 252,100 m in the supplied configuration.
- H2O molecular mass: 2.991507624607566e-26 kg.

Latitude is signed planetocentric latitude. The published tables often give
positive magnitudes for south latitude; those values are stored as negative
latitudes in the configuration. Longitude is positive westward. Internally,
the spherical transformation uses east longitude equal to the negative of the
configured west longitude.

For latitude `lambda` and east longitude `phi`, the outward spherical normal
is

```text
n = (cos(lambda) cos(phi), cos(lambda) sin(phi), sin(lambda)).
```

The local north and east unit vectors form the tangent basis. Jet tilt is
measured from the outward radial normal. Azimuth is clockwise from local north.
The jet axis is

```text
a = cos(tilt) n
    + sin(tilt) [cos(azimuth) north + sin(azimuth) east].
```

Configuration validation rejects inward-pointing axes, nonunit directions,
invalid temperatures or speeds, and incompatible frames or units.

## Neutral source physics

This section gives the mathematical model actually evaluated by
`plume_model.py`. Bold symbols below are vectors. Each fitted component is
normalized to a physical H2O production rate; no arbitrary kernel amplitude
is left in the forward model.

### Point-source geometry and shifted Maxwellian

Consider source (i) at position \(\mathbf{x}_i\), with outward surface
normal \(\mathbf{n}_i\), unit jet axis \(\mathbf{a}_i\), temperature \(T_i\),
bulk speed \(u_i\), and molecular mass \(m\). At a requested point
\(\mathbf{x}\), define

\[
\begin{aligned}
\mathbf{d}_i &= \mathbf{x}-\mathbf{x}_i, &
R_i &= |\mathbf{d}_i|,\\
\widehat{\mathbf{r}}_i &= \mathbf{d}_i/R_i, &
c_i &= \widehat{\mathbf{r}}_i\mathbin{\cdot}\mathbf{a}_i,\\
\mu_i &= \widehat{\mathbf{r}}_i\mathbin{\cdot}\mathbf{n}_i, &
\gamma_i &= \mathbf{a}_i\mathbin{\cdot}\mathbf{n}_i .
\end{aligned} \tag{1}
\]

Here \(c_i\) is the cosine of the angle between the source-to-observer ray and
the bulk-flow axis, while \(\mu_i\) determines whether that ray leaves the
surface. Only the outward hemisphere is populated: the kernel is set to zero
for \(\mu_i<0\). The analytic point-source solution is singular at \(R_i=0\),
so evaluation exactly on a source is prohibited.

The most-probable thermal speed and dimensionless drift speed are

\[
v_{\mathrm{th},i}=\sqrt{\frac{2k_{\mathrm B}T_i}{m}},
\qquad s_i=\frac{u_i}{v_{\mathrm{th},i}} . \tag{2}
\]

Along a particular ray, molecular velocity is \(\mathbf{v}=v
\widehat{\mathbf r}_i\), with \(v\ge 0\). The implemented velocity kernel is a
shifted Maxwellian,

\[
\exp\!\left[-\frac{|v\widehat{\mathbf r}_i-u_i\mathbf a_i|^2}
{v_{\mathrm{th},i}^2}\right]
=\exp[-t^2+2s_i c_i t-s_i^2],
\qquad t=\frac{v}{v_{\mathrm{th},i}} . \tag{3}
\]

This is a collisionless free-streaming kernel, not a locally re-thermalized
Maxwellian at \(\mathbf{x}\). The velocity direction at a field point is fixed
by the ray from the source; only the speed is integrated. Equivalently, one
may write the ray distribution as

\[
f_i(v,\widehat{\mathbf r}_i;R_i)=
\frac{C_i}{R_i^2}L(\mu_i)
\exp[-t^2+2s_i c_i t-s_i^2]H(\mu_i), \tag{4}
\]

where \(H\) is the Heaviside outward-hemisphere restriction and \(L\) is the
chosen angular source law. With the angular delta function implied by radial
free streaming, number density and source-centered radial number flux are

\[
n_i=\int_0^\infty f_i v^2\,dv,
\qquad
F_{R,i}=\int_0^\infty f_i v^3\,dv. \tag{5}
\]

Equation (5) explains why the second velocity moment produces density and the
third moment produces flux.

### Radial velocity moments and their numerical evaluation

Define the dimensionless moments

\[
K_p(c,s)=\int_0^\infty t^p
\exp[-t^2+2sct-s^2]\,dt, \qquad p\in\{2,3\}. \tag{6}
\]

Putting \(z=sc\) and completing the square gives

\[
K_p(c,s)=e^{-s^2(1-c^2)}J_p(z),\qquad
J_p(z)=\int_0^\infty t^p e^{-(t-z)^2}\,dt. \tag{7}
\]

For \(z\ge-2\), the code evaluates the following recurrence rather than
performing a new numerical integral for every point:

\[
\begin{aligned}
J_0(z)&=\frac{\sqrt{\pi}}{2}\operatorname{erfc}(-z),\\
J_1(z)&=\frac{1}{2}e^{-z^2}+zJ_0(z),\\
J_2(z)&=\frac{1}{2}J_0(z)+zJ_1(z),\\
J_3(z)&=J_1(z)+zJ_2(z).
\end{aligned} \tag{8}
\]

For \(z<-2\), direct use of (8) can lose precision through cancellation. The
equivalent confluent-hypergeometric form is used instead,

\[
K_p(c,s)=e^{-s^2}\frac{\Gamma(p+1)}{2^{p+1}}
U\!\left(\frac{p+1}{2},\frac12,z^2\right), \tag{9}
\]

where \(U\) is the Tricomi confluent hypergeometric function. Tests compare
both branches with independent adaptive quadrature over representative
\((c,s)\) values.

### Selectable source laws

Two laws are implemented:

1. `tenishev2010_point` (default): the outward angular factor is
   `L(mu) = 1`.
2. `surface_flux_maxwellian`: the projected-area factor is
   `L(mu) = max(mu,0)`.

The first is the 2010 point-source emission law used by the reference model.
The second is selectable for controlled sensitivity studies. The difference
is physical: a tangent ray can have nonzero density for the point-source law,
whereas the surface-flux law vanishes at `mu = 0`.

In the Tenishev et al. (2010) point-source derivation, the injection
distribution contains a directional cosine and the conversion from particles
crossing a surface to volume density contains the reciprocal directional
cosine. Those factors cancel in the far-field density, leaving no additional
\(\mu\) multiplier. That cancellation is represented here by \(L=1\).
The surface-flux option deliberately retains a projected-area factor as a
different boundary-law sensitivity test; it should not be described as the
same 2010 emission law.

### Exact production-rate normalization

Let the polar axis of the angular integration be the surface normal. For
\(\mu=\cos\theta\in[0,1]\), azimuth \(\phi\in[0,2\pi)\), and
\(\gamma=\mathbf a\cdot\mathbf n\), the ray-axis cosine is

\[
c(\mu,\phi)=\gamma\mu+
\sqrt{1-\gamma^2}\sqrt{1-\mu^2}\cos\phi . \tag{10}
\]

The dimensionless angular normalization is therefore

\[
\mathcal N(s,\gamma)=
\int_0^1\int_0^{2\pi}
L(\mu)K_3[c(\mu,\phi),s]\,d\phi\,d\mu . \tag{11}
\]

Gauss-Legendre quadrature is evaluated at increasing orders 24, 48, 96, 192,
and 384 until the relative change is at most 1e-9. The result is cached by
speed ratio, axis-normal cosine, and source law.

Substitution of (6) into (5) gives a factor \(v_{\mathrm{th}}^3 K_2\)
for density and \(v_{\mathrm{th}}^4K_3\) for flux. Requiring the integral of
radial flux through any source-centered hemisphere to equal the physical
production rate \(Q_i\),

\[
Q_i=R_i^2\int_{\mu\ge0}F_{R,i}\,d\Omega, \tag{12}
\]

sets \(C_i=Q_i/[v_{\mathrm{th},i}^4\mathcal N_i]\). The final implemented
point-source expressions are

\[
\boxed{
n_i(\mathbf{x})=
\frac{Q_iL(\mu_i)K_2(c_i,s_i)}
{\mathcal N_i v_{\mathrm{th},i}R_i^2}H(\mu_i)
} \tag{13}
\]

and

\[
\boxed{
F_{R,i}(\mathbf{x})=
\frac{Q_iL(\mu_i)K_3(c_i,s_i)}
{\mathcal N_iR_i^2}H(\mu_i)
}. \tag{14}
\]

Equation (12) then follows identically from (11) and (14). The units are
\([n_i]=\mathrm{m}^{-3}\) and \([F_{R,i}]=\mathrm{m}^{-2}\mathrm{s}^{-1}\).
The \(R_i^{-2}\) factor is geometric dilution. Angular anisotropy arises from
the shifted Maxwellian through \(K_p(c_i,s_i)\), from the outward surface
horizon, and, for the surface-flux option, from \(L(\mu_i)\).

Thus every fitted source coefficient is a physical total emission rate rather
than an arbitrary angular amplitude. Tests verify flux conservation for both
source laws and exact inverse-square density scaling along fixed rays.

### Density from all sources at an arbitrary point

For \(N_s\) point or composite stripe sources, the localized-plume density is
linear in their rates:

\[
n_{\mathrm{plume}}(\mathbf{x})=
\sum_{i=1}^{N_s}Q_i\,G_i(\mathbf{x}),\qquad
G_i(\mathbf{x})\equiv\frac{n_i(\mathbf{x})}{Q_i}, \tag{15}
\]

where \(G_i\) is equation (13) evaluated at unit production. Including the
empirical background gives the value compared with an INMS density datum,

\[
n_{\mathrm{model}}(\mathbf{x})=
\sum_iQ_iG_i(\mathbf{x})+c_0+b
\left(\frac{|\mathbf{x}|}{R_E}\right)^\alpha . \tag{16}
\]

Terms not selected in the configuration are omitted. In the reference model
\(R_E=252.1\) km and \(\alpha=-2\). Equation (16), evaluated independently at
every spacecraft position, is the complete local-density forward operator;
there is no streamline tracing, spatial mesh, interpolation, or time stepping.

### Tiger-stripe line sources

Each stripe component is a polyline of surface nodes. A radially directed
point source is created at every node with the stripe's common temperature,
speed, and source law. Trapezoidal arc-length weights distribute a single
component rate uniformly along physical line length. If node positions are
\(\mathbf{x}_0,\ldots,\mathbf{x}_{N-1}\) and segment lengths are
\(\Delta\ell_j=|\mathbf{x}_{j+1}-\mathbf{x}_j|\), the unnormalized nodal
trapezoid weights are

\[
\widetilde w_0=\frac{\Delta\ell_0}{2},\qquad
\widetilde w_j=\frac{\Delta\ell_{j-1}+\Delta\ell_j}{2},\qquad
\widetilde w_{N-1}=\frac{\Delta\ell_{N-2}}{2}, \tag{17}
\]

and \(w_j=\widetilde w_j/\sum_k\widetilde w_k\). For one fitted stripe rate
\(Q_\ell\), the stripe density is

\[
n_\ell(\mathbf{x})=Q_\ell\sum_{j=0}^{N-1}w_jG_j(\mathbf{x}). \tag{18}
\]

Because \(\sum_jw_j=1\) and every nodal kernel is separately flux-normalized,
the surface-integrated production remains exactly \(Q_\ell\). This is a
piecewise-linear quadrature of a uniform source per unit mapped arc length;
it is not a claim that the physical vent output is uniform.

The density of a stripe component is the weighted sum of its nodes. A stripe
with many nodes still has only one fitted total rate. The reference hybrid
case contains seven stripe components and 102 nodes; it does not fit 102
independent vent rates.

### Empirical background

The local-density operator supports `none`, `constant`, `power_law`, or `both`.
The reference configuration uses

```text
n_background(r) = c0 + b (r/R_E)^alpha,
alpha = -2,
R_E = 252.1 km.
```

`b` is the fitted background density at the reference radius. The equivalent
coefficient for `c0 + c1/r^2` is `c1 = b R_E^2` when `alpha = -2`.

Occultation integrations support only the power-law background. A spatially
constant density has a divergent infinite-chord column and is therefore
deliberately excluded from the UVIS operator. Background terms are nuisance
parameters and are not included in reported localized plume production.

## Source geometry cases

The source case selects a spatial basis. It does not import fixed production
rates: all selected physical components receive newly fitted nonnegative
amplitudes.

| Source case | Physical components | Direction and kinematics | Columns with reference background |
|---|---:|---|---:|
| `tenishev2014_hybrid_8` | 8 jets + 7 stripe components | Published per-source 2014 values; stripe families at their configured speeds | 17 |
| `porco2014_98` | 98 discrete jets | Published tilts/azimuths; controlled `T=180 K`, `u=1500 m/s` | 100 |
| `portyankina2022_solar_41` | 41 selected jets | Local radial axes; `T=180 K`, `u=1500 m/s` | 43 |

### Hybrid eight-jet plus stripe case

This is the default and is self-contained in the code tree. Its eight point
sources use the locations, directions, temperatures, and speeds from Tenishev
et al. (2014). It includes four normal line components on Alexandria,
Baghdad, Cairo, and Damascus at 180 K and 350 m/s, plus three high-speed line
components on Baghdad, Cairo, and Damascus at 180 K and 950 m/s.

### Porco 98-jet case

This case uses all jets in Porco et al. (2014), Table 2, through the external
catalog. There are no additional continuous stripe components. Porco et al.
constrained particulate-jet geometry, not unique H2O gas kinematics for every
jet; 180 K and 1.5 km/s are explicit controlled-comparison assumptions.

With two 79-sample INMS profiles and 100 fitted columns, this case is
underdetermined even before geometrically similar jets are considered. Its
individual amplitudes, and potentially its total rate, can be extremely
sensitive to geometry and gas kinematics.

### Portyankina 41-jet case

This case uses the 41 Porco sources having strictly positive coefficients in
the Portyankina et al. (2022) 2010 solar-occultation solution. The catalog's
published coefficients are retained only as provenance; this application
refits all 41 rates. Directions are local outward normals following that
calculation. The other 57 Porco jets are not proven inactive; the 41-source
list is one selected nonnegative basis.

The 98- and 41-jet cases require:

```text
/data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json
```

The path can be supplied with `--jet-catalog`, by `jet_catalog_file` in a
copied configuration, or by `ENCELADUS_JET_CATALOG`.

## Observation operators

### INMS local density

At spacecraft sample \(j\), equation (16) is evaluated at the SPICE-derived
position \(\mathbf{x}_j\). For each localized source, the stored matrix column
is its density response to a reference production
\(Q_\star=10^{28}\ \mathrm{s}^{-1}\):

\[
A_{ji}=Q_\star G_i(\mathbf{x}_j),
\qquad i=1,\ldots,N_s . \tag{19}
\]

If enabled, the background columns are

\[
A_{j,c_0}=1,
\qquad
A_{j,b}=\left(\frac{|\mathbf{x}_j|}{R_E}\right)^\alpha . \tag{20}
\]

The coefficient vector is mixed-unit by construction:

\[
q_i=Q_i/Q_\star\quad\hbox{for localized sources},\qquad
q_{c_0}=c_0,\qquad q_b=b. \tag{21}
\]

Thus the predicted density vector is simply

\[
\widehat{\mathbf y}=\mathbf A\mathbf q, \tag{22}
\]

and the reported total localized production is
\(Q_{\mathrm{tot}}=Q_\star\sum_{i=1}^{N_s}q_i\). Empirical-background
coefficients are excluded from this sum.

`plume_fit.py` accepts only qualified number-density CSV files accompanied by
metadata declaring:

- observable `number_density`;
- frame `IAU_ENCELADUS`;
- position units `m`;
- density and uncertainty units `m^-3`;
- uncertainty model `independent_gaussian`;
- instantaneous sampling; and
- either `synthetic` or `calibrated_observation` data kind.

Raw detector counts are rejected. A mass-44 count cannot be relabeled as H2O
or CO2 density without a separate inlet/surface response and calibration
operator.

The supplied E3/E5 observations are reconstructions of the published
sticking-corrected profiles, not new inversions of L1A counts. Each has 79
outbound samples from +10 to +400 s at 5 s spacing. Their exact construction,
absolute normalization, digitization uncertainty, common 23% scale systematic,
and trajectory provenance are recorded in the adjacent metadata JSON files.

### UVIS line-of-sight columns

For UVIS occultation sample \(j\), the fundamental model observable is the
H2O line-of-sight column,

\[
N_j=\int_{s_{\min,j}}^{s_{\max,j}}
n_{\mathrm{model}}(\mathbf{x}_j(s))\,ds,
\qquad
\mathbf{x}_j(s)=\mathbf b_j+s\widehat{\mathbf d}_j, \tag{23}
\]

where \(\mathbf b_j\cdot\widehat{\mathbf d}_j=0\), \(\mathbf b_j\) is the
closest-approach vector, and \(\widehat{\mathbf d}_j\) is a unit line-of-sight
direction. Because density is linear in the source rates, the column is also
linear:

\[
N_j=\sum_iQ_i\mathcal G_{ji}+N_{\mathrm{bg},j},
\qquad
\mathcal G_{ji}=\int G_i(\mathbf b_j+s\widehat{\mathbf d}_j)\,ds. \tag{24}
\]

The source columns stored in the UVIS design matrix are
\(B_{ji}=Q_\star\mathcal G_{ji}\), with units m\(^{-2}\) per fitted unit of
\(10^{28}\ \mathrm{s}^{-1}\).

For a distant star, the ray is an infinite chord
\((s_{\min},s_{\max})=(-\infty,+\infty)\). Defining
\(b_j=|\mathbf b_j|\), the transformation

\[
s=b_j\tan\theta,\qquad
ds=b_j\sec^2\theta\,d\theta,
\qquad -\frac{\pi}{2}<\theta<\frac{\pi}{2} \tag{25}
\]

maps the unbounded integral to a finite Gauss-Legendre interval. If the
observer occupies signed coordinate \(s_{\mathrm{obs},j}\) on that same line,
the lower angular limit becomes
\(\theta_{\min}=\tan^{-1}(s_{\mathrm{obs},j}/b_j)\). The upper limit remains
\(\pi/2\).

If the tangent of the infinite line lies behind the finite observer, the
physical ray is integrated directly from observer position \(\mathbf x_{0,j}\)
to infinity:

\[
\mathbf{x}_j(s)=\mathbf{x}_{0,j}+s\widehat{\mathbf d}_j,
\qquad 0\le s<\infty . \tag{26}
\]

Using a Gauss-Legendre coordinate \(\xi\in(-1,1)\), the mapping is

\[
s=|\mathbf{x}_{0,j}|\frac{1+\xi}{1-\xi},
\qquad
\frac{ds}{d\xi}=\frac{2|\mathbf{x}_{0,j}|}{(1-\xi)^2}. \tag{27}
\]

Equations (25) and (27) retain the asymptotic \(1/r^2\) tails without imposing
an arbitrary outer boundary. The code normalizes every ray direction and
rejects nonorthogonal tangent geometry, a tangent at or below the reference
sphere, nonfinite geometry, and inward rays passed to the outward-ray
operator.

For a power-law background \(n_{\mathrm{bg}}=b(r/R_E)^\alpha\), the numerical
integral exists only for \(\alpha<-1\). In the reference \(\alpha=-2\) case,
the infinite-chord result provides a useful analytic verification:

\[
N_{\mathrm{bg}}(b_j)=
bR_E^2\int_{-\infty}^{\infty}\frac{ds}{b_j^2+s^2}
=\frac{\pi bR_E^2}{b_j}. \tag{28}
\]

A constant background is not admitted in the UVIS operator because its
infinite column diverges.

The production studies use order 192 and compare a deterministic 12-ray subset
with order 224. The matrix columns have units m^-2 per `1e28 molecules/s`.

The distant 2007 zeta Orionis event is treated as an infinite chord. Its line
direction is constructed as the normalized cross product of tangent radius and
tangent-track velocity, with signs made continuous. The 2010 solar event uses
SPICE-derived finite Cassini-to-Sun directions and observer coordinates. Rays
for which the infinite-line tangent lies behind the observer are integrated
outward from the finite observer instead of through a fictitious chord.

**Current 2010 geometry limitation.** The generic outward-ray operator,
equations (26)--(27), is correct for an actual finite origin. However, the
present study driver passes the tabulated tangent vector as that origin for
the 308 of 1,133 E10 rows whose observer coordinate is nonnegative. The
generated geometry table stores spacecraft radius and its scalar coordinate
along the Sun ray, but not the complete three-dimensional Cassini position
needed to reconstruct the true origin. Thus the current outward-branch
columns begin at the tangent point, not at Cassini, and are an approximation.
Before publication-grade use, the SPICE product must retain the full
IAU_ENCELADUS spacecraft vector and the caller must pass that vector to the
outward-ray operator. The 825 rows with negative observer coordinate use the
finite lower limit on the tangent-line parameter as documented above. This
limitation must remain disclosed in results produced by the current code.

### UVIS calibration boundary

UVIS measures attenuation of a stellar or solar source, not local plume
emission. In an ideal spectrally resolved absorption model, Beer--Lambert
attenuation would give

\[
I_j(\lambda)=I_{0,j}(\lambda)e^{-\tau_j(\lambda)},
\qquad
\tau_j(\lambda)=\sum_k\sigma_k(\lambda)N_{k,j}, \tag{29}
\]

where \(\sigma_k\) is the absorption cross section of species \(k\). For H2O
alone, a monochromatic calibrated measurement would satisfy
\(N_{\mathrm{H_2O}}=\tau/\sigma_{\mathrm{H_2O}}\). A broadband detector instead
measures a response-weighted integral of \(I_0(\lambda)e^{-\tau(\lambda)}\), so
one scalar cross section cannot in general convert every broadband count into
an absolute column.

The prepared data use a fixed off-plume reference count \(C_{0}\), observed
count \(C_j\), transmission \(T_j=C_j/C_0\), and raw optical-depth proxy

\[
\tau_j^{\mathrm{raw}}=-\ln T_j=-\ln(C_j/C_0),
\qquad
\sigma_{\tau,j}\simeq C_j^{-1/2}. \tag{30}
\]

The uncertainty in (30) is the observed-count Poisson term only; uncertainty
in \(C_0\), spectral calibration, cross sections, and the published absolute
column are separate systematics. Samples with \(C_j>C_0\) legitimately have
negative raw \(\tau\) and are retained in the shape fit rather than clipped.

Because the wavelength-dependent calibration operator is unavailable, the
current model fits only resolved optical-depth *shape*. For one event, let
\(\mathbf B\) be the column-density design matrix from equation (24),
\(\mathbf p\ge0\) a set of nonnegative shape coefficients that absorb the
unknown common optical-depth/column conversion, and \(\beta\) a signed
reference-level nuisance offset. The fitted proxy is

\[
\widehat{\boldsymbol\tau}=\mathbf B\mathbf p+\beta\mathbf 1. \tag{31}
\]

The code minimizes

\[
\chi^2_{\tau}=\left\|\mathbf W_\tau
(\mathbf B\mathbf p+\beta\mathbf1-
\boldsymbol\tau^{\mathrm{raw}})\right\|_2^2,
\qquad p_i\ge0, \tag{32}
\]

with \(\mathbf W_\tau=\operatorname{diag}(1/\sigma_{\tau,j})\). To allow a
signed \(\beta\) while retaining a nonnegative active-set solve for
\(\mathbf p\), define \(\mathbf b=\mathbf W_\tau\mathbf1\),
\(\mathbf D=\mathbf W_\tau\mathbf B\), \(\mathbf z=\mathbf W_\tau
\boldsymbol\tau^{\mathrm{raw}}\), and the weighted baseline-removal projector

\[
\mathbf P=\mathbf I-\frac{\mathbf b\mathbf b^{\mathsf T}}
{\mathbf b^{\mathsf T}\mathbf b}. \tag{33}
\]

The Lawson--Hanson NNLS problem is

\[
\min_{\mathbf p\ge0}\|\mathbf P(\mathbf D\mathbf p-\mathbf z)\|_2^2, \tag{34}
\]

after which the analytically optimal baseline is

\[
\beta=\frac{\mathbf b^{\mathsf T}(\mathbf z-\mathbf D\mathbf p)}
{\mathbf b^{\mathsf T}\mathbf b}. \tag{35}
\]

This is the exact operation performed by `bounded_weighted_fit(...,
signed_last=True)`.

For each event the procedure then:

1. augments the nonnegative column matrix with a signed constant optical-depth
   baseline;
2. eliminates the signed nuisance baseline analytically;
3. solves the projected source problem with nonnegative least squares;
4. reconstructs the best optical-depth shape; and
5. applies one published peak or representative H2O column to set the common
   physical amplitude.

The supplied constraints are:

- 2007 zeta Orionis: `(1.40 +/- 0.14)e20 m^-2`;
- 2010 E10 solar: `(0.90 +/- 0.23)e20 m^-2`.

Let \(S_j=(\mathbf B\mathbf p)_j\) and let \(N_{\mathrm{pub}}\) be the one
published representative/peak column for that event. The physical coefficient
vector is defined by one common scale factor,

\[
g=\frac{N_{\mathrm{pub}}}{\max_j S_j},\qquad
\mathbf q_{\mathrm{UVIS}}=g\mathbf p,
\qquad
\widehat N_j=(\mathbf B\mathbf q_{\mathrm{UVIS}})_j. \tag{36}
\]

Equation (36) enforces the published amplitude once; it does not turn the
individual broadband samples into independent absolute-column retrievals.

For UVIS-combined and all-data fits, the independently fitted baseline is
subtracted and negative values are clipped to zero. The remaining optical-depth
profile is converted to a physical-column proxy by

\[
N_j^{\mathrm{proxy}}=
\max(\tau_j^{\mathrm{raw}}-\beta,0)
\frac{N_{\mathrm{pub}}}
{\max_k\max(\tau_k^{\mathrm{raw}}-\beta,0)}, \tag{37}
\]

and \(\sigma_{\tau,j}\) is multiplied by the same scale. Equation (37) is
explicitly a proxy, not a new wavelength-resolved UVIS calibration. In a
future calibrated implementation, equations (31)--(37) should be replaced by
forward convolution of equation (29) with the UVIS spectral response.

## Inference and diagnostics

### Weighted nonnegative amplitude fit

For fixed source geometry and kinematics, amplitudes are linear and
nonnegative. The minimized objective is

```text
chi2 = sum_j [(A q - y)_j / sigma_j]^2,
q_i >= 0.
```

Before solving, every weighted column `A[:,i]/sigma` is divided by its Euclidean
norm. This reduces numerical scaling differences between source and background
columns. SciPy `lsq_linear` then solves the scaled bounded problem with a
relative tolerance of 1e-12 and at most 5,000 iterations. Unobservable zero
columns are assigned zero amplitude.

In matrix notation, for measured vector \(\mathbf y\), standard deviations
\(\boldsymbol\sigma\), and design matrix \(\mathbf A\), define

\[
\mathbf W=\operatorname{diag}(\sigma_1^{-1},\ldots,\sigma_N^{-1}),
\qquad
\widetilde{\mathbf A}=\mathbf W\mathbf A,
\qquad
\widetilde{\mathbf y}=\mathbf W\mathbf y. \tag{38}
\]

The fitted vector is therefore

\[
\widehat{\mathbf q}=
\arg\min_{\mathbf q\ge0}
\|\widetilde{\mathbf A}\mathbf q-\widetilde{\mathbf y}\|_2^2,
\qquad
\chi^2=\sum_{j=1}^{N}
\left(\frac{(\mathbf A\widehat{\mathbf q})_j-y_j}{\sigma_j}\right)^2 .
\tag{39}
\]

Nonnegativity is imposed on plume production and on empirical density
backgrounds. The only signed fitted quantity in the complete workflow is the
UVIS optical-depth reference offset \(\beta\), handled separately by equations
(33)--(35).

The columns of \(\widetilde{\mathbf A}\) can differ greatly in magnitude
because the vector mixes production-rate and density-background parameters.
For column \(k\), the solver calculates

\[
d_k=\|\widetilde{\mathbf A}_{:k}\|_2,
\qquad
\mathbf C_{:k}=\widetilde{\mathbf A}_{:k}/d_k, \tag{40}
\]

solves \(\min_{\mathbf z\ge0}\|\mathbf C\mathbf z-\widetilde{\mathbf y}\|^2\),
and transforms back with \(q_k=z_k/d_k\). A column with \(d_k=0\) is
unobservable and assigned zero.

Goodness of fit and matrix identifiability are kept conceptually separate. If
\(s_1\ge\cdots\ge s_r\) are singular values of \(\mathbf C\), the reported
condition number and threshold-dependent effective rank are

\[
\kappa=\frac{s_1}{s_r},
\qquad
r_{\mathrm{eff}}(\epsilon)=
\#\{k:s_k/s_1\ge\epsilon\},
\quad \epsilon\in\{10^{-2},10^{-3},10^{-4}\}. \tag{41}
\]

The result includes:

- physical coefficients and predictions;
- normalized residuals and chi-square;
- formal matrix rank;
- the full singular spectrum of the scaled matrix;
- effective ranks at relative singular-value thresholds 1e-2, 1e-3, and 1e-4;
- scaled condition number;
- unobserved column indices; and
- explicit rank-deficiency and ill-conditioning warnings.

Formal floating-point rank is not sufficient for physical interpretation.
Nearby jets can generate nearly parallel signatures even when the matrix is
technically full rank. Individual source rates must not be reported as unique
measurements when the condition number and effective rank show otherwise.
The nominal reduced statistic is
\(\chi^2_\nu=\chi^2/[N-r-p_{\mathrm{nl}}]\), where \(r\) is formal matrix rank
and \(p_{\mathrm{nl}}\) is the number of optimized nonlinear shape parameters.
It is labeled *nominal* because active inequality constraints, digitization
errors, and correlated calibration systematics violate the assumptions of a
textbook unconstrained chi-square distribution. A covariance matrix is
deliberately not reported.

### Optional variable projection for source shape parameters

`shape_parameters` can list bounded per-jet changes to:

- `temperature_K`;
- `bulk_speed_m_s`;
- `tilt_deg`; or
- `azimuth_deg`.

For every trial nonlinear parameter vector, the code rebuilds the physical
matrix and solves all linear nonnegative amplitudes. SciPy `least_squares`
optimizes the outer residual with explicit bounds and tolerances of 1e-10.
This is variable projection: linear amplitudes are never promoted to outer
nonlinear parameters.

If \(\boldsymbol\theta\) contains the selected nonlinear kinematic parameters
and \(\mathbf A(\boldsymbol\theta)\) is the regenerated forward matrix, the
outer problem is

\[
\widehat{\boldsymbol\theta}=
\arg\min_{\boldsymbol\theta_{\min}\le\boldsymbol\theta\le
\boldsymbol\theta_{\max}}
\left\|\mathbf W\left[\mathbf A(\boldsymbol\theta)
\widehat{\mathbf q}(\boldsymbol\theta)-\mathbf y\right]\right\|_2^2, \tag{42}
\]

where \(\widehat{\mathbf q}(\boldsymbol\theta)\) is the inner NNLS solution of
equation (39). The outer solver uses Jacobian scaling, tolerances
\(f_{\mathrm{tol}}=x_{\mathrm{tol}}=g_{\mathrm{tol}}=10^{-10}\), and at most
1,000 forward evaluations. The final output records the outer Jacobian rank
and parameters that terminate at a bound.

The reference configuration has an empty `shape_parameters` list. E3/E5 alone
cannot support free kinematics for many sources without strong external priors.

### Combined fits and dataset balancing

`solve_combined_blocks()` fits one nonnegative vector to multiple datasets.
When balancing is enabled, every point uncertainty in dataset `d` is multiplied
by `sqrt(N_d)`. The minimized objective becomes

```text
chi2_balanced = sum_d (1/N_d)
                sum_j-in-d [(A_d q - y_d)_j / sigma_dj]^2.
```

Equivalently,

\[
\chi^2_{\mathrm{bal}}=
\sum_{d=1}^{D}\frac{1}{N_d}
\sum_{j=1}^{N_d}
\left(\frac{(\mathbf A_d\mathbf q)_j-y_{dj}}
{\sigma_{dj}}\right)^2 . \tag{43}
\]

This prevents the 1,133-sample solar occultation from dominating a 97-sample
stellar occultation solely because it has more rows. The result records the
raw chi-square, chi-square per sample, and balanced contribution for each
dataset. This balanced objective is a diagnostic, not a formal joint
likelihood across heterogeneous instruments.

The study calculates:

- independent E3, E5, 2007 UVIS, and 2010 UVIS fits;
- one E3+E5 INMS vector;
- one 2007+2010 UVIS vector; and
- one balanced vector using all four observations.

### Fixed-shape, independently refitted rates

Orbital variability should not be estimated by comparing a single combined
rate across epochs. The code therefore separates spatial shape from total
normalization. Given a fitted vector, only jet and stripe coefficients are
retained and normalized to unit sum:

```text
f_i = q_i / sum_source q_i,
sum_source f_i = 1.
```

In symbols, with \(\mathcal S\) denoting the jet and stripe components,

\[
f_i=\frac{q_i}{\sum_{k\in\mathcal S}q_k},
\qquad i\in\mathcal S,\qquad
\sum_{i\in\mathcal S}f_i=1. \tag{44}
\]

For a target observation, the physical source columns are collapsed to
`A_shape = sum_i f_i A_i`. One nonnegative scalar multiplying this column is
then fitted together with target-specific background nuisance amplitudes. The
scalar is the target's total localized source rate in units of 1e28/s.

The collapsed response and target-specific inverse problem are

\[
\mathbf a_{\mathrm{shape}}=
\sum_{i\in\mathcal S}f_i\mathbf A_{:i}, \tag{45}
\]

\[
(\widehat q_{\mathrm{tot}},\widehat{\mathbf b})=
\arg\min_{q_{\mathrm{tot}}\ge0,\,\mathbf b\ge0}
\left\|\mathbf W\left[
q_{\mathrm{tot}}\mathbf a_{\mathrm{shape}}+
\mathbf A_{\mathrm{bg}}\mathbf b-\mathbf y\right]\right\|_2^2 . \tag{46}
\]

The physical rate is
\(Q_{\mathrm{tot}}=10^{28}\widehat q_{\mathrm{tot}}\ \mathrm{s}^{-1}\).
Background coefficients \(\mathbf b\) are independently refitted at every
epoch and never enter \(\mathbf f\) or the reported plume rate.

Four strategies are produced:

1. `independent_each_observation`: every observation independently determines
   its own shape and rate.
2. `shared_shape_all_observations`: all four datasets train one shape, and
   every observation refits its total rate.
3. `shared_shape_inms_only`: E3 and E5 train the shape, which is transferred to
   all INMS and UVIS observations.
4. `shared_shape_uvis_only`: the two UVIS events train the shape, which is
   transferred to all observations.

Backgrounds never enter the normalized plume shape or reported source rate.

### Component ablation and cross-prediction

Every study also fits:

- `full`;
- `no_stripes`;
- `no_jets`; and
- `background_only`.

For the 98- and 41-jet cases, `no_stripes` equals `full` because those cases
contain no separate line components. Cross-prediction evaluates every fitted
vector against every observation. A poor transfer can indicate temporal
variability, missing transport physics, source geometry error, or instrument
response error; it does not identify which explanation is correct.

## Model assumptions and validity limits

### Explicit physical assumptions

- The plume is steady during each observation.
- Molecules travel ballistically along straight lines after emission.
- Neutral-neutral collisions are negligible.
- Enceladus gravity is neglected.
- Sources lie on a sphere of radius 252.1 km.
- Every stripe node uses the spherical radial normal rather than local terrain
  topography.
- Source temperature, bulk speed, direction, and law are fixed unless an
  explicit bounded `shape_parameters` entry is supplied.
- The solution is H2O only.
- Empirical backgrounds are nonnegative and time-independent within one fit.
- INMS processed density samples are treated as instantaneous.
- UVIS extinction shape is proportional to the modeled H2O column after a
  signed reference baseline; one published column provides absolute scale.

### Omitted physics

The model omits:

- neutral-neutral collisions and collisional acceleration;
- Enceladus gravity, surface interception after launch, and ballistic return;
- surface topography and terrain-normal source directions;
- condensation, sublimation, and phase change;
- photodissociation, photoionization, electron-impact chemistry, and charge
  exchange;
- pickup ions, plasma feedback, fields, and magnetospheric flow;
- dust and gas-dust momentum or thermal coupling;
- source intermittency, tidal stress response, and propagation delay;
- spacecraft exposure-time convolution;
- the full INMS 1,293-element adsorption/re-emission response; and
- a wavelength-resolved UVIS spectral response and retrieval likelihood.

### Interpretation rules

- E3 and E5 are development/fitting observations, not independent validation
  of a source tuned to them.
- A smaller residual from a larger jet catalog does not prove that its
  individual amplitudes are physical.
- Source totals are conditional on source geometry, kinematics, emission law,
  and background.
- Pointwise fit sigma and common absolute calibration systematics have
  different meanings and must not be combined as independent errors.
- The UVIS proxy is suitable for controlled shape and published-scale
  comparisons, not for claiming a new absolute column retrieval.
- Held-out observations lose their validation role if they are used for
  parameter tuning; any such change must be declared before inspecting results.

## Configuration reference

The reference file is `config/tenishev2014_stage0.json`.

### Top-level physics fields

| Field | Meaning |
|---|---|
| `schema_version` | Input schema revision. |
| `physics` | Must be `gravity_free_collisionless`. |
| `frame` | Must be `IAU_ENCELADUS`. |
| `reference_radius_m` | Spherical source/reference radius. |
| `molecular_species` | Documented species; `H2O` in the reference case. |
| `molecular_mass_kg` | Molecular mass used in thermal speed. |
| `source_law` | `tenishev2010_point` or `surface_flux_maxwellian`. |
| `source_case` | Default source geometry. |
| `jet_catalog_file` | Optional external catalog path. |
| `stripe_geometry_file` | Path relative to the configuration file. |
| `shape_parameters` | Optional bounded nonlinear parameters. |

### Discrete jet fields

Each entry requires `id`, `latitude_deg`, `west_longitude_deg`, `tilt_deg`,
`azimuth_deg`, `temperature_K`, and `bulk_speed_m_s`. A per-source `source_law`
may override the global law.

### Stripe fields

Each entry specifies `id`, `geometry_id`, `temperature_K`, and
`bulk_speed_m_s`. `geometry_id` is resolved to a checked polyline at load time.

### Background fields

```json
"background": {
  "model": "both",
  "power_law_exponent": -2.0
}
```

For infinite ray integration, a power-law exponent must be below -1 so the
column converges.

### Example nonlinear parameter

The reference fit does not enable this, but a copied configuration can contain:

```json
"shape_parameters": [
  {
    "source_id": "I_baghdad",
    "field": "bulk_speed_m_s",
    "lower": 400.0,
    "upper": 1400.0
  }
]
```

Use such freedom only with a declared physical prior and sufficient data.

## Installation

The application requires Python 3, NumPy, SciPy, Matplotlib, and Pillow. It was
tested with Python 3.8.10, NumPy 1.24.4, SciPy 1.10.1, and Pillow 10.4.0.

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus/standalone
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

The SPICE C++ tools additionally require a C++11 compiler and CSPICE. The
default `CSPICE_ROOT` is `/home/vtenishe/SPICE/cspice`.

## Running one INMS fit

The simplest command fits the default hybrid model to E3:

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus/standalone

python3 plume_fit.py \
  --config config/tenishev2014_stage0.json \
  --observations /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.csv \
  --metadata /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.metadata.json \
  --output /tmp/e3_hybrid_fit.json
```

The E5 command substitutes `e5` in the two input names and output name.
`plume_fit.py` refuses to overwrite an existing output file.

Select another source case explicitly:

```bash
python3 plume_fit.py \
  --config config/tenishev2014_stage0.json \
  --source-case porco2014_98 \
  --jet-catalog /data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json \
  --observations /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.csv \
  --metadata /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.metadata.json \
  --output /tmp/e3_porco98_fit.json
```

Replace `porco2014_98` with `portyankina2022_solar_41` for the 41-jet case.

### Single-fit JSON output

The result includes:

- resolved configuration and source labels;
- component amplitudes and source rates in molecules/s;
- background coefficients;
- prediction and normalized residual arrays;
- chi-square and nominal reduced chi-square;
- formal/effective ranks and singular values;
- condition number and fit warnings;
- optional shape-optimization diagnostics;
- observation metadata; and
- SHA-256 hashes of configuration, data, metadata, source geometry, catalog,
  and relevant software.

## Running the complete E3/E5/UVIS study

Run one source case into a new directory:

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus/standalone

python3 joint_fitting_study.py \
  --config config/tenishev2014_stage0.json \
  --source-case tenishev2014_hybrid_8 \
  --jet-catalog /data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json \
  --data-root /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1 \
  --output-dir /data/vtenishe/cassini_validation_data/studies/my_new_study/tenishev2014_hybrid_8 \
  --quadrature-order 192
```

The output directory must not already exist. The catalog option is unnecessary
for the hybrid case but can be left in a common scripted command.

Run all three cases:

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus/standalone

MODEL_CONFIG=config/tenishev2014_stage0.json
FIT_DATA=/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1
JET_CATALOG=/data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json
STUDY_ROOT=/data/vtenishe/cassini_validation_data/studies/my_new_study
mkdir -p "$STUDY_ROOT"

for SOURCE_CASE in \
  tenishev2014_hybrid_8 \
  porco2014_98 \
  portyankina2022_solar_41
do
  python3 joint_fitting_study.py \
    --config "$MODEL_CONFIG" \
    --source-case "$SOURCE_CASE" \
    --jet-catalog "$JET_CATALOG" \
    --data-root "$FIT_DATA" \
    --output-dir "$STUDY_ROOT/$SOURCE_CASE" \
    --quadrature-order 192
done

python3 compare_source_case_studies.py --study-root "$STUDY_ROOT"
```

Each case directory contains:

| Product | Contents |
|---|---|
| `study_results.json` | Complete vectors, independent and combined fits, fixed-shape strategies, input/software hashes, convergence, and limitations. |
| `component_amplitudes.csv` | Every fitted component for independent and combined vectors. |
| `source_rate_shape_strategies.csv` | Four strategies x four observations, with total rates and fit metrics. |
| `model_comparison.csv` | Full/no-stripe/no-jet/background-only ablations. |
| `cross_prediction_metrics.csv` | Transfer of each fitted vector to every observation. |
| `fit_profiles.png`, `.pdf` | INMS and UVIS profile comparisons. |
| `fitting_study.md` | Generated human-readable interpretation. |

The comparison program creates cross-case CSV, JSON, Markdown, PNG, and PDF
files in `STUDY_ROOT`. It checks case identities and refuses to overwrite an
existing synthesis.

## Orbital phase, tidal forcing, and plume-rate correlation

### What is and is not calculated

The current tidal analysis does **not** calculate a physical stress tensor or
an instantaneous thermodynamic heating rate. It calculates three explicitly
named orbit-phase proxies:

1. the amplitude of degree-2 tide-raising forcing, proportional to
   \((a/r)^3\);
2. an apoapsis-centered sinusoidal crack-opening proxy; and
3. a squared forcing-rate proxy motivated by rate-dependent dissipation.

This distinction is essential in a publication. Saturn's tide-raising
potential is an external load. Love numbers and an interior/shell model convert
that load to displacement and strain; rheology converts strain to stress and
phase lag; fracture orientation converts the tensor to fault-normal and shear
traction; and a hydraulic/vent model converts traction to escaping gas. Heat
production additionally requires the out-of-phase part of the stress--strain
response. None of those mappings is uniquely determined by orbital distance
alone.

### Keplerian orbital coordinate

The independent orbital coordinate is mean anomaly \(M\), measured from
periapsis and increasing uniformly in time:

\[
M(t)=n(t-t_{\mathrm p})\pmod{2\pi},\qquad
n=\frac{2\pi}{P}. \tag{47}
\]

The default eccentricity and period are \(e=0.0047\) and
\(P=1.370218\) days. For each requested \(M\), the eccentric anomaly \(E\) is
the solution of Kepler's equation,

\[
M=E-e\sin E. \tag{48}
\]

The implementation initializes \(E_0=M\) and applies Newton iteration

\[
E_{k+1}=E_k-
\frac{E_k-e\sin E_k-M}{1-e\cos E_k}, \tag{49}
\]

for at most 20 iterations, stopping when the largest correction is below
\(10^{-14}\) rad. Orbital radius and true anomaly \(f\) then follow from

\[
\frac{r}{a}=1-e\cos E, \tag{50}
\]

\[
f=2\operatorname{atan2}\!\left(
\sqrt{1+e}\sin\frac{E}{2},
\sqrt{1-e}\cos\frac{E}{2}\right). \tag{51}
\]

Mean anomaly is used on phase plots because it is linear in time; true anomaly
is retained for geometric transformations. The event phases currently read
from metadata/tidal_event_phases.csv are \(274^\circ\) for E3,
\(200^\circ\) for E5, \(236^\circ\) for the 2007 zeta Orionis occultation,
and \(98^\circ\) for the 2010 E10 solar occultation. They are literature table
values. A final expanded analysis should recompute osculating phases from
SPICE at every observation epoch and use the table only as an independent
check.

The alternative second-focus coordinate is the unsigned angle at Saturn's
occupied focus \(F_1\), between the direction to the empty focus \(F_2\) and
the direction to Enceladus:

\[
\alpha_{F_2}=\arccos[-\cos f]. \tag{52}
\]

It is \(0^\circ\) at apoapsis and \(180^\circ\) at periapsis. Because equation
(52) folds the pre- and post-apoapsis halves of the orbit onto one interval,
the plotting code shows those branches separately. This coordinate is useful
for geometric visualization but discards causal timing information.

### Tide-raising amplitude

The degree-2 tide-raising potential at a surface location can be written,
apart from a convention-dependent sign and the removal of degree-zero and
degree-one terms, as

\[
U_2(\psi,t)=
\frac{GM_{\mathrm S}R_E^2}{r(t)^3}P_2(\cos\psi),\qquad
P_2(x)=\frac{1}{2}(3x^2-1), \tag{53}
\]

where \(M_{\mathrm S}\) is Saturn's mass and \(\psi\) is angular distance from
the Saturn-facing direction. The implemented orbit-only amplitude removes the
constant prefactor and surface angular dependence:

\[
\boxed{
\mathcal F(M)=\left(\frac{a}{r(M)}\right)^3
=\left(1-e\cos E\right)^{-3}
}. \tag{54}
\]

\(\mathcal F\) is largest at periapsis and smallest at apoapsis. It measures
the magnitude of equilibrium tidal forcing, not the time-varying differential
stress relative to a reference figure and not heating. On plots, this proxy is
min--max normalized over a dense 0--360 degree grid:

\[
\mathcal F_{\mathrm{norm}}=
\frac{\mathcal F-\min\mathcal F}
{\max\mathcal F-\min\mathcal F}. \tag{55}
\]

### Crack-opening proxy

The default quantity correlated with plume production is an intentionally
simple one-harmonic opening response,

\[
\boxed{
\mathcal O(M;\delta)=
\frac{1-\cos(M-\delta)}{2}
}. \tag{56}
\]

At zero lag, \(\mathcal O=0\) at periapsis and \(\mathcal O=1\) at apoapsis.
A positive \(\delta\) shifts maximum opening to
\(M=180^\circ+\delta\), representing a delayed response. The default
\(\delta=0^\circ\) is fixed rather than optimized on the four available
epochs. A literature-motivated approximately \(56^\circ\) lag can be evaluated
as a declared sensitivity case, but choosing the lag that maximizes the
present four-point correlation would be circular.

Equation (56) is not obtained by projecting an elastic stress tensor onto a
tiger-stripe fault. It represents the hypothesis that orbitally modulated
extension opens a hydraulically connected fracture most strongly near
apoapsis. It is useful for a baseline empirical test precisely because it has
only one optional phase parameter.

### Flexing-rate proxy

A heuristic rate-sensitive proxy starts from equation (54). Differentiating
Kepler's equation gives

\[
\frac{dE}{dM}=\frac{1}{1-e\cos E}. \tag{57}
\]

Consequently,

\[
\frac{d\mathcal F}{dM}=
-\frac{3e\sin E}{(1-e\cos E)^5}, \tag{58}
\]

and the implemented flexing-power proxy is

\[
\boxed{
\mathcal H_{\mathrm{flex}}(M)=
\left[\frac{d\mathcal F}{dM}\right]^2
=\frac{9e^2\sin^2E}{(1-e\cos E)^{10}}
}. \tag{59}
\]

Because \(dM/dt=n\), a dimensional squared-rate quantity would be
\(n^2\mathcal H_{\mathrm{flex}}\). The constant \(n^2\) cancels under the
normalization used for correlation, so the script does not multiply by it.
This proxy vanishes at both apses and has two maxima per orbit. It resembles
dissipation proportional to squared deformation rate in a simple viscous
element, but it is not a Maxwell, Burgers, Andrade, ocean, porous-core, or
turbulent-slot heating solution.

### What a physical stress and heating calculation would require

For a specified layered interior, complex degree-2 Love numbers convert
equation (53) to surface displacement. Schematically,

\[
u_r=\frac{h_2}{g}U_2,\qquad
\mathbf u_{\mathrm h}=\frac{l_2}{g}\nabla_{\!s}U_2, \tag{60}
\]

where \(h_2\) and \(l_2\) depend on forcing frequency, shell thickness,
viscosity, ocean structure, and core properties. Surface or three-dimensional
strain follows from spatial derivatives of \(\mathbf u\); a constitutive law
then gives stress. For linear isotropic elasticity,

\[
\sigma_{ij}=2\mu_{\mathrm e}\epsilon_{ij}
+\lambda_{\mathrm e}\epsilon_{kk}\delta_{ij}, \tag{61}
\]

while a viscoelastic calculation replaces the elastic moduli by complex,
frequency-dependent moduli or evolves the constitutive law in time.

For a mapped fracture segment with horizontal unit normal \(\mathbf m\) and
along-strike unit vector \(\mathbf t\), the relevant resolved tractions would
include

\[
\sigma_{\mathrm n}=\mathbf m^{\mathsf T}\boldsymbol\sigma\mathbf m,
\qquad
\tau_{\mathrm{strike}}=\mathbf t^{\mathsf T}\boldsymbol\sigma\mathbf m. \tag{62}
\]

Normal opening, Coulomb failure, and strike-slip hypotheses use different
functions of these quantities. Optical and physical libration, the changing
Saturn direction in the body-fixed frame, nonsynchronous components, and the
reference equilibrium figure must also be specified. The current script
implements none of equations (60)--(62); therefore its output must be called a
forcing or opening *proxy*, not calculated tidal stress.

The commonly quoted orbit-averaged eccentricity-tide dissipation is

\[
\left\langle\dot E\right\rangle=
\frac{21}{2}\frac{k_2}{Q}
\frac{n^5R_E^5e^2}{G}
=\frac{21}{2}\frac{k_2}{Q}
\frac{GM_{\mathrm S}^2R_E^5ne^2}{a^6}, \tag{63}
\]

to the usual low-eccentricity synchronous approximation. Equation (63) is a
cycle average and cannot supply the phase curve used in the present
correlation. Locally, irreversible heating depends on the out-of-phase
stress--strain work,

\[
\langle q\rangle=
\left\langle\sigma_{ij}\dot\epsilon_{ij}\right\rangle, \tag{64}
\]

so stress amplitude without rheological phase information is insufficient.
Volumetric shell heating, ocean/core dissipation, and turbulent dissipation in
water-filled slots should be calculated and validated separately.

### Correlation with fitted source rate

For each source geometry and each of the four spatial-shape strategies,
equation (46) supplies four epoch-specific total localized rates \(Q_j\).
Empirical backgrounds are excluded. Given a selected proxy \(X_j\), the
reported Pearson coefficient is

\[
r_{\mathrm P}=
\frac{\sum_j(X_j-\overline X)(Q_j-\overline Q)}
{\sqrt{\sum_j(X_j-\overline X)^2
\sum_j(Q_j-\overline Q)^2}}, \tag{65}
\]

and Spearman \(\rho\) is the same correlation applied to ranks. Two-sided
p-values are calculated by SciPy. The plotted guide line is the ordinary
unweighted least-squares relation

\[
Q_j=\beta_0+\beta_1X_j+\varepsilon_j. \tag{66}
\]

The correlation and guide line do not use the plotted vertical uncertainties
as statistical weights. Those bars contain only common absolute-scale
systematics: 23% for each INMS reconstruction and the published fractional
column uncertainty for each UVIS event. They omit source-map degeneracy,
transport-model discrepancy, cross-instrument calibration, proxy uncertainty,
and epoch-to-epoch stochastic jet variability.

With \(N=4\) heterogeneous epochs, equations (65)--(66) are exploratory
diagnostics. The three source geometries and four shape-training strategies
are retained separately to expose structural sensitivity. No lag, proxy, or
source geometry should be selected by maximizing this same four-point
correlation and then reported as independently significant.

## Running the tidal-proxy analysis

The tidal analysis is stored with the external Cassini data because it
combines fit products, orbital phases, and literature metadata:

```bash
cd /data/vtenishe/cassini_validation_data

python3 scripts/plot_tidal_proxy_correlation.py \
  --fit-study-root /data/vtenishe/cassini_validation_data/studies/source_shape_tidal_20261014 \
  --output-dir /data/vtenishe/cassini_validation_data/studies/my_tidal_analysis \
  --proxy opening \
  --phase-lag-deg 0
```

Available proxies are `opening`, `forcing`, and `flexing_power`. Do not tune a
phase lag on the present four epochs and then interpret the optimized
correlation as an independent test.

## Rebuilding SPICE geometry

Build both C++ tools:

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus/standalone
make CSPICE_ROOT=/home/vtenishe/SPICE/cspice
```

Build only one:

```bash
make spice-ground-tracks CSPICE_ROOT=/home/vtenishe/SPICE/cspice
make spice-observer-ranges CSPICE_ROOT=/home/vtenishe/SPICE/cspice
```

Ground-track interface:

```text
build/spice_ground_tracks LSK PCK SPK TAG UTC_REFERENCE T0 T1 DT OUTPUT.csv
```

Rebuild the 2010 UVIS observer geometry using its metakernel:

```bash
cd /data/vtenishe/cassini_validation_data/data/raw/spice/metakernels

/home/vtenishe/Mars2/AMPS/srcEnceladus/standalone/build/spice_observer_ranges \
  cas_uvis_2010_e10_solar.tm \
  /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/uvis_2010_e10_solar_optical_depth.csv \
  /tmp/uvis_2010_e10_solar_observer_range.csv
```

Generate into a new file and compare it with the checked external product;
do not silently replace provenance evidence.

## Tests and verification

Run the complete Python test suite with external-data checks enabled:

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus/standalone

ENCELADUS_PLUME_FIT_DATA=/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1 \
ENCELADUS_JET_CATALOG=/data/vtenishe/cassini_validation_data/data/derived/porco2014_jet_catalog.json \
  python3 -m unittest discover -s tests -v
```

Perform syntax checks:

```bash
python3 -m py_compile \
  plume_model.py fitting.py plume_fit.py joint_fitting_study.py \
  compare_source_case_studies.py
```

The tests cover:

- analytic velocity moments against independent quadrature;
- total flux normalization for both source laws;
- inverse-square behavior and the outward-horizon boundary;
- coordinate convention checks;
- uniform line-source weighting;
- analytic inverse-square background columns;
- finite and infinite UVIS ray operators;
- all three source-case memberships and directions;
- synthetic nonnegative recovery and unobservable columns;
- input metadata enforcement and raw-count rejection;
- published INMS profile maxima;
- combined-fit dataset balancing;
- fixed-shape rate recovery with separately fitted backgrounds;
- source/stripe reconstruction and figure fixtures; and
- SPICE geometry regression products when the required local resources exist.

These tests verify implementation consistency and mathematical conservation.
They do not validate the collisionless physical approximation or prove that a
published source geometry is unique.

## Adding a new observation

1. Store raw and processed data under the external Cassini validation-data
   tree, not in this source directory.
2. Record original source, units, frame, calibration, selection, uncertainty
   semantics, transformations, and SHA-256 values.
3. Add a dedicated observation operator if the observable is not instantaneous
   H2O number density or the currently supported UVIS optical-depth shape.
4. Use SPICE to compute geometry in `IAU_ENCELADUS`; do not copy a plotted
   ground track when kernels are available.
5. Add regression and physical-limit tests before fitting.
6. Decide and document whether the observation is for development, shape
   training, calibration, or held-out validation before viewing its residuals.
7. Extend the study's dataset registry, balancing, plots, and provenance
   output rather than silently concatenating rows.

## Adding a new source model

1. Preserve the `IAU_ENCELADUS` signed-latitude and west-longitude convention.
2. Store published source catalogs outside the code when they are validation
   data products; include an exact selection and provenance record.
3. Declare temperature, speed, direction, and source-law assumptions for every
   source family.
4. Ensure each fitted source is normalized to a physical total `Q`.
5. Add the case to `SOURCE_CASES` and materialize it in `load_config()`.
6. Test exact source membership, component count, direction mode, and
   kinematics.
7. Inspect rank, singular values, and condition number before interpreting
   individual amplitudes.

## Handoff to later AMPS stages

Stage 0 should hand off an ensemble rather than a single supposedly unique
vent map:

- alternative published source geometries;
- total-rate ranges under the four spatial-shape strategies;
- common calibration-scale uncertainties;
- background alternatives;
- source-law and kinematic sensitivity cases; and
- explicit residual/identifiability diagnostics.

Stage 1 can use these as boundary-source alternatives in the neutral AMPS
calculation and compare pickup-ion/plasma consequences with BATS-R-US. Stage 2
can replace the analytic neutral and pickup-ion treatment with fully coupled
AMPS components. Stage 3 adds charged dust. None of those stages should treat
the Stage-0 NNLS partition among nearly dependent jets as ground truth.

## Primary references

- Tenishev et al. (2010), *An approach to numerical simulation of the gas
  distribution in the atmosphere of Enceladus*, DOI
  `10.1029/2009JA015223`.
- Teolis et al. (2010), *Detection and measurement of ice grains and gas
  distribution in the Enceladus plume by Cassini's Ion Neutral Mass
  Spectrometer*, DOI `10.1029/2009JA015192`.
- Tenishev et al. (2014), *Effect of the Tiger Stripes on the water vapor
  distribution in Enceladus' exosphere*, DOI `10.1002/2014JE004700`.
- Porco, DiNino, and Nimmo (2014), *How the Geysers, Tidal Stresses, and
  Thermal Emission across the South Polar Terrain of Enceladus Are Related*,
  DOI `10.1088/0004-6256/148/3/45`.
- Portyankina et al. (2022), *Modeling the complete set of Cassini's UVIS
  occultation observations of Enceladus' plume*, DOI
  `10.1016/j.icarus.2022.114918`.
- Teolis et al. (2015), INMS sensitivity and calibration review, DOI
  `10.1007/s11214-014-0133-8`.
- Hansen et al. (2011), 2010 E10 UVIS solar-occultation H2O column, DOI
  `10.1029/2011GL047415`.
- Hansen et al. (2020), complete Cassini UVIS occultation assessment, DOI
  `10.1016/j.icarus.2019.113461`.
- Hurford et al. (2007), *Eruptions arising from tidally controlled periodic
  openings of rifts on Enceladus*, DOI `10.1038/nature05821`.
- Olgin, Smith-Konter, and Pappalardo (2011), *Limits of Enceladus's ice shell
  thickness from tidally driven tiger stripe shear failure*, DOI
  `10.1029/2010GL044950`.
- Hedman et al. (2013), *An observed correlation between plume activity and
  tidal stresses on Enceladus*, DOI `10.1038/nature12371`.
- Nimmo, Porco, and Mitchell (2014), *Tidally modulated eruptions on
  Enceladus: Cassini ISS observations and models*, DOI
  `10.1088/0004-6256/148/3/46`.
- Běhounková et al. (2015), *Timing of water plume eruptions on Enceladus
  explained by interior viscosity structure*, DOI `10.1038/ngeo2475`.
- Kite and Rubin (2016), *Sustained eruptions on Enceladus explained by
  turbulent dissipation in tiger stripes*, DOI
  `10.1073/pnas.1520507113`.
- Nimmo et al. (2023), *Origin and Evolution of Enceladus's Tidal
  Dissipation*, DOI `10.1007/s11214-023-01007-4`.
- Berne et al. (2024), *Jet activity on Enceladus linked to tidally driven
  strike-slip motion along tiger stripes*, DOI
  `10.1038/s41561-024-01418-0`.
