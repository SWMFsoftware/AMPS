# Comet 1P/Halley parameters from Rubin et al. (2011)

This document records the provenance and limits of
`input/halley1d.input`. The source is Rubin et al.,
“Monte Carlo modeling of neutral gas and dust in the coma of Comet 1P/Halley,”
Icarus 213 (2011) 655–677, doi:10.1016/j.icarus.2011.04.006. The local PDF
`~/Downloads/Rubin-2011-icarus.pdf` was checked directly on 2026-10-27.

The file is named for **Halley**. The Hartley 2 cases remain
`input/hartley2.input` and `input/hartley2.radial.input`. Both use the same
generic application code in the historically named `srcHartley2` directory.

## Geometry and boundaries

| Quantity | Paper value | Application representation |
| --- | ---: | --- |
| Symmetry | 1D spherical | `Dimension=1`, `Geometry=Spherical` |
| Nucleus radius | 6 km | `NucleusRadius=6000 m` |
| Inner boundary | top of Knudsen layer, a few metres above the surface | represented at the stated 6 km coordinate; no exact offset is published |
| Outer boundary | `10^6 km` cometocentric | `DomainOuterRadius=1.0e9 m` |
| Heliocentric distance | 0.90 AU | `HeliocentricDistanceAU=0.90` |
| Parent injection temperature | 180 K | every parent `SurfaceSource` |
| Initial parent bulk speed | 0 m/s | every parent `SurfaceSource` |
| Surface distribution | uniform | `UNIFORM` |
| Returning gas | fully thermalized and re-emitted | DIM1 inner-boundary callback at 180 K |
| Outer crossing | particle removed | DIM1 outer-boundary callback |

The paper says the near-surface collisional mean free path is “on the order of
a fraction of a meter.” It does not publish an exact first-cell width, radial
cell count, AMR hierarchy, or inner-boundary offset. The input uses the
user-supplied 2,000 ideal cells and the reproducible mean-free-path estimate
below; these are numerical case choices rather than paper values.

## Parent composition and rates

Table 2 gives abundances relative to water. Multiplying them by the stated
`Q(H2O)=7.0e29 s^-1` produces the rates below. Their sum is
`8.54e29 s^-1`, matching the paper.

| Parent | Relative abundance | Rate, s^-1 |
| --- | ---: | ---: |
| H2O | 100% | `7.00e29` |
| CO | 13% | `9.10e28` |
| CO2 | 2.5% | `1.75e28` |
| CH3OH | 1.7% | `1.19e28` |
| H2CO | 1.5% | `1.05e28` |
| NH3 | 1.5% | `1.05e28` |
| CH4 | 0.6% | `4.20e27` |
| C2H6 | 0.4% | `2.80e27` |
| C2H4 | 0.3% | `2.10e27` |
| C2H2 | 0.3% | `2.10e27` |
| HCN | 0.2% | `1.40e27` |

The neutral network tracks 27 species:

```text
H H2 C CH CH2 N CH3 NH CH4 NH2 O NH3 OH H2O C2 C2H C2H2 CN HCN
C2H4 CO C2H5 HCO C2H6 H2CO CH3OH CO2
```

## Surface mean-free-path estimate

An estimate can be calculated without attributing a nonexistent exact value
to the paper. For each zero-drift half-Maxwellian parent,

```text
Gamma_i = Q_i/(4 pi R^2)
<v_r,+> = sqrt(2 k T/(pi m_i))
n_i = Gamma_i/<v_r,+>.
```

At `R=6000 m` and `T=180 K`, the 11 parent densities sum to
`8.598305522e18 m^-3`. Weighting the Table 3 H2O-parent cross sections by
those densities gives `2.496830337e-19 m^2`; the H2O-H2O component uses the
published Crifo law at 180 K. The conventional hard-sphere expression

```text
lambda = 1/(sqrt(2) n sigma)
```

then gives `lambda ~= 0.329 m`. Omitting the equilibrium `sqrt(2)` factor gives
`0.466 m`. The outward multi-species gas is not an equilibrium single-species
Maxwellian, so these are documented resolution estimates, consistent with the
paper’s “fraction of a meter,” rather than an exact mesh prescription.

The AMPS boundary routine samples a flux-weighted Maxwellian. Residence-time
weighting of those particles in the first cell recovers the outward spatial
half-Maxwellian used in the density formula above; an extra factor of two is
therefore not applied. `halley1d.input` uses
`NucleusSurfaceResolution=0.3293692849882 m`.

## Derived numerical controls

The derivation can be repeated with:

```text
python3 srcHartley2/halley2011/derive_halley1d_controls.py
```

For the 2,000-cell geometric profile, the result is:

| Control | Value | Basis |
| --- | ---: | --- |
| first ideal cell | `0.3293692849882 m` | surface H2O mean free path above |
| geometric width ratio | `1.008573952117` | 2,000 widths spanning 6,000 m to `1.0e9 m` |
| final ideal cell | `8.501013718e6 m` | geometric profile |
| maximum AMR refinement | `29` | first dyadic AMR cell no wider than the target |
| allocated AMR blocks | `373` | clean one-process mesh construction |
| allocated AMR cells | `2,984` | dyadic approximation to the 2,000-cell ideal profile |
| smallest allocated cell | `0.2328292 m` | level-29 AMR cell width |
| Courant number | `0.5` | supplied numerical control |
| direct-parent baseline speed | `400 m/s` | supplied numerical control |
| reference injection count | `16` | at least 50 parent particles in the surface cell |
| particle-buffer capacity | `5,000,000` | initial capacity consistent with the paper’s “several 10^6”; monitor during convergence |
| sampling-window length | `500 iterations` | supplied numerical control |
| hard iteration ceiling | `10,000 iterations` | supplied numerical control |
| convergence tolerance | `0.05` | supplied control; volume-weighted relative L1 density-profile change |
| consecutive converged windows | `3` | supplied control; applies to every compiled species |

For a zero-drift parent, the expected surface-cell population is

```text
N_cell = N_ref V_char / (C sqrt(2 k T/(pi m))).
```

The lightest direct parent, CH4, is limiting: `N_ref=16` gives `52.454`
particles. The other direct parents range from `54.068` for NH3 to `86.985`
for CO2. This estimate applies to directly injected parents. Daughter
populations are set by their production and loss rates and cannot be forced to
50 particles per cell through the surface injection count.

Across the 2,984 allocated cells, the no-loss occupancy estimate sums to
approximately `2.21e6` direct-parent model particles. The five-million-slot
buffer leaves about `2.79e6` slots for transient excess and daughters, but its
high-water mark must be measured in the convergence run.

The input uses species-local particle weights proportional to local timestep.
This prevents model-particle occupancy from growing with the width of the
logarithmic cells. Individual AMPS weight corrections preserve the physical
weight when a reaction changes species.

The reaction-product timestep speeds start from `400 m/s` for every directly
injected parent. For each Table 4 channel, the maximum center-of-mass recoil
of product `i` is estimated from energy and momentum conservation as

```text
u_i = sqrt(2 E M_other/(m_i M_total)).
```

The script propagates `V_parent + u_i` through successive reactions and rounds
up to a whole metre per second. A species that is injected at 400 m/s but is
also produced by a faster dissociation receives the larger propagated value.
This is a conservative timestep scale for the implemented energy-sharing
algorithm. A Maxwellian has no finite strict maximum, so runtime displacement
statistics still need to be checked.

## Collisions and chemistry

Table 3 supplies 66 hard-sphere cross sections for all unordered pairs of the
11 parent species. The water self-collision law is

```text
sigma_H2O,H2O = 1.66e-19 (T/300 K)^-0.6 m^2
```

and the paper limits it to 50–300 K. The generated database enforces that
interval. The other listed parent-pair values are constant. Molecules and
atoms scatter over `4*pi`.

Table 3 does **not** contain the full collision matrix for the 16 daughter-only
species. The application deliberately reports a missing pair rather than
assigning zero or inventing a cross section. A chemically active physical run
therefore needs a cited daughter collision dataset or an explicitly approved
mixing prescription.

Table 4 has been transcribed into the named `RUBIN2011` chemistry database.
Rates are referenced to 1 AU and scale as inverse heliocentric distance
squared. Electronic states are collapsed into the corresponding neutral
species while retaining the channel excess energy. Ionization channels are
neutral sinks because ions are outside this gas model.

Rubin Table 4 gives four neutral OH photodissociation branches at 1 AU:

| Excess energy | Rate, s^-1 |
| ---: | ---: |
| `0.36 eV` | `3e-6` to `6e-6` |
| `0.67 eV` | `5.0e-7` |
| `2.69 eV` | `1.4e-6` |
| `3.84 eV` | `3.0e-7` |

The three fixed neutral branches sum to `2.2e-6 s^-1`, so the neutral OH
photodissociation total is `5.2e-6` to `8.2e-6 s^-1`. Table 4 separately gives
OH photoionization as `1.2e-5 s^-1`; including that neutral sink makes the
total modeled OH loss `1.72e-5` to `2.02e-5 s^-1` at 1 AU.

Rubin cites Combi (1996) for the four neutral branches. Combi’s tabulation
resolves the low-energy predissociation rate as a heliocentric-velocity and
cross-section-dependent interval, with published alternatives of
`3.0e-6`–`6.1e-6 s^-1` and `2.5e-6`–`5.6e-6 s^-1`. Schleicher and A'Hearn
(1988, ApJ 331, 1058--1077, doi:10.1086/166622) calculate an A-state lifetime
of `1.8e5`–`4.0e5 s` at 1 AU, equivalent to `2.5e-6`–`5.56e-6 s^-1` as the
heliocentric radial velocity changes. They also quote the separate
van Dishoeck--Dalgarno combined higher-state contribution as `1.85e-6 s^-1`
at solar minimum and `3.43e-6 s^-1` at solar maximum, with about 25%
uncertainty. These values describe different adopted solar spectra/channel
groupings and do not select a unique number from Rubin’s input interval.

The Halley input selects `3.0e-6 s^-1`, the lower endpoint of Rubin's printed
range, for the 0.36 eV ground-state OH branch. This is an explicit run choice,
not a value inferred from the excess energy. The methanol-to-formaldehyde
branch is disabled as requested; this is the alternate case described by
Table 4 footnote (a).

## Cooling and omitted physics

Rubin et al. include collisional relaxation and the Crovisier (1984)/Huebner
(1985) water radiative-cooling approximation presented by Gombosi et al.
(1986). The legacy DIM3 Hartley 2 source contains an implementation, but the
generic DIM1 execution path has not yet been connected to it and the current
paper does not reproduce the derivation. A Halley production run must not be
described as complete until this generic cooling path and internal degrees of
freedom are validated.

This Halley case is gas-only. Dust size groups, gas-dust momentum exchange,
gravity-limited grains, and dust instrument comparisons are not enabled.

## Numerical values absent from the paper

The paper reports several million representative particles but does not give
enough information to reconstruct all numerical controls. It does not state:

- radial cell count, exact cell boundaries, or first-cell width;
- timestep/Courant policy or per-species transport-speed bounds;
- statistical particle weights or reference injection-particle count;
- convergence metric, tolerance, or required consecutive-window count; or
- a complete daughter collision matrix.

The input records the user-supplied or derived controls above. Convergence is
evaluated only after a complete 500-iteration sample. For species `s` and
successive windows `k-1` and `k`, the implemented metric is

```text
sum_cells V |n_s,k - n_s,k-1| /
max(sum_cells V |n_s,k|, sum_cells V |n_s,k-1|).
```

The maximum metric over all compiled species must not exceed `0.05` for three
successive comparisons. A species absent from both windows has zero change;
appearance or disappearance has relative change one. The first complete
window establishes the reference and cannot count toward the three passes.

## Figure-comparison scope

Figures 1–21 compare modeled neutral densities, NMS mass channels, water bulk
velocity, remote column densities, and OH aperture counts. Figures 22–23 are
dust results. Tables 5 and 6 describe the NMS instrument response and belong
in a Halley validation case rather than the generic chemistry database.

The PDF does not provide numerical source tables for the plotted curves or
measurements. Reproducing those comparisons will require traceable digitized
data (with digitization uncertainty) or the original datasets, plus the
unresolved physical items above.
