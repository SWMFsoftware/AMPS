# Neutral-collision data review for the Halley species set

## Scope

The Halley input compiles 27 neutral species:

```text
H, H2, C, CH, CH2, N, CH3, NH, CH4, NH2, O, NH3, OH, H2O,
C2, C2H, C2H2, CN, HCN, C2H4, CO, C2H5, HCO, C2H6, H2CO,
CH3OH, CO2
```

AMPS needs a cross section that controls neutral-neutral DSMC collision
frequency and momentum transfer. Electron-impact, ion-neutral charge-exchange,
photodissociation, and purely reactive rate-coefficient databases do not
provide that quantity. State-to-state rotational de-excitation data are useful
for an internal-energy model, but they are not a replacement for a total or
momentum-transfer collision cross section.

There are 378 unordered pairs, including self pairs, for 27 species. No online
source found in this review provides measured or ab-initio total
momentum-transfer cross sections for the complete 378-pair set at cometary
temperatures.

## Directly applicable cometary hard-sphere data

[Rubin et al. (2011)](https://doi.org/10.1016/j.icarus.2011.04.006),
Table 3, provides all 66 unordered pairs among the 11 Halley parent species.
Those values are already transcribed in `input/CometPhysicsDatabase.pm`.
Rubin states that the CO/H2O/daughter values come from Combi (1996), that the
water self-collision law comes from Crifo (1989), and that the remaining parent
values are estimates. The paper does not provide the daughter-daughter matrix.

[Fougere's dissertation](https://deepblue.lib.umich.edu/bitstream/2027.42/110338/1/fougere_1.pdf),
Table 2.1, reproduces the 21-pair Combi (1996) hard-sphere set for
H2O/OH/H2/H/O/CO. Cross sections below are in `cm^2`:

| Pair | Cross section | Pair | Cross section |
|---|---:|---|---:|
| H2O-H2O | `1.66e-15 (T/300 K)^-0.6` | H2O-OH | `3.2e-15` |
| H2O-H2 | `3.2e-15` | H2O-H | `1.8e-15` |
| H2O-O | `1.8e-15` | H2O-CO | `3.2e-15` |
| OH-OH | `3.0e-15` | OH-H2 | `3.0e-15` |
| OH-H | `1.5e-15` | OH-O | `1.5e-15` |
| OH-CO | `3.0e-15` | H2-H2 | `3.0e-15` |
| H2-H | `1.5e-15` | H2-O | `1.5e-15` |
| H2-CO | `3.0e-15` | H-H | `1.2e-15` |
| H-O | `1.2e-15` | H-CO | `1.5e-15` |
| O-O | `1.2e-15` | O-CO | `1.5e-15` |
| CO-CO | `3.2e-15` |  |  |

At the Halley collision-temperature floor of 180 K, the Crifo water law gives
`2.2554e-15 cm^2`. [Huang et al. (2016)](https://doi.org/10.3847/1538-4357/833/2/160),
Table 2, independently prints ten of the same H2O/OH/H2/H/O values; it does not
extend the species coverage. Combining Rubin and Combi covers 84 distinct
pairs because H2O-H2O, H2O-CO, and CO-CO occur in both sets. Thus 294 of 378
pairs still lack a directly published cometary hard-sphere value.

## Complete Lennard-Jones parameter coverage

The [Ansys CHEMKIN transport database](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/pdf/Ansys_Chemkin_Input_Manual.pdf),
Table 5.1, contains Lennard-Jones parameters for every compiled species. The
manual defines `epsilon/k` as the potential well depth and `d` as the
Lennard-Jones collision diameter. It also states that this historical database
mixes best fits to experimental macroscopic transport measurements, estimates
following Svehla, and computational results. It does not identify which method
was used for each row, so none of these individual entries can be labeled
"measured" without another species-specific source.

The table below transcribes those parameters. `pi*d^2` is a derived geometric
self-collision estimate, not a value printed by CHEMKIN. One square angstrom is
`1e-16 cm^2`.

| Species | epsilon/k (K) | d (angstrom) | pi*d^2 (`1e-15 cm^2`) |
|---|---:|---:|---:|
| H | 145.0 | 2.050 | 1.320 |
| H2 | 38.0 | 2.920 | 2.679 |
| C | 71.4 | 3.298 | 3.417 |
| CH | 80.0 | 2.750 | 2.376 |
| CH2 | 144.0 | 3.800 | 4.536 |
| N | 71.4 | 3.298 | 3.417 |
| CH3 | 144.0 | 3.800 | 4.536 |
| NH | 80.0 | 2.650 | 2.206 |
| CH4 | 141.4 | 3.746 | 4.408 |
| NH2 | 80.0 | 2.650 | 2.206 |
| O | 80.0 | 2.750 | 2.376 |
| NH3 | 481.0 | 2.920 | 2.679 |
| OH | 80.0 | 2.750 | 2.376 |
| H2O | 572.4 | 2.605 | 2.132 |
| C2 | 97.53 | 3.621 | 4.119 |
| C2H | 209.0 | 4.100 | 5.281 |
| C2H2 | 209.0 | 4.100 | 5.281 |
| CN | 75.0 | 3.856 | 4.671 |
| HCN | 569.0 | 3.630 | 4.140 |
| C2H4 | 280.8 | 3.971 | 4.954 |
| CO | 98.1 | 3.650 | 4.185 |
| C2H5 | 252.3 | 4.302 | 5.814 |
| HCO | 498.0 | 3.590 | 4.049 |
| C2H6 | 252.3 | 4.302 | 5.814 |
| H2CO (`CH2O` in CHEMKIN) | 498.0 | 3.590 | 4.049 |
| CH3OH | 481.8 | 3.626 | 4.131 |
| CO2 | 244.0 | 3.763 | 4.449 |

The open [Stanford FFCM-2 transport table](https://web.stanford.edu/group/haiwanglab/FFCM2/docs/TrialModel/Transport/)
and [Cantera GRI-Mech data](https://github.com/Cantera/cantera/blob/main/data/gri30.yaml)
provide an independent machine-readable transcription of the same Sandia
transport family for 26 species. GRI-Mech omits C2; the CHEMKIN table supplies
its `epsilon/k=97.53 K` and `d=3.621 angstrom` entry.

The repeated parameters show where the database uses surrogate estimates:
C and N are identical; CH and OH are identical; CH2 and CH3 are identical;
NH and NH2 are identical; C2H and C2H2 are identical; C2H5 and C2H6 are
identical; and HCO and H2CO are identical. They should not be presented as 27
independent measurements.

## Transparent all-pair estimates

A hard-sphere matrix compatible with the current AMPS collision interface can
be generated from the CHEMKIN diameters using the Lorentz arithmetic rule:

```text
d_ij = (d_i + d_j)/2
sigma_HS,ij = pi*d_ij^2
```

This fills all 378 pairs and is traceable, but it is an estimate. For example,
the derived O-O value is `2.376e-15 cm^2`, whereas Combi uses
`1.2e-15 cm^2`; the H2O-H2O geometric value is `2.132e-15 cm^2`, close to the
Crifo value `2.255e-15 cm^2` at 180 K. The agreement is therefore
species-dependent and the estimates must not overwrite published values.

A temperature-dependent Lennard-Jones transport model would additionally use
the Berthelot well-depth rule

```text
epsilon_ij/k = sqrt((epsilon_i/k)*(epsilon_j/k))
Tstar_ij = T/(epsilon_ij/k)
```

and a selected transport collision integral. [Neufeld, Janzen, and Aziz
(1972)](https://doi.org/10.1063/1.1678363) provide fitted Lennard-Jones
collision integrals for `0.3 <= Tstar <= 100`. At 180 K the CHEMKIN parameters
put every Halley pair inside that stated range (`Tstar=0.315` to `4.74`). A
collision integral calibrates diffusion or viscosity; it is not by itself the
total DSMC encounter cross section. Implementing it requires an explicit
choice of diffusion-equivalent, viscosity-equivalent, variable-hard-sphere,
or variable-soft-sphere mapping and treatment of polar species.

[Svehla's NASA TR R-132](https://ntrs.nasa.gov/api/citations/19630012982/downloads/19630012982.pdf)
is the historical estimation basis named by CHEMKIN. It covers 100--5000 K,
includes radicals, derives force constants from experimental viscosity when
available, and otherwise estimates them. That temperature coverage supports
using the parameters at 180 K, but it does not remove the pair-model choice.

## Detailed inelastic data that must remain separate

- [BASECOL](https://basecol.org/) contains state-to-state collisional
  excitation data for many relevant targets, including C, O, H2, CH, OH, CN,
  NH, CO, CO2, HCN, H2O, NH3, CH3OH, H2CO, C2H, CH2, HCO, and NH2. Most
  colliders are H, H2, or He. These rates can support a future internal-state
  model but do not form the required 27-by-27 total cross-section matrix.
- [NASA/GISS molecular collision data](https://data.giss.nasa.gov/mcrates/)
  includes H2O-H2 and comet-specific CO-H2O rotational excitation rates. These
  are likewise internal-energy-transfer data.
- [Loreau et al. (2020)](https://academic.oup.com/mnras/article/498/4/5489/5902862)
  calculate H2O-H2O rotational de-excitation cross sections from 100 to 800 K.
  Individual values depend on the initial and final rotational states and are
  not a total elastic cross section.
- Crossed molecular-beam measurements exist for state-resolved H2O-H2
  rotationally inelastic scattering at 44.8 and 71.2 meV
  ([Yang et al.](https://doi.org/10.1063/1.3589360)). Those collision energies
  are above the 180 K thermal scale and the measurements do not supply the
  missing all-pair momentum-transfer matrix.

## Selected and implemented model

The literature supports a source-priority database:

1. retain Rubin Table 3 for all parent-parent pairs;
2. add the remaining non-overlapping Combi pairs for H2O/OH/H2/H/O/CO; and
3. use a clearly labeled estimate only for the remaining 294 pairs.

The Halley input explicitly selects step 3 as the CHEMKIN/Lorentz geometric
hard-sphere rule:

```text
NeutralCollisionCompletion=COMBI1996_CHEMKIN_LORENTZ_HS
```

`input/CometPhysicsDatabase.pm` implements the priority order without
overwriting an existing pair. Its generated 27-species matrix contains 66
Rubin rows, 18 additional Combi rows, and 294 CHEMKIN/Lorentz estimates. Each
row retains its source tier. The temperature-dependent Lennard-Jones
transport-equivalent alternative remains unimplemented because it would
require a different explicit transport mapping, especially for polar
molecules. The completion is case-selected rather than an implicit fallback.
