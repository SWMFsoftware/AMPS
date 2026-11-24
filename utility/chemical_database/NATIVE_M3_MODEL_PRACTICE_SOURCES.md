# How published models treat ion creation, breakup and surface return

The closest observational comparison for the Moon is the ARTEMIS pickup-ion
study. The original authors deposited ionization frequencies for 25 species
at 1 au. Other investigations resolve different parts of a reaction depending
on what their instrument measures. The following choices can be reproduced
and compared separately in an exosphere model.

## Neutral sodium and pickup ions

| Investigation and observable | Implemented choice | Use in a comparison |
| --- | --- | --- |
| Tenishev et al., lunar neutral Na density and tail | Section 3.3.3 uses an illuminated loss probability of the form `1 − exp(−Δt/τ)` with a 36 h active-Sun lifetime. The paper also discusses 47 h for the quiet Sun. It tracks the Moon and Earth's shadows. | Reproduce neutral survival and shadow switching under the same illumination and lifetime. Neutral Na observations have no sensitivity to how the photoelectron and ion divide the excess energy after removal. [Article](https://doi.org/10.1016/j.icarus.2013.08.021) |
| Poppe, Halekas & Harada, lunar ARTEMIS high-energy pickup-ion flux | The model uses one ion-production frequency per mechanism and species. Newborn singly charged ions start at rest in the lunar frame and follow cycloids in uniform fields. The observations integrate over ion mass. | Compare total source rates at 1 au and large-scale pickup flux under their assumptions. The archived minimum, median and maximum describe changes in solar or solar-wind conditions. They supply no uncertainty on a single reaction and no electron or recoil spectrum. [Study](https://doi.org/10.1029/2022JE007422) [Deposited rates](https://doi.org/10.5281/zenodo.6994648) |
| Paral et al., Mercury hybrid and particle-tracing model of MESSENGER flyby ion populations | Neutral photoionization weight is accumulated in a mesh cell. New ions sample the local neutral velocity distribution. | Provides a published alternative birth-velocity assumption for testing sensitivity of heavy-ion transport. A draw from the local distribution need not reproduce the velocity of the particular ionized neutral. [Model, paragraph 12](https://doi.org/10.1029/2010GL044413) |

Shen, Halekas & Poppe also adopted zero birth velocity to infer source
locations from 115 ARTEMIS events. They used measured fields for selected
steady intervals and tested multiple ion masses because the instrument lacks
mass discrimination. Poppe et al.'s baseline also underpredicted the filtered
high-energy flux. Agreement in trajectory geometry alone would therefore leave
composition and source strength unresolved. [Shen et al., sections 2–3](https://doi.org/10.3847/1538-4357/ad40a4)
[Poppe et al., section 4.2](https://doi.org/10.1029/2022JE007422)

The archived Na median photoionization frequency, 7.59 × 10⁻⁶ s⁻¹,
corresponds to a 36.6 h lifetime. Its proximity to the earlier AMPS 36 h
choice is a check on units and illumination conditions. Both values trace
to the Huebner compilations. [Native file](https://zenodo.org/api/records/6994648/files/exospheric_ionization_rates.txt/content)
[Tenishev et al., section 3.3.3](https://doi.org/10.1016/j.icarus.2013.08.021)

For a lunar sensitivity study, use both the stationary birth and a birth
carrying the moving parent's velocity, with identical ionization rate,
fields, parent population, and sampling. Report the change in pickup flux
or energy spectrum. The native rate file and five-species projection can
also be used by other applications. [Data and units](NATIVE_M3_PICKUP_PRECEDENT_SOURCES.md)

## Electron impact

The heliospheric model of Sokół et al. averages a total electron-impact
cross section over core and halo electron populations and calculates losses
of interstellar neutrals. This is a useful independent check on how the chosen
electron population changes survival. Its two-temperature population and
solar-wind dependence differ from the controlled single Maxwellian used by
AMPS. [Model and inputs, section II.3](https://arxiv.org/abs/1809.09823)

For He at 1 au, electron impact accounts for 0.043% of the sum of Poppe's
three median frequencies. Sokół et al. obtain a mean electron-impact share
of 15% in the ecliptic over 1985–2018. A ratio of three medians and a mean
of time-dependent fractions are different statistics. Publish the electron
distribution, solar-wind conditions and averaging convention with any such
comparison. [Poppe archive](https://doi.org/10.5281/zenodo.6994648)
[Sokół et al., Table 1](https://arxiv.org/html/1809.09823v3#S3.T1)

WarpX 26.06 provides an example for a different observable. Its binary
ionization sampler subtracts a threshold, samples directions for two
electrons, imposes a common electron momentum magnitude in its construction,
and determines the recoil from momentum and energy conservation. This gives
a reproducible generic comparison for event-ledger verification. Knowing
that a sampler conserves energy leaves the physical partition and residual
He⁺ electronic state to be checked independently.
[WarpX ionization implementation](https://warpx.readthedocs.io/en/26.06/_static/doxyhtml/_split_and_scatter_func_8_h_source.html#l00504)

## Water photodissociation and dissociative ionization

Three published uses resolve different observables.

| Model | How it creates products | What can be checked |
| --- | --- | --- |
| Shinnaka et al., comet 67P H Lyman-α images | Table 2 assigns H and OH speeds of 17.5 and 1.05 km/s to neutral breakup over 1357–1860 Å. Its row for ionization below 984 Å has a fraction and leaves the product speeds unspecified. | Neutral H coma morphology and the neutral-fragment speed assumption for this wavelength interval. [Published Table 2](https://doi.org/10.3847/1538-3881/153/2/76) |
| Lindkvist et al., hybrid plasma model of 67P | Intact H₂O⁺ production at 1 au uses photoionization frequencies of 3.31 × 10⁻⁷ and 8.28 × 10⁻⁷ s⁻¹ for weak and strong solar activity. The model separately accounts for destruction of the neutral coma. | The bulk ion source and plasma response for the selected activity regime. [Model and Table 1](https://doi.org/10.1051/0004-6361/201732353) |
| EXOSPHID, published exosphere software | In its 0–662 Å band it assigns fractions 0.630, 0.310, 0.043 and 0.017 to H₂O⁺, OH⁺, O⁺ and H⁺ channels, respectively. The dissociative routine passes the photon energy through an ionization call and a heavy-pair breakup call, returning two heavy-product velocities. | Compare branching conventions within that band. A single-photon event with an explicit outgoing electron still needs a joint budget and a state-appropriate partition. [Pinned species table](https://github.com/charligolea/EXOSPHID/blob/62e670f9d0344ced1c1414a86b40fbba43b373ef/src/database/species/H2O.jl) [Pinned dispatcher](https://github.com/charligolea/EXOSPHID/blob/62e670f9d0344ced1c1414a86b40fbba43b373ef/src/photoreactions/photodestruction_logic.jl#L229-L249) |

In AMPS the Hartley 2 neutral-coma reference also records loss-only charged
channels. Their omitted charged products cannot be reconstructed from the
neutral-fragment speeds. The material to reproduce the published model
choices is described in the [water source record](NATIVE_M3_WATER_PRECEDENT_SOURCES.md).

The Leiden photodissociation database also provides wavelength-dependent
partial cross sections for the neutral OH-producing and O-producing channels.
Its fractions marginalize over excitation and combine O + H₂ with O + 2H.
Its radiation-field normalization must be accounted for before using a
frequency in a lunar calculation. [Leiden database and definitions](https://home.strw.leidenuniv.nl/~ewine/photo/H2O_photodissociation_branching.html)

## Returning OH on lunar regolith

Jones et al. assume all returning OH chemisorbs without an activation
barrier. In their model half of the produced radicals escape before return.
The surface rule acts on the returning half. Smolka et al. assume full
thermal accommodation and sample a thermal **flux** upon launch, with speed
weight proportional to `v³ exp[−mv²/(2kTₛ)]`. Laferriere et al.'s
archived lunar hopping code instead samples a thermal speed weight
proportional to `v² exp[−mv²/(2kTₛ)]`. Its smooth-surface residence
function for OH and water uses `1 / [f A exp(−E_d/(kTₛ))]` with
`f = 0.01`, `A = 10¹³ s⁻¹`, and `E_d = 0.6 eV`.
An implementation that copies only the barrier and attempt frequency would
shorten the residence time by a factor of 100 at the same temperature.
[Jones et al.](https://doi.org/10.1029/2018GL080008)
[Smolka et al.](https://www.hou.usra.edu/meetings/lpsc2023/pdf/2089.pdf)
[Laferriere et al., archived implementation](https://doi.org/10.6084/m9.figshare.28548857.v1)

The three choices govern whether an incoming radical
binds chemically, how long a bound species remains, and how a desorbing
particle leaves. In a comparison, report the returning flux, retained
inventory, residence times, and outgoing energy distribution separately.
Measurements of incident-energy-dependent sticking and accommodation on the
specified regolith remain a distinct input. [Surface source record](NATIVE_M3_LUNAR_SURFACE_PRECEDENT_SOURCES.md)

## Proton–sodium charge exchange at 4 keV/u

The lunar pickup-ion compilation gives an integrated Na ion-production
frequency for a solar-wind population. The I24 case instead selects a
particular collision energy, initial Na(3s) and outgoing H(1s) plus
ground-state Na⁺. The earlier archive review found integral theory benchmarks
and state-resolved angular results at other energies. This search did not
establish a differential distribution or recoil error bound for that exact
case. [Lunar rate archive](https://doi.org/10.5281/zenodo.6994648)
[Ground-channel source and units](NATIVE_M3_SODIUM_BENCHMARK_SOURCES.md)
[Angular-data review](NATIVE_CEX_SCATTERING_SOURCES.md)

Swaczyna et al. compared angle-preserving and differential-scattering
calculations for resonant H⁺ + H in the outer heliosheath. Their particular
20 km/s, 7500 K example shifted the secondary neutral bulk speed by about
2 km/s. This shows how to test sensitivity to an angular law. Both the target
species and impact energy differ from I24, so the reported shift cannot set
a lunar Na recoil tolerance. Run the approved isotropic solid-angle case
(uniform in cosine of the angle in the incoming pair center-of-mass frame)
as an experimental comparison. Keep its imposed distribution visible wherever
product energies or angles are reported. [Angular-scattering study](https://arxiv.org/abs/1911.10230)
