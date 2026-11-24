# Molecular water in Schörghofer's lunar models

The Planetary Code Collection contains two separate models for H₂O. One
predicts the flux leaving a grain as a function of water coverage and temperature.
The other propagates ballistic hops and removes water exposed to solar photons.
The two calculations have different inputs and state variables. Combining them
would require explicit coupling of the adsorbed and exospheric inventories.
The source packet identifies tag v1.2.3, commit
`ab8baaecb0c8f4414f651b7d31cad610339be0cb`. [Archived release](https://doi.org/10.5281/zenodo.15161597)

## Desorption of retained H₂O

In the efficient-diffusion model, adsorbed water preferentially occupies the
higher binding energy sites. The grain desorption-energy distribution begins at
0.65 eV and has an exponential width of 0.22 eV. The pure ice limit uses
0.529 eV. The monolayer population is 10¹⁹ molecules m⁻². At positive surface
temperature T in kelvin, the attempt frequency in inverse seconds is
`2.2 × 10¹⁷ / √(T/K)`. These are model inputs for activated grains. For the
cited 2025 subsurface calculations the authors also examine the opposite
limit, in which each molecule stays on its original site.
[Tagged model](https://github.com/nschorgh/Planetary-Code-Collection/blob/ab8baaecb0c8f4414f651b7d31cad610339be0cb/Lunar/desorptionrates.f90#L4-L70)
[Adsorption model](https://doi.org/10.3847/PSJ/acf19b)
[Diffusion endmembers and Eq. 4](https://doi.org/10.3847/PSJ/ade5b2)

The positive input coverage θ is molecules per square meter of grain surface.
Let v denote θ divided by the monolayer population and b denote 1/(k_BT) in
inverse electronvolts. For v ≤ 1 the emitted flux in molecules m⁻² s⁻¹ is
`ν θ exp(−Ec b) v^(W b) / (1 + W b)`. This gives longer individual residence
times at lower coverage. For v > 1 the code returns a second expression with
γ = 1/v², tending to sublimation from pure ice at high coverage. The archived
routine evaluates a candidate labeled S2 but assigns S1 to its output. Zero
coverage gives zero flux. The packet reports these active branches and labels
the applicable species. A temperature interval or an uncertainty on OH sticking
cannot be inferred from the fitted energies. [Source routine](https://github.com/nschorgh/Planetary-Code-Collection/blob/ab8baaecb0c8f4414f651b7d31cad610339be0cb/Lunar/desorptionrates.f90#L36-L72)

## Ballistic H₂O and incident OH

The associated example for the lunar exosphere selects molecular water, draws
two tangential Gaussian velocities and an outward half-normal velocity at the
surface temperature, and uses an illuminated destruction rate of 1.26 × 10⁻⁵
s⁻¹ at 1 au (mean lifetime about 22 h). Inspection of the tagged implementation
matters. Its Monte Carlo removal probability is the integrated hazard h up to
0.2, and 1 − exp(−h) above 0.2. It records a removed parent and supplies no
electron, ion or neutral daughter trajectory. The exact exponential survival
formula is the physical probability for a constant Poisson hazard. [Ballistic source](https://github.com/nschorgh/Planetary-Code-Collection/blob/ab8baaecb0c8f4414f651b7d31cad610339be0cb/Exospheres/montecarlo.f90#L26-L139)
[Lunar migration calculation](https://doi.org/10.1002/2014GL060820)

Farrell et al. modeled polar ice sputtering and impact vaporization, followed
by H₂O migration and long-lived retention of the OH photoproduct. Their
argument uses a desorption temperature above 400 K. The underlying Hibbitts
et al. experiments dosed initially dehydrated analog materials with molecular
water and interpreted the strongly bound hydroxyl fraction using recombinative
desorption evidence. In related measurements Poston et al. observed possible
water release from recombination of native hydroxyls above 425 ± 25 K. This
supports a chemically retained surface reservoir. Determining the probability
of gas-phase OH sticking and its scattered energy at a specified impact still
requires a measurement or bounded model for that distinct process.
[Farrell et al.](https://doi.org/10.1016/j.pss.2013.05.009)
[Hibbitts et al.](https://www.sciencedirect.com/science/article/abs/pii/S001910351100073X)
[Poston et al.](https://doi.org/10.1029/2012JE004283)

## Replay and database checks

From the AMPS root, use a fresh directory. The first invocation downloads two
immutable tagged files, checks their SHA-256, and records an HTTP receipt.
The second invocation reads those same originals offline. The shared catalogue
exposes the model parameter metadata for reuse and retains the source URLs.
This release uses the existing metadata-only policy for the Fortran originals.

```bash
work=$(mktemp -d)
source=utility/chemical_database/native_m3_lunar_water_coverage.py
python3 -B "$source" --fetch --acquisition "$work/raw" --output-dir "$work/packet"
python3 -B "$source" --acquisition "$work/raw" --output-dir "$work/replay"
cmp "$work/packet/source_packet.json" "$work/replay/source_packet.json"
```

Given a lock selecting this optional source, query
`db.describe("schorghofer-lunar-water-models-v123")["files"]`, then the
`quantity_metadata` on each returned file. The checks are enabled with
`AMPS_CHECK_LUNAR_WATER_COVERAGE_RELEASE=1` and the selected lock on the
ordinary shared database test runner. They compare the active code branch,
the physical scope, the coverage and energy parameters, the photo-loss
condition and the metadata access policy. A failure to fetch or reverify an
original prevents creation of a publishable input packet. No lunar physics
switch or observational status follows from the catalogue entry.

| Tagged source | SHA-256 |
| --- | --- |
| Lunar/desorptionrates.f90 | `0c15a23cc629a096913085e7df12f5039bb9e90003696a0efbb0f54b29e33e9f` |
| Exospheres/montecarlo.f90 | `ffe2012fae505c7cd57241215f1155fdba22c53f3c310555ba53d38a73eddef0` |
