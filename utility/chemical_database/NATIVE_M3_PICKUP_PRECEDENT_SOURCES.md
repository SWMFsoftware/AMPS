# Ionization frequencies used in the Poppe lunar pickup model

Poppe, Halekas and Harada (2022) computed lunar pickup fluxes from exospheric neutral distributions. Their public data supplement provides the three input frequencies per neutral at 1 au: photoionization, solar-wind proton charge exchange and solar-wind electron impact. Each has a minimum, median and maximum. The native file contains 25 rows; the shared CSV exposes all nine values for He, Ne, Na and both argon isotopes. The authors also deposited distinct files for reduced ARTEMIS observations, neutral distributions and volume ion production. The per-neutral frequencies come from the model input file. [Author dataset](https://doi.org/10.5281/zenodo.6994648) · [Paper, Table 1 and §2.2](https://doi.org/10.1029/2022JE007422)

## Interpreting the data

The native file specifies inverse seconds, so a frequency multiplied by the local neutral number density gives an ion production rate per unit volume under the stated conditions. Photoionization primarily follows the Huebner–Mukherjee compilation. For each species, the authors add the median frequencies of all three processes. Figure 2 spans activity levels for photons and a two-standard-deviation range of solar-wind proton flux or electron density for the collision processes. Neither those endpoints nor the extra printed digits of the ASCII source provide laboratory measurement errors. A per-species experimental citation for each individual collision frequency is unavailable from the native file. [Poppe et al., §2.2 and Figure 2](https://doi.org/10.1029/2022JE007422)

| Deposited label | Published species | One-based line of the deposited file |
|---|---|---:|
| He | He | 9 |
| Ne | Ne | 15 |
| Na | Na | 16 |
| Ar36 | ³⁶Ar | 23 |
| Ar | ⁴⁰Ar | 26 |

The source-file row order, species labels and decimal strings survive the projection. Table 1 identifies the two argon species and gives the same median rates to its printed precision for each isotope. Rates containing zeros for other species elsewhere in the full file require individual investigation; the projected subset contains only positive frequencies. [Native frequency file](https://zenodo.org/api/records/6994648/files/exospheric_ionization_rates.txt/content) · [Printed species and units](https://doi.org/10.1029/2022JE007422)

The model traces analytic trajectories for singly charged test ions in constant perpendicular solar-wind flow and magnetic field. It neglects the initial velocity of the parent neutral, lunar gravity, and electromagnetic perturbations from the wake and crustal fields. The paper compares predicted flux with distributed ARTEMIS measurements lacking ion composition; it reports a deficit in the filtered high-energy pickup flux. These assumptions give a precedent for a bulk-flux sensitivity study. The frequencies leave the product ion's electronic state, electron or photon energy sharing, and birth recoil undetermined. [Poppe et al., §§2.2 and 4](https://doi.org/10.1029/2022JE007422)

## How other models handle the missing quantities

Choudhary et al. (2016) used a different lunar observable, near-surface electron density inferred from a Chandrayaan-1 radio occultation. Their fluid model included 16 ion species, among them He⁺, Ne⁺ and Ar⁺. They calculated photoproduction from EUVAC solar fluxes and Schunk–Nagy and Huebner et al. cross sections. The density comparison constrains the combined production, loss and transport assumptions. It supplies no measurement of the birth velocity or electronic state of an individual noble-gas ion. [Choudhary et al., §3](https://doi.org/10.1002/2016GL070612)

Sokół, Bzowski and Tokumaru (2019) model helium and neon survival in the heliosphere using time-, distance- and latitude-dependent ionization. Their electron-impact term integrates a Lotz cross section over cool-core and hot-halo electron distributions. Over 1985–2018 at 1 au in the ecliptic, their mean modeled electron-impact shares are 15% for He and 9.6% for Ne. By comparison, electron impact contributes about 0.043% to the sum of the three Poppe helium medians. These numbers illustrate the sensitivity to the chosen electron environment and method. They represent different temporal statistics and cannot by themselves select a local lunar rate. [Sokół et al., §§II–III and Table 1](https://arxiv.org/html/1809.09823v3) · [Poppe data](https://doi.org/10.5281/zenodo.6994648)

WarpX 26.06 provides an explicit electron-impact event approximation in the binary DSMC operator. It subtracts the ionization energy from the incoming center-of-mass kinetic energy, draws two independent isotropic electron directions, assigns equal outgoing electron momentum magnitudes and calculates ion recoil to conserve momentum and kinetic energy. It demonstrates how an assumed angular and energy-sharing law can close the kinematic ledger. He-specific excited-state fractions and an experimental joint product distribution require separate evidence. The AMPS packet contains the Poppe source only. [WarpX 26.06, source lines 459–542](https://warpx.readthedocs.io/en/26.06/_static/doxyhtml/_split_and_scatter_func_8_h_source.html#l00459)

## Reproduction and access

Run from the AMPS root using Python 3 and curl. Each destination must be new. Only the first command that runs the recipe needs network access.

```sh
work=$(mktemp -d)
recipe=utility/chemical_database/native_m3_pickup_precedent.py
python3 -B "$recipe" fetch "$work/acquisition"
python3 -B "$recipe" prepare "$work/acquisition" "$work/packet"
python3 -B "$recipe" prepare "$work/acquisition" "$work/replay"
cmp "$work/packet/source_packet.json" "$work/replay/source_packet.json"
python3 -B "$work/packet/native_m3_pickup_precedent.py" project \
  "$work/packet/exospheric_ionization_rates.txt" "$work/reproduced.csv"
cmp "$work/packet/poppe2022-atomic-1au-rates.csv" "$work/reproduced.csv"
```

The admitted native file has SHA-256 `a208a94f318bbe2c63c59e6e28cfd15bf78901b7f7560480c8298bb235b7826e`; the Zenodo member has published MD5 `149afa8c4f7ee6d8d276f097ffc59ec1`. The projection SHA-256 is `76fb5cdb5181a5a4a6ceb0e9fff40da28fb3069affe107e440644863ed083b61`. The acquisition records retrieval times, effective HTTPS URLs, and the exact hash of the Zenodo rights record used. A fresh fetch receives new retrieval times; replay against retained originals reproduces the same manifest.

The deposit grants CC BY 4.0 to its data and identifies Poppe, Halekas and Harada as creators. Redistribution of the native numerical data or the derived CSV requires credit, a link to the licence and identification of modifications. The original paper carries a separate AGU copyright notice. Its packet entry exposes bibliographic and method metadata while the original PDF stays with the acquisition. [Versioned deposit](https://doi.org/10.5281/zenodo.6994648) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) · [Article](https://doi.org/10.1029/2022JE007422)

For a release selected with its exact database lock, `db.describe("poppe2022-lunar-pickup-rates")` lists the source records and limitations, `db.source_path("poppe2022:ionization-rates")` returns the complete licensed file, and `db.table("poppe2022:atomic-1au-rates").rows` returns the five projected rows. To exercise the selected release through both shared database readers, supply its lock:

```sh
AMPS_CHECK_M3_PICKUP_PRECEDENT=1 python3 -B test/chemical_database/run_all_tests.py \
  --database-lock "$release_lock" --output-dir "$work/test-results"
```

The three checks compare 15 median values against the independently printed table, replay all 45 values from the native lines using the production Python and C++ table readers, and enforce rights and scientific scope. These are database access and transcription checks. No AMPS ionization law, daughter sampler or observational validation is enabled by the source packet.
