# Neutral noble-gas loss in a lunar exosphere model

Alexander Peschel, Cesare Grava, Jakob Breternitz and Philipp Reiss published
*Probabilistic Modeling of Lunar Topography for Exosphere Simulations and
Implications for Helium, Neon, and Argon*, Planetary Science Journal 7, 69
(2026). The packet preserves their original PDF with attribution under the
CC BY 4.0 notice on page 1 and a projection of the photoionization row in
Table 1 on page 3. The table uses Huebner and Mukherjee (2015) as its source.
It supplies a reproducible published model input for testing a neutral-loss
law under the stated conditions. [Primary article](https://doi.org/10.3847/PSJ/ae4901)
[Licence](https://creativecommons.org/licenses/by/4.0/)

| Species | Paper’s quiet-Sun rate (s⁻¹) | Local D04 quiet (s⁻¹) | Local D04 active (s⁻¹) |
| --- | ---: | ---: | ---: |
| He | 5.25 × 10⁻⁸ | 5.64 × 10⁻⁸ | 1.68 × 10⁻⁷ |
| Ne | 5.9 × 10⁻⁷ | 1.84 × 10⁻⁷ | 5.87 × 10⁻⁷ |
| ⁴⁰Ar (paper), elemental Ar (D04) | 3.06 × 10⁻⁷ | 3.42 × 10⁻⁷ | 7.74 × 10⁻⁷ |

D04 values come from rows 2–3, 6–7 and 8–9 respectively of
`photochemistry/processed/phidrates/photo_rates_1au.csv` under the external
validation-data root. These rows have optical depth zero, distance 1 au and
solar activity 0 or 1. The data file SHA-256 is
`3aa54915a08258371dbe975949e6afb2ab74618164fadca8bceebe876c926847`.
The published Ne number is 3.21 times the D04 quiet value and lies within
0.6% of the D04 active value. The cause has not been established. Retain
the authors’ quiet-Sun label and cite the selected source in each run.
The three rates belong to different species. They carry no published
uncertainties or residual-ion state probabilities. [PHIDRATES](https://phidrates.space.swri.edu/)

## Conditions and interpretation

The calculation assumes a continuous source and a stationary distribution
in the rotating lunar frame with zero obliquity. He and Ne enter on the
dayside with cosine solar-incidence weighting. All returning He and Ne
relaunch with a local-temperature thermal flux distribution. Radiogenic
⁴⁰Ar enters uniformly over the globe and all incident Ar is adsorbed.
The authors use a 0.31 eV Ar adsorption energy in equation (5).
Table 1 gives source fluxes of 4.80 × 10¹⁰, 1.65 × 10⁸ and 1.12 × 10⁷
m⁻² s⁻¹ respectively. The ⁴⁰Ar value conflicts with 1.12 × 10⁻⁷ m⁻² s⁻¹
in the left column of the same printed page. Both exponents were checked
on the rendered page. Source fluxes are excluded from the projection.

The neutral survives a directly illuminated time step with probability
`exp(−k Δt)`. Entry into shadow suppresses the photon loss. A sampled
ionization removes the particle, and gravitational escape terminates
trajectories separately. The paper discusses extending the calculation to OH
and H₂O. Its three simulated species provide no OH sticking coefficient,
charged water products or gas-phase recoil law. [Model §§2.1, 4 and Table 1](https://iopscience.iop.org/article/10.3847/PSJ/ae4901/pdf)

The linked Zenodo archive is Norbert Schörghofer’s *Planetary-Code-Collection*
v1.2.3, cited by the authors in their thermal-model discussion. The TUM team
also publishes its ExESS project. The exact run configuration or independent
supplement reconciling the Ne and Ar discrepancies was not established.
[Code archive](https://doi.org/10.5281/zenodo.15161597)
[TUM ExESS project](https://github.com/Smolkaa/ExESS.jl)

## Replay and query

Run from the AMPS checkout with Python 3, Poppler and curl. Only the fetch
uses the network. The archived PDF SHA-256 is
`2aa3fc35bf93027e5879c904250c840e0a3db2ee37fd4330f5bdd5feffd978ae`.
The recipe reads the PDF text layer, checks the table header and units,
and preserves the original mantissas and signed exponents. The output CSV
SHA-256 is `ba6baa024ab3437f60c9b399cb19023d2ecd2f6a4e37ceb79c819d56bb5b5748`.
Changes to the PDF edition require review before reuse.

```bash
w=$(mktemp -d)
p=utility/chemical_database/native_m3_lunar_noble_gas.py
python3 -B "$p" fetch "$w/acquired"
python3 -B "$p" prepare "$w/acquired" "$w/packet"
python3 -B "$p" project "$w/acquired/peschel-et-al-2026.pdf" "$w/replayed.csv"
cmp "$w/packet/peschel2026-table1-photoionization.csv" "$w/replayed.csv"
python3 -B "$p" compare-d04 /data/vtenishe/moon_validation_data/photochemistry/processed/phidrates/photo_rates_1au.csv
```

With the chosen release lock, `db.source_path("peschel2026:article")` gives the
original and `db.table("peschel2026:quiet-sun-neutral-loss").rows` the three
rows. `db.describe("peschel2026-lunar-noble-gas-loss")` exposes the solar
label, physical scope, Ar source warning and attribution. Select the release
checks with `AMPS_CHECK_M3_LUNAR_NOBLE_GAS=1` and an explicit lock when
running `test/chemical_database/run_all_tests.py`. The tests use the primary
printed values, inspect the PDF identity, and check both language interfaces.
