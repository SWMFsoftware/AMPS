# Thermodynamic and mapping input audit

The new context assets are measurements at the Ulysses observer. They do not
provide a solar-to-Ulysses background solution. The original 61 PHA flux records,
comparison criteria and `NOT_RUN` verdict are preserved.

| Product | Retained interval | Available quantities |
|---|---|---|
| NASA CDAWeb `UY_COHO1HR_MERGED_MAG_PLASMA` | 1 October–30 November 2003; 1464 hourly records | HGI position, signed RTN magnetic field, bulk speed/angles, proton/alpha densities, separate large/small proton-temperature estimates |
| NASA CDAWeb `UY_PROTON-MOMENTS_SWOOPS` | 28–29 October 2003; 401 observation records | Signed RTN field, actual velocity components, proton density, parallel/perpendicular proton temperatures |

The observer radius is 5.18–5.27 AU. Seven hourly plasma records are archive
fill. The shorter tensor interval has 401 available proton-temperature pairs;
that availability says nothing about electron temperature or magnetic linkage.
The original CSV bytes and HAPI metadata are retained. The quality summary in
`reference/thermal_context_provenance.json` counts fill without modifying it.
No resampling, temporal joins, field rotations, scalar-temperature replacement,
constant mass-loading inference or event forcing is performed. COHO's
large/small alternatives are not the separate product's parallel/perpendicular
components. RTN and HGI metadata are not silently interpreted as HCI.

A complete observationally constrained wind still requires a same-population
electron pressure/temperature definition, composition and uncertainties,
independent stream/transient selection, and dated magnetic mapping/harmonics.
A scalar electron/proton EOS cannot be inferred from proton measurements alone.
An outer observation does not extend the approved 0.30–1.20 AU profile.

## Concrete future early-radius route

[Telloni et al. (2022)](https://doi.org/10.3847/1538-4357/ac8103) documents
Metis observations on 17 January 2021 at 16:30 UTC over 3.5–6.3 solar radii,
linked to two separate PSP streams. Start with sample 1 alone, 18 January
17:49–20:12 UTC near 0.11 AU; sample 2 is a different source and sector.
The coronal Doppler-dimming inference uses scaled PSP proton properties and
assumed coronal electron temperature/composition, so it is not an independent
measurement of every thermodynamic variable.

An obtainable data route is ESA's [Solar Orbiter Archive](https://soar.esac.esa.int/soar/)
for the dated Metis pB/Ly-alpha observations, plus NASA CDAWeb/SPDF for the
matching SWEAP proton moments, FIELDS magnetic data and
`PSP_FLD_L3_SQTN_RFS_V1V2` electron density/core-temperature records. The exact
sample-1 SQTN interval was retrieved in this investigation: 2432 records,
2359 positive density entries and 1449 positive core-temperature entries.
Historical density-quality flags are unavailable before September 2023 in
that product. These records need calibration and population checks: their
positive electron densities span approximately 1308–2347 cm^-3, whereas the
paper's sample-1 proton-density summary is 4540 cm^-3. No reconciliation or
charge-neutral merge is invented here. Core temperature also does not supply
an undocumented total core/halo/strahl pressure.

Recover the actual dated field boundary and mapping controls, rather than
using the general textbook source-surface height: the paper's Figure 3 uses
2.0 solar radii and maps down to 1.2 solar radii. Coefficients, order,
normalization, frame, sector, map uncertainty and observation operators must
be identified. This would be a separate modern stream case. It cannot serve
as unmeasured forcing for the 2003 Ulysses interval and does not cover 1.05
solar radii or establish a connection to 5 AU.

## Other primary constraints and their limits

[Dakeyo et al. (2022)](https://lesia.obspm.fr/perso/pascal-demoulin/22/Dakeyo22_SW_isopoly.pdf)
combines statistical PSP/Helios wind families over approximately 0.1–1 AU.
Its inner isothermal/polytropic continuation is a fitted model, not a measured
coronal temperature profile; the study explicitly identifies the unconstrained
region below 0.1 AU. Different speed families must retain their separate
proton and electron thermodynamics.

[Hahn and Savin (2013)](https://arxiv.org/abs/1302.5403) supplies coronal-hole
ion/wave constraints. Separating ion thermal broadening from wave broadening
uses stated assumptions; those observations are not a co-temporal proton
measurement for any historical case in this campaign. A density/electron
constraint or a heavy-ion line width must not be relabelled as a proton
boundary temperature.

No Stage13 case is newly qualified by this audit. The existing OV01 continuity
verdict remains its own result. Actual same-domain context improves input
readiness while full coronal-to-outer thermodynamics, topology, response and
production-comparison evidence remain open.

The early-route availability check used this exact public query, retained only
as a retrieval recipe (not as Ulysses forcing):

```text
https://cdaweb.gsfc.nasa.gov/hapi/data?id=PSP_FLD_L3_SQTN_RFS_V1V2&parameters=Time,electron_density,electron_core_temperature&time.min=2021-01-18T17:49:00Z&time.max=2021-01-18T20:12:00Z&format=csv
```

Retrieved CSV SHA256: `feaa6af42ddd55af9cd669b31cb122b46fd900c967e3e4a4c8126f7a4f0d919a`.
Its parameter metadata are available through the same HAPI server's
`info?id=PSP_FLD_L3_SQTN_RFS_V1V2` endpoint. These numbers are a documented
readiness check, not a successful cross-instrument or wind-model comparison.

## Direct electron coverage check

The opaque dataset name `UY_M0_BAE` was checked, including the official
[SWOOPS electron guide](https://spdf.gsfc.nasa.gov/pub/data/ulysses/plasma/swoops/electron/swoops_electron_users_guide_update_20030213.txt).
For the October 28–30 request, HAPI returned one record dated May 31, 2003.
It was rejected. The actual 2003 CDF directory lists files only through
May 31. `reference/electron_coverage_audit.json` records the URLs, hashes and
coverage discrepancy. No out-of-period temperature was installed.

The local reproducer checks HAPI record times as well as hashes. A dataset's
nominal start/stop interval cannot turn an event gap into a measurement. The
guide also documents use of interpolated ion charge density in the electron
core density, so even an available electron-density column is not automatically
an independent mass-loading or charge-neutrality test.
