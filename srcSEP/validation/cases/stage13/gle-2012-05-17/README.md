# 2012 May 17: low-corona constraints and GLE71/near-Earth observations

`case.json` describes a proposed observational validation case; it does not
supply unmeasured shock state, acceleration history or a new tolerance. The
low-corona inference, GOES integral flux and neutron-monitor count rates are
separate observables. Available references alone do not establish a successful
particle-transport or ground-level-enhancement comparison.

The reference package contains:

- Two unmodified event rows selected from the existing EV01 CCMC table. They
  preserve the >10 MeV/10 pfu and >100 MeV/1 pfu definitions and rounded official
  event quantities. EV01 and its complete source table remain unchanged.
- The complete 8096-byte [official GOES-13 processed JSON](https://ccmc.gsfc.nasa.gov/assessment/topics/SEP/campaign2020/docs/json/GOES-13_integral.2012-05-17T000000Z.json),
  with five integral-channel event summaries and original processing metadata.
  Its `sep_profile` strings name time-series files; those series are not embedded
  and have not been acquired. Onset peak and maximum peak remain distinct.
- `gopalswamy2013_low_corona_constraints.csv` and
  `low_corona_inference_basis.json`, transcribed from the [primary paper](https://cdaw.gsfc.nasa.gov/publications/gopal/gopal2013ApJL.pdf),
  DOI `10.1088/2041-8205/765/2/L30`, retaining method, clock and uncertainty limits.

The reported 1.38 solar-radii value is a heliocentric CME/shock proxy at the
corrected metric-type-II onset. It uses an earlier EUV disturbance image and
short kinematic extrapolation with the paper's spherical construction. The
2.32 solar-radii value is a COR1 CME height associated with the inferred
particle-release time, not a directly observed particle-release position.
The stored 01:40 release clock is normalized to electromagnetic arrival timing;
its inference assumes a scatter-free first-arrival path and particle energy.
The 01:40:26 COR1 exposure epoch is stored separately. Neither heliocentric
radius is a height above the photosphere.

The [Rouillard et al. alternative](https://arxiv.org/pdf/1605.05208v2), DOI
`10.3847/1538-4357/833/1/45`, uses multi-species velocity dispersion and a longer
inferred path. Its physical solar-release clock and electromagnetic comparison
clock remain distinct. The paper's release uncertainty is written differently
in Section 2, Section 8 and Appendix A; all statements are preserved in the
inference JSON. No error convention is silently selected or averaged with the
Gopalswamy result. The two reported heights have no numeric sigma supplied by
this package, so their comparison cannot receive a fabricated tolerance.

The [official CCMC event page](https://ccmc.gsfc.nasa.gov/assessment/topics/SEP/campaign2020/20120517.php)
was independently checked against the existing rows. Ground-level counts must
come from the [Oulu GLE database](https://gle.oulu.fi/) or [NMDB](https://www.nmdb.eu/nest/),
with station corrections, original-file provenance, yield functions and viewing
geometry. A GOES proton flux is not a substitute for a neutron-monitor response.
No GLE station time series or response conversion is claimed here.

Reproduce the pinned small official JSON into fresh evidence:

```sh
python3 -B reference/fetch_reference.py --output-dir /absolute/fresh/evidence
```

A restricted runner can first download the exact JSON URL using top-level
`curl -fsSL`, then use `--source-cache /cache`. The receipt explicitly records
cache replay, verifies bytes/SHA256 and leaves source references untouched.
The selected CSV is reproduced by retaining exactly `event_id=20120517` from
`../../EV01/reference/ccmc_2021_goes_observations.csv`; the parent checksum is
recorded in `reference/provenance.json`. No interpolation, background removal,
normalization or synthetic sample is introduced.

The scientific gate remains incomplete until the event-specific background,
front/field-line history, source and boundary contract, instrument response and
uncertainty/scoring prescriptions are reviewed. See `missing_required_inputs`;
none is replaced by controlled-model defaults.
