# D02 illumination research and operator requirements

Research date: 2026-11-12 (America/Chicago). **D02 remains PARTIAL; I15/I20
are NOT VALIDATED.** The user authorized research and verification, not an
assumed integration epoch or replacement of a frozen reference product.

## Primary evidence inspected

| Evidence | What it establishes | Limit |
|---|---|---|
| Native `AVGVISIB_65{N,S}_240M_201608.LBL` and matching IMG | Product identity, signed 16-bit little-endian samples, scale/offset, projection/frame and observable definition. | The labels do not identify their date fields as illumination-simulation integration bounds. |
| Native `LPSR_65{N,S}_240M_201608.LBL` and IMG | Independent external PSR classifications and native projection/scaling. | A spatial classification does not contain a timestamped illumination history. |
| [NASA PGDA product 69](https://pgda.gsfc.nasa.gov/products/69) | General description of hourly sampling over approximately one 18.6-year cycle; finite solar-disc geometry; links to the frozen products. | No exact simulation start/end epochs or endpoint/sample-count convention for these files. |
| [LOLA archive README](https://imbrium.mit.edu/EXTRAS/ILLUMINATION/README_ILLUMINATION_2016.TXT) | Product families, supported map resolutions/extent, IMG as primary format and a warning about terrain-related artifacts. | No simulation interval/cadence specification. |
| [NASA NTRS record 20120010094](https://ntrs.nasa.gov/citations/20120010094), DOI 10.1016/j.icarus.2010.10.030 | The 2011 paper describes long-period calculations with six-hour sampling. | This older publication does not establish the exact 201608 generation settings; NTRS offers no full-paper download at that record. |
| [Author's LEAG 2011 abstract 2030](https://www.lpi.usra.edu/meetings/leag2011/pdf/2030.pdf) | Horizon-method example with 720 azimuths; a figure describes one 18.6-year cycle at six-hour sampling. | This is an older calculation and cannot supply the newer product's missing controls. |

The cadence statements need not describe the same calculation/version. Treat
them as separate claims, not an inconsistency to resolve by choosing the value
that makes a model comparison convenient. Product-specific metadata remain
required. No contact message was sent to the instrument team.

Exact downloaded documentation, access dates, resolved URLs and SHA-256 values
are preserved under
`test_output/srcMoon/m2-next-20261113T042634Z/research/source_provenance.json`.
They were downloaded into campaign artifacts, never into canonical `raw/`.

| Documentation snapshot | SHA-256 |
|---|---|
| PGDA HTML | `aca046555f29da9fcbe926c0b90450ce69e78a44cc406bc2a60d9a010eb3fe42` |
| Archive README | `f806f9c7fb410b17dd75c2f44f91755b5ca890ac047004a1e8fbb819c24b9ec6` |
| NTRS metadata JSON | `aac53b01e26e3c16be87ea55d2ecc1855ae168b69fddafaddc8c53f0d973b2fc` |
| LEAG abstract PDF | `93d66fd0c0a7f9121a9504e73cdefdca18ec11e5d92e02eb15724bf8191e8def` |

## Correct native observable

The frozen AVGVISIB labels define an occurrence frequency: the fraction of
sampled times for which any nonzero portion of the solar disc is visible.
For a declared set of N epochs this would be

`V = sum(visible_disc_fraction(t_j) > 0) / N`.

This definition comes from the native label; the epochs t_j and N do not.
Retain these separate model diagnostics:

- any-disc visibility at an epoch (Boolean);
- visible solar-disc area fraction at an epoch;
- projected incident radiative flux in W/m2;
- Earth-occultation fraction, separate from terrain blocking.

The thermal boundary must use the appropriate emitted/absorbed energy flux;
an occurrence Boolean cannot supply the partially obscured radiative flux.
Likewise, temporal mean visible area is a different statistic from AVGVISIB.
Do not compare one against the other and call it validation.

Native DN converts with the label's 0.00004 scale and zero offset for the
frozen AVGVISIB products. One DN is quantization, not a documented complete
measurement/model uncertainty. Do not invent a fixed error bar from it.
Geometry is planetocentric/east-positive, polar stereographic, in the native
DE421 mean-Earth/polar-axis frame; native map scale is 240 m/pixel and the
projection sphere is 1737.4 km. The production reference sphere/datum remains
distinct and must be converted explicitly, as in the existing geometry tests.

## What the current files cannot supply

The AVGVISIB north label contains START_TIME 2009-07-13T17:33:17.246 and
STOP_TIME 2013-07-18T14:12:13. That span is not an 18.6-year integration and
the label does not state that it is the averaging interval. Retain these as
native metadata with an unresolved role; do not use them as simulation bounds.
The file's 201608 suffix and creation date also do not specify an integration.

A product-specific contract must resolve exact start/end epochs, time system,
cadence, sample count/endpoint convention, solar-disc/ephemeris correction
conventions, observer height and Earth-occultation policy. The current sources
do not establish this full set. The read-only audit reports explicit nulls and
SKIPPED instead of supplying defaults.

Even a fully documented average map cannot reconstruct individual illumination
epochs or sunrise/sunset times. I15's unchanged epoch-classification/timing
requirements additionally need an independent time-resolved reference (for
example, producer-supplied states or identified/calibrated timestamped imagery
with documented geometry). Such a reference/selection must be acquired and
frozen before examining AMPS agreement. No averaged DN threshold is substituted
for that information.

## Changes and authorized verification

The older spatial processor placed general 18.6-year/one-hour/0.5-degree claims
into its QA as though they were exact-product controls. The correction leaves
every CSV science row, scaling, selected pixel and edge flag unchanged; its
QA v2 reports `spatial_processing_status=READY`, overall `status=PARTIAL`, and
an unresolved temporal operator. Old/new CSV identity is independently checked
in the campaign artifacts. Canonical external metadata/raw products are not
rewritten or promoted.

`audit_d02_temporal_contract.py` verifies native contracts and frozen hashes.
`audit_m2_packages.py` invokes its ten protocol self-checks automatically and
includes its detailed report in the existing global M2 supporting campaign.
No new U/I ID, physics law or acceptance threshold is introduced. Formal
I14-I20 retain same-campaign I13 and their separate gates. Neither D02 spatial
decoding nor the metadata audit produces an observational comparison figure.

Remaining factual input: authoritative generation settings/time-resolved
benchmark for the exact frozen D02 products. Research has narrowed the gap;
it has not resolved it. The separate I19/I20 numerical proposal is in
M2_I19_I20_PROTOCOL_PROPOSAL.md/JSON and remains unapproved.
