# Burgess 2016 HxOy table views

`native_nist_hxoy.py` exposes literal cells from Tables 3 and 13 of Donald R.
Burgess, Jr., *An Evaluation of Gas Phase Enthalpies of Formation for
Hydrogen-Oxygen (HxOy) Species*, Journal of Research of NIST **121**, 108–138
(2016), [DOI 10.6028/jres.121.005](https://doi.org/10.6028/jres.121.005).
The original sources are the [official PDF](https://nvlpubs.nist.gov/nistpubs/jres/121/jres.121.005.pdf)
and [PMC7339710 native article HTML](https://pmc.ncbi.nlm.nih.gov/articles/PMC7339710/).

The adapter requires the exact SHA-256 supplied by the caller. It retains each
cell once, its row/cell indices, literal row/column spans, text, and original
inner HTML. Superscripts and footnote links remain available in `native_html`.
It never fills spanning cells, selects an evaluated value, converts units,
assigns a confidence level, or maps a generic species to an isotope/state.
Table 3 concerns OH formation enthalpy at 298.15 K and includes historical and
computed entries; Table 13 contains several different ion-energy quantities.
Neither table is automatically a zero-K reaction-state profile.

The narrow `nist-jres-121005-public-use` policy applies to this NIST-employee
article and its lossless views, separately from NIST SRD/WebBook/ASD products.
Retain the exact [NIST notice](https://www.nist.gov/open/license), its disclaimer,
the citation above, and “Republished courtesy of the National Institute of
Standards and Technology.” CSV extraction is a change in representation only.
An article's older source constants and quoted uncertainties remain native.

Generate new views (an existing output directory is refused):

```sh
python3 -B utility/chemical_database/native_nist_hxoy.py \
  --source /path/to/PMC7339710.html --expected-sha256 SHA256 \
  --output-dir /new/path/to/views
```

Run local administrative parser controls and actual-source checks:

```sh
AMPS_NIST_HXOY_SOURCE=/path/to/PMC7339710.html \
AMPS_NIST_HXOY_SHA256=SHA256 \
python3 -B test/chemical_database/test_native_nist_hxoy.py
```

Without both explicit source variables, the authoritative-source class is
SKIPPED. The class declares `requires_authoritative_source = True` for the
shared protocol. A source-integrity PASS does not qualify molecular physics.
