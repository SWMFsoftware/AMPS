# Calculated spectrum exchange format

## 1. Status and purpose

This document is the normative Stage-11 contract for transferring a calculated
energetic-particle spectrum from `srcEarth` to `shieldSim`. The filename retains
the requested `SPECTUM_INPUT_FORMAT.md` spelling.

The implemented schema identifier is:

```text
shieldSim.spectrum_exchange.v1
```

The `shieldSim` reader, case binding, unit conversion, exact spectrum sampler,
schema tests, and end-to-end Geant4 test are implemented. No file under
`/data/vtenishe/GEO/AMPS/srcEarth/` outside `shieldSim/` was modified. A future
change to the `srcEarth` flux/spectrum writers must emit this document directly
or through a separately tested exporter. Until that producer change exists,
Stage 11 has a tested consumer boundary but not a completed native AMPS producer.

This is a data-transfer contract, not a statement that a spectrum is validated
for a particular mission, anomaly, material, geometry, or radiation quantity.
The complete JSON bytes are content-hashed when read and are archived in the
resolved shieldSim case.

## 2. Read-only audit of current `srcEarth` products

The format below was derived from the code that actually creates and tests the
current products. The relevant implementation paths are:

- `boundary/spectrum.h` and `boundary/spectrum.cpp`: boundary-spectrum energy
  basis, intensity units, time selection, uncertainty, and table interpolation;
- `util/FluxNumerics.h`: canonical SI quantity definitions and angular factors;
- `gridless/DensityGridless.cpp`: `gridless_points_spectrum.dat`,
  `gridless_points_flux.dat`, and `gridless_points_density.dat`;
- `3d/DensityMode3D.cpp`: the corresponding `mode3d_*` point and shell products;
- `test/F1`--`F5`, `F11`, `F12`, `F15`, and `F16`: analytic and product-
  consistency checks; and
- `test/SA01`--`SA04`, `test/GOES01`, and
  `test/USWMFCoupledProducts`: time/location binding, observational comparison,
  coupled-product provenance, and strict Tecplot parsing.

### 2.1 Current spectrum storage

The point products are Tecplot ASCII files with one zone per point or epoch.
The common spectrum columns include:

```text
E_MeV
T, T_lower, T_upper
unresolved_fraction
N_sampled, N_resolved
J_boundary_perMeV
J_local_perMeV, J_local_lower_perMeV, J_local_upper_perMeV
J_boundary_lower_perMeV, J_boundary_upper_perMeV
J_omni_perMeV, J_omni_lower_perMeV, J_omni_upper_perMeV
J_planar_perMeV, J_planar_lower_perMeV, J_planar_upper_perMeV
```

Mode3D also writes `N_allowed`. The gridless point spectrum does not currently
write that column. Metadata in `AUXDATA` states the energy basis, mass number,
intensity unit, relative boundary-spectrum uncertainty, temporal selection,
and planar convention. Mode3D adds field-snapshot and coupled-product
fingerprints.

`J_local_perMeV` is the quantity required by shieldSim: local differential
directional intensity in `m^-2 s^-1 sr^-1` per declared energy coordinate.
For an isotropic local population, the other two columns are derived values:

```text
J_omni   = 4 pi J_local
J_planar = pi J_local
```

Neither `J_omni_perMeV` nor `J_planar_perMeV` may be exported as the exchange
`value`. ShieldSim applies the inward-plane angular factor `pi` itself. Passing
an already angular-integrated column would double-count angular normalization.

The current AMPS calculations evaluate spectra at energy nodes. The test
reconstruction in `test/SA01/sa01_common.py` uses a power law between two
positive nodes (linear in log energy and log intensity) and a linear segment
when either endpoint is zero. Physical-zero nodes are retained; they are not
deleted and interpolation never bridges over them.

### 2.2 Current integrated flux storage

`gridless_points_flux.dat` and `mode3d_points_flux.dat` store integrated scalar
products at each point:

- total omnidirectional flux, `F_tot_m2s1`;
- conservative lower and upper values;
- requested energy-channel fluxes;
- one-way planar flux; and
- configured detector-response rates.

Those products are useful cross-checks, but they cannot reconstruct a unique
differential spectrum. They therefore must not be used as the shieldSim source.
A future producer should independently integrate the exchange rows using the
declared interpolation and compare them with the AMPS integrated products
before marking `quality.status` as `accepted`.

### 2.3 Why the Tecplot files are not the interface

The current products are effective diagnostics and plotting artifacts, but not
a sufficiently strict software boundary:

- point index, position, coordinate frame, and epoch are split across zones,
  companion products, run input, and `AUXDATA`;
- gridless and Mode3D have related but non-identical schemas;
- species identity is supplied by run configuration rather than each numeric
  spectrum component;
- `E_MeV` alone does not say total particle energy or energy per nucleon;
- point samples and bin averages are not distinguished by the numeric rows;
- unsupported or stale schema variants cannot be rejected by version; and
- angularly differential and angularly integrated columns are adjacent and
  could be confused by a column-selection script.

The exchange JSON makes each of these semantics machine-enforced.

## 3. One context per file

One exchange file describes one time context and one spatial location. It may
contain multiple explicitly identified particle components.

An AMPS product with multiple zones must be split into one JSON document per
location/epoch. A producer must not select the first zone implicitly. The
location and time written to `context` must be those used for the selected
zone. A manifest may group and hash many exchange documents, but that manifest
is outside schema v1.

## 4. Complete top-level structure

```json
{
  "schema_version": "shieldSim.spectrum_exchange.v1",
  "dataset_id": "unique_dataset_identifier",
  "science_status": "calculated_unvalidated",
  "producer": {
    "software": "AMPS/srcEarth",
    "software_version": "project release or build identity",
    "code_revision": "source revision",
    "run_id": "unique run/snapshot/point identity",
    "created_utc": "2026-10-18T00:00:00Z",
    "calculation_mode": "mode3d",
    "source_product": "mode3d_points_spectrum",
    "input_fingerprint": "0000000000000000000000000000000000000000000000000000000000000000",
    "field_snapshot_id": "snapshot identity required for mode3d",
    "field_snapshot_fingerprint": "1111111111111111111111111111111111111111111111111111111111111111",
    "boundary_spectrum_fingerprint": "2222222222222222222222222222222222222222222222222222222222222222",
    "product_control_fingerprint": "3333333333333333333333333333333333333333333333333333333333333333",
    "description": "optional nonempty description"
  },
  "context": {
    "time": {
      "interpretation": "instantaneous_quasi_static",
      "reference_utc": "2026-10-18T00:00:00Z"
    },
    "location": {
      "location_id": "point_000000",
      "frame": "GSM",
      "unit": "km",
      "position": [6.6, 0.0, 0.0]
    }
  },
  "quality": {
    "status": "accepted",
    "unresolved_semantics": "lower_forbidden_upper_allowed",
    "maximum_unresolved_fraction": 0.0,
    "acceptance_limit": 0.01
  },
  "components": []
}
```

Unknown fields are fatal. JSON duplicate keys, non-finite numbers, malformed
Unicode, and trailing content are also fatal.

### 4.1 Top-level fields

| Field | Requirement |
|---|---|
| `schema_version` | Exactly `shieldSim.spectrum_exchange.v1`. |
| `dataset_id` | Required, nonempty, stable scientific identity. It is not a filename. |
| `science_status` | One of `calculated_unvalidated`, `calculated_validation_qualified`, or `verification_fixture`. |
| `producer` | Required provenance object described below. |
| `context` | Required time and point-location object. |
| `quality` | Required aggregate trajectory-resolution decision. |
| `components` | Nonempty array with unique `component_id` values. |

`calculated_validation_qualified` may be used only after the project has
documented the validation evidence supporting that claim. The reader validates
the spelling, not the scientific truth of the claim. Both
`calculated_unvalidated` and `verification_fixture` are accepted only when the
shieldSim case purpose is `verification_test`; a science-facing anomaly or
preliminary-mission case requires `calculated_validation_qualified`.

### 4.2 Producer provenance

`software`, `software_version`, `code_revision`, `run_id`, `created_utc`,
`calculation_mode`, and `source_product` are required nonempty strings. The
calculation mode is `gridless`, `mode3d`, or `verification_fixture`.
`source_product` names the producer product from which the component was
formed; it is provenance, not a path that shieldSim opens. UTC strings must have the form
`YYYY-MM-DDTHH:MM:SS[.fraction]Z`; local-time offsets are not accepted.

All fingerprint fields are exactly 64 lowercase hexadecimal SHA-256 characters.
`input_fingerprint` and `boundary_spectrum_fingerprint` are required for
`gridless` and `mode3d`; only a hand-calculated `verification_fixture` may omit
them. The input fingerprint covers the canonical run controls needed to
reproduce the spectrum: species, boundary spectrum, field snapshot or field
drivers, domain, location, epoch, energy grid, angular grid, retry/extension
policy, and unresolved tolerance. Mode3D additionally requires
`field_snapshot_id`, `field_snapshot_fingerprint`, and
`product_control_fingerprint`. Field-snapshot fields are forbidden for other
modes. Mode3D producers should also retain current channel, detector, and
observation-state fingerprints in their normal run manifest. Schema v1 does
not invent empty snapshot fields when a gridless calculation has no
corresponding Mode3D object.

### 4.3 Time

Two explicit forms are supported:

```json
{"interpretation":"instantaneous_quasi_static","reference_utc":"2026-10-18T00:00:00Z"}
```

or

```json
{
  "interpretation":"interval_average",
  "reference_utc":"2026-10-18T00:30:00Z",
  "start_utc":"2026-10-18T00:00:00Z",
  "end_utc":"2026-10-18T01:00:00Z"
}
```

An instantaneous file must not contain `start_utc` or `end_utc`; an interval
average must contain both. Current Mode3D snapshot products are
`instantaneous_quasi_static`. A producer must not relabel a snapshot as an
interval average merely because downstream exposure integration will later
assign it a duration.

### 4.4 Location

`location_id`, `frame`, and `unit` are required. `position` is exactly three
finite numbers in x/y/z order. Version 1 requires `unit: "km"`. The frame name
must describe the coordinates actually exported; the producer must not infer
GSM from a filename or silently convert SM/GEO coordinates. ShieldSim archives
the context but does not use spacecraft position to alter local slab transport.

### 4.5 Aggregate quality

The required unresolved semantics are exactly:

```text
lower_forbidden_upper_allowed
```

This matches the AMPS trajectory classification bounds: the lower curve treats
unresolved directions as forbidden and the upper curve treats them as allowed.
These are conservative numerical-classification bounds, not Gaussian standard
uncertainties or confidence intervals.

`maximum_unresolved_fraction` must equal the maximum of all component rows.
Both it and `acceptance_limit` are in `[0,1]`. `status` is `accepted` exactly
when the maximum is less than or equal to the limit. ShieldSim can parse a
`rejected` document for diagnosis, but refuses it as a transport source.

## 5. Component structure

Each component has this structure:

```json
{
  "component_id": "local_proton",
  "particle": {"name": "proton", "Z": 1, "A": 1},
  "energy": {
    "coordinate": "total_kinetic_energy",
    "unit": "MeV",
    "representation": "point_samples",
    "interpolation": "log_log_positive_linear_zero",
    "support": "closed_interval_no_extrapolation"
  },
  "quantity": {
    "type": "differential_directional_intensity",
    "unit": "m^-2 s^-1 sr^-1 MeV^-1",
    "angular_representation": "isotropic"
  },
  "uncertainty": {
    "semantics": "conservative_bounds",
    "source": "AMPS local bounds from boundary and trajectory classification"
  },
  "data": []
}
```

### 5.1 Particle identity

Accepted names and identities are:

| `name` | Required `Z`, `A` |
|---|---|
| `proton` | 1, 1 |
| `alpha` | 2, 4 |
| `neutron` | 0, 1 |
| `gamma` | 0, 0 |
| `ion` | explicit positive integers with `A >= Z` |

The file reader preserves all valid component identities. The current
continuous shieldSim environment sampler admits only explicit proton and alpha
components. Neutrons, photons, and other ions are rejected at case admission;
they are never relabeled. A separate case is required for each selected
component in schema v1.

The producer must map the actual AMPS run species into this identity. It must
not guess a species from mass number alone or use `alpha` for an inclusive
helium isotope mixture.

### 5.2 Energy coordinate and units

The allowed pairs are:

| `coordinate` | `unit` | Required intensity unit |
|---|---|---|
| `total_kinetic_energy` | `MeV` | `m^-2 s^-1 sr^-1 MeV^-1` |
| `kinetic_energy_per_nucleon` | `MeV/nucleon` | `m^-2 s^-1 sr^-1 (MeV/nucleon)^-1` |

Per-nucleon energy requires a massive particle with explicit positive `A` and
is prohibited for gamma components.

ShieldSim converts a selected component to total particle kinetic energy. If
`x = E/A`, then

```text
E_total = A x
dJ/dE_total = (1/A) dJ/dx
```

It also converts SI area normalization once:

```text
1 m^-2 = 10^-4 cm^-2
```

Thus the integrated physical intensity is conserved while matching
shieldSim's `cm2` source-plane normalization. No other unit inference or
conversion is performed.

### 5.3 Quantity and angular representation

Version 1 accepts only:

```text
type                   = differential_directional_intensity
angular_representation = isotropic
```

This maps the AMPS `J_local_perMeV` scalar to an isotropic directional
intensity. The shieldSim case must set `source.mode` to `isotropic`; the source
plane samples directions with `p(mu)=2 mu` and uses the inward-plane integral
`pi J`. A directional PAD, pitch-angle grid, gyrophase distribution, or
anisotropic arrival map requires a later schema version and a corresponding
directional source implementation. Such data must not be reduced to
`angular_representation: "isotropic"` merely to pass the reader.

### 5.4 Uncertainty semantics

Two values are allowed:

- `conservative_bounds`: every data row must contain `lower` and `upper`, with
  `0 <= lower <= value <= upper`; or
- `not_supplied`: rows must not contain `lower` or `upper`.

`source` is always required and explains how the uncertainty statement was
formed. Conservative bounds remain absolute lower/upper intensity curves in
memory. They are not stored as statistical standard deviations. The current
shieldSim transport samples the nominal `value`; it preserves bounds and
provenance but does not propagate the lower and upper AMPS curves as separate
transport ensembles automatically.

## 6. Energy representations

### 6.1 Point samples

Required declarations are:

```json
{
  "representation": "point_samples",
  "interpolation": "log_log_positive_linear_zero",
  "support": "closed_interval_no_extrapolation"
}
```

At least two rows are required. Energies must be finite, positive, unique, and
strictly increasing. A row with conservative bounds is:

```json
{
  "energy": 100.0,
  "value": 2.5,
  "lower": 2.0,
  "upper": 3.1,
  "valid": true,
  "unresolved_fraction": 0.01,
  "n_sampled": 100,
  "n_resolved": 99,
  "n_allowed": 61
}
```

`n_allowed` is optional because the current gridless spectrum file does not
contain it, but it must be present on every row or no row of a component.
`n_sampled` is a positive integer, `0 <= n_resolved <= n_sampled`, and

```text
unresolved_fraction = (n_sampled - n_resolved) / n_sampled
```

must agree within `1e-12`. If supplied, `0 <= n_allowed <= n_resolved`.

For adjacent positive nodes `(E0,J0)` and `(E1,J1)`, shieldSim uses

```text
k = ln(J1/J0) / ln(E1/E0)
J(E) = J0 (E/E0)^k
```

and integrates/samples that power law analytically. If either endpoint is
zero, it uses the straight line through the two nodes. A zero-zero interval has
zero probability. The first and last node are the exact closed support; there
is no extrapolation. Case energy limits may clip the support, but clipping does
not change the declared interpolation of the retained segment.

### 6.2 Bin averages

Required declarations are:

```json
{
  "representation": "bin_averages",
  "interpolation": "piecewise_constant",
  "support": "closed_interval_no_extrapolation"
}
```

At least one row is required:

```json
{
  "energy_lower": 10.0,
  "energy_upper": 20.0,
  "value": 40.0,
  "valid": true,
  "unresolved_fraction": 0.0,
  "n_sampled": 100,
  "n_resolved": 100
}
```

Every bin requires `energy_upper > energy_lower`. Bins must be ordered and must
not overlap; gaps are allowed and remain zero-support gaps. `value` is the
differential intensity averaged over that exact bin and is sampled uniformly
in energy inside the bin. No effective center or edge is inferred.

## 7. Validity and missing data

Every row requires an explicit Boolean `valid`. The file reader preserves
invalid rows, but shieldSim refuses to transport a selected component
containing any invalid row. It does not drop the row, sort rows, replace a fill
value, interpolate across a gap, clamp a negative intensity, or manufacture an
uncertainty.

All energies and intensities must be finite. Intensities may be zero but not
negative. A selected component must have positive integrated weight somewhere
inside the requested case energy limits.

## 8. Required producer mapping from current AMPS products

For each selected AMPS location/epoch, a conforming producer must:

1. Bind the exact point/snapshot and write its actual position, frame, and UTC.
2. Obtain the explicit run species, charge, mass number, energy basis, and
   intensity unit from configuration/state; do not infer them from column names.
3. Copy the exact `E_MeV` nodes in emitted order and verify strict increase.
4. Use `J_local_perMeV` as `value`.
5. Copy `J_local_lower_perMeV` and `J_local_upper_perMeV` as conservative
   bounds, not statistical errors.
6. Copy `unresolved_fraction`, `N_sampled`, `N_resolved`, and Mode3D
   `N_allowed` when available.
7. Compute the document-wide maximum unresolved fraction and compare it with
   the configured AMPS unresolved tolerance.
8. Reintegrate the emitted curve with the interpolation in Section 6 and
   compare total/channel integrals against the corresponding AMPS flux product.
   Any difference above a preregistered floating-point tolerance is fatal.
9. Record stable producer/run identity and the input fingerprint.
10. Write to a temporary file, flush/close it, validate it with an independent
    reader, and atomically rename it to the final name. A partially written
    JSON document must never look accepted.

For a time-dependent gridless trajectory, the zone UTC is authoritative only
after it has been checked against the per-sample field and boundary-spectrum
selection. For a Mode3D snapshot, `SNAPSHOT_EPOCH_UTC` and the spectrum
evaluation epoch must be checked for the intended coupling. A producer must
fail on conflicting time identities rather than choose one silently.

## 9. Selecting the file in a shieldSim case

The file and component are selected under `source.spectrum`:

```json
"source": {
  "mode": "isotropic",
  "spectrum": {
    "format": "spectrum_exchange_v1",
    "path": "../calculated_spectra/run42_point7_proton.json",
    "component_id": "local_proton"
  },
  "energy_limits": {
    "proton": {
      "min": {"value": 10.0, "unit": "MeV"},
      "max": {"value": 10.0, "unit": "GeV"}
    },
    "alpha": {
      "min": {"value": 40.0, "unit": "MeV"},
      "max": {"value": 40.0, "unit": "GeV"}
    }
  }
}
```

`path` is resolved relative to the case file unless it is absolute. Both
proton and alpha energy-limit objects remain required by the current case
schema, even when only one component is selected. Limits are total kinetic
energy per particle after any per-nucleon conversion.

A complete runnable example is
`examples/stage11_spectrum_exchange_case.json`; its two-component input is
`tests/data/stage11_spectrum_exchange.v1.json`.

From the repository root:

```bash
cmake -S . -B test_output/build -DBUILD_TESTING=ON \
  -DGeant4_DIR=/path/to/Geant4/lib/cmake/Geant4
cmake --build test_output/build -j 8 --target shieldSim

test_output/build/shieldSim \
  --validate-case=examples/stage11_spectrum_exchange_case.json

test_output/build/shieldSim \
  --case=examples/stage11_spectrum_exchange_case.json
```

## 10. Failure conditions

The v1 reader or case admission fails on, among other conditions:

- a different or absent schema version;
- unknown/misspelled fields;
- missing time, location, species, units, interpolation, support, or quality;
- invalid UTC spelling or SHA-256 spelling;
- unsupported quantity or angular representation;
- contradictory particle name, `Z`, or `A`;
- a per-nucleon gamma or missing mass number;
- duplicate component IDs;
- fewer than two point nodes, unordered/duplicate nodes, or overlapping bins;
- negative/non-finite values or inconsistent uncertainty bounds;
- inconsistent trajectory counts/unresolved fraction;
- an aggregate quality maximum that differs from row data;
- `quality.status` inconsistent with the acceptance limit;
- a rejected quality status at transport admission;
- an invalid selected row;
- a selected component absent from the file;
- a `verification_fixture` or `calculated_unvalidated` source used for a
  mission/anomaly purpose;
- a non-isotropic v1 source or a case source mode other than `isotropic`;
- a selected species other than explicit proton or alpha; or
- zero selected weight after energy-limit clipping.

## 11. Verification evidence and limitations

The CTest gates are:

```bash
cd test_output/build
ctest -R '^stage11_spectrum_exchange(_transport)?$' --output-on-failure
```

`stage11_spectrum_exchange` checks ingestion, provenance/quality retention,
unknown-version rejection, missing-interpolation rejection, invalid-node
rejection, explicit component selection, and a hand-calculated alpha
per-nucleon plus SI-area Jacobian.

`stage11_spectrum_exchange_transport` sends the point-spectrum fixture through
the normal case parser, exact continuous sampler, Geant4 transport, and proxy
writers. Its source normalization is compared with an independent analytic
integral: one positive power-law segment plus one linear-to-zero segment. It
also checks the `pi` isotropic factor, source-plane rate, source hash, context,
quality, component identity, and rejection boundaries.

Evidence is written under:

```text
test_output/build/validation_evidence/
test_output/build/validation_figures/stage11_spectrum_exchange_analytic.png
test_output/build/validation_figures/stage11_spectrum_exchange_analytic.svg
```

The PNG/SVG pair plots the analytic source integral and isotropic plane-rate
comparisons. It is a review artifact; the numeric tests remain authoritative.

The following are deliberately not claimed by v1:

- the native AMPS writer has not yet been changed to emit the format;
- no frozen native AMPS export has yet been compared byte-for-byte or
  quantity-for-quantity with an independently prepared exchange fixture;
- only isotropic scalar local intensity is transportable;
- only one explicit proton or alpha component is selected per case;
- lower/upper AMPS curves are retained but not automatically propagated as
  separate shieldSim ensembles; and
- passing the interface tests does not validate AMPS environmental physics or
  any downstream proxy outside its separately documented validation domain.

These boundaries must be removed by new versioned implementation and tests,
not by relaxing schema-v1 checks.
