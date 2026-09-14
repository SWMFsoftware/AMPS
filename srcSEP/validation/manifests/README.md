# External evidence manifests

These files are completion templates, not reference data and not passing
evidence. Keep `status=INCOMPLETE` until a campaign has actually run and every
required field is populated.

For each entry in `inputs`, `path` is resolved relative to the manifest. Record
the exact archived bytes used by the analysis, then calculate lowercase SHA-256
with a platform tool such as `sha256sum`. The campaign runner recomputes the
digest and rejects missing or changed bytes. `source_url_or_pid`, `access_utc`,
and `license_or_acknowledgment` make the archive traceable beyond a local path.

The SWMF record must describe a native coupled execution and its field-line
selection/cadence, and must include `background_replay` and
`transport_response` metric families. The observational record must identify event time bounds,
mission, instrument, processing level/version, variables, cadence, coordinate
system, quality flags, held-out policy, uncertainty method, and the model-to-
instrument forward operator. Across all event records, `metric_family` must
cover:

- `onset`;
- `anisotropy`;
- `spectra`;
- `fluence`;
- `decay`;
- `multi_spacecraft_longitude`.

Every metric also carries its value, threshold, comparison operator, and units.
The runner will never infer a pass from a plot or free-form narrative.

WP39 replaces the former free-form forward-operator label with the structured
`configuration.forward_operator` object in
`observational_event.template.json`. It must identify the response matrix and
angular model, cadence and species, detector dead time and saturation policy,
background subtraction, uncertainty propagation, and operator version. Archive
response files through `inputs` and record their SHA-256 like observations.
