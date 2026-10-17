# Stage-5 silicon damage-response data

This directory contains the only damage-response package admitted by the
current versioned calculation path. It is a testable, checksum-pinned numerical
input, not a general material or device-damage library.

## Validated calculation boundary

| Target | Incident contribution | Response | Domain | Evaluation |
|---|---|---|---|---|
| natural silicon (`Si_natural_G4_11_4`) | proton | DDD | 0.0002–10000 MeV total kinetic energy | log-log between positive SR-NIEL nodes |
| natural silicon | helium-4 | DDD | 0.0001–4000 MeV total kinetic energy | log-log between positive SR-NIEL nodes |
| natural silicon | neutron | DDD and ASTM `n_eq` | 1e-10–20 MeV total kinetic energy | constant within each explicit E722 group |

Electrons, ions heavier than helium-4, and non-silicon absorbers have no
enabled Stage-5 damage response. There is no extrapolation. If any nonzero
absorber-entry contribution required by a requested result lies outside the
table domain or belongs to an unsupported DDD species, the result is written as
`not_assessable`; it is not clamped, dropped, or replaced by zero.

## Files and physical conventions

- `damage_response_library.v1.json` is the strict runtime manifest. It records
  the package identity/version, target identity and atomic weight, interpolation
  rules, domains, source URLs/versions/retrieval times, raw-response hashes,
  canonical-table hashes, and the E722 normalization.
- `sr_niel_v11_1_si_proton.csv` and
  `sr_niel_v11_1_si_alpha.csv` contain the positive response nodes returned by
  SR-NIEL-7 v11.1 for natural Si, displacement threshold `Ed = 21 eV`, Coulomb
  plus hadronic contributions, no hadronic-threshold scaling, and unit fluence.
  Zero-valued leading proton nodes are excluded because log-log interpolation
  is undefined there; the first positive node is the lower supported boundary.
- `astm_e722_19_si_neutron_groups.csv` contains 640 explicit group intervals,
  the SR-NIEL response value, and its hardness ratio to the prescribed silicon
  1-MeV-neutron reference. Group edges are the E192 boundaries reproduced in
  SAND2019-15194 Appendix B; the response values are the ASTM E722-19 option
  from SR-NIEL-7 v11.1.

The neutron reference denominator is calculated, and checked at load time, as

```text
S_ref = 95 MeV mb × 1e-27 cm2/mb × N_A / 28.0855 g/mol
      = 0.002037006185398159 MeV cm2/g.
```

Only neutron damage is divided by this value. Proton and helium-4 DDD are never
converted to E722 `n_eq`.

## Runtime integrity and interpretation

The loader rejects unknown manifest fields, missing provenance fields, wrong
units or energy coordinates, duplicate/missing species, malformed or
nonmonotonic nodes/groups, group gaps, a response/reference mismatch, a
95-MeV-mb conversion mismatch, and every CSV hash mismatch. The manifest hash
written to resolved cases and results is SHA-256 of its canonical JSON; each
table hash is SHA-256 of the exact committed CSV bytes.

The authoritative production observable is the absorber-entry fluence fold.
For a slab of entry area `A`, inward normal `n`, and track direction `u`, each
accepted crossing contributes `S_NIEL(E)/(A |n dot u|)` to DDD; a neutron
crossing contributes its E722 hardness factor with the same projected-cosine
weight to `n_eq`. Geant4 nonionizing energy deposited in the finite slab
divided by slab mass is emitted separately as a diagnostic and is never added
to the table fold. DDD is an environmental exposure proxy, not a prediction of
device degradation.

The manifest records hashes of the retrieved calculator result pages and the
SAND report. Those upstream bytes are not runtime inputs. Production and CI use
only these committed canonical files and never contact a mutable web service.
No redistribution license was identified during retrieval; the manifest records
that fact rather than inventing one, and the project must verify upstream terms
before redistributing the data package.
An upstream update requires a new package version, new filenames/hashes,
rerunning the full Stage-5 suite, and review; do not replace these files in
place.

`scripts/generate_stage5_niel_tables.py` makes the transformation reproducible
from saved upstream result pages and the saved SAND PDF. It strictly extracts
the numeric tables and E192 block, checks the expected 143/133/640 response
counts and every neutron midpoint, writes stable CSV formatting, and can compare
the generated hashes with this manifest. For the pinned source bytes:

```bash
python3 scripts/generate_stage5_niel_tables.py \
  --proton-html /path/to/proton-result.html \
  --alpha-html /path/to/alpha-result.html \
  --neutron-html /path/to/neutron-result.html \
  --sandia-pdf /path/to/SAND2019-15194.pdf \
  --output /path/to/rebuilt-niel \
  --expected-manifest data/niel/damage_response_library.v1.json
```

Formal physical validation is still pending a project-selected accelerator or
device-degradation data set with enough metadata to reconstruct particle,
energy, fluence, silicon target/device, temperature/bias/annealing conditions,
and the measured degradation observable. CRaTER LET data do not validate
semiconductor displacement damage.
