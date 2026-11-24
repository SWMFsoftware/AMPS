# I23 — Electron-impact ionization rate coefficients

The complete specification is [reference/acceptance.json](reference/acceptance.json).
Original acceptance: Coefficient relative error <=2%; ensemble lifetime error <=3%; explicit zero/invalid density/temperature handling. Existing U16 <=1e-12 criterion remains.

The [approved2026-11-18 addendum](../../M3_QUALIFICATION_ADDENDUM_20261118.md)
allocates relative0.003 numerical bias and fixes endpoint-survival lifetime
as primary. Existing direct/local coefficient gates are unchanged. The owner
records approval; missing physical electron-product/bath/source profiles and
actual native/predecessor evidence still block formal qualification.

This individual owner is currently a **formal capability/predecessor gate**.
It preserves source/input/executable and same-campaign predecessor evidence,
but full particle execution/scoring is unavailable while the recorded physical
kinetic policies/production event dispatch are missing. It reports SKIPPED,
never a scalar-rate PASS relabeled as chemistry qualification. The capability
register and [completion proposal](../../M3_COMPLETION_PROPOSAL.md) identify
the required implementation/data; new source recommendations are separate.

From the repository root (an existing dedicated campaign directory):

```sh
python3 -B srcMoon/test/integration-tests/I23/test.py \
  --amps "$PWD/amps" --input srcMoon/examples/m3_reaction_rates.in \
  --mpi-ranks 8 --campaign-mode regression \
  --output-dir test_output/srcMoon/M3-formal-NEW
```

The runner creates `<output>/I23/` with its real result/manifest,
stdout/stderr and an exact input snapshot when available. No `receipt.txt`
or comparison figure exists for a run that did not execute. Exit codes:
`0=PASS`, `1=FAIL`, `2=ERROR`, `77=SKIPPED`. Do not overwrite prior evidence.

The owner now retains `protocol_preparation.json`, binding the concrete
[physical model](../../M3_PHYSICAL_MODEL_SPECIFICATION.md) and
[qualification protocol](../../M3_QUALIFICATION_PROTOCOL.md). Matching
identities and a structural preparation PASS do not launch particle chemistry
or change this formal gate. Missing inputs remain SKIPPED; corrupt identities
are ERROR. Full source-family, product and predecessor requirements remain.

## Stationary helium pilot

Use a successful [atomic-family preparation](../../../../utility/chemical_database/state_profiles/ATOMIC_FAMILY_CONTRACT.md)
and its output `profile.json`.

```bash
/usr/bin/python3 -B srcMoon/test/integration-tests/I23/run_helium_pilot.py \
  --profile /path/to/profile.json --profile-sha256 PROFILE_SHA256 \
  --output-dir test_output/srcMoon/i23-helium-NEW
```

The runner compiles the production helium provider, authenticates D05 and
the physical state records, and compares the coefficient, density response,
fit slope, incident mean and thresholds at 10 eV. It writes
`comparison.csv`, the profile snapshot, logs, `run_manifest.json` and
`result.json`. The guarded composer is included in the verification.
Linked I23 still needs authenticated model-error evidence before a
helium event can commit. [Pilot contract](HELIUM_PILOT_CONTRACT.md)

Add `--kim-1992-pdf /path/to/jresv97n6p689_A1b.pdf` to run the optional
independent fit comparison. The runner verifies the exact PDF digest and
extracts its Table 1 neutral-helium row using `pdftotext`. SciPy integrates
the published fit over its finite range and records the result separately
in `kim_model_comparison.json`. It compares the predicted rate and incident
mean, with numerical quadrature errors distinct from physical data and
fit uncertainty. [Kim, 1992](https://doi.org/10.6028/jres.097.032)
