# I22 — Photoionization for He, Ne, Ar, and future H2O/OH

The [atomic-family preparer](../../../../utility/chemical_database/state_profiles/ATOMIC_FAMILY_CONTRACT.md)
emits six explicit source-effective photon candidates, one for each
isotope and source endpoint. Each retains the D04 printed threshold and
mean, physical ground-ion energy and native row identity. The worksheet
supports source/state review before recoil or ground-ion yields are selected.
The U15/U24 and full linked I22 gates remain applicable.

The complete specification is [reference/acceptance.json](reference/acceptance.json).
Original acceptance: Ensemble rate error <=2%; linked branching statistically consistent at95% confidence; complete declared family sums to1 within numerical precision. U15 separately retains99% intervals and source-precision transcription.

The [approved2026-11-18 addendum](../../M3_QUALIFICATION_ADDENDUM_20261118.md)
fixes total-removal endpoint rate `-log(S/N)/T` as primary, relative0.002
numerical bias, simultaneous95% exact branch intervals and10% relative
half-width for every positive source branch. Partial rates/event-time MLE are
separately named diagnostics. The complete family and full U24 dependency stay
required. Protocol approval does not supply missing physical/native evidence.

This individual owner is currently a **formal capability/predecessor gate**.
It preserves source/input/executable and same-campaign predecessor evidence,
but full particle execution/scoring is unavailable while the recorded physical
kinetic policies/production event dispatch are missing. It reports SKIPPED,
never a scalar-rate PASS relabeled as chemistry qualification. The capability
register and [completion proposal](../../M3_COMPLETION_PROPOSAL.md) identify
the required implementation/data; new source recommendations are separate.

From the repository root (an existing dedicated campaign directory):

```sh
python3 -B srcMoon/test/integration-tests/I22/test.py \
  --amps "$PWD/amps" --input srcMoon/examples/m3_reaction_rates.in \
  --mpi-ranks 8 --campaign-mode regression \
  --output-dir test_output/srcMoon/M3-formal-NEW
```

The runner creates `<output>/I22/` with its real result/manifest,
stdout/stderr and an exact input snapshot when available. No `receipt.txt`
or comparison figure exists for a run that did not execute. Exit codes:
`0=PASS`, `1=FAIL`, `2=ERROR`, `77=SKIPPED`. Do not overwrite prior evidence.

The owner now retains `protocol_preparation.json`, binding the concrete
[physical model](../../M3_PHYSICAL_MODEL_SPECIFICATION.md) and
[qualification protocol](../../M3_QUALIFICATION_PROTOCOL.md). Matching
identities and a structural preparation PASS do not launch particle chemistry
or change this formal gate. Missing inputs remain SKIPPED; corrupt identities
are ERROR. Full source-family, product and predecessor requirements remain.
