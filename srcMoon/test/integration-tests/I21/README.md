# I21 — Sodium photoionization survival and Na-to-Na+ conversion

The complete specification is [reference/acceptance.json](reference/acceptance.json).
Original acceptance: Ensemble lifetime relative error <=2%; inverse-square rate error <=0.1%; Na-to-Na+ species budget closes under the reviewed physical convention.

The [approved2026-11-18 addendum](../../M3_QUALIFICATION_ADDENDUM_20261118.md)
allocates relative0.002 numerical bias and fixes endpoint-survival lifetime
as the homogeneous primary. The owner snapshots that approved protocol;
actual source/state providers, native execution and formal predecessors remain
required. Numerical/statistical approval is no longer a missing permission.

This individual owner is currently a **formal capability/predecessor gate**.
It preserves source/input/executable and same-campaign predecessor evidence,
but full particle execution/scoring is unavailable while the recorded physical
kinetic policies/production event dispatch are missing. It reports SKIPPED,
never a scalar-rate PASS relabeled as chemistry qualification. The capability
register and [completion proposal](../../M3_COMPLETION_PROPOSAL.md) identify
the required implementation/data; new source recommendations are separate.

From the repository root (an existing dedicated campaign directory):

```sh
python3 -B srcMoon/test/integration-tests/I21/test.py \
  --amps "$PWD/amps" --input srcMoon/examples/m3_reaction_rates.in \
  --mpi-ranks 8 --campaign-mode regression \
  --output-dir test_output/srcMoon/M3-formal-NEW
```

The runner creates `<output>/I21/` with its real result/manifest,
stdout/stderr and an exact input snapshot when available. No `receipt.txt`
or comparison figure exists for a run that did not execute. Exit codes:
`0=PASS`, `1=FAIL`, `2=ERROR`, `77=SKIPPED`. Do not overwrite prior evidence.

The owner now retains `protocol_preparation.json`, binding the concrete
[physical model](../../M3_PHYSICAL_MODEL_SPECIFICATION.md) and
[qualification protocol](../../M3_QUALIFICATION_PROTOCOL.md). Matching
identities and a structural preparation PASS do not launch particle chemistry
or change this formal gate. Missing inputs remain SKIPPED; corrupt identities
are ERROR. Full source-family, product and predecessor requirements remain.
