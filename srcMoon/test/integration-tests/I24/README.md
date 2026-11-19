# I24 — Charge-exchange reactions with solar-wind ions

The complete specification is [reference/acceptance.json](reference/acceptance.json).
Original acceptance: Direct rate error <=1%; ensemble lifetime error <=3%; correct products and exact charge/mass conservation under the reviewed physical convention; unsupported channels OFF. U17 retains >=100000 probabilistic-branch events.

The [approved2026-11-18 addendum](../../M3_QUALIFICATION_ADDENDUM_20261118.md)
allocates relative0.003 numerical bias and fixes endpoint-survival lifetime
as primary. Complete source/state products, U17's existing event minimum and
actual native/numerical/physical-ledger/predecessor evidence remain required.
An experimental partial channel cannot replace this complete formal scope.

This individual owner is currently a **formal capability/predecessor gate**.
It preserves source/input/executable and same-campaign predecessor evidence,
but full particle execution/scoring is unavailable while the recorded physical
kinetic policies/production event dispatch are missing. It reports SKIPPED,
never a scalar-rate PASS relabeled as chemistry qualification. The capability
register and [completion proposal](../../M3_COMPLETION_PROPOSAL.md) identify
the required implementation/data; new source recommendations are separate.

From the repository root (an existing dedicated campaign directory):

```sh
python3 -B srcMoon/test/integration-tests/I24/test.py \
  --amps "$PWD/amps" --input srcMoon/examples/m3_reaction_rates.in \
  --mpi-ranks 8 --campaign-mode regression \
  --output-dir test_output/srcMoon/M3-formal-NEW
```

The runner creates `<output>/I24/` with its real result/manifest,
stdout/stderr and an exact input snapshot when available. No `receipt.txt`
or comparison figure exists for a run that did not execute. Exit codes:
`0=PASS`, `1=FAIL`, `2=ERROR`, `77=SKIPPED`. Do not overwrite prior evidence.

The owner now retains `protocol_preparation.json`, binding the concrete
[physical model](../../M3_PHYSICAL_MODEL_SPECIFICATION.md) and
[qualification protocol](../../M3_QUALIFICATION_PROTOCOL.md). Matching
identities and a structural preparation PASS do not launch particle chemistry
or change this formal gate. Missing inputs remain SKIPPED; corrupt identities
are ERROR. Full source-family, product and predecessor requirements remain.
