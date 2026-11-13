# U14 — Na photoionization lifetime and survival

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U14_na_photoionization/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

U14 measures the lifetime/shadow kernel only. The November 12 production audit
found that the previously configured generic Exosphere reaction processor had
no Na entry, despite the local lifetime passing. The canonical input now uses
`Moon::Chemistry::ProcessParticle`, which applies the same 130500 s legacy
lifetime and shadow gate through PIC's actual particle-list dispatcher.
The linked 4096-particle sunlit/shadow survival and event-ledger checks are
registered as an I11 prerequisite; see
[BUDGET_CONTRACT.md](../integration-tests/I11/BUDGET_CONTRACT.md).
Neither that wiring check nor U14 is a D04 rate validation or Na+ product
qualification. A charged Na production configuration fails explicitly until
product weights/kinetics are qualified.
