# U02 — Legacy-target isolation and source-level guard

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U02_target_isolation/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

The source guard uses whitespace-normalized complete token sequences for the
two callback assignments. This permits ordinary C++ line wrapping while still
requiring the legacy sphere interaction callback and the overload-resolved
production injection callback. It remains a source/configuration guard, not a
linked-runtime result.
