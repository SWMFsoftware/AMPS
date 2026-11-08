# Enceladus forward-model global adapter

This suite executes `srcEnceladus/test/run_tests.py` and imports all 11 local
software checks plus all 31 independent numeric reference comparisons. The
application-local receipt remains authoritative; the adapter rejects missing,
extra, reordered, duplicated, or unsupported results and checks that the local
exit status agrees with the receipt.

The run is dependency-light. It compiles `ForwardModel.cpp` with the compiler
argv selected by the global machine configuration, but it does not configure
or execute a native AMPS particle simulation. PASS therefore verifies the
forward-model parser/law/reference contracts only, not a flyby-density
qualification or a native Enceladus application run.

Run it from the repository root:

```bash
python3 src/test/run_tests.py --test APPLICATION-ENCELADUS-FORWARD-001
```

Evidence is retained below the global run's
`artifacts/APPLICATION-ENCELADUS-FORWARD-001/` directory, including the local
JSON receipt, transcript, and SVG reference comparison.
