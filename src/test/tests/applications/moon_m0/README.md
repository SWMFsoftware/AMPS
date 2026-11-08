# Moon M0 global boundary

This suite gives every `srcMoon/test/run_integration_tests.py` I01–I10 case a
stable global identity and verifies that its individual runner and README are
present. It currently reports every scientific case as `SKIP`: the local M0
orchestrator performs multiple checkout-wide builds but has not yet been
updated to the repository's mandatory process check, exact `build/` removal,
configuration, and `MAKEFLAGS=-j16 make amps` preflight for every rebuild.

The SKIPs are intentional evidence of non-execution, not failures and not
passes. Run the boundary audit with:

```bash
python3 src/test/run_tests.py --test APPLICATION-MOON-M0-001
```

The authoritative local manual campaign remains:

```bash
python3 srcMoon/test/run_integration_tests.py --output-dir /new/output/path
```

Do not interpret that manual command as globally qualified until its native
build workflow is brought into compliance and this adapter imports its exact
result registry.
