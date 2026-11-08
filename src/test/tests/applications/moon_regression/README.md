# Moon regression global adapter

This exclusive suite integrates the complete registry in
`srcMoon/test/run_all_test.py`: U01–U26 and linked I11–I20. It performs a
clean Moon/SPICE build using the repository-mandated native preflight, then
invokes the application aggregate with `--skip-build`. Every local
PASS/FAIL/SKIPPED/ERROR becomes an exact global PASS/FAIL/SKIP/ERROR case;
the adapter never treats process launch as scientific success.

The default local mode is the bounded regression campaign. I11 and I13 remain
intentional non-execution SKIPs there, and dependent M2 cases may also SKIP.
Those statuses are preserved. This suite does not claim M1/M2 qualification.

Machine-specific SPICE locations are configuration data. Create an overlay
like this (outside the repository is allowed):

```json
{
  "schema_version": 1,
  "paths": {
    "spice_path": "/site/SPICE/cspice",
    "spice_kernels": "/site/SPICE/Kernels"
  },
  "capabilities": {
    "moon-spice": {"available": true, "policy": "optional"}
  }
}
```

Then run:

```bash
python3 src/test/run_tests.py \
  --test APPLICATION-MOON-REGRESSION-001 \
  --config /path/to/moon-machine.json
```

The adapter deletes only the generated repository-root `build/` directory.
It leaves the final Moon configuration in place and retains the executable,
build identity, build transcript, local aggregate summary, and all local
per-test artifacts beneath the global artifact directory.
