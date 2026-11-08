# Enceladus E3 native campaign adapter

This exclusive extended suite owns all 13 checks from the 10,000-step E3
H2O+CO2 native runner: process/configuration, time step, completion/injection,
species loss and conservation ledgers, output/finalization, and flyby Tecplot,
mesh-resolution, and PNG products.

Routine machines retain every identity as `SKIP`. Enable the costly execution
only in a machine overlay approved for MPI and the required runtime:

```json
{
  "schema_version": 1,
  "capabilities": {
    "enceladus-e3": {"available": true, "policy": "optional"}
  }
}
```

```bash
python3 src/test/run_tests.py \
  --test APPLICATION-ENCELADUS-E3-001 \
  --config /path/to/e3-machine.json
```

The enabled path performs a clean native Enceladus build, uses the configured
MPI launcher argv/rank flag, and preserves the complete local run directory.
A routine SKIP is not E3 validation; compiler/startup evidence also cannot
replace this evolving simulation.
