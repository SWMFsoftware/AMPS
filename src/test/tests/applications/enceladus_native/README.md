# Native Enceladus startup adapter

This exclusive suite performs a clean native `enceladus` build through the
mandatory repository preflight and imports all six gates from
`srcEnceladus/test/run_native_startup.py`. The gates cover CLI/default-deck
validation, rejected inward stripe geometry, rejected inaccessible plasma
cache, and generated pointer-mover/particle-processing hooks. They stop before
MPI initialization and do not claim an evolving E3 simulation.

```bash
python3 src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --test APPLICATION-ENCELADUS-NATIVE-STARTUP-001
```

The suite retains the configured executable and its hash, build transcript,
local receipt, and startup transcript. It changes the checkout-wide generated
tree and therefore always runs alone.
