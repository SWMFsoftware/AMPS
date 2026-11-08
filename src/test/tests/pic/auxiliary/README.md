# PIC auxiliary numerical kernels

This dependency-light T8 suite exercises production kernels shared by AMPS ray
tracing, column integration, tabulated stopping-power providers, and background
field ownership.  It uses exact slab-intersection geometry, analytic integrals
of constant/affine fields, and hand-computed table values.  Out-of-range or
malformed inputs are checked through typed status codes rather than accepting
silent extrapolation.

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-PIC-AUXILIARY-001 \
  --output-dir /tmp/amps-core-t8-auxiliary
```

This suite qualifies the dependency-light kernels and the production call from
ray/column code.  It does not by itself qualify distributed mesh traversal,
application-specific physical stopping-power tables, or optional coupling
libraries.
