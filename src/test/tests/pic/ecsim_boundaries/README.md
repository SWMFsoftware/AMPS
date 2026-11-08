# ECSIM boundary-row policies

The ECSIM solve advances `delta E`.  This suite checks the production algebraic
row policy independently of mesh traversal: periodic faces create no physical
row, Dirichlet gives `delta E=0`, and homogeneous Neumann gives the boundary
value minus the average of one to three inward face-normal donors.  The current
face policy lacks the component/normal data needed for a physical perfect-
conductor closure, so that request is explicitly unsupported rather than being
silently aliased to Dirichlet.

```sh
python3 src/test/run_tests.py --test CORE-ECSIM-BOUNDARIES-001 \
  --output-dir /tmp/amps-core-t8-boundaries
```
