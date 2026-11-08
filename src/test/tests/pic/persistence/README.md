# T9 restart, data-file, and output integrity

The three manifest owners share one dependency-light executable but own
disjoint case IDs. Restart cases exercise the production AMPSPF1 envelope,
canonical state round-trip, deterministic continuation, empty/large payloads,
corruption/version/layout rejection, native endian/width rejection, canonical
cross-build width independence, atomic publication, and explicit low-level
legacy opt-in. The production restart compatibility switch defaults on to
preserve existing same-build AMPS files and can be disabled before recovery.

Data-file cases use a committed analytic affine background, concatenate writer
rank chunks, and consume them under another partition. Output cases pin the
production trajectory variable order/units and reject nonfinite records.

The suite does not claim that inventoried raw legacy AMR, trajectory scratch,
or stack files are cross-build portable. Their unsupported boundary is part of
the acceptance evidence rather than being hidden behind a successful same-host
read.

The actual production call path has a separate configured receipt:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-PIC-RESTART-NATIVE-001
```

That one-rank case saves and reads AMPS sampling and particle restart files,
then checks completed sample bytes/counters, AMR cell ownership, species,
position, velocity, global weight/time step, internal energy, opaque extension
bytes, and the scalar MPI RNG seed. Hybrid per-thread RNG continuation and
changed-rank ownership remain outside this T9 boundary.
