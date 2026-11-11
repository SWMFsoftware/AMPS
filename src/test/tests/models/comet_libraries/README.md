# Imported reusable comet libraries

`CORE-COMET-LIBRARIES-001` executes the maintained reaction-network,
collision-database, radiative-cooling and resonance-fluorescence C++ tests,
and authenticates the reviewed reference tables and generated include. The
adapter builds in the global runner's per-suite scratch directory and requires
only the configured C++17 compiler and Python. It does not configure AMPS or
read production archives.

Run from the repository root:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-COMET-LIBRARIES-001
```

These checks validate the reusable equations, units, lookup and data contracts.
Native DIM1 transport/chemistry and observational science validation are
separate Hartley2 campaigns. New registered cases must also be added to
`src/test/TEST_CATALOG.md` through its generator.

## Exact case inventory and reference boundaries

| Case | Source and independent contract |
| --- | --- |
| `COMET-REACTION-NETWORK-001` | `src/models/chemical_solver/test_reaction_network.cpp`: explicit unimolecular/binary/productless reactions, multiplicities, rates, indexed parent-channel order and validation diagnostics. Protects the generic network representation used by comet chemistry. |
| `COMET-COLLISION-DATABASE-001` | `src/models/cometary_physics/test_collision_database.cpp`: known constant/temperature-power laws, symmetric unordered lookup, validity range, missing pair and duplicate rejection. Prevents an absent/out-of-range record being silently accepted as physical zero collision. |
| `COMET-RADIATIVE-COOLING-001` | `src/models/cometary_physics/test_radiative_cooling.cpp`: independently evaluated temperature branches, cgs-to-SI energy conversion, known optical attenuation and disabled/zero limits. Protects units and branch selection. |
| `COMET-FLUORESCENCE-001` | `src/models/resonance_fluorescence/resonance_fluorescence.cpp` with `tests/test_resonance_fluorescence.cpp`: maintained independent line-profile/normalization and metadata checks. Protects reusable spectroscopy without invoking a comet transport model. |
| `COMET-REFERENCE-DATA-001` | Read-only `tools/generate_reference_tables.py --check` plus `tests/test_reference_data.py`: reviewed data hashes and generated-include agreement. Does not download, regenerate or replace authoritative scientific tables. |

The adapter consumes `tools.cxx` and `tools.python` from the required machine
profile. C++ builds use C++17, -O2, -Wall, -Wextra and -Wpedantic and write only
to the suite work directory. The maintained assert-based fixtures must not be
compiled with assertions disabled. It executes only requested owned cases;
missing/setup failures retain ERROR instead of deleting their IDs. A nonzero
compiled test or failed data check publishes FAIL for that case. Native AMPS
configuration/archives, SPICE and particle allocation are not involved.

The libraries' physics/data/API descriptions remain authoritative in their
own READMEs: [reaction network](../../../../models/chemical_solver/README.md),
[collision/cooling](../../../../models/cometary_physics/README.md), and
[fluorescence](../../../../models/resonance_fluorescence/README.md).
The five cases pass in the completed 2026-11-10 core campaign. That evidence
is distinct from Hartley2 Stage-1--10 application qualification and does not
close the separate legacy 3-D geometry/runtime issue.

Together with generic INP-05 and the native-configuration governance case,
these are the earlier seven additions from the 1028-case preintegration
baseline to the 1035-case repair baseline. Subsequent numerical/cross-ray/
coverage development adds 32 more; the complete current core registry has
1067 cases. Exact final counts, commands and identity comparisons are in
[`src/test/STATUS.md`](../../../STATUS.md).
