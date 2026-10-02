# SWCME/SSE test-fixture correction

This update fixes the eight failures reported by `test/run_tests.py --all`:
`SWBG3D01–07` and `SSE3D05`. It is a source overlay for the existing finite-SSE
checkout, including its preceding runtime/weak-shock corrections. It contains
directly edited C++ and Python sources, both input examples, detailed comments,
README updates and captured test evidence. There is no installer script.

## Why the tests failed

| Tests | Cause | Correction |
| --- | --- | --- |
| `SWBG3D01–07` | The C++ fixtures still opened the retired `examples/sep3d_swcme_mesh_background_20rs_1au.in`. Renaming the installed input therefore broke both fixture loading and parser-mutation setup. | Both accesses use one constant naming `examples/sep3d_swcme_sphere_mesh_background_20rs_1au.in`; the configuration gate explicitly verifies canonical `Sphere` geometry. |
| `SSE3D05` | The source test imported a fixed Earth observer from `sep3d_analytic_parker.in`. That observer does not intersect the finite active corridor selected in the SSE propagation example. Injection requires a valid observer, so configuration validation rejected the fixture before source preparation. | Construct a named fixed test probe on the example's normalized Parker curve, at half its finite arc length, and serialize it at 17-digit precision. Then exercise actual canonical finite-source preparation and the existing physical patch checks. |

The configuration factory continues to reject outside-corridor observers.
The default narrow-corridor fixture explicitly checks that its mirrored probe
is rejected. No production physics, observer position, spatial guard or
coupling mode is changed by the test correction.

The generated probe is a test location, not an Earth ephemeris. The source-free
propagation input needs no observers. A real injection run should select the
Parker source angles and corridor width for its actual fixed observer, or use
full-domain allocation. Comments in the SSE input explain the two allocation
modes, the required inactive sentinels, the box settings and the independent
refinement tube. The SSE input included here retains `parker-tube`, a 30-solar-
radius active solar neighborhood and a 0.05-AU corridor radius at 1 AU.

## Install and rerun

Extract the `srcSEP3D` entries over the existing AMPS checkout. For example,
from its root, with the downloaded archive beside it:

```sh
tar -xzf SEP3D_test_fixture_fix_20261002.tar.gz srcSEP3D
cd srcSEP3D
env MAKEFLAGS="-j16" test/run_tests.py --all --amps-source .. --make-config ../Makefile.conf --output-dir test_output/all --rebuild
```

`--rebuild` is necessary: changing the Python runner alone leaves the old C++
fixture filenames in a previously compiled `test/stage1`. The canonical
spherical input is included under its new name. The retired input can remain
locally, but these tests do not read it. This overlay assumes the finite-SSE
implementation and shared model libraries are already installed; it is not
a replacement for the entire AMPS checkout.

For a quick focused rerun from `AMPS/srcSEP3D`:

```sh
test/run_tests.py --group SSE3D --group SWBG3D --rebuild --output-dir test_output/swcme-fixtures
```

## Verification

The exact reported failures were reproduced with the renamed spherical input
absent under its old name and the current Parker-tube SSE input: **9 PASS,
8 FAIL** in the two affected groups. After correction, the same 17 tests gave
**17 PASS, 0 FAIL, 0 SKIP, 0 ERROR**. `SSE3D05` also passed with the four documented
full-domain allocation assignments; the original Parker-tube input bytes were
restored afterward.

The complete 198-case runner catalog was rebuilt and executed locally:

```sh
env MAKEFLAGS="-j16" test/run_tests.py --all --amps-source .. --output-dir ../../runner-fixture-verification/all --rebuild
```

Result: **168 PASS, 0 FAIL, 30 SKIP, 0 ERROR**. This environment has no configured
AMPS `Makefile.conf`, so `BLDL3D01` was skipped. The other skips require a linked
native executable, explicitly supplied comparison sources, reviewed external
evidence or a deferred live-coupling capability. The local command therefore
omitted `--make-config`; it does not claim a successful NASA AMPS/MPI build or
native run. Rerun the original configured command above on the target machine.

The archive's optional `verification/SEP3D_test_fixture_fix_20261002` directory
contains before/after group JSON and logs, the full-domain check and the full
runner JSON/JUnit summary and console output. Extracting only `srcSEP3D`, as
above, leaves that evidence out of the source tree.
