# SCCM validation runner

The test system mirrors the staged release gates in
`model/testing_validation.md`. Every canonical identifier can run alone, while
a stage gate is cumulative: Stage 10 always re-runs Stages 0 through 9.

```sh
make test
make test-stage0
make test-stage1
make test-stage2
make test-stage3
make test-stage4
make test-stage5
make test-stage6
make test-stage7
make test-stage8
make test-stage9
make test-stage10
python3 test/run_tests.py --test PFSS3D01
python3 test/run_tests.py --all
python3 test/run_tests.py --list
```

The cumulative release-gate sizes are fixed by contract: Stages 0 through 10
select 13, 22, 53, 82, 96, 106, 117, 128, 153, 173, and 196 tests,
respectively. Every Make
target passes the corresponding value through `--expect-count`, while the
runner independently validates its registry before executing anything.  Thus
a partial update cannot misreport the old 53-test Stage 0--2 registry as a
successful later-stage run. A correct Stage 10 log starts with:

```text
SELECTION: stage=10; tests=196; registry-total=196
```

To diagnose a copied source tree, `python3 test/run_tests.py --list | wc -l`
must print `196`, and the last listed identifier must be `XM3D02`. If it does
not, replace the complete model package (at minimum both `makefile` and
`test/run_tests.py`) instead of mixing files from different release stages.

The runner writes JSON and JUnit evidence below `build/test-results/`. Each
`test/individual/ID/` directory contains a stand-alone launcher, the roadmap
description, and reference metadata. The launchers delegate to the global
registry, avoiding a second implementation that could drift.

`ARCHSCCM01` runs five source-audit regressions before the actual built
archive/public-header audit. Compiled `sep_*.o` files beside neutral sources
are not UTF-8 text: only regular C/C++ source/header files enter the text scan.
Invalid encoding in a genuine source remains a named failure. Nested source
checks, forbidden includes, upward neutral dependencies, archive symbol checks
and the external C++17 consumer remain enforced. See
[`individual/ARCHSCCM01/README.md`](individual/ARCHSCCM01/README.md) for the
failure cause and a one-case rerun command. These regressions are internal to
the existing canonical gate, so the Stage-12 aggregate still has 206 shared IDs.

The Python entry points support Python 3.7 and newer. In particular, type
annotations are postponed so the runner imports correctly on Python 3.8-based
NASA HEC software stacks, where evaluating `list[str]` directly would raise
`TypeError: 'type' object is not subscriptable`.

Generated and mutated text is written through `Path.open(..., newline="\n")`
rather than the newer `Path.write_text(..., newline="\n")` convenience API.
That preserves byte-deterministic LF output while remaining compatible with
Python 3.7 and 3.8.

Documentation tests are launched with `unittest discover -s test` rather than
by converting a file path to a dotted module name. This avoids both the missing
`test/__init__.py` failure and collisions with Python's own `test` package on
Python 3.8 installations.
