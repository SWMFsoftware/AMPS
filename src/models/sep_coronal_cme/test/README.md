# SCCM validation runner

The test system mirrors the staged release gates in
`model/testing_validation.md`. Every canonical identifier can run alone, while
a stage gate is cumulative: Stage 6 always re-runs Stages 0 through 5.

```sh
make test
make test-stage0
make test-stage1
make test-stage2
make test-stage3
make test-stage4
make test-stage5
make test-stage6
python3 test/run_tests.py --test PFSS3D01
python3 test/run_tests.py --all
python3 test/run_tests.py --list
```

The cumulative release-gate sizes are fixed by contract: Stages 0 through 6
select 13, 22, 53, 82, 96, 106, and 117 tests, respectively.  Every Make
target passes the corresponding value through `--expect-count`, while the
runner independently validates its registry before executing anything.  Thus
a partial update cannot misreport the old 53-test Stage 0--2 registry as a
successful Stage 6 run.  A correct Stage 6 log starts with:

```text
SELECTION: stage=6; tests=117; registry-total=117
```

To diagnose a copied source tree, `python3 test/run_tests.py --list | wc -l`
must print `117`, and the last listed identifier must be `RH3D11`.  If it does
not, replace the complete model package (at minimum both `makefile` and
`test/run_tests.py`) instead of mixing files from different release stages.

The runner writes JSON and JUnit evidence below `build/test-results/`. Each
`test/individual/ID/` directory contains a stand-alone launcher, the roadmap
description, and reference metadata. The launchers delegate to the global
registry, avoiding a second implementation that could drift.

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
