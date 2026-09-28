# SCCM validation runner

The test system mirrors the staged release gates in
`model/testing_validation.md`. Every canonical identifier can run alone, while
a stage gate is cumulative: Stage 2 always re-runs Stages 0 and 1.

```sh
make test-stage0
make test-stage1
make test-stage2
python3 test/run_tests.py --test PFSS3D01
python3 test/run_tests.py --list
```

The runner writes JSON and JUnit evidence below `build/test-results/`. Each
`test/individual/ID/` directory contains a stand-alone launcher, the roadmap
description, and reference metadata. The launchers delegate to the global
registry, avoiding a second implementation that could drift.

The Python entry points support Python 3.7 and newer. In particular, type
annotations are postponed so the runner imports correctly on Python 3.8-based
NASA HEC software stacks, where evaluating `list[str]` directly would raise
`TypeError: 'type' object is not subscriptable`.

Documentation tests are launched with `unittest discover -s test` rather than
by converting a file path to a dotted module name. This avoids both the missing
`test/__init__.py` failure and collisions with Python's own `test` package on
Python 3.8 installations.
