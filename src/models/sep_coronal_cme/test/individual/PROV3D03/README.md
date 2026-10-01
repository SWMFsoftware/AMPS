# PROV3D03

missing metadata, mixed epochs, ambiguous units, response reuse, and calibration/validation data leakage are rejected.

Run from any directory with `python3 test.py`. The launcher delegates to the global registry, so individual and cumulative gates execute the identical implementation.

Implementation: the identically named unittest class in `test/test_stage12.py` exercises `tools/preprocessing/`. It uses synthetic maps, geometry, covariance and candidate products; no event observation is hardcoded in C++. The [Stage-12 guide](../../../docs/STAGE12_OBSERVATION_PREPROCESSING.md) explains inference assumptions, immutable metadata and role ownership. The [synthetic CLI example](../../../examples/stage12/README.md) provides a complete runnable asset/selection workflow.
