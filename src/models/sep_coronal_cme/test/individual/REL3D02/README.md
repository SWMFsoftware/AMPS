# REL3D02

the cross-application driver receives independent executable and bundle paths; neither application searches a sibling source tree.

Run from any directory with `python3 test.py`. The launcher delegates to the global registry, so individual and cumulative gates execute the identical implementation.

Implementation: `test/test_stage13.py` verifies the release machinery in [the Stage-13 guide](../../../docs/STAGE13_RELEASE_QUALIFICATION.md). Actual production qualification additionally requires clean owning-app builds, MPI/convergence/campaign evidence and D1--D10.
