# SHEATH3D03

Stage-11A and Stage-11B capability gates remain independent; enabling either cannot silently enable the other.

Run from any directory with `python3 test.py`. The launcher delegates to the global registry, so individual and cumulative gates execute the identical implementation.

Implementation: `test/tests_stage11.cpp` calls the public `discontinuity_transport.h` owning providers. It uses independent manufactured field/orbit/finite-volume references. Supplied-family domains, frame invariants, passive-wave policy and integration limits are documented in the [Stage-11 guide](../../../docs/STAGE11_DISCONTINUITY_TRANSPORT.md). Passing these planar-family checks does not qualify a curved CME sheath or finite composite transition.
