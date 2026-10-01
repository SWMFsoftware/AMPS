# ELL3D11

reserved, non-release test for future CME deflection and rotation. Analytic pure-deflection, pure-rotation, and combined spheroid histories reproduce center velocity, angular velocity, surface velocity, and normal speed, including both `d_c*de_r/dt` and `Q_dot` terms. The rotation remains in `SO(3)`, fixed direction/attitude is recovered exactly, and componentwise quaternion interpolation or omitted history derivatives fail.

Run from any directory with `python3 test.py`. The launcher delegates to the global registry, so individual and cumulative gates execute the identical implementation.

Implementation: `test/test_stage14.py` calls public research kernels and immutable offline producers, with synthetic independent references and negative fixtures. [The Stage-14 guide](../../../docs/STAGE14_RESEARCH_EXTENSIONS.md) records the implemented domains and outstanding host/campaign gates. A software PASS is not an observational-campaign or production-adapter qualification.
