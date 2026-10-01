# CPL3D10

reserved, non-release test for a future generalized nonradial steady-spatial Piola-pushforward winding provider. It must verify the Section 6.5 characteristic/boundary equations, orientation/invertibility, divergence and surface-flux preservation, pushed HCS/separatrix/interface events, transformed gradients/focusing length, nonradial wind mass-per-flux, 3-D/1-D parity, fold rejection, and exact reduction to the schema-5 radial map at `R_w=R_scs`; it is not a schema-5 release dependency.

Run from any directory with `python3 test.py`. The launcher delegates to the global registry, so individual and cumulative gates execute the identical implementation.

Implementation: `test/test_stage14.py` calls public research kernels and immutable offline producers, with synthetic independent references and negative fixtures. [The Stage-14 guide](../../../docs/STAGE14_RESEARCH_EXTENSIONS.md) records the implemented domains and outstanding host/campaign gates. A software PASS is not an observational-campaign or production-adapter qualification.
