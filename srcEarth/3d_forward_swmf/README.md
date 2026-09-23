# SWMF-coupled Step 3 field snapshot and Step 4 trajectories

Before each requested backward-product callback,
`PrepareGlobalSWMFCoupledMagneticFieldForCutoff()` gathers authoritative owner-cell
magnetic field and plasma velocity into the compact Mode3D arrays. Electric field is
derived cell-by-cell as

```text
E = -v × B.
```

The arrays are published with one `FieldProvider.h` metadata record containing source
`PIC::CPLR:SWMF`, the configured reference epoch plus authoritative PT simulation-time
offset, GSM/SI units, domain, derived-E interpolation mode, and a deterministic
generation ID. Cutoff and density/flux retain that ID and fail if a replacement is
published between products.

Assembly rejects unavailable buffers, incomplete/duplicate ownership, and non-finite
owner or reduced values. Those errors are not mapped to forbidden particle access.
Step 3 does not change the coupled field values, AMR interpolation, trajectory mover,
termination policy, or test thresholds.

Roadmap Step 4 uses that published generation through the same backend-neutral
`TrajectoryRequest`/`TrajectoryResult` contract as standalone Mode3D and gridless.
Every request fingerprint is checked against the active coupled snapshot before a
trajectory starts, and every result records the verified fingerprint, mover, backward
convention, exact termination, retry/extension provenance, and optional complete outer-
boundary phase-space state. A coupled snapshot replacement therefore cannot be mixed
silently into an in-progress cutoff/flux batch.

Although the coupled snapshot contains the derived electric field, the released Step 4
backtracer remains frozen magnetic-only. Electric or explicitly time-dependent
characteristics and `PhysicalBackwardTime` requests fail fast; the code does not treat
the presence of an E array as permission to use the static antiparticle shortcut.
