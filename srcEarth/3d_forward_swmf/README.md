# SWMF-coupled Step 3 field snapshot

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
