# Charge-exchange scattering evidence and isotropic model

The user requested cataloguing the
[I24 source review](../../srcMoon/test/integration-tests/I24/ANGULAR_SOURCE_FOLLOWUP_20261124.md)
and its conditional uniform-angle model. The shared database uses
an explicit release lock. Its existing catalogue and CSV readers
can expose both entries without changing the runtime ABI or adding
a reaction or species selection.

## Source inventory and publication contract

The shared r8 catalogue contains the byte-verified IAEA D119990
integral table, the live CollisionDB record, and metadata for the
Knoop and Courbin publications. Preserve all 27 r8 source packets.
Expose an authored source review for the D119990 crosswalk and add a
reference from the integral table's catalogue metadata. Keep the
integral table's 12 rows and original label verbatim. Publish the
isotropic direction law as a species-independent mathematical
model through `--include-model-definitions`. The new definitions
carry `physical_activation=false`.

The locally retained PDF Chapter 5
(`7daca10c4ecf37bd79f0a5d353cfc9c6d6adc9c2b0277d74aa01d5381c5c1785`)
is identified by `chemistry-archive-20261124:knoop-h-na`.
Keep it metadata-only. The native IAEA corroboration
`iaea-d119990-20261124:live-html` has SHA-256
`1015ec7a843e02ef7fc1ddf97f3589ac4dc4823f12bb28e9434313876dd90fd6`.
The exact existing integral CSV SHA-256 is
`2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93`.
The release must introduce zero numerical angular data rows. The
reported figures do not establish an absolute angular law or a
bound outside their measured range.

At 4 keV/u the IAEA compilation lists `6.34e-15 cm²`. In Knoop
Fig. 5.6 the prompt outer-shell `H(n=1)` symbol lies between
`1e-17` and `1e-16 cm²`. This is a figure-scale consistency
check. It does not establish the IAEA record's original channel.
Knoop §5.2.3 reports the corresponding small-angle laboratory
measurements for energies from 7 keV/u upward, with a 25 keV/u
peak-normalized trace in Fig. 5.10. The 4 keV/u panel omits that
branch. Chapter 4 Eq. (4.10) defines the small-angle measurement
convention. [IAEA CollisionDB](https://db-amdis.org/collisiondb/datasets/119990/)
[Knoop Chapter 5](https://pure.rug.nl/ws/portalfiles/portal/2855778/c5.pdf)
[Knoop Chapter 4](https://pure.rug.nl/ws/portalfiles/portal/2855777/c4.pdf)

## Model and verification scope

`isotropic-com-two-body-scattering-v1` is an experimental fallback.
Draw two independent uniform variates and sample `mu = cos(theta)`
uniformly in `[-1,1]` and the azimuth uniformly in `[0,2π)`. The
resulting directions have density `1/(4π)` per sr. One sampled
direction sets the orientation of both opposite final momenta in
the pair center-of-mass frame. The existing
`two-body-conservative-v1` definition supplies their magnitude and
frame transformation from caller-supplied masses, available energy,
reaction Q and total incoming momentum. The angular descriptor
supplies no integral rate, physical state mapping or applicability
bound. Any sensitivity run requires a separately declared rate and
explicitly selected collision model.

The publication oracle is the existing r8 table digest and source
record identities. Check both descriptors against the reviewed
source findings and the elementary isotropic solid-angle measure.
The CLI and Python must return the new entries. C++ must read their
rows using the existing model-definition table and see the physical
scope warning on the IAEA table. The source permission for the
Knoop PDF must remain held. These are software and catalogue checks.
Formal I24 remains SKIPPED and no lunar test ID changes.

The AMPS roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
The canonical validation-data roadmap SHA-256 is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The AMPS copy controls the approved early development and resource
amendments. The shared scientific rules control source admission and
the formal I24 gate.
