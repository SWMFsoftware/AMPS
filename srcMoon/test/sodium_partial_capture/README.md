# Approved P07 Na ground-product partial capture controls

This supporting runner compiles normal production `SodiumPartialCapture`,
`ReducedKinematics`, and the source-identity service. It does not replace the
existing all-state He/Na charge-exchange registry, enable particle reactions,
or complete U17/I24. The measured source channel is H+ + Na -> H(1s) +
Na+(2p6); its twelve native4–25keV/u points are read without interpolation.
The approved forward-COM direction is explicitly an experimental closure,
not an angular measurement. Actual source/beam state and Q are preserved.

```sh
python3 -B srcMoon/test/test_m3_sodium_partial_capture.py \
  --output-dir test_output/srcMoon/m3-na-partial-NEW \
  --reference-dir test_output/srcMoon/m3-atomic-online-databases-20261118T083400Z \
  --mass-inventory test_output/srcMoon/m3-approved-products-20261118T163142Z/source-acquisition-v2/reference_mass_inventory.json
```

Inputs are explicit: the reference directory contains the prior native IAEA
extraction/manifests; the mass inventory comes from the maintained
[NIST preparation script](../physical_state_sources/README.md). The runner
verifies source-response hashes and numerical provenance before using values.
An unavailable required source yields SKIPPED, altered bytes or malformed
metadata yield ERROR. A genuine evaluated numerical mismatch yields FAIL.
Exits0/1/2/77 distinguish those outcomes. Outputs are new-directory logs,
result/manifest, native node/product CSV, and PNG/EPS comparison figures.

The declared source-target-rest proton beam uses energy E[eV/u], so
`g=sqrt(2*E*e/u)` and `lambda=n*sigma*g`. Here u[kg] and e[J/eV] come from
the acquired CODATA source. Reference density2.5e6m^-3 is a mathematical
input, not a measured lunar wind. Atomic state energies use the acquired
NIST ground-state ionization records; separate neutral-atom reference zeros
are energy definitions. Unavailable isotope shifts/uncertainties stay explicit.
No microscopic accuracy claim follows from an exact conservation result.

The independent numerical oracle uses Decimal arithmetic with sixty digits
for unit conversion, speed, rate and state-energy difference, plus direct
COM/momentum/energy identities. Comparisons use a predeclared1e-12 relative
representation-level allowance, not an observational or physical tolerance.
Controls also cover all native nodes, zero density, nonnode/out-of-range
energies, invalid constants, forged receipts, wrong states/isotopes, and
transactional preservation after refusal. The source SHA256 fixes the partial
channel; a different file cannot be authorized by a caller-supplied hash.

These are supporting extensions within existing U17 meaning. No new permanent
ID is assigned and native PIC/mover transport remains a separately qualified
adapter. Normal M0 and shared-core behavior remain unchanged by this pure
library; the current application build still links these same objects.

The native D119990 label specifies the outgoing states but does not supply
the initial Na state, underlying bibliography or an uncertainty table. The
3s initial state is the user-approved reference cohort, not recovered
experimental metadata. `physicalTargetStateQualified()` therefore stays false;
this component cannot activate physical lunar CEX until the state crosswalk
is qualified. That limit is retained in product provenance and the manifest.

Reproduce the exact twelve native records with the maintained parser:

```sh
python3 -B srcMoon/test/sodium_partial_capture/prepare_source.py \
  --source-html test_output/srcMoon/m3-atomic-online-databases-20261118T083400Z/collisiondb_Na_H1s_119990.html \
  --output-dir test_output/srcMoon/na-partial-source-NEW
```

Use that output as `--reference-dir`. The parser binds the original response,
Scatter renderer, native units and exact reproduced CSV hash. The runner
independently repeats extraction and rejects an altered derived NIST inventory
even when every original source hash is intact.
