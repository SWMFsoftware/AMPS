# Retired BON2020 v2.3 intake proposal

BON2020 is not used in production and is not required to complete Stage 2.
Stage-2 validation instead uses checksum-pinned downloaded published spectra,
including a direct AMS-02 proton case that drives full Geant4 transport and
the proxy writers. This directory is retained only to document why nearby
unverified BON2020-named files remain prohibited.

If a future, separately approved study chooses to evaluate BON2020, place the
delivered archive and documentation in a controlled intake location (the
archive need not be committed if its license does not permit that). Before
implementing such an optional adapter, record:

- the original delivered filename, byte length, SHA-256 digest, delivery date,
  license/redistribution terms, and NASA catalog identifier `MSC-26835-1`;
- every included document and its digest;
- the package's own version identifier and supported execution platform;
- documented input fields, units, energy coordinate, valid energy domain,
  species/isotope convention, radial coordinate, and output normalization;
- the package-defined controls for solar-minimum and solar-maximum conditions,
  including epoch/averaging or modulation inputs exactly as documented; and
- one untouched vendor example with its documented output, for an adapter
  acceptance oracle.

Only after that inspection could the project owner approve exact 1 AU scenario
records for that future study. The implementation would need to preserve nuclei `Z=1..28`,
energy-per-nucleon semantics, valid-domain limits, package/version identity,
all generation inputs, and output hashes in the canonical dataset.

Under the current approved plan:

- there is no `BON2020Adapter` and no approved BON2020 transport input;
- Stage 2 is complete without BON2020;
- the nearby `gcr_BON2020_LIS_*` files are unverified legacy artifacts and are
  prohibited for physical/science-facing runs; and
- AMS-02 may test GCR ingestion mechanics but may not be relabeled as either
  approved deep-space BON2020 scenario or production environment.

Catalog record:
<https://software.nasa.gov/software/MSC-26835-1>
