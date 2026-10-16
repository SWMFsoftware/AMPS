# BON2020 v2.3 delivery gate

Stage 2 requires the official NASA BON2020 v2.3 package, but those bytes are
not present in this repository. NASA's software catalog distributes this
release through a repository request followed by secure email. The adapter
must not be inferred from the publication, older files, filenames, or another
release.

The project owner must place the delivered archive and documentation in a
controlled intake location (the archive need not be committed if its license
does not permit that). Before implementation resumes, record:

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

Only after that inspection may the project owner approve the two exact 1 AU
scenario records. The implementation must then preserve nuclei `Z=1..28`,
energy-per-nucleon semantics, valid-domain limits, package/version identity,
all generation inputs, and output hashes in the canonical dataset.

Until this gate is satisfied:

- there is no `BON2020Adapter` and no approved BON2020 transport input;
- Stage 2 is incomplete;
- the nearby `gcr_BON2020_LIS_*` files are unverified legacy artifacts and are
  prohibited for physical/science-facing runs; and
- AMS-02 may test GCR ingestion mechanics but may not be relabeled as either
  approved deep-space BON2020 scenario.

Catalog record:
<https://software.nasa.gov/software/MSC-26835-1>
