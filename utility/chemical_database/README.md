# Shared chemical catalogue tools

The user-requested persistent installation is
`/data/vtenishe/CHEMICAL_DATABASE`. Its detailed API/usage manual is generated
from [DATABASE_README.md](DATABASE_README.md). The independent C++ runtime and
scalar services are in [src/models/chemical_database](../../src/models/chemical_database/README.md).

These tools provide an offline, permission-aware source inventory, immutable
publisher, Python API and query CLI. They preserve native source bytes and
true scientific gaps. They do not select a particle reaction network or turn
source/integrity checks into physical validation. Read the preimplementation
[contract](IMPLEMENTATION_CONTRACT.md) and [source policy](RIGHTS.md).

Publish into a **new, nonexistent** staging directory:

```bash
python3 -B utility/chemical_database/build_database.py \
  --amps-root "$PWD" \
  --validation-root /data/vtenishe/moon_validation_data \
  --output-root /tmp/NEW_CHEMISTRY_PUBLICATION \
  --release-id REVIEWED_RELEASE_ID
```

The native source inventory currently uses the explicitly retained research
packets documented in the database proposal. Missing packets are represented
explicitly, not fetched or replaced. `publication.json` retains an incomplete
state on ERROR. A finished convenience lock is written only after archive
verification, SDK and README creation succeed.

Query any completed staging or installed snapshot:

```bash
python3 -B utility/chemical_database/chemdb.py \
  --lock /path/to/publication/database-lock.json verify
python3 -B utility/chemical_database/chemdb.py \
  --lock /path/to/publication/database-lock.json list --kind tables
```

`python/amps_chemical_database` is the source package copied into the standalone
installation. Native source authentication, units/domains and scalar formulas
are shared C++ production code. Python table access retains source strings and
explicit unit guards; it does not reimplement the rate formulas independently.

The shared test runner is
[test/chemical_database/run_all_tests.py](../../test/chemical_database/README.md).
Its source-independent protocol is registered with the normal core test
routine. Data-dependent comparisons use an explicit real installation lock;
missing data remain SKIPPED. Existing Moon namespace wrappers link the same
shared numerical implementations, preserving their original source contracts.

No runtime command downloads data, mutates canonical originals, selects a
latest source, interpolates a gap or grants missing physical/permission scope.
External installation and its verification are recorded with the dated
development evidence and `srcMoon/STATUS.md`.

The November 19 acquisition extends the catalogue through explicit, offline
imports. Read [SOURCE_PACKET_CONTRACT.md](SOURCE_PACKET_CONTRACT.md) before
adding a provider. A reproducible acquisition plan names each HTTPS URL,
relative output path and maximum byte count; acquisition writes a fresh
directory with original responses and per-request outcomes:

```bash
/usr/bin/python3 -B utility/chemical_database/acquire_sources.py \
  --plan /path/to/explicit-acquisition-plan.json \
  --output-dir /path/to/NEW_ACQUISITION --timeout-seconds 25
/usr/bin/python3 -B utility/chemical_database/build_database.py \
  --amps-root "$PWD" --validation-root /data/vtenishe/moon_validation_data \
  --output-root /path/to/NEW_PUBLICATION --release-id NEW_RELEASE_ID \
  --source-packet /path/to/acquisition/source_packet.json
```

Repeat `--source-packet` for separate packets. The default source inventory is
unchanged when this argument is absent. Failed HTTP requests are acquisition
evidence, not scientific source objects. Each admitted packet pins original
bytes, primary rights evidence and any native-to-CSV recipe/parents; the release
also archives the exact packet JSON. CSV projections preserve numeric tokens,
signs, order and gaps. Parser tests use the real originals when explicitly
available; unavailable optional originals produce SKIPPED.

The new installation is an independent immutable snapshot at
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r2`.
Select its own `database-lock.json` explicitly. The original root lock and r1
remain available. The installed README contains real NASA CEA, IAEA, ExoMol and
TU Wien table examples and the physical-input gap matrix. This source-access
release enables no reaction model and changes no Moon production configuration.

The later November19
[M3 literature follow-up](../../srcMoon/test/M3_PHYSICS_LITERATURE_20261119.md)
identifies further measured inputs and documented community models. The
approved follow-up publishes CHIANTI O I, Apollo TPD, water-fragment workbooks,
JWST and ARTEMIS author products in the separate immutable r3 snapshot. Source
roles, native values, units, gaps and provider-specific rights remain explicit.
r1/r2 and application selections remain unchanged. The source atomic note also
corrects an earlier truncated-threshold interpretation using authenticated D04.

Read [PHYSICAL_MODELS_INTEGRATION_CONTRACT.md](PHYSICAL_MODELS_INTEGRATION_CONTRACT.md)
for the preimplementation contract. The optional publication flag
`--include-model-definitions` adds application-independent model descriptions
and their ordinary CSV table, without adding physical capabilities to the
profile. Python `models()`/`model(id)` and CLI `list --kind models` inspect this
metadata; C++ reads the same table through the existing Catalogue. New
`src/models/chemical_database/ReducedModels.h` helpers require explicit caller
parameters and preserve the old runtime and rate services. They implement
general mathematical/reduced kernels rather than selecting a Moon network.
Installed README section9 gives executable CLI/Python and shared C++ guidance.

The additional database-first/publication search is published separately as
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r4`.
Read [PHYSICAL_INPUT_STATUS_20261119_R4.md](PHYSICAL_INPUT_STATUS_20261119_R4.md)
for the new molecular radiative matrices, IAEA excited-Na composites, Apollo
15221 TPD, soil adsorption/migration and OH(A)+H2 source-data scope. Native
count/label gaps, measured/model distinctions and source-specific permissions
remain explicit. The new release corrects the MYTHOS transition-completeness
interpretation with matching-version native evidence; all earlier release
bytes and source/CSV values are preserved.

Select the r4 lock explicitly. Its installed README section10 gives executable
queries and examples; the same Python/CLI/C++ Catalogue paths serve every AMPS
application. Source readers `new_molecular_exomol.py`, `new_surface_*.py` and
their exact helper dependencies ship with release provenance and the installed
software. `amend_exomol_metadata.py` reproduces the version-bound metadata
correction without changing science payloads. No provider or reaction is
automatically selected, and no source-access PASS closes M3.

The [targeted follow-up](PHYSICAL_INPUT_STATUS_20261120.md) distinguishes
available physical primitives and derivable reference models from unresolved
source metadata or new model choices. Its separate r5 snapshot adds NIST H/O
evaluation cell views and the original Figshare sticking file, with exact
source-specific permission binding; additional literature/model originals with
uncleared payload rights remain metadata-only. Earlier releases and application
selections remain unchanged. Installed README section11 gives source/table
queries. The new `native_nist_hxoy.py` reader performs lossless cell extraction,
not thermochemical state selection or unit/uncertainty conversion.
