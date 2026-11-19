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
