# AMPS shared chemical data runtime

This C++17 HOST library has no PIC, MPI, SPICE, application, network, JSON, SQL,
OpenSSL or HDF5 dependency. It owns strict source I/O, the immutable detached
runtime catalogue and the existing PHIDRATES/Voronov/total-capture evaluators.
It does not integrate particles, choose a bath/geometry/spectrum, activate
chemistry, infer daughter states or create permission/scientific approval.

`IMPLEMENTATION_CONTRACT.md` records the pre-change inventory and guards.
`RUNTIME_SCHEMA.md` specifies the exact block/CSV export consumed by C++.
The generic grammar header is reused through a relative include inside the
portable `src/models/chemical_database` + `src/general/input` source tree.

## Files and public interfaces

- `DataIO`: `ReactionData` SHA-256, strict CSV and explicit bounded reads;
  source-compatible with the historical Moon helper names.
- `Runtime`: `BundleRequest`, typed `Status`, `Catalogue`,
  `OpenRuntimeBundle`, source/table metadata and guarded payload access.
- `Photochemistry`, `ElectronImpact`, `ChargeExchange`: maintained production
  services mechanically migrated into `AMPS::ChemicalDatabase`, retaining
  public names, reviewed historical schema strings, hashes, formulas, source
  knots, endpoint conventions and invalid-domain behavior.
- `RateServices`: explicit-role, catalogue/permission-guarded preparation of
  those same services; no filename guessing or automatic source selection.
- `example.cpp`: a standalone catalogue-inspection executable, not AMPS.
- `ReducedModels.h`: explicit-input mass balance, radiative branching, Opal
  energy-marginal quantiles, two-body kinematics, unresolved-reservoir bookkeeping
  and finite-inventory thermal release. These new helpers do not alter the
  historical rate services. See `REDUCED_MODELS.md` for scope/units/guards;
  source selection, species constants and physical activation belong to callers.

The old `srcMoon-D04-phidrates-v1` and `srcMoon-total-capture-v1` schema strings
are preserved historical **data identities**, not application dependencies.
Changing a namespace must not silently authorize different bytes under them.
The parent application can keep old headers as namespace aliases and empty
compatibility translation units; only this shared source owns implementations.

## Open a pinned immutable bundle

```cpp
#include "Runtime.h"
namespace DB = AMPS::ChemicalDatabase;
DB::BundleRequest request;
request.root = configured_absolute_archive_root;
request.exportPath = configured_root_relative_export;
request.expectedExportSha256 = configured_export_digest;
request.expectedReleaseId = configured_release_id;
request.expectedReleaseSha256 = configured_final_release_digest;
request.expectedProfileId = configured_profile_id;
request.expectedProfileSha256 = configured_profile_digest;
request.requiredCapabilities = {"table-access"};
DB::Catalogue catalogue;
DB::Status status = DB::OpenRuntimeBundle(request, &catalogue);
if (!status.ok()) { /* owner records status.code/message and refuses use */ }
```

All pins above must come from the actual reviewed input/release. A JSON payload
is hash-authenticated but not parsed in C++. The detached export describes the
runtime meaning, and its own digest binds that description. The Python
exporter's independent reproducibility/schema checks establish its relation
to canonical data. Avoid mutually self-referential release/export digests.

`GetSource`/`GetTable` inspect metadata, including missing and restricted
records. `SourcePath`, `TablePath` and `TableRows` additionally enforce source,
table and profile guards. Unknown/pending/prohibited/not-applicable permission
is never treated as allowed. A profile scope only permits a lookup attempt;
it does not certify that each record is usable. Raw missing/negative data may
be inventoried with unavailable/quarantined capabilities, without exposing a
hazard evaluator. No physical zero stands for refusal.

## Prepare an explicitly selected rate service

The optional `scalar-rates` scope is separate from metadata/table inspection.
Pass exact table IDs in `PhotochemistryTables` or `TotalCaptureTables`, or an
explicit Voronov table ID. `RateServices.h` lists their schema/format roles.
The factories check all selected table and source permissions before invoking
the SAME original service kernels. PHIDRATES retains all17 catalogue rows,
including its two NON-reactive absorption totals, and15 reactive channels;
quiet and active endpoints are separate. Voronov remains an explicitly
Maxwellian coefficient service for neutral He/Na/Ne/Ar, with406 source records
validated and the same documented Te domain. Total capture retains its named
He source fit/Na interpolation and supplies no daughter kinematics.

Legacy direct service APIs remain source-compatible for applications that
already own frozen-source provenance. Such direct calls do not establish that
a new shared profile has permission/qualification; a release-selected runtime
must use the guarded preparation path. Neither path changes the original
rate formulas or installs a PIC reaction processor. Current metadata-only
sources cannot be promoted by selecting their species name or editing an
input hash. Numerical evaluation, kinetics and observational qualification
remain separate capabilities.

## Compile the standalone example

From the AMPS repository root, without configuring AMPS:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -Werror -pedantic \
  src/models/chemical_database/DataIO.cpp \
  src/models/chemical_database/Runtime.cpp \
  src/models/chemical_database/example.cpp -o /tmp/chemdb-inspect
/tmp/chemdb-inspect ROOT EXPORT_PATH EXPORT_SHA RELEASE_ID RELEASE_SHA \
  PROFILE_ID PROFILE_SHA catalogue
```

Supply real paths/pins; symbolic names in this command are not data/reference
values. The tool prints the typed open outcome and source/table inventory.
Other services add their corresponding `.cpp` plus `DataIO.cpp`; the guarded
factories also link `Runtime.cpp` and `RateServices.cpp`. Normal AMPS build
integration belongs to the parent-owned `src/general/makefile`: use unique
`ChemDb*.o` archive members and explicit source/dependency rules. There is no
separate root CMake project and generated copies are not maintained sources.
GPU device queries are not supplied; host/CUDA compiler integration must be
verified separately before claiming that configuration supported.

## Verification scope

The parent shared runner exercises published SHA vectors, independent hashlib,
strict CSV/schema/refusal fixtures, source/permission laundering, path/symlink
containment, transactional opens, immutable lifetime/concurrent access, and
old/new numerical services. Analytical infrastructure fixtures are never
observations. Application adapters must retain actual compiler/executable/
source/data identities and run the same previous core/Moon/comet gates.
A source-only smoke cannot establish linked or observational PASS. There is
no new permanent test ID, threshold relaxation, reference update or canonical
source mutation in this component.

## Explicit reservoir products

`MaxwellianVoronovTable::LogTemperatureDerivative` returns the analytic
dimensionless slope of the loaded fit on its native temperature domain.
The slope can supply a fit-derived incident mean in an explicitly admitted
isotropic Maxwellian reservoir. Source state, physical cross-section and
model-error evidence still control that admission.

[ReservoirProducts.h](ReservoirProducts.h) adds reusable, header-only C++17
services for the approved limited heavy-ion/reservoir models; its full API,
units, evidence and refusal contract is in [ReservoirProducts.md](ReservoirProducts.md).
It supports the collision-biased Maxwellian mean from an explicitly supplied
rate derivative, inherited-velocity ionization with signed electron-reservoir
moments, and Q-corrected forward partial charge exchange with separately tracked
neutral/plasma isotope carriers. It selects no source, bath, energy frame,
differential law or physical error budget. Transactional refusals preserve all
outputs. The registered core chemical-database suite runs its independent
analytic/conservation/refusal checks. These source-tree additions do not modify
the immutable installed r3 snapshot or activate application transport.
