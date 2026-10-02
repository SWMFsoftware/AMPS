# Runtime background selection and publication

The runtime freezes application configuration, owns the integer clock and
controls which field generation may be consumed. Model implementations belong
to [`../background`](../background/README.md); native AMPS storage and MPI
operations belong to `main_lib.cpp`. A new model joins this boundary without
adding its own center-node write loop or halo protocol.

## Configuration and factory responsibilities

| Entry point | Contract |
|---|---|
| `ParseConfigurationText` / `LoadConfigurationFile` | Parse explicit model selections and canonically resolve SWCME assignments |
| `RunConfiguration3D::Create` | Freeze normalized options, storage layout and physics identity; check authority/ID pairing |
| `CreateBackgroundProvider` | Construct a built-in or registered candidate; validate before replacing the caller's output |
| `RegisterBackgroundModel` | Register a nonempty unique extension ID on each rank before acquisition |
| `BackgroundAuthorityMatches` | Check metadata provenance against the configured authority |
| `StandaloneAdapter` / `SwmfAdapter` | Publish a validated descriptor through the correct acquisition ownership path |

Schema 4 accepts `background.provider=swcme` for the canonical runtime, or
`background.provider=runtime-model` with a registered `background.model_id`.
Schema 3 retains its original standalone Parker contract. Prescribed turbulence
remains required for these standalone inputs. SWMF import still uses its typed
host. The reserved Python provider remains unimplemented.

SWCME construction re-resolves the retained canonical assignments and compares
both manifest and fingerprint with the frozen shock configuration. It never
creates a second CME from application defaults. A FULL_ICME mesh control must
use RESOLVED_COMPRESSION and source-free shock propagation, with Sphere or
finite SSE geometry. SSE must keep its tangent-flank ejecta outside the inner
ambient handoff. The native mover receives shape, normalized axis and half width
with the same epoch's apex radius/speed; see [SSE.md](../SSE.md).
Combining that compression with the existing DSA particle source is rejected.

Standalone request construction preflights runtime provider preparation at
launch before AMR allocation. This catches missing registrations, bad coverage
and incompatible metadata early. The actual owner-cell provider is still
created or installed during acquisition and validated against its real mesh
samples; preflight is not a field-publication test.

## Cadence transaction in AMPS

`RefreshBackgroundAtBoundary` runs only when the Runtime background event is
due, at `CurrentTimeS()` derived from the integer tick and configured step.

1. Prepare the runtime provider and build a temporary immutable owner-cell
   snapshot, or obtain the complete pair staged by a SWMF host.
2. Prepare/evaluate matching turbulence in scratch. Validate cell order,
   sample completeness, native buffer layout and temporal coverage.
3. Join readiness on every rank, including failures and empty owner ranks.
   Compare exact generation and a deterministic metadata hash, and require
   positive global physical-cell coverage.
4. Check continuity with the installed provider identity/configuration and an
   advancing epoch/generation. Request, begin and stage the Runtime descriptor;
   join stage validation before any live write.
5. Copy each validated sample and its waves to the application slice and
   allocated native DATAFILE fields. Both allocated native time slots receive
   the same epoch: fields are frozen between cadences.
6. Exchange associated data and finish enabled guiding-center auxiliary work,
   then commit the Runtime descriptor, swap the immutable active snapshot and
   increment the completed-refresh counter before the next particle phase.

Initial fill uses the same candidate validation and field mappings. Its initial
descriptor is published before the final halo operation; initialization output
and movers remain blocked until that operation completes and readiness is set.
The staged refresh transaction above is specifically the later-update order.

Candidate and staging failures abort with a typed diagnostic before live bytes
change. The application does not recover and continue after a write or MPI
failure. `Runtime` exposes failure acknowledgement for lifecycle tests/hosts;
that API does not imply an implemented native AMPS recovery campaign.

## Add a future model

Implement `BackgroundProvider` with frozen typed parameters, transactional
preparation, complete SI samples, deterministic identity and truthful capability
flags. Its resolved manifest must include source data, coordinate conversion,
derivative policies, coverage and any heating/interpolation assumptions. Publish
metadata with `ProviderKind::RuntimeModel` for a named extension.

Register a factory on **each MPI rank** before parsing/initialization; see the
[complete registration example](../../SWCME_MESH_BACKGROUND_CHANGES.md#add-another-source).
The registry is process-local and protected by a mutex. Factory invocation runs
outside that mutex, allowing dependency construction. Empty, duplicate,
reserved built-in names and null callbacks are rejected. A failed factory,
null successful result or failed validation leaves the caller's existing
provider unchanged.

Select the factory with:

```ini
[background]
provider = runtime-model
model_id = my-model
```

The model ID alone does not add custom text keys or a subprocess protocol.
Provide a typed parser/host configuration for extra model parameters. A typed
host may alternatively call `SEP3D::InstallBackgroundProvider` after immutable
configuration and before acquisition; the API retains shared ownership and
rejects a second owner. Acquisition still checks metadata against configured
authority/frame. The public `SEP3D.h` only forward-declares provider types so
canonical model headers do not enter unrelated PIC translation units.

## Restart and evidence

The default provider `RestorePreparedGeneration` hook rejects unsupported
restoration. Parker reconstructs the checkpoint epoch and restores its live
counter so the next update exceeds the active snapshot's generation. This hook
does not add SWCME checkpoint support: the source-free propagation control
currently requires a fresh run.

The completed runtime manifest records the installed mesh authority, cadence,
epoch and generation independently of the shock history. They can differ when
the front is sampled every tick but the mesh cadence is longer. Native capture
reads owner bytes and one physical representative per received remote block
without preparation, writes or an extra exchange. Its committed-update counter
is compared with due cadence events, rather than inferred from generation.
Readback covers the mapped primitive/transport/E fields when allocated; native
current and electron pressure are filled but not separately compared by these
gates. Optional application tensor slots follow the frozen `[storage]` layout.

See [`../test/README.md`](../test/README.md) for portable and native test
selection, prerequisite SKIPs and aggregate log paths. A configured enclosing
AMPS checkout is required for native MPI evidence.
