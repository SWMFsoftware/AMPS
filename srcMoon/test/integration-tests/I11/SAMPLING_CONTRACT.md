# I11 sampling-window regression contract

This extends I11's existing sampled-observable wiring coverage. It assigns no
new U/I identifier and cannot establish spatial convergence on its own.

| Field | Contract |
|---|---|
| Purpose / capability | Verify independent completed cell-sampling windows when species visualization files are disabled. |
| Production entry points | Normal `amps_init`, `PIC::Sampling::Sampling`, `PIC::Mesh::switchSamplingBuffers`, `PIC::Mesh::flushCollectingSamplingBuffer`, completed-sample callback. |
| Configuration | SPICE Moon, approved 480 km child input; any recorded positive MPI rank count, OpenMP=1. No particle advances occur in this deterministic storage fixture. |
| Inputs / oracle | Add exactly representable, index-dependent markers (the eight-value pattern repeats) to one interior production cell per rank before each window, selected with GetCenterNode(i,j,k). The single-index zero entry is a ghost cell and must not be used. Three equal windows must each return the same markers, and the new collecting buffer must be exactly zero. The reference is the identity sum(constant)/one step, not a frozen model output. |
| Units / frame / time | Raw sampling-storage markers are dimensionless; no physical observation or force is supplied. The production mesh retains LSO coordinates and its configured frozen TDB epoch. |
| Controls / procedure | Record all compiled species; exercise all species file flags off, first species on, and all on. Reset the fixture's output index below the first visualization index to avoid expensive files while retaining callbacks. Runtime particle collection is disabled only in this storage fixture. |
| Metrics / acceptance | Exact zero marker error in each completed window; exact zero collecting-buffer magnitude; exactly one completed callback per window on each rank. Nine windows cover three flag layouts and both alternating buffers. |
| Expected PASS | Equal input windows stay equal for every compiled species and MPI owner, independently of output flags. |
| Failure modes | Stale alternating buffers, reset conditional on species output, erased completed data, missing/duplicate callback, absent allocated cell. |
| Data / hashes | Only the normal configured SPICE kernels; executable, configuration, source and input identities in the manifest. No Dxx observational package. |
| Predecessors | Build/configuration and approved input preflight; this fixture is an I11 infrastructure prerequisite. |
| Artifacts | Input snapshot, stdout/stderr, receipt, result, manifest and process-resource CSV. |
| Status | PASS means this exact storage identity was evaluated and met; FAIL means an evaluated identity failed; build/runtime/receipt failure is ERROR; unavailable configuration is SKIPPED. None is a formal I11 PASS. |

The live I11 path additionally records instantaneous weighted Na particle
inventory alongside the completed-window density integral. Their distinct
time meanings must be preserved: instantaneous inventory need not exactly
equal a window mean during spin-up. These diagnostics do not change the scoring
vector or its thresholds. Existing source/chemistry/boundary event counters
must be audited before claiming a complete physical budget; I19 still owns
production-wide event-budget qualification.

`particle_diagnostics.csv` is indexed in the per-level result and hashed by
the campaign manifest for both qualified and capped outcomes. A missing file
is infrastructure ERROR. This provenance requirement does not add the direct
inventory to the accepted scoring vector or change a tolerance.

For compiled Na+Ar storage coverage, `build_m2_qualification.py --purpose
sampling` retains the exact derived compile-time species/source settings and
executable checksum. `run_sampling_check.py --build-receipt <receipt>` uses
the already-documented M2 Ar source parameters so AMPS initializes a positive
Ar weight normally; it selects spherical geometry, frozen epoch and the
simple-cosine temperature service. No particle advances occur. The ordinary
M2 builder's default I19 chemistry-off check is unchanged and a sampling-only
build is explicitly not an I19 configuration qualification.

The four formal mesh levels may run concurrently in isolated working
directories, using the same rank count and ordinary I11 scorer. A shared I12
receipt may be captured at initialization and independently scored, avoiding a
second mesh allocation. The approved 500/5000-step contract and I13 review
requirements remain binding. Resource termination is ERROR and preserves
completed children; it never changes a science threshold.
