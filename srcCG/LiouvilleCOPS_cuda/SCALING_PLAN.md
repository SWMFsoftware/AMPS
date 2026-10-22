# Scaling plan: SHAP8 20M and 44M surface models on a PCIe GPU cluster

Status: **decisions of 2026-10-21 in §0a; P1, P2 (incl. P2b: restart/fine-tuning kernels) and P3 COMPLETE (STATUS.md C15); P3 deferrals in §9 (R4); P4 next** (2026-10-21, user request: "prepare the model to run with a 20M surface
model, and then 40M ... implement pipelining (non-blocking communication with the host + concurrent
processing of what is in GPU memory); the model will be pushed to a GPU cluster with multiple CPUs per
node (no NVLink, only PCIe); make the necessary modifications in data management and load balancing;
prepare a detailed implementation plan and suggestions before starting").

Everything here stays inside LiouvilleCOPS_cuda. `../LiouvilleCOPS` (fit_cops.py, cops_io.py, the C++
model, the CPU file formats and workflow) is not changed.

## 0. Summary

* **What breaks at 20M/44M.** The kernel records (about 350 per facet) no longer fit in GPU memory even as
  FP32: 140 GB (20M) and 308 GB (44M), against 87 GB for 12M on two 48 GB GPUs. Five other parts also stop
  scaling:
  * the per-face arrays, 25+ basis columns per GPU (10 GB at 44M);
  * the host-side activity field `Y a`, recomputed on one CPU thread every evaluation;
  * the single-threaded loader, about 3e6 records/s per GPU;
  * the single-GPU tools: masks, `StateEvaluator`, the base builder;
  * disk space: the 44M FP64 kernel files take about 617 GB.
* **Design in one sentence.** Each GPU owns a contiguous range of 4096-record *blocks*. It keeps as many of
  them resident as fit, and streams the rest every pass through pinned triple buffers on a copy stream,
  while the compute stream works on what is already on the device. Ranks exchange only parameters and
  per-point sums, which are O(1 MB) per evaluation and independent of the mesh. A cost model calibrated
  per GPU (compute and measured PCIe bandwidth) balances the blocks.
* **Exactness is kept by construction.** Every per-block partial sum depends only on the block's records
  and is reduced in block order. Streamed, resident, one GPU, N GPUs and N nodes therefore give
  **bitwise identical** pressures and Jacobians for the same block table. This extends the current
  "two GPUs = one GPU" guarantee and is testable on the local machine.
* **Resident distribution is the fastest; streaming is GPU-hour-efficient.** 44M takes ≈ 1.9 h on 8
  resident GPUs and ≈ 3.9 h on 4 streaming GPUs (PCIe gen4; section 6). That is about the same
  GPU-hours, because the compute-bound Jacobian passes hide most of the transfer. Both run on the same
  code path: a rank with enough memory simply streams nothing.
* **Order (after the decisions of §0a):**
  * P1–P4 (store with face compaction, streaming engine, MPI distribution, load balancing) are built and
    validated on the local 2 × A6000.
  * Then P5a (several GPUs in one base-build process) and **20M on this machine**.
  * Voyager bring-up (P0) whenever an account exists. One 8-A100 node holds 20M fully resident (≈ 17 min
    per fit) and 44M with little or no streaming.
  * **44M waits for the 20M results.**

## 0a. Decisions (user, 2026-10-21)

| # | Question | Decision |
|---|---|---|
| Q1 | cluster hardware | **8 A100 per node, PCIe, no NVLink**; the cluster is UAH's Voyager (RCHAU); details in §1a |
| Q3 | activity field | **split the calculation between the GPUs and combine**: each GPU forms the activity of its own faces from its basis columns (§3.5); rank 0 combines the per-point sums |
| Q4 | FP32 record store; validate once | **yes** (§3.1, §3.6) |
| Q5 | 44M | **wait for the 20M results** |
| Q6 | where to start | **on this machine** (2 × A6000); disk: §1b |
| Q7 | base build | **start with one node** (several GPUs in one process); **plan multi-node with MPI** between nodes (§3.7) |
| new | load only the faces that can contribute | **mandatory face compaction** (§3.5): a GPU holds per-face data only for faces that occur in its records |

## 1. Sizes (measured for 12M, projected for 20M/44M)

Records per facet measured: 788k 311.7, 1M 298.5, 3M 272.4, SHAP8 12M 350.6. The projections use 350
(±15%). Throughputs are measured on one RTX A6000 in the 12M fits: pressures 2.5e9 records/s, Jacobian
0.71e9 records/s (three tile passes), mixed gravity 7.2e10 pairs/s. Pinned H2D is 12.3 GB/s per GPU here
(PCIe gen3 x16, no contention with two GPUs at once).

| | 12M (measured) | 20M | 44M |
|---|---|---|---|
| facets / genus | 12,364,290 / 4 | 19,999,504 / 4 | 43,966,474 / 6 (repaired BDF; the original has 9 zero-area faces) |
| kernel records | 4.33e9 | ≈ 7.0e9 | ≈ 15.4e9 |
| CPU-format kernel files (40 B/record) | 173 GB | ≈ 280 GB | ≈ 617 GB |
| FP32 device records (20 B/record) | 87 GB | ≈ 140 GB | ≈ 308 GB |
| per-face device arrays per GPU, uncompacted (≈ 235 B/face) | 2.9 GB | 4.7 GB | 10.3 GB |
| 48 GB GPUs for fully resident records | 2 (headroom 0.06) | 4 | 9 (8 with face compaction, §3.5) |
| 80 GB GPUs for fully resident records | 2 | 2 | 5 |
| mixed gravity (N² pairs, GPU time) | 1.5e14, 2,130 GPU-s | 4.0e14, ≈ 5,600 GPU-s | 1.9e15, ≈ 27,000 GPU-s |
| gravity tables per GPU (240 + 80 B/face) | 4.0 GB | 6.4 GB | 14.1 GB |
| GPU base build (∝ facets) | 3,906 s, 1 GPU | ≈ 6,300 GPU-s | ≈ 13,900 GPU-s |
| one fit, 2004 P + ≈1900 J evaluations (3.48e-6 GPU-s/record) | 15,100 GPU-s | ≈ 24,400 GPU-s | ≈ 53,600 GPU-s |
| fit outputs written by fit_cops (text, maps) | ≈ 7 min | ≈ 12 min | ≈ 25 min |

On this machine, /data has 239 GB free. A 20M base (280 GB) does not fit here without freeing space, so it
has to be built on the cluster (question Q6).

## 1a. Target cluster: Voyager (UAH, RCHAU)

From public sources:
* Voyager is the main system of the Regional Computing Center for Alabama Universities (RCHAU), hosted at
  UAH and funded by an NSF grant.
* About 3000 CPU cores and 24 NVIDIA Ampere GPUs (UAH announcement, January 2024). The grant description
  gives 64–128 cores per node and 240–360 TFLOPS theoretical FP64 for the GPU subsystem.
* Batch jobs run under Slurm; MPI tools, CUDA tools and Python environments are installed (rchau.uah.edu).

From the user: 8 A100 per node, connected by PCIe only. With 24 GPUs in total, that means **three
8-GPU nodes**: up to 24 A100 for one job (multi-node), 8 within one node.

Not public; P0 measures these with `tools/cluster/probe_node.sh` (`nvidia-smi topo -m`, `lscpu`,
`numactl -H`, `free`, `lsblk`, `ibstat`, a concurrent H2D bandwidth test):
* A100 memory (40 or 80 GB);
* CPU sockets per node, and which GPUs hang off which socket and PCIe switch;
* host RAM; node-local NVMe; inter-node network (likely InfiniBand); Slurm partitions, walltime, quotas.

What A100 changes:
* **FP64:** 9.7 TFLOPS against 1.25 on the A6000 (7.8×). Memory bandwidth: 1.6–1.9 TB/s against 0.77 TB/s
  (2–2.5×). The Jacobian kernels gather basis columns per record, so they are probably partly
  memory-bound. The projections below assume **3× (range 2–6×)** the A6000 rates; P0 measures the real
  factor.
* **PCIe:** eight GPUs per node share the host's PCIe lanes, often two to four per switch. The
  concurrent H2D bandwidth per GPU can be well below the 25 GB/s of a single gen4 link. The calibration
  (§3.4) measures it with all GPUs copying at once.

## 1b. Disk on this machine (2026-10-21)

| Filesystem | Free | Use |
|---|---|---|
| `/data` | 220 GB | bases and products so far |
| `/` (root) | 322 GB | system |
| `/home` | **2.7 GB (100% used)** | the code and its `.venv`; writes there may fail |

A local 20M run needs about 280 GB of CPU-format kernel files plus about 140 GB of FP32 store, more than
either filesystem holds alone (question R1, §9).

## 2. Current architecture: what limits scaling

| Component | Now | Limit at 20M/44M |
|---|---|---|
| `device_store.DeviceDataset.load` | all records of the assigned rows resident; one Python thread per GPU reads, validates, hashes and converts | memory; ≈ 3e6 records/s per GPU (12M bind 715 s; 44M on 8 GPUs ≈ 11 min per bind, on 4 ≈ 22 min) |
| `evaluator.CudaBackend` | per-block partial sums, reduced per point in block order (`first_block`); Jacobian in tiles of ≤ 32 sums, **each tile a full pass over the records** (NC31: 3 passes) | tile-major order would re-stream every tile |
| `multi_gpu.MultiGpuBackend` | spawned worker per GPU, pipes, one host `act_lin` in shared memory; contiguous **row** blocks balanced by records | one node only; row granularity; `Y a` on one host thread (44M: 8.8 GB read per evaluation, ≈ 1 s, about 1 h per fit with all GPUs idle) |
| `gpu_model._Registry` | a new basis (template fit) creates a new backend and **reloads all records** | 12M: 2 × 12 min binds per template fit |
| `base_validation` | rank 0 reads and hashes all kernel files (≈ 3 GB/s locally) | 617 GB from a parallel file system on every fit |
| `masks.py`, `restart.StateEvaluator`, `restart_tool`, `c05_write_bundle` | one GPU, all records resident | cannot hold 20M (`fine_tuning.SplitEvaluator` is row-split but single-node) |
| `native/liouville_cops_cuda` | one GPU per run | 44M ≈ 3.9 h; the C++ main supports `FirstPoint`/`LastPoint` |
| `c04_gravity_equivalence --devices` | target blocks per GPU, one node | fine for memory; multi-node and an FP64 spot check needed |

## 3. Design

### 3.1 Data: a device-ready record store (derived product)

* New tool `tools/c15_record_store.py`. It converts a validated base, once, into an FP32
  structure-of-arrays store with the layout of the device: `face` int32, `kn`, `kphi`, `kv`, `cos`
  float32, shadow folded into cos.
* The conversion applies the checks the FP32 loader already makes: zero-ness, finiteness and sign kept,
  values FP32 cannot hold refused, and the one-ulp COS75 adjustment.
* The store is split into **chunks aligned with blocks and points**. Target 256–512 MB per chunk; one
  point is about 270 MB at 44M, and a large point is split at block boundaries.
* `store_manifest.json` records:
  * the kernel-set digest of the source base, from `run_manifest.json`;
  * the conversion audit;
  * one SHA-256 per chunk;
  * the block table.
* A fit checks the manifest against the base it validates, so a store can never be paired with another
  base.
* **Effect:**
  * half the I/O (20 instead of 40 B/record);
  * no per-bind Python conversion or re-validation; chunk digests are verified as the chunks are read;
  * any rank can read exactly its block range;
  * a job can stage the store to node-local NVMe.
* The CPU-format kernel files stay the reference. The store is a cache derived from them: it can be
  deleted and rebuilt, and it is never read by the CPU tools.

### 3.2 Streaming engine (one GPU): resident set plus a pipelined stream

* New module `cops_cuda/streaming.py` with a `ChunkEngine`, used by `CudaBackend` when its memory plan
  does not fit.
  * **Resident set:** the first chunks of the rank's range, up to the budget left after the per-face
    arrays and the stream buffers.
  * **Stream set:** the remaining chunks, held in **pinned host memory** (`cudaHostAlloc`, allocated on
    the GPU's NUMA node, §3.6), or read from node-local disk into a pinned ring when host RAM is short.
  * **Pipeline:** K = 3 device slots of one chunk each; a copy stream (`cudaMemcpyAsync` H2D) and a
    compute stream, ordered by CUDA events:

    ```
    copy  : [H2D c0][H2D c1][H2D c2]  wait(slot0 free) [H2D c3] ...
    comp  : [resident chunks ........][c0 all tiles][c1 all tiles][c2 ...]
    ```

    The resident chunks are processed first, so the first copies are hidden behind them.
* **Chunk-major Jacobian:** every requested tile is computed for a chunk while it is on the device, so a
  Jacobian costs one transfer, not three.
* **Results:** per-block partials go into the block-indexed `partial` array exactly as now. The per-point
  reduction in block order is unchanged, so **streamed equals resident bitwise**.
* **Overlap rules:**
  * a pass takes max(compute of all records, streamed bytes / bandwidth);
  * pressure passes are transfer-bound, about 2.5e9 records/s compute against 0.6e9 (gen3) or 1.2e9
    (gen4) records/s of transfer;
  * Jacobian passes are compute-bound at 0.71e9 records/s, so their streaming is largely hidden.
* **Suggested option S5:** when the rank is transfer-bound, compute the Jacobian columns during each
  pressure pass and cache them under the parameter vector. TRF asks for J at points where it has
  already evaluated f, which saves one transfer per accepted step. It is optional and measured first.
* **Masks:** `seen_positive`, `visibility_union` and the restart kernels (`facet_point_shares`,
  `point_change_partials`, `first_*`) run through the same engine. They are OR, min and max reductions or
  block partials, all independent of order.

### 3.3 Distribution: one rank per GPU, across nodes

* **Transport abstraction** (`cops_cuda/transport.py`):
  * `LocalTransport`: the current spawn, pipes and shared memory, kept for single-node runs and the
    tests;
  * `MpiTransport`: mpi4py, one rank per GPU, launched with `srun` or `mpirun`.
* Rank 0 runs the unchanged `fit_cops.main()` through `gpu_fit.py`; ranks ≥ 1 run the worker loop.
  `MultiGpuBackend` keeps its interface: `evaluate`, `seen_faces`, `identity`, `close`.
* **Per evaluation:**
  * rank 0 broadcasts the parameter state (≈ 300 B) and the activity field: option A or B, §3.5;
  * each rank returns its per-point sums, at most nsums × npt × 8 B ≈ 0.6 MB.
  * No GPU-to-GPU or device-to-device traffic, so NVLink is irrelevant; PCIe carries only the
    record stream and small arrays.
* **Block-granular partition with exact boundaries:**
  * Ranks own block ranges, not rows. At most one point per rank boundary is split.
  * For those boundary points (≤ nranks − 1), the ranks return their **block partials**, a few hundred
    blocks. Rank 0 reduces them in block order.
  * Result: bitwise identical to one GPU for any partition, and much finer balance than rows. A 44M
    point is about 13e6 records, about 1/1150 of the work, which is too coarse for 16 ranks.
* **Identity and failure handling** as now:
  * every rank reports its block range, store chunk digests and device identity; rank 0 records them
    in `fit_manifest.json` (`compute_backend`);
  * a rank failure or timeout raises `BackendError` on rank 0, which ends the job with
    `MPI_Abort`; nothing is retried silently;
  * timeouts scale with the records (as the bind timeout does now).
* **Template fits load once:** the `_Registry` keeps one record store per row set. A basis change
  re-uploads only the basis columns, not the records (suggestion S2 below), which saves about 23 min per
  template fit at 12M.

### 3.4 Load balancing

* **Calibration at start-up, about 10 s per GPU.** Each rank measures on its own data:
  * records/s for the pressure and the full Jacobian kernels;
  * pinned H2D bandwidth with **all GPUs of the node copying at once**, which captures shared PCIe
    switches and root complexes;
  * free device memory.
* **Cost model:**
  * resident records cost r/C (C = compute rate);
  * streamed records cost max(r/C, 20 B·r / BW).
  * The planner chooses each rank's block range and resident share to minimize the maximum predicted
    pass time, weighted 2004:1900 for pressure:Jacobian passes, subject to device memory and host pinned
    memory per node. Heterogeneous GPUs (for example a mix of A100 and A40) are handled by their measured
    C and BW.
* **One rebalance:** after the first about 20 evaluations, the planner compares the measured pass times
  with the model. If the imbalance is above 10%, the ranges are moved once (records reloaded from the
  store) before the fit proceeds. The final block table is recorded.
* **Within a node, the streamed chunks can be dynamic** (optional, phase P4b):
  * the node's streamed chunks sit in one portable pinned pool (`cudaHostAllocPortable`);
  * any GPU of the node takes the next chunk from a node-local queue;
  * results land in the block-indexed partials, so the assignment does not change the bits;
  * this absorbs PCIe contention that varies with time. The resident parts stay static.
* **NUMA and PCIe placement:**
  * each rank is bound to the cores of the socket its GPU hangs off (`--gpu-bind=closest`, `--cpu-bind`,
    hwloc), and allocates its pinned buffers there by first touch from a bound thread;
  * crossing UPI or Infinity Fabric for H2D can halve the bandwidth;
  * `nvidia-smi topo -m` and the measured bandwidths go into the manifest.

### 3.5 Per-face data: only the faces that can contribute; the activity field split over the GPUs

**What can contribute.**
* A facet f contributes to point p exactly when (p, f) has a kernel record. The base builder already
  stores only those: `k_compact` keeps the facets with a nonzero Kn, KPhi or KV.
* The *records* therefore cannot be thinned further without changing the pressures.
* The *per-face* arrays can. Today every GPU holds them for all faces:
  * the clipped activity `af` and `active`;
  * ACO2;
  * 25 harmonic columns and the template columns.
* A face that occurs in none of a GPU's records is never read by that GPU.

**Face compaction (mandatory, user decision).**
1. **Global:** the faces that occur in any record of the base (the "seen" set).
2. **Per rank:** the faces that occur in the rank's own blocks.
   * The store (§3.1) writes, per chunk, the sorted unique face list.
   * At bind the rank forms the union for its blocks and a map `global → local`: a bitmap over nfaces,
     a prefix sum, and an int32 map that is used once and dropped.
   * It rewrites its records' `face` field to local indices (int32 on the device). Locality is kept,
     because faces stay in ascending order within a point.
3. **Uploads per rank:**
   * the basis columns (Y and the templates) for its local faces only, read by the rank from
     `elements.bin` and the bundle;
   * ACO2 for its local faces only.
   * Per-face device memory scales with the rank's local face count, not with nfaces.
4. **Outputs keep global indices:**
   * `seen_faces`: each rank returns its local list, mapped back; rank 0 forms the union;
   * the restart/mask kernels' per-face results are scattered back to global indices on the host.
5. **Exactness:** compaction renames faces and does not reorder records or sums, so every invariant of
   §4 holds bitwise.

**Measured on 12M** (all 4.33e9 records; contiguous row blocks balanced by records, as now):

| Faces held | Faces | Share of 12,364,290 |
|---|---|---|
| any record of the base (global compaction) | 9,212,140 | 74.5% |
| largest rank, 2 ranks | 8,896,998 | 72% |
| largest rank, 4 ranks | 7,352,083 | 59% |
| largest rank, 8 ranks | 6,387,718 | 52% |
| largest rank, 16 ranks | 5,907,581 | 48% |

* With 8 ranks the per-face device memory falls to about half: 44M ≈ 5.4 GB instead of 10.3 GB per GPU.
* The records stay the bulk (308 GB at 44M). They are already restricted to contributing (point, face)
  pairs, so compaction frees per-face memory but does not shrink them.
* Row blocks ordered by time are already spatially coherent along the trajectory. A spatial ordering of
  rows could compact further, but it changes the block table and costs balance; it is not planned.

**Activity field: split between the GPUs, then combined (user decision Q3).**
* Each rank forms `act_lin = Y_local a + A_scale B_local b` for its local faces with a fixed-order device
  kernel. The 25 + nb terms are added in the order of the host expression, so the result is deterministic
  and independent of the partition.
* `af = max(act_lin, 0)` and `active = act_lin > 0` follow on the device.
* No per-face data crosses the network; rank 0 broadcasts only the parameters (about 300 B).
* The per-point sums are combined on rank 0 as now (§3.3).
* A face seen by several ranks is computed identically by each, since the inputs and the order are
  equal.
* **Check:** at the first evaluation, and every 100th, rank 0 compares the device `act_lin` of a sample
  of faces with the host NumPy expression of `Model.evaluate`. It records the maximum ulp difference and
  the count of active-set flips (expected: none, unless a face sits exactly at zero) in
  `fit_manifest.json`.
* The old exact host path stays available (`--activity host`) for verification runs. In it, rank 0
  computes `act_lin` and scatters each rank its compacted faces.
* `fit_cops` itself still computes the final activity on the host for its outputs, unchanged.

**Option S7 (later), Y on the fly:** degree-4 harmonics from the facet direction (25 values by
recurrence) would remove the harmonic columns from the GPUs. That changes rounding against
`elements.bin` Y, so it waits until the compacted columns are measured.

### 3.6 Host side: threads, NUMA, memory

* **Loader with several CPU threads per GPU** ("multiple CPUs per node"):
  * per rank, a pool of reader threads (NVMe or parallel file system to pinned memory) and a hash/verify
    thread;
  * with the FP32 store the bind becomes I/O-bound; on NVMe, 4–8 GB/s per node replaces 120 MB/s per GPU.
  * Without the store (fallback), the current conversion runs in a thread pool; NumPy releases the GIL.
* **Host RAM budget per node:** pinned stream pool, plus about 2 × 512 MB ring per GPU, plus rank 0's
  host model. Rank 0 holds `elements.bin` (12 GB at 44M), Y (8.8 GB), the per-face outputs and fit_cops
  state: about 60–80 GB. Put rank 0 on a node with enough RAM; the planner knows its share.
* **Base validation (the reference audit at fit start):**
  * each rank validates and hashes its own kernel files in parallel;
  * rank 0 runs the reference `validate_fit_base` with workers that return those precomputed results.
    This is the mechanism of `base_validation.py`, now distributed, and the result is identical.
  * **Validate once (decision Q4):**
    * the first fit on a base runs the full distributed audit and writes `base_audit.json` beside the
      store: the reference audit result, per-file sizes and SHA-256, and the kernel-set digest;
    * later fits check, per file, size and mtime plus `run_manifest.json` against it, then reuse the
      audit;
    * any mismatch triggers a full re-validation;
    * `--revalidate` forces one.

### 3.7 The other stages at scale

* **Base build (C11), decision Q7: one node first, multi-node with MPI planned.**
  * **P5a, one node:** several GPUs inside one `liouville_cops_cuda` process.
    * One host thread per GPU takes the next point from a shared counter.
    * Each thread writes its kernel file through the unchanged C++ writer; the table row is kept in
      memory and written in point order.
    * Same files and manifest as a one-GPU run.
    * First, a check that per-point results do not depend on earlier points (RNG state, caches). Gate:
      the 788k build is bitwise the one-GPU build.
    * 20M: ≈ 6,300 s on one A6000, so ≈ 53 min on the two local A6000. On 8 A100, ≈ 3–7 min if an
      A100 is 2–6× an A6000 for this ray casting.
  * **P5b, several nodes:** an MPI driver (`tools/c11_build_mpi.py`).
    * Rank r runs the native builder with `FirstPoint`/`LastPoint` for its point range (a generated
      input differing only in those keys) into a rank directory.
    * Rank 0 then merges into one base:
      * kernel files moved, not rewritten;
      * `COPS_model.dat` rows concatenated in point order;
      * `run_manifest.json` written with the reference's own digest definitions (kernel set, table,
        elements).
    * The merge is accepted only if the CPU's `validate_fit_base` passes on it and, on 788k, the merged
      base is bitwise the one-process base (kernel files and table rows).
* **Gravity (C04):**
  * the target split already exists; extend it to several nodes with MPI (one target range per rank, then
    the existing concatenation);
  * the field is needed at every facet (cliff templates and maps cover the whole surface), so gravity
    targets are not compacted; the sources are all facets by definition;
  * memory per GPU is 14 GB at 44M;
  * add an **FP64 spot check** on a target subset, so mixed precision is validated on each new mesh:
    random targets, plus late-window, cliff-boundary and near-threshold facets. 1e4 targets at 44M
    ≈ 4e11 FP64 pairs ≈ 1 min.
  * Admission: `--allow-genus 4` (20M), `--allow-genus 6` (44M repaired).
* **Bundles, masks, fine tuning, restart tool:**
  * all move to the distributed and streaming backend (masks through the ChunkEngine, §3.2);
  * `SplitEvaluator` and `StateEvaluator` become thin wrappers over `MultiGpuBackend`;
  * `restart_tool enhance|translate` gain `--devices` and MPI.
* **Fit outputs:** `fit_cops` writes per-face text files and maps on rank 0, unchanged. The restart file
  grows to about 1.6 GB at 44M.

## 4. Correctness and validation

Invariants, each with a test:

| # | Invariant | Test (local 2×A6000 unless noted) |
|---|---|---|
| I1 | streamed = resident, bitwise (same device, same blocks) | synthetic and 788k with a forced tiny budget (`COPS_CUDA_MAX_RESIDENT_BYTES`); every pressure/Jacobian column |
| I2 | any block partition over 1..N ranks = one GPU, bitwise, including split boundary points | `mpirun -n 2/3/4` on the 2 local GPUs (several ranks per GPU allowed in tests); random block ranges |
| I3 | store-loaded = kernel-file-loaded device arrays, bitwise | device digests (`verify_roundtrip`) |
| I4 | dynamic intra-node chunk assignment = static, bitwise | randomized chunk order |
| I5 | distributed base validation = reference audit, identical | as `test_base_validation.py`, with ranks |
| I6 | multi-GPU native base build = one-GPU build, bitwise | 788k (5 min) |
| I7 | option B activity field within a few ulp of the host expression; active set identical except exact-zero cases | synthetic and 12M |
| I8 | failure handling: killed rank or truncated store chunk gives a clean `BackendError`, nothing partial | fault injection |

Production gates, in order:
1. **On the cluster:** the test suite, then 788k B29 FP64 bitwise equal to the local result. CUDA FP64
   with `--fmad=false`, correctly rounded division and sqrt is expected to be bitwise portable across
   architectures; this gate confirms it.
2. **12M mixed:** B29 rerun distributed and streamed, bitwise equal to today's fit (same block table).
3. **20M:** base, gravity (with the FP64 spot check), bundles, B29/N30/C30/NC31; comparison package with
   12M, 3M and 1M.
4. **44M:** after decision Q5.

## 5. Phases (order after the decisions of §0a)

| Phase | Content | Where | Relative size |
|---|---|---|---|
| P1 | FP32 record store with per-chunk face lists (`c15_record_store.py`); store loader with reader threads; **face compaction**; registry fix (one record load per row set, basis-only rebind); validate-once audit | local (788k, 3M, 12M) | M |
| P2 | `ChunkEngine`: resident plus streamed, pinned triple buffer, copy/compute streams, chunk-major tiles; masks and restart kernels on it; I1, I3 | local | L |
| P3 | `transport.py`: Local and MPI (mpi4py on the local Open MPI 4.0.3 for tests; the cluster's MPI); block-granular partition with exact boundary reduction; activity split over the GPUs (§3.5); distributed base validation; I2, I5, I7, I8 | local (`mpirun -n 2..4`) | L |
| P4 | Calibration, cost model, planner, one rebalance, NUMA/PCIe binding; P4b dynamic streamed chunks (I4) | local, then cluster | M |
| P5a | Several GPUs in one native base-build process (one node); I6 | local (788k gate), then 20M | M |
| P6 | **20M on this machine:** base (P5a, 2 GPUs), mixed gravity (`--allow-genus 4`, FP64 spot check), bundles, store, B29/N30/C30/NC31 streamed on 2 A6000 (≈ 4 h per fit); comparison package with 12M, 3M, 1M | local | — |
| P0 | Voyager bring-up: probe script, uv venv, CuPy for the cluster CUDA, mpi4py on the cluster MPI, `native` (ARCH=sm_80), Slurm templates (`tools/cluster/`), A100 benchmark; then the gates of §4 | cluster | S |
| P5b | Multi-node base build (MPI point ranges and merge); multi-node gravity (MPI) | cluster | M |
| P7 | Decide on 44M from the 20M results (Q5) | — | — |

Each phase ends with a full test suite, a STATUS.md section and an unchanged CPU workflow. P0 can run as
soon as an account on Voyager exists. It needs nothing beyond P1–P3 to be useful, and runs the 20M fits
in minutes instead of hours.

## 6. Projected fit times (one fit, 2004 P + ≈ 1900 J evaluations)

**On Voyager (A100, 3× the A6000 rates, range 2–6×; PCIe gen4 shared by 8 GPUs, ≈ 12 GB/s each when all
copy at once — to be measured):**

| Case | GPUs | Streamed per GPU per pass | P pass | J pass | Fit |
|---|---|---|---|---|---|
| 20M, one node | 8 × A100 (40 or 80 GB) | 0 | 0.12 s | 0.42 s | ≈ 17 min (8–25) |
| 44M, one node | 8 × A100 80 GB | 0 | 0.26 s | 0.92 s | ≈ 38 min |
| 44M, one node | 8 × A100 40 GB | ≈ 0.3e9 (5.5 GB) | 0.46 s | 0.92 s | ≈ 45 min |
| 44M, three nodes | 24 × A100 | 0 | 0.09 s | 0.31 s | ≈ 13 min |

Add per fit: bind from the store (minutes), and `fit_cops` writing its outputs on rank 0 (12–25 min). At
these speeds the outputs become a large share of a run.

**On the local 2 × A6000 and generic 48 GB GPUs (measured A6000 rates):**

From the measured A6000 rates. A100/H100 have 8–25× the FP64 rate of an A6000, which would make the
streamed cases transfer-bound and the resident cases several times faster (question Q1).

| Case | GPUs | Records streamed per GPU per pass | P pass | J pass | Fit |
|---|---|---|---|---|---|
| 20M resident | 4 × 48 GB | 0 | 0.70 s | 2.5 s | ≈ 1.7 h |
| 20M, local 2 × A6000, streamed (gen3 12.3 GB/s) | 2 | 1.6e9 (32 GB) | 2.6 s | 4.9 s | ≈ 4.0 h |
| 44M resident | 8–9 × 48 GB | 0 | 0.75 s | 2.7 s | ≈ 1.9 h |
| 44M resident | 16 × 48 GB | 0 | 0.39 s | 1.4 s | ≈ 1.0 h |
| 44M, one node 4 × 48 GB, streamed, gen4 (≈ 25 GB/s) | 4 | 2.2e9 (44 GB) | 1.8 s | 5.4 s | ≈ 3.9 h |
| 44M, same, gen3 | 4 | 2.2e9 (44 GB) | 3.6 s | 5.4 s | ≈ 4.9 h |

The streamed cases assume no face compaction (§3.5); compaction leaves more room for resident records.
The 44M single-node case needs about 176 GB of pinned host memory on the node (4 × 44 GB).

Add per fit:
* bind: with the store on NVMe, minutes;
* fit_cops outputs: 12–25 min;
* the start-up validation: §3.6.

## 7. Risks

* **The 44M source** has 9 zero-area triangles. Only the repaired BDF passes the C++ check. That repair
  (18 faces changed, no vertices moved; `../LiouvilleCOPS/tools/repair_collinear_ply_bdf.py`) is a
  scientific choice (Q5).
* **FP32 underflow.** On smaller facets, Kn values could fall below what FP32 holds. The loader and
  store refuse such values loudly; if that happens, a per-point scale factor is the remedy.
* **Records per facet may grow** with roughness (SHAP8 12M has 350 against SHAP5's 272–312). The sizes
  in §1 carry ±15%.
* **Determinism on the cluster:** the cluster's driver, CUDA and CuPy versions may differ from the local
  ones; gate 1 checks this. If not bitwise, the tolerance policy of VALIDATION.md applies.
* **Pinned memory limits** (ulimit, cgroup, `vm.max_map_count`) on cluster nodes; P0 checks them.
* **Walltime:** a 44M fit is 1–5 h. A crash loses the run, because `fit_cops` cannot resume TRF
  exactly. Suggestion S8 writes the state at every accepted step, so `--restart-from` can warm-start.

## 8. Suggestions (beyond the request)

* **S1. FP32 record store (§3.1)** even for 12M: removes the 12-min bind and halves I/O.
* **S2. Registry: one record load per row set,** basis-only rebind. Saves about 23 min of every 12M
  template fit today.
* **S3. Single-pass Jacobian.** All sums of a request in one launch with more accumulators, or a
  two-level reduction. In mixed mode each tile recomputes the static factors, so the Jacobian is 3.5× a
  pressure pass; a single pass could bring it close to 1.5×. Measure the register pressure first.
* **S4. GPU type matters more than count.** A6000, A40 and L40-class GPUs have weak FP64 (1/32–1/64); an
  A100 or H100 node is likely several times faster per GPU for these FP64-arithmetic kernels.
* **S5. Speculative Jacobian during pressure passes** when transfer-bound (§3.2).
* **S6. Validate-once policy** for very large bases (§3.6, Q4).
* **S7. Harmonics on the fly** instead of 25 basis columns per GPU (§3.5).
* **S8. Accepted-step checkpoints** for warm restarts after a lost job.
* **S9. Science first, size second:**
  * all four 12M fits stopped at max_nfev 2000;
  * 12M NC31 found a new branch that needs a converged check;
  * the cliff templates changed with the SHAP8 roughness (≥ 60°: 14.4% of the area against 4.1% on
    3M).
  * A 20M run tells whether 12M is converged in the quantities that matter. 44M (3.5× the cost) is
    worth it only if 12M and 20M differ. Raising max_nfev or the optimizer scaling (a CPU-branch
    decision) may matter more for the result than resolution.

## 9. Still open

* **R1. Disk for the local 20M run: RESOLVED** (user, 2026-10-21: base in `/tmp`, which survives until a reboot;
  store and everything else on `/data`). Measured: base 194 GB (5.06e9 records), store 101.1 GB.
* **R2.** `/home` is 100% full (2.7 GB left); the code, its `.venv` and this repository live there.
* **R3, measured in P0.** A100 memory (40/80 GB), sockets and GPU/PCIe-switch attachment, host RAM,
  node-local NVMe, network, Slurm limits and quota on Voyager.
* **R4. Deferred from P3 (2026-10-21):**
  * block-granular partition with split boundary points: row blocks are within 2.1% of ideal at 24 ranks on 20M
    (largest point 3.4% of a rank's share); revisit in P4 with the cost model, where resident/streamed shares
    matter more than the record count;
  * distributed base validation (I5): validate-once makes it a one-time cost on rank 0;
  * per-rank element subsets: every rank reads `elements.bin` (Y alone is 4 GB at 20M, 8.8 GB at 44M per rank);
  * the cluster gates of §4 (P0).
