# srcGeospace status

As of **2026-10-27 20:05 CDT**.  Repository `/data/vtenishe/GEO/AMPS`; builds and test
results outside this tree in `/data/vtenishe/GEO/AMPS_runs/srcGeospace` (below `$RUNS`).

This file describes the current state.  The GPU work (parked 2026-10-27) is documented in
[GPU_STATUS.md](GPU_STATUS.md).  The history, the reasoning behind each decision and the
session log are in [handoff.md](handoff.md) (sections 10 and 11 cover 2026-10-26/27); the
design and build references are in [doc/](doc/README.md).

---------------------------------------------------------------------------------------

## 1. Summary

* **CPU validation of the current code is complete**, in double and float, for the full
  test list of that time (62 tests plus SA01-SA03; srcGeospace's former copy of srcEarth's
  list, see section 2):
  * double: 64 of 65 as listed; only C10 fails (it also fails in the 2026-10-22 baseline
    and in srcEarth);
  * float: 59 of 65 as listed; the six failures are the known float failures, with
    outputs identical to the 2026-10-23 float pass.
* **AMR allocator fix in the AMPS core** (`src/meshAMR/meshAMRdef.h`, committed by the
  user as `b06e7e7`): the mesh setup on one MPI rank was quadratic in the number of
  mesh elements.  The C2 3d MESH test takes 100 s instead of 672 s on 4 ranks and now
  completes on 1 rank (about 4 min; it never finished before).  Cutoff results are
  byte-identical (section 3).
* **GPU work parked** (user decision 2026-10-27; all GPU runs stopped at 12:20).  Status,
  results, open issues and how to resume: [GPU_STATUS.md](GPU_STATUS.md).  In short: the
  54 double and 55 float tests that finished on one GPU give the CPU's outputs and
  pass/fail results except C3 3d MESH (section 4.3); one GPU took 7.6 (double) and 5.3
  (float) times the CPU time over 53 tests, 5-83 times in the tracing of the mesh decks.
  The engine `gpu/` remains the CPU path of srcGeospace; only the nvc++ device pool is
  parked.
* **Open issue**: the merged GPU build gives different C3 3d MESH cutoffs than the CPU
  (double and float); the GPU build of 2026-10-22 matched the CPU (section 4.3).
* **Performance** (quiet machine, 2026-10-27, `list_timing`, two rounds, same core and
  flags; section 5.1): srcEarth 327.9 s, srcGeospace 186.4 s; srcEarth with the numerics
  tuning of 2026-10-27 180.0 s (0.55 of the original, 0.97 of srcGeospace).  Float is not
  faster than double on the CPU (ratio 1.00).
* **Roles** (user decision 2026-10-27): srcEarth stays the main tool for SEP tracing in
  geospace and keeps everything it has; srcGeospace may replace it once the GPU slowdown is
  solved.  srcGeospace compiles on its own (every file the build needs is in it) and keeps
  only what srcEarth does not have; its end-to-end tests are srcEarth's tests, run with
  srcEarth's runner (`test/tools/run_srcearth_tests.sh`).
* **For the user to commit** (section 6): the reduced srcGeospace tree, by amending
  `465b15b` before pushing; the srcEarth tuning (4 files; full list validated except the
  comparison of C19, whose rerun runs: same pass/fail results and outputs as the
  unmodified srcEarth so far), copied into the repository srcEarth once C19 is compared.

## 2. Code state

| item | state |
|---|---|
| git `HEAD` | `465b15b` (2026-10-27 02:20, srcGeospace committed by the user), on top of `bf68a2a` (srcEarth shieldSim script) and `b06e7e7` (the allocator fix).  `55ecae0` and `f9c2c00` (between `47fccd5` and `b06e7e7`) do not touch `src/` or `srcEarth/` |
| core interpolation work | committed by the user: `b6cb6cb` (stencil-cell cache), `9698bcd` (`Stencil.Length` fix), `47fccd5` (row-only construction, row-structure reuse, dropped-weight normalization) |
| `src/meshAMR/meshAMRdef.h` | committed by the user as `b06e7e7`: doubling growth of the AMR element stacks' chunk lists (section 3).  `src/` has no uncommitted changes |
| `srcGeospace/` | committed by the user as `465b15b` (596 files; not pushed).  Working tree reduced on 2026-10-27 (user decision, section 2.1) to 152 files: 446 deleted (`test/` except srcGeospace's own unit test and tools, `studies/`, 56 non-build copies of srcEarth files), 2 new (`GPU_STATUS.md`, `test/tools/run_srcearth_tests.sh`), documents and tools updated; to be committed by amending `465b15b`.  `shieldSim/` removed earlier the same day (kept only in `srcEarth/shieldSim`) |
| `srcEarth/` | unchanged in the repository so far.  The numerics tuning of 2026-10-27 (4 files: `gridless/CutoffRigidityGridless.cpp`, `3d/GlobalMagneticField.cpp`, `3d/FieldInitializationPlan.h`, `3d/Mode3D.cpp`; work copy `$RUNS/work/srcearth_tune/srcEarth`) is copied in once its validation is complete (section 6), to be committed by the user (message `$RUNS/work/COMMIT_MSG_srcEarth_tuning.txt`) |
| `Makefile` | modified before this work (user's state), not touched |

Decisions in force (details in handoff.md sections 1, 5 and 10):

* Normalization rule of the row stencil: **option (a)**, the new rule in core and engine
  (2026-10-26).  The GRIDDED round-off differences against the 2026-10-22 baseline are
  accepted; the 2026-10-23 and 2026-10-26 merged passes are the reference.
* CPU validation covers the full list including SA01-SA03.
* Heavy tests (C2 3d MESH, C19, SA01-SA03) never overlap each other; other tests run next
  to them while about 80 GiB of memory stays free (memory-based lock, section 8).
* Builds are clean builds in private stages (`scripts/build_amps.sh`); results stay outside
  `srcGeospace`; no commits by this work.
* 2026-10-27: GPU development parked (GPU_STATUS.md).  srcEarth: no engine; its existing
  numerics are tuned with small changes (this lifts "srcEarth is never modified" for that
  tuning).  Roles and layout of srcEarth and srcGeospace: section 2.1.

### 2.1 Layout of srcGeospace against srcEarth (2026-10-27)

* **Builds on its own**: all sources the AMPS build compiles or includes (srcEarth's,
  namespaces renamed), the engine `gpu/`, `makefile`, `input/`, `scripts/`.
* **Its own, not in srcEarth**: `gpu/`, the engine glue in 10 files, `scripts/`, `doc/`,
  `tools/`, `BUILD.md`, `STATUS.md`, `GPU_STATUS.md`, `handoff.md`, the unit test
  `test/UEngineEquivalence` (list `test/list`) and the tools `test/tools/` (runs of
  srcEarth's tests, list subsets with `-np/-nt` rewriting, output comparison, cross-build
  report, timing).
* **Taken from srcEarth, not duplicated**: the tests and their runner
  (`srcEarth/test/list`, `srcEarth/test/test_runner.py`), studies, examples, release
  validation, standalone-campaign tools, plotting/post-processing scripts and the
  application documentation (`srcEarth/README.md`).  srcEarth's tests run against a
  srcGeospace executable from a run root (`test/tools/run_srcearth_tests.sh`); its unit and
  source-text tests check srcEarth's sources, srcGeospace's code is covered by the
  end-to-end tests, `UEngineEquivalence` and `AMPS_EARTH_ENGINE_VERIFY=1`.
* Files: 152; 84 have a counterpart in srcEarth (4 identical, 65 differing only by the
  renaming, 15 with other changes), 68 exist only in srcGeospace (doc/INVENTORY.md).

## 3. AMR allocator fix (`src/meshAMR/meshAMRdef.h`, commit `b06e7e7`)

**Problem.**  `cAMRstack`, `cAssociatedDataAMRstack` and `cAMRheap` allocate mesh elements
in chunks of 10 (`_STACK_DEFAULT_BUFFER_BUNK_SIZE_`) and keep lists of pointers to the
chunks.  `initMemoryBlock()` grew those lists by a fixed 100 entries
(`_STACK_DEFAULT_BUFFER_LIST_SIZE_`) and copied the whole list each time: one full copy per
1,000 elements, quadratic in the elements per rank.  One rank holding the whole mesh is the
worst case; the 1-GPU C2 3d MESH spent 47 min in this phase without finishing.

**Fix.**  Macro `_STACK_NEXT_LIST_SIZE_(s)` (after line 175, with a comment): the first
growth adds 100 entries, later growths double the list.  It replaces the fixed step in the
allocation, the memory accounting, the clearing loop and the size update of the three
`initMemoryBlock()` functions.  `cAssociatedDataAMRstack` sizes its lists from
`BaseElementStack.dataBufferListSize`, so both use the same macro.

**Verification** (C2 3d MESH, CPU double; `$RUNS/results/2026-10-26_allocfix/`):

| run | time | peak memory | outputs |
|---|---|---|---|
| 4 ranks, without the fix (CPU double pass) | 672.5 s | 36.9 GiB/rank, 143.0 GiB total | reference |
| 4 ranks, with the fix | 100 s | 31.4 GiB/rank, 118.5 GiB total | cutoff files, summaries, plot, manifest byte-identical |
| 1 rank, with the fix | 246 s | 99.4 GiB | byte-identical to the 4-rank reference |
| 1 GPU float, with the fix (gpu1_float #020) | 266 s | 99.4 GiB | byte-identical to the CPU float run |

Without the fix the 1-rank setup did not finish in 47 min (GPU, 2026-10-26) or was killed
(exit 137 at 792 s, 2026-10-22).

Side effects:

* `cAssociatedDataAMRstack::getAllocatedMemory()` counts list capacity, so the "Memory
  usage" lines of the AMPS log can show up to 2 times larger values.  Nothing else reads it.
* The saved mesh file stores the tree-node list size (50200 -> 51200 in C2) and unset
  padding bytes (`cAMRnodeID`; the tree node's char refinement levels written with
  `sizeof(int)`), so it is not byte-identical; every differing byte was accounted for.
  Loading an old mesh file reads the stored size and should work, but was not tested.

Binaries with the fix: `builds/merged_gpu_float`, `builds/cpu_double_2026-10-26_allocfix`,
`builds/c9tail_gpu_{double,float}`; `builds/merged_gpu_double` is rebuilt with it after the
1-GPU double pass (section 6).  The CPU passes and the 1-GPU double pass used binaries
without the fix.

## 4. Test results

### 4.1 CPU (2026-10-26/27)

Builds `builds/cpu_{double,float}_2026-10-26` (clean, `HEAD` 47fccd5, without the allocator
fix).  `test_runner.py -j 6` (six tests at a time, like srcEarth) next to the GPU stream.
Results: `$RUNS/results/2026-10-26_merged_cpu/{double,float}`.

| pass | as listed | failed | wall time |
|---|---|---|---|
| double | 64 / 65 | #054 C10 | 17:02-21:29 (4.5 h) |
| float | 59 / 65 | #016 C1 3d MESH, #020 C2 3d MESH, #052 C9 FULL_SCAN, #053 C9 DIRECT_ACCESS, #054 C10, #057 C17 | 21:29-01:17 (3.8 h) |

* The six float failures are those of the 2026-10-23 float pass.  C1/C2 3d MESH summaries
  are byte-identical to it; the C9 x2, C10 and C17 result files are identical apart from
  paths and times.
* SA01-SA03 pass as workflows in double (2528 / 2027 / 1490 s) and, for the first time, in
  float (2620 / 2246 / 1494 s).  Their scientific status in `qualification_summary.json`
  is **FAIL** in both precisions, as in the 2026-10-23 srcGeospace runs and the 2026-10-16
  srcEarth runs (the list uses the WORKFLOW exit-code policy, which records but does not
  gate the scientific result).  Not a regression.

### 4.2 GPU (stopped 2026-10-27 12:20; details in GPU_STATUS.md section 3.2)

GPU list `list_cpu` = the full list without SA01-SA03 (62 tests), runner options
`--np 1 --nt 32` (one GPU; C13/C18 keep their own settings).
Results: `$RUNS/results/2026-10-26_merged_gpu/<pass>`.

| pass | GPU | binary | final state |
|---|---|---|---|
| `gpu1_double` + `gpu1_double_rest` | 0 | `merged_gpu_double` (2026-10-23, without the fix) | #001-#055 done as on the CPU (C10 fails, identical scores); #020 C2 3d MESH killed in single-rank setup; #056 C19 killed at snapshot 2 of 3 |
| `gpu1_float` | 1 | `merged_gpu_float` (with the fix) | #001-#055 done; the failures (C1 and C2 3d MESH, C9 FULL_SCAN and DIRECT_ACCESS, C10) are exactly the CPU float failures; #056 C19 killed 2.5 min into tracing |
| not run | | | C17, C16, C15, C13, C12, C8 (both precisions); the 2-GPU passes; the rebuild of `merged_gpu_double` with the fix; the C9 tail study |

No run root has `test_runner_report.*` files (runners killed); results are in
`<pass>_runner_stdout.log`, the per-test logs and `test_output/`.

### 4.3 GPU outputs against the CPU

`test/tools/compare_outputs.py` over the finished GPU tests
(`$RUNS/results/2026-10-26_merged_gpu/gpu1_{double,float}_vs_cpu_*_partial.{txt,csv}`):

| | identical | numeric differences | different |
|---|---|---|---|
| double (#001-#051) | 575 + 8 after paths, 44 known nondeterministic | 104 | 1 |
| float (#001-#051) | 473 + 2 after paths, 40 known nondeterministic | 104 | 1 |
| double and float (#052-#055: C9 FULL_SCAN, C9 DIRECT_ACCESS, C10, C18; GPU_STATUS.md section 3.3) | 775 + 15 after paths, 136 known nondeterministic (each) | 0 | 0 |

* Most numeric differences are input decks written with other `np`/`nt`, plus summation
  order at 1e-14 relative in F5 and F11 (one rank with 32 threads against 4 x 16).
* **C3 3d MESH** (list line 64) is the exception: BINARY-search cutoffs land at other
  locations and UPPER_SCAN cutoffs differ (median 1.0 %, maximum 17 % in double), in
  double and float.  The size is that of a rounding-level change (CPU float against CPU
  double differs as much); GPU_STATUS.md section 3.4 has the full comparison.
  * CPU C3 is reproducible: the 2026-10-22, 10-23 and 10-26 CPU runs are identical.
  * The 2026-10-22 GPU build (before the merge of the experimental engine) matched it, but
    that run used 4 ranks x 16 threads on one GPU, the 2026-10-26 runs 1 x 32.
  * Excluded as causes **on the CPU** (`$RUNS/work/c3_spec_check/`): speculative search
    evaluation (CPU with `AMPS_EARTH_SPECULATION=1` matches the CPU reference) and the rank
    count (CPU with `--np 1 --nt 32` matches too).
  * Likely, but not established, a GPU-only part of the merged engine; the isolation plan
    is GPU_STATUS.md section 5, item 2.
* Gridded C6/C7 decks (72 files identical, 8 known nondeterministic), C1 3d MESH, C2 3d MESH
  (float; the double run was killed) and C11 agree with the CPU byte for byte; their only
  differing files are the input decks (`np`/`nt`).

### 4.4 Failures and their status

| test | double CPU | float CPU | GPU | status |
|---|---|---|---|---|
| C10 (GRIDDED, DIRECT_ACCESS) | F | F | F in both precisions (identical scores and outputs) | fails since the 2026-10-22 baseline and in srcEarth (2026-10-16); not analysed |
| C1 3d MESH | P | F | float F (same output) | known float failure (`doc/PRECISION.md`) |
| C2 3d MESH | P | F | float F (same output); double killed | known float failure |
| C17 | P | F | not run | known float failure |
| C9 FULL_SCAN, C9 DIRECT_ACCESS | P | F | double P, float F (outputs identical to the CPU) | float failures known since 2026-10-21 (time of the minimum cutoff latitude; also with the FMA rounding baseline) |
| C3 3d MESH | P | P | P, but cutoffs differ from the CPU | open, section 4.3 |
| SA01-SA03 scientific status | FAIL | FAIL | not run | same in srcEarth; workflow passes |

## 5. Performance

### 5.1 srcGeospace against srcEarth (CPU)

* **Same core, 12 tests, two rounds** (2026-10-22, `doc/VALIDATION_CPU_2026-10-21.md`
  section 5; machine loaded, ratios good to 10-30 %): total 0.69 of srcEarth's time;
  gridless decks with many directions 0.20-0.41; C8, F11 0.61-0.80; Mode3D MESH decks
  about 1.0.
* **Full list, today against srcEarth's last full pass** (srcEarth: repository-root
  `test_runner_report.json`, 2026-10-16, up to 5 tests at a time, old core; srcGeospace:
  2026-10-26 double, six at a time next to the GPU stream, lock waits subtracted):

  | | srcEarth | srcGeospace | ratio |
  |---|---|---|---|
  | all 64 tests | 31,107 s | 17,167 s | 0.55 |
  | 8 tests over 1000 s | | | 0.50 |
  | SA01 / SA02 / SA03 | 6172 / 4757 / 3680 s | 2527 / 2027 / 1490 s | 0.41 / 0.43 / 0.40 |
  | C19 SMOKE / C9 FULL_SCAN | 4900 / 3807 s | 2059 / 1818 s | 0.42 / 0.48 |
  | 18 tests of 30-1000 s | | | median 0.77 |
  | 38 tests under 30 s | | | 0.95 |

  Not a clean comparison: srcEarth was built before the core caching commits and the
  allocator fix (srcEarth gets those too when rebuilt), and the load differed.  The gridless
  C6 subtests looked 1.7-2.1 times slower in double, but in the float pass, run one after
  another instead of overlapping, they took about srcEarth's time (35 / 60 / 63 s against
  33 / 54 / 64 s): most likely load, to be confirmed.
* **Quiet machine, same core `b06e7e7` and flags, `list_timing`, two rounds of opposite
  order, each test alone** (2026-10-27; `$RUNS/results/2026-10-27_timing/`,
  `timing_final.log`; round totals within 0.6 %).  "tuned" is srcEarth with the numerics
  tuning of 2026-10-27 (section 6):

  | test | srcEarth | srcEarth tuned | srcGeospace | tuned/srcEarth | tuned/srcGeospace |
  |---|---|---|---|---|---|
  | C1 3d MESH | 30.3 | 28.8 | 31.3 | 0.95 | 0.92 |
  | C1 gridless | 1.9 | 1.6 | 1.6 | 0.82 | 0.97 |
  | C6 INITIAL GRIDDED | 58.9 | 46.5 | 49.9 | 0.79 | 0.93 |
  | F1 | 2.8 | 2.6 | 1.8 | 0.93 | 1.49 |
  | F3 | 38.5 | 10.4 | 14.9 | 0.27 | 0.70 |
  | F4 | 10.7 | 2.4 | 2.2 | 0.23 | 1.09 |
  | F11 | 20.9 | 15.2 | 14.1 | 0.73 | 1.08 |
  | C17 | 15.6 | 3.9 | 4.2 | 0.25 | 0.93 |
  | C13 | 43.7 | 12.4 | 13.8 | 0.28 | 0.90 |
  | C12 | 28.3 | 10.0 | 11.8 | 0.35 | 0.85 |
  | C8 | 76.2 | 46.0 | 40.7 | 0.60 | 1.13 |
  | **total (12, incl. F4 self-tests 0.1 s)** | **327.9 s** | **180.0 s** | **186.4 s** | **0.55** | **0.97** |

  srcGeospace with `AMPS_EARTH_ENGINE=LEGACY` took 0.97 of the original srcEarth's time
  (baseline timing of the morning, `timing_baseline.log`), so the original difference was
  the batched engine.  The gridless C6 subtests (`timing_c6g.log`): srcEarth 33.5 / 54.5 /
  64.0 s, srcGeospace 37.1 / 54.2 / 62.9 s; the earlier "1.7-2.1 times slower" was load.

### 5.2 Double against float (CPU, 2026-10-26)

60 tests with the same outcome in both passes, lock waits subtracted: double 13,709 s,
float 13,673 s (ratio 1.00).  Tests over 1000 s 0.94-1.11; 100-1000 s 0.54-0.85;
C8 1.10 (float slower, as on 2026-10-22).  Float is meant for the GPU, not for the CPU.

### 5.3 GPU against CPU

Final numbers in GPU_STATUS.md section 4 (they supersede the partial ones given here
earlier, 6.6 and 3.2 on other test sets).  Over the 53 tests finished in all four passes,
lock waits subtracted: one GPU took 7.6 times (double, binary without the fix) and 5.3
times (float, with the fix) the CPU time; the tracing phases of the Mode3D mesh decks
5-83 times; C9 FULL_SCAN 25,804 s (double) and 18,389 s (float) against 1,818 / 1,781 s.
Not a clean comparison (shared machine, 1 rank x 32 threads against 4 or 8 x 16).  The
C9 tail observation and the prepared study: GPU_STATUS.md section 4.6.

### 5.4 Peak memory of the heavy tests

Total resident memory of all processes of a test, 30 s samples (`memlog_proc.txt`):

| test | double | float |
|---|---|---|
| SA01-SA03 | 37.3-37.4 GiB | 37.7 GiB |
| C19 SMOKE | 115.4 GiB | 123.6 GiB |
| C2 3d MESH, 4 ranks, without the fix | 143.0 GiB | 159.9 GiB |
| C2 3d MESH, 4 ranks, with the fix | 118.5 GiB | -- |
| C2 3d MESH, 1 rank, with the fix | 99.4 GiB | 99.4 GiB (GPU) |

## 6. srcEarth tuning, its validation, and what follows

The GPU runs and `memmon.sh` were stopped on 2026-10-27 at 12:20 (user decision; what had
not run is listed in GPU_STATUS.md section 3.2).  User decisions of the afternoon: no
engine in srcEarth, instead tune srcEarth's existing numerics with small changes; quiet
timing of srcEarth and srcGeospace (done, section 5.1); roles and layout of the two
applications (section 2.1); run tests concurrently where memory allows.

1. **srcEarth numerics tuning**, 4 files (work copy `$RUNS/work/srcearth_tune/srcEarth`;
   copied into the repository srcEarth for the user to commit once the validation is
   complete, message `$RUNS/work/COMMIT_MSG_srcEarth_tuning.txt`):
   * `gridless/CutoffRigidityGridless.cpp`: the field evaluator normalized and compared the
     model name on every field evaluation (`Model()`, 69 % of F4's CPU); the canonical name
     and the branch are now fixed at construction.  The main thread of the gridless
     cutoff driver polled the end of each chunk with a 200 ms sleep while all workers
     idled; the last worker now wakes it.
   * `3d/GlobalMagneticField.cpp`: `CheckedGlobalCellIndex_` inlinable and without the
     redundant per-cell presence read; interior trilinear fast path of the Mode3D field
     interpolation (the AMPS stencil's own arithmetic for the case "trilinear, all 8
     cells inside the leaf"; `AMPS_MODE3D_INTERP_FASTPATH=0` off, `=verify` compares every
     result with the AMPS stencil bit for bit and aborts on a difference).  The switch was
     first named `AMPS_EARTH_INTERP_FASTPATH` (builds `tune3`-`tune5`); renamed at 19:10 to
     srcEarth's naming of Mode3D switches (`AMPS_MODE3D_*`).
   * `3d/FieldInitializationPlan.h`, `3d/Mode3D.cpp`: the field-writer plan reserves its
     hash map.
   * FMA (`-march=native` or `-mavx2 -mfma`, measured on top): 6-7 % on gridded Mode3D,
     about 1.6 % over `list_timing`, results change at round-off; not adopted.
2. **Validation** (`$RUNS/results/2026-10-27_srcearth_tune/`; both builds core `b06e7e7`,
   identical flags):
   * **Full list, tuned** (`validate.sh`, run root `tuned/`, build `srcearth_tune5`,
     `srcEarth` -> work copy so that the unit and source-text tests see the tuned sources;
     srcEarth's runner `-j 6`, 14:53-19:00): 62 of 64 as listed; `UStep12Release` (list
     digest) and C10 fail, as in srcEarth's own pass of 2026-10-16.
   * **Full list, unmodified** (`run_base_j12.sh`, run root `base_j12/`, build
     `srcearth_timing_base`, `-j 12` under the memory lock, 16:43-19:56): 60 of 64 as
     listed; besides `UStep12Release` and C10, C19 (refused to start, item 3; rerun
     running since 19:56) and C14 line 136 (collision, item 3; rerun P).
   * **Outputs** (`compare_all.sh`; `tuned_vs_base*.csv/txt`; all but C19, 20:00): 1,711
     files identical, 25 identical after replacing run-root paths, 198 AMR mesh-cache files
     (they differ between runs of the unmodified build too); numeric differences only in
     C13's timing and rank files and in F5/F11 last-bit sums (relative <= 4.2e-16 in the
     data files; the DYNAMIC MPI scheduler, they differ between runs of the unmodified
     build too, `doc/VALIDATION_CPU_2026-10-21.md`).  SA01-SA03, compared per run directory:
     identical apart from the mesh cache and the time stamps and pids in their validation
     reports.  C19: comparison after its rerun.
   * **Fast path in verify mode** (`verify_fastpath.sh`, `AMPS_EARTH_INTERP_FASTPATH=verify`
     with `tune5`): the 9 light Mode3D mesh entries (C1 and C3 3d MESH, C6 GRIDDED x3, C11,
     C14 CROSS_SOLVER MESH, C7 GRIDDED, C15) pass, no verify abort, outputs identical to
     `base_j12`; the variable reaches every `amps` rank (`verify_env_check/`, read from
     `/proc/<pid>/environ` during a run).
   * **Renamed build** (`check_tune6.sh`, `builds/srcearth_tune6` = the committed sources):
     the executable names only `AMPS_MODE3D_INTERP_FASTPATH`; srcEarth's 13 unit and
     source-text entries (`srcEarth` -> work copy) give the full-pass results
     (`UStep12Release` with the same digests), C1 3d MESH, C6 INITIAL GRIDDED, F4, C17 pass
     with outputs identical to `base_j12`.
3. **Incidents of the concurrent passes** (details in handoff.md section 12):
   * the unmodified pass ran from a copy of the list without its `! runner: exclusive`
     lines (the lock serializes the heavy tests), so C19 refused to start (it requires the
     runner's `AMPS_TEST_RUNNER_EXCLUSIVE=1`); rerun with that variable (`rerun_c19.sh`);
   * the two C14 entries (list lines 136 and 139) write the same files and collided at
     `-j 12` ("Cannot write standalone_mode3d_manifest.json"); both rerun one after the
     other (`rerun_base.sh`): P;
   * the tuned pass (srcEarth's runner, own memory gate) ran outside the lock, so its SA01
     and SA02 overlapped the unmodified pass's C2 3d MESH and SA01 and the lock admitted
     light tests without counting them: MemAvailable fell to 16.8 GiB at 17:22 (no OOM).
4. **srcGeospace layout** (section 2.1): working tree reduced and checked (clean build
   `builds/srcgeospace_reduced`; `UEngineEquivalence` through `run_srcearth_tests.sh`: P;
   C1 gridless, F4, C17 of srcEarth's list: P, outputs identical to the full-tree runs,
   `$RUNS/results/2026-10-27_srcgeospace_layout/`); the user amends `465b15b` before
   pushing.

## 7. Incidents of 2026-10-26

* **17:02-17:49, stall.**  Under the first (strict) lock the 1-GPU C2 3d MESH held the
  machine for 47 min in single-rank setup; killed at the user's request (`gpu1_double`
  #020, exit 143).  Root cause: the allocator (section 3).  The lock became memory-based.
* **18:07-18:13, investigation build wrote into the repository.**  A private build of the
  investigation agent used a link root whose `srcGeospace` was a symbolic link, so the AMPS
  configure step processed the repository `srcGeospace` in place (five files changed, all
  902 mtimes touched).  The five files were restored from the 2026-10-22 copy and the whole
  tree was verified against the 16:50 build-stage copy; `src/`, `srcEarth/` and the
  repository build files were unchanged.  Lost: the original mtimes.
  Do not build from `$RUNS/work/single_rank_setup/AMPS_repo_fix` (links to repository files).
* **23:50, GPU runner killed.**  To keep the original GPU stream script from also starting
  `gpu1_float` (the user asked to start it on the idle GPU 1), its shell was stopped with
  SIGSTOP; the Claude Code task manager then ended that task and killed its process group,
  including the `gpu1_double` runner after test #051.  #052 survived (its output goes to its
  own log), the remaining tests run in `gpu1_double_rest`.  Lesson: never signal the shell
  of a background task whose children must keep running.

## 8. Test infrastructure of these runs

All in `$RUNS/results/2026-10-26_parallel/`:

| file | purpose |
|---|---|
| `lockrun.sh` | cross-runner lock: `lockrun.sh heavy|light CMD`.  Heavy tests take `heavy.lock` (never overlap) and start when MemAvailable >= reserve + 30 GiB; light tests start one at a time while MemAvailable minus (reserve - current memory of the running heavy test) >= 80 GiB.  Reserve per test type from the measured peaks: SA 50 GiB, C19 135 GiB, C2 3d MESH and other heavy tests 170 GiB (`HEAVY_RESERVE_GIB` overrides).  Events in `lock_events.log`.  Earlier versions: `lockrun_strict_1648.sh`, `lockrun_res120_1749.sh`, `lockrun_res150_2035.sh` |
| `make_locked_list.py` | prefixes every list command with `lockrun.sh heavy|light` (heavy: after `! runner: exclusive`, or the C2 3d MESH command) |
| `cpu_stream.sh` | the CPU stream (finished): clean CPU builds, then double and float with `-j 6` |
| `gpu_stream.sh`, `gpu_stream2.sh` | the GPU stream (the first one died at 23:50; the second was stopped on 2026-10-27 at 12:20) |
| `c19_overlap_helper.sh` | 2026-10-27 12:14: held `heavy.lock` for the float C19 that the user allowed to overlap the double C19 (both stopped minutes later) |
| `memmon.sh` | `memlog.txt` (MemAvailable and amps total, 60 s), `memlog_proc.txt` (each amps process, 30 s) |

The runner reports' `command` fields carry the `lockrun.sh` prefix: strip it before
comparing with older run roots (`compare_builds.py` matches tests by command).  Timings of
these passes include machine sharing; subtract the waits in `lock_events.log`.

## 9. Pending decisions and tasks

User decisions or actions:

* amend `465b15b` with the reduced srcGeospace tree (section 2.1) before pushing (the user
  commits); commit the srcEarth tuning (section 6; 4 files, message in
  `$RUNS/work/COMMIT_MSG_srcEarth_tuning.txt`);
* progress bar for srcGeospace (proposal in handoff.md section 5.2) -- awaiting go-ahead;
* rename `gpu/` -> `engine/` -- offered, not requested.

Option B of 2026-10-27 (one copy of the tests and studies, in srcEarth) is done by the
layout of section 2.1; srcEarth's runner is not changed (the run root does what the
former `--run-root/--amps` options of srcGeospace's runner did, `sublist.py` rewrites
`-np/-nt`).

Tasks (suggested order):

1. Analyse C10 (fails in double and float, in srcEarth and in the tuned srcEarth).
2. `UStep12Release` (srcEarth): its approved list digest no longer matches
   `srcEarth/test/list` (fails in every pass since at least 2026-10-16); the list's
   approval belongs to srcEarth's maintainers.
3. Engine: skip `BuildMeshCache` when only B changes (same geometry; C10 uploads 10
   snapshots).
4. Documentation: `doc/VALIDATION_CPU_2026-10-21.md` section 7; a validation report of the
   2026-10-26 passes; `doc/BUILD.md` note on the allocator fix.
5. Cleanup: obsolete builds under `$RUNS/builds` (handoff.md section 8 lists the
   deletable ones; also `srcearth_tune1`-`tune4`, `srcearth_tune5_native`,
   `srcearth_tune5_avx2fma`, and `srcearth_tune5` once the tuning is committed),
   `$RUNS/work/srcearth_port/` (the abandoned engine port; its `rename_diff.py` is
   reusable); the still unanalysed
   `$RUNS/results/2026-10-22_srcGeospace_cpu/double_vs_earlier_double.txt`.
6. GPU tasks, when the GPU work resumes: GPU_STATUS.md section 5 (C3 3d MESH, verify on the
   GPU, rebuilt binaries, the rest of the validation, the C9 tail study, launch
   configuration, GPU cleanup).

## 10. Locations

Builds (`$RUNS/builds/`):

| stage | contents | allocator fix |
|---|---|---|
| `cpu_double_2026-10-26`, `cpu_float_2026-10-26` | CPU passes of 2026-10-26 | no |
| `cpu_double_2026-10-26_allocfix` | C2 3d MESH check of the fix | yes |
| `srcearth_timing_base`, `srcgeospace_timing_base` | 2026-10-27 timing baseline: unmodified srcEarth and srcGeospace, core `b06e7e7` | yes |
| `srcearth_tune6` | srcEarth with the numerics tuning as committed (section 6) | yes |
| `srcearth_tune5` | the same with the fast-path switch's first name `AMPS_EARTH_INTERP_FASTPATH`: the build of the full-list validation and the quiet timing; `tune1`-`tune4` intermediate, `tune5_native`, `tune5_avx2fma` FMA trials | yes |
| `srcgeospace_reduced` | the reduced srcGeospace tree (section 2.1), clean-build check | yes |
| `merged_gpu_double` | 1-GPU double pass (2026-10-23 build; never rebuilt) | no |
| `merged_gpu_float` | 1-GPU float pass (rebuilt 2026-10-26 20:05) | yes |
| `c9tail_gpu_double`, `c9tail_gpu_float` | C9 tail study (instrumented engine copy) | yes |
| `merged_cpu_double`, `merged_cpu_float` | 2026-10-23 merged CPU passes | no |
| `cpu_double` | 2026-10-22 baseline (pre-merge) | no |

Results (`$RUNS/results/`):

| directory | contents |
|---|---|
| `2026-10-27_timing/` | quiet-machine timing (scripts `build_srcearth.sh`, `run_timing_ab.sh`; logs `timing_baseline.log`, `timing_c6g.log`, `timing_final.log`; run roots `timing_<label>[_tag]_<round>`) |
| `2026-10-27_srcearth_tune/` | validation of the srcEarth tuning (section 6): run roots `tuned/` (`validate.sh`), `base_j12/` (`run_base_j12.sh`, reruns `rerun_base.sh`, `rerun_c19.sh`), `verify_fastpath/`, `verify_env_check/`, `tune6_check/`; comparisons `compare_all.sh` -> `tuned_vs_base*.csv/txt`, `verify_vs_base.*`, `tune6_vs_base.*`.  `base` is a guard file (a regular file, so that `validate.sh` did not start a second unmodified pass; `base_j12/` replaced it) |
| `2026-10-27_srcgeospace_layout/` | checks of the reduced srcGeospace tree (`unit/`, `e2e/`) |
| `2026-10-26_merged_cpu/{double,float}` | CPU passes (`double_j1_aborted`: the first, stopped -j 1 attempt) |
| `2026-10-26_merged_gpu/<pass>` | GPU passes, `stream.log`, `<pass>_runner_stdout.log`, GPU-vs-CPU comparisons |
| `2026-10-26_allocfix/` | builds and C2 3d MESH checks of the allocator fix |
| `2026-10-26_parallel/` | stream scripts, lock, memory logs (section 8) |
| `2026-10-23_merged_cpu/` | merged CPU passes of 2026-10-23 (reference since decision (a)) |
| `2026-10-22_srcGeospace_cpu/` | pre-merge baseline |

Work directories (`$RUNS/work/`): `srcearth_tune/` (work copy of srcEarth with the tuning,
`srcEarth_tuning.diff`, profiler `prof/`: `libsampler.so` LD_PRELOAD sampler with per-thread
CPU or wall timers, `report.py`, `callers.py`, `prof_run.sh`), `srcearth_port/` (abandoned
engine port; `rename_diff.py` = srcGeospace against srcEarth after reversing the renames),
`single_rank_setup/` (allocator investigation, diff, measurements, incident record),
`c9_gpu_tail/` (C9 study set-up), `c3_spec_check/` (C3 checks), `gpu_opt/` (experimental
engine copy, identical to the merged code).

srcEarth's last full CPU pass: `/data/vtenishe/GEO/AMPS/test_runner_report.json` and
`test_runner_logs/` (2026-10-16).

## 11. Reproduce

From `/data/vtenishe/GEO/AMPS`:

```
RUNS=/data/vtenishe/GEO/AMPS_runs/srcGeospace

# clean builds (each call wipes its stage first)
srcGeospace/scripts/build_amps.sh --stage $RUNS/builds/<name> -j 12             # CPU double
srcGeospace/scripts/build_amps.sh --stage $RUNS/builds/<name> --float -j 12     # CPU float
srcGeospace/scripts/build_amps.sh --stage $RUNS/builds/<name> --gpu [--float] -j 12
# srcEarth with the same configuration steps and flags (APP= a work copy; default the repository)
[APP=<srcEarth-dir>] $RUNS/results/2026-10-27_timing/build_srcearth.sh $RUNS/builds/srcearth_<name>

# srcEarth's tests against a srcGeospace build (six at a time), and srcGeospace's own unit test
T=srcGeospace/test/tools
$T/run_srcearth_tests.sh --run-root $RUNS/results/<run> --amps $RUNS/builds/<name>/amps -- -j 6
$T/run_srcearth_tests.sh --run-root $RUNS/results/<run>_unit --amps $RUNS/builds/<name>/amps \
    --list srcGeospace/test/list

# GPU pass on one GPU (GPU index in GEOSPACE_GPUS; list without SA01-SA03)
GEOSPACE_AMPS=$RUNS/builds/<gpu-stage>/amps GEOSPACE_GPUS=0 \
$T/run_srcearth_tests.sh --run-root $RUNS/results/<run> --amps srcGeospace/scripts/amps_gpu_select.sh \
    --exclude 'test/SA0[123]/' --np 1 --nt 32 -- -j 1

# compare two run roots
python3 $T/compare_outputs.py <ref-root> <new-root> --out <file>.csv

# a srcEarth build in a private stage (APP= a work copy, EXTRA_CXX= extra C++ flags)
$RUNS/results/2026-10-27_timing/build_srcearth.sh $RUNS/builds/srcearth_<name>
```

When several runners share the machine, use lists prefixed by `make_locked_list.py` so
that heavy tests do not overlap.  Never pass `--rebuild` to the SA runners (it rebuilds
the repository root).
