# srcGeospace

srcGeospace is the AMPS geospace energetic-particle application with a batched trajectory
engine (`gpu/`) that runs on CPU threads or on an NVIDIA GPU (OpenACC) from one source;
the original per-trajectory tracer remains available as the legacy path
(`AMPS_EARTH_ENGINE=LEGACY`).  It was created from `srcEarth` at commit `a29ec499`
(`SOURCE_REVISION`).

srcEarth remains the main application for SEP tracing in geospace.  srcGeospace is the
development line that may replace it once the GPU performance problem is solved (the GPU
work is parked; see [`GPU_STATUS.md`](GPU_STATUS.md)).  Therefore:

* **srcGeospace compiles on its own.**  Every file the build needs is in this tree: the
  application sources (srcEarth's, with the C++ namespaces `Earth`, `EarthUtil`,
  `Earth3DForward` renamed to `Geospace`, `GeospaceUtil`, `Geospace3DForward`), the
  engine `gpu/`, the input deck `input/` and the build scripts `scripts/`.
* **Everything else that srcEarth already provides is not duplicated here:**
  * application documentation (modes, products, input reference, test contracts):
    `srcEarth/README.md`, which applies to srcGeospace with `srcEarth` read as
    `srcGeospace` in source paths and `Earth` as `Geospace` in namespaces;
  * end-to-end tests: srcEarth's tests (`srcEarth/test/list`), run against a srcGeospace
    executable with srcEarth's runner; srcGeospace keeps only its own unit test and its
    tools ([`test/README.md`](test/README.md));
  * examples, release validation, standalone campaign tools, plotting and
    post-processing scripts, studies: `srcEarth/examples`, `srcEarth/release_validation`,
    `srcEarth/standalone_campaign`, `srcEarth/plot`, `srcEarth/studies` and the scripts
    at the top of `srcEarth`.

Where to look:

* Build: `srcGeospace/scripts/build_amps.sh [--float] [--gpu]` (clean builds in stages
  under `$GEOSPACE_RUNS/builds`, outside this tree); step by step in
  [`BUILD.md`](BUILD.md), technical reference in [`doc/BUILD.md`](doc/BUILD.md).
* Tests: `srcGeospace/test/tools/run_srcearth_tests.sh --run-root <dir> --amps <exe>`
  ([`test/README.md`](test/README.md)).
* Engine design, precision, validation and inventory: [`doc/README.md`](doc/README.md).
* State of the work: [`STATUS.md`](STATUS.md); GPU: [`GPU_STATUS.md`](GPU_STATUS.md).
