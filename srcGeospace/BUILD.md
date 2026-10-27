# Building and running srcGeospace (CPU and GPU)

Step-by-step guide from a fresh AMPS checkout to running the CPU and the GPU build.
Technical details (build internals, all engine environment variables, recorded build
outcomes): [`doc/BUILD.md`](doc/BUILD.md).  Engine design: [`doc/DESIGN.md`](doc/DESIGN.md).

Reference machine (all commands below were used there): Ubuntu 20.04, g++/gfortran 9.4,
Open MPI 4.0.3, Python 3.8, 2 x NVIDIA RTX A6000 (driver 550.163.01), NVIDIA HPC SDK
24.5 (nvc++, CUDA 12.4).

## 1. Prerequisites

| needed for | software |
|---|---|
| CPU build | C++17 compiler and Fortran compiler (g++/gfortran), MPI with `mpicxx`, `mpif90` (Open MPI), GNU make, Perl |
| field models, driver files | SPICE toolkit (`cspice`, NAIF) and kernels (section 2.3) |
| test runners, plots | Python 3 with `numpy`, `scipy`, `pandas`, `matplotlib`, `aacgmv2` (a few observational validation scripts also import `cdflib`) |
| GPU build only | NVIDIA HPC SDK (`nvc++` with OpenACC; 24.5 was used) and an NVIDIA GPU with a matching driver |

Python packages in a virtual environment (as in `test/list`):
```bash
python3 -m venv --system-site-packages .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy scipy pandas matplotlib 'aacgmv2>=2.6,<3'
```

## 2. Install AMPS (once)

### 2.1 Get the code
```bash
git clone git@github.com:SWMFsoftware/AMPS.git
cd AMPS                       # the repository root; srcGeospace/ is inside it
```

### 2.2 Configure AMPS
```bash
./Config.pl -install -compiler=gfortran,gcc_mpicc     # Makefile.conf: gfortran + gcc/mpicxx
./Config.pl -spice-path=/path/to/cspice               # SPICE library (cspice/lib/cspice.a)
./Config.pl -spice-kernels=/path/to/Kernels           # kernel root (section 2.3)
./Config.pl -f-link-option=-lmpi++                    # see below
```

* `-f-link-option=-lmpi++` sets `EXTRALINKEROPTIONS_F=-lmpi++` in `Makefile`.  srcGeospace
  is linked by `mpif90` (the Fortran field-model interface is on), which then also needs
  Open MPI's C++ library.  Use the C++ MPI library name of your MPI installation.
* Do **not** run `./Config.pl -openmp=on`: it also switches AMPS to its hybrid
  MPI+OpenMP compilation mode, which srcGeospace builds do not use.  The srcGeospace
  build script sets `OPENMP=on` (the `-fopenmp` flag) itself.
* `./Config.pl -application=...` is not needed: the srcGeospace build script uses the
  deck in `srcGeospace/input/` and sets the application entries in its build stage.
* `./Config.pl -show` prints the resulting settings.

### 2.3 SPICE kernels
`srcGeospace/input/srcGeospace.input` loads these files relative to the kernel root
given with `-spice-kernels`:
```
NAIF/NAIF0011.TLS  NAIF/pck00010.tpc  NAIF/DE405.BSP
OTHER/GSE.tf  OTHER/GSM.tf  OTHER/SM.tf  OTHER/earth_assoc_itrf93.tf  OTHER/earth_latest_high_prec.bpc
```
Without SPICE the code builds, but runs that read Tsyganenko driver files or use SPICE
frames stop with an error.

### 2.4 Where builds and results go
Builds are made in private stages, never in the repository: by default
`$GEOSPACE_RUNS/builds/<cpu|gpu|multicore>_<double|float>` with
`GEOSPACE_RUNS=<repository>/../AMPS_runs/srcGeospace`.  Test results also belong outside
the source tree (section 5).  Set `GEOSPACE_RUNS` to choose another location.

## 3. Build the CPU path

From the repository root:
```bash
srcGeospace/scripts/build_amps.sh -j 12                 # double -> $GEOSPACE_RUNS/builds/cpu_double/amps
srcGeospace/scripts/build_amps.sh -j 12 --float         # float  -> .../builds/cpu_float/amps
srcGeospace/scripts/build_amps.sh -j 12 --stage DIR     # any stage directory outside the repository
```
* Every call is a clean build: the stage is removed and recreated (only stages created
  by the script, or empty directories, are removed).
* The script links the read-only parts of the repository, copies `share/`, `input/`,
  `Makefile.local` and the `.conf` files into the stage, adds the srcGeospace deck, and
  sets in the stage `InputFileAMPS=geospace.input`, `APPLICATION=geospace`,
  `INTERFACE=on`, `OPENMP=on`.  SPICE paths and compilers come from the repository
  configuration of section 2 (it warns if SPICE is not configured).
* The trajectory engine runs on CPU threads in this build.
* Check: `$GEOSPACE_RUNS/builds/cpu_double/amps -h` prints the usage.

## 4. Build the GPU path

```bash
source /path/to/nvhpc/env.sh          # put nvc++ on PATH (here: /data/vtenishe/openacc/env.sh)
srcGeospace/scripts/build_amps.sh -j 12 --gpu           # double -> .../builds/gpu_double/amps
srcGeospace/scripts/build_amps.sh -j 12 --gpu --float   # float  -> .../builds/gpu_float/amps
```
* If `nvc++` is not on `PATH`, the scripts source the environment file given by
  `NVHPC_ENV=<file>` (or `--gpu-opts "--nvhpc-env <file>"`).
* The script first compiles only `gpu/EngineDevice.cpp` with nvc++/OpenACC into
  `<stage>/gpu_engine/libGeospaceEngine.so` (default `-gpu=cc86,cuda12.4,nordc,nofma`),
  then builds AMPS with g++ as for the CPU and links that library (RUNPATH set; no
  `LD_LIBRARY_PATH` needed).  It verifies that `amps` does not contain the CPU device pool.
* Another GPU or CUDA version: `--gpu-opts "--arch ccNN --cuda X.Y"` (e.g. `cc80` for
  A100, `cc90` for H100).  `--gpu-opts "--fma"` allows FMA (faster; results then differ
  from the CPU build at round-off level).  `--gpu-target multicore` builds the OpenACC
  multicore variant.  All options: `doc/BUILD.md` section 4.
* Double GPU results were bitwise identical to the CPU double build (`-gpu=nofma`) in the
  2026-10-22 checks; the merged build of 2026-10-23 differs on C3 3d MESH
  (`GPU_STATUS.md` section 3.4).  The GPU work is parked: `GPU_STATUS.md`.

## 5. Run

### 5.1 Command line (CPU and GPU builds alike)
```bash
mpirun -np <ranks> <stage>/amps -mode gridless    -i AMPS_PARAM.in [options]
mpirun -np <ranks> <stage>/amps -mode 3d          -i AMPS_PARAM.in [options]
mpirun -np <ranks> <stage>/amps -mode 3d_forward  -i AMPS_PARAM.in [options]
<stage>/amps -h                                     # all options and input-file sections
```
* `gridless`: direct IGRF/Tsyganenko (or dipole) field evaluation and backtracing, no
  mesh.  `3d`: AMR mesh with the field interpolated from the mesh (or analytic with
  `-mode3d-field-eval ANALYTIC`).  `3d_forward`: forward transport on the mesh.
* Threads per rank: `-mode3d-parallel THREADS -mode3d-threads N` (Mode3D cutoff),
  `-density-parallel THREADS -density-threads N` (gridless density/spectra); the test
  runners show typical combinations (e.g. `-np 4 -nt 16` in `srcEarth/test/list`).
* Common options: `-epoch YYYY-MM-DDTHH:MM:SS`, `-mover BORIS|HC4|RK4|...`,
  `-cutoff-search UPPER_SCAN|...`, `-mode3d-field-eval MESH|ANALYTIC`.
* Run in an empty working directory outside the repository; outputs are written there.
* Example inputs: `srcEarth/examples/` (e.g. `standalone_step7_combined.in`, gridless,
  dipole field; not part of the regression list).  Inputs of the regression tests are
  generated by their runners (`srcEarth/test/C*/run_*.py`), which are the best starting
  points.

### 5.2 CPU specifics
The batched engine is the default (`AMPS_EARTH_ENGINE=LEGACY` selects the original
per-trajectory code).  `AMPS_EARTH_ENGINE_VERBOSE=1` prints engine statistics, including
the backend (`CPU-threads/double` on the CPU build).

### 5.3 GPU specifics
* One GPU per MPI rank.  One rank: choose the device with `CUDA_VISIBLE_DEVICES=0`
  (or `ACC_DEVICE_NUM=0`).  Several ranks: start `amps` through
  `scripts/amps_gpu_select.sh`, which gives each local rank its own GPU:
  ```bash
  GEOSPACE_AMPS=$GEOSPACE_RUNS/builds/gpu_double/amps GEOSPACE_GPUS=0,1 \
      mpirun -np 2 srcGeospace/scripts/amps_gpu_select.sh -mode 3d -i AMPS_PARAM.in -mode3d-threads 16
  ```
* Prefer one rank per GPU with more host threads.  The GPU needs many trajectories in
  flight; small decks (hundreds of locations) run faster on the CPU build.
* `AMPS_EARTH_ENGINE_VERBOSE=1` shows the OpenACC backend; speculative cutoff searches
  are on by default on the GPU (`AMPS_EARTH_SPECULATION=0` turns them off; results are
  identical either way).  `AMPS_EARTH_STEPS_PER_LAUNCH` (default 2048 on GPU) sets the
  integration steps per kernel launch.  Other switches: `doc/BUILD.md` section 7.1.

### 5.4 Test suite
srcGeospace's end-to-end tests are srcEarth's tests (`srcEarth/test/list`), run against a
srcGeospace executable with srcEarth's runner; srcGeospace keeps only its own unit test
(`srcGeospace/test/list`: `UEngineEquivalence`) and tools (`srcGeospace/test/README.md`).
From the repository root (results go to the run root, outside the source tree):
```bash
R=$GEOSPACE_RUNS/results/<run name>
T=srcGeospace/test/tools
$T/run_srcearth_tests.sh --run-root $R/cpu_double --amps $GEOSPACE_RUNS/builds/cpu_double/amps -- -j 6
$T/run_srcearth_tests.sh --run-root $R/unit --amps $GEOSPACE_RUNS/builds/cpu_double/amps \
    --list srcGeospace/test/list

# GPU: one rank per GPU (C13/C18 keep their own counts)
GEOSPACE_AMPS=$GEOSPACE_RUNS/builds/gpu_double/amps GEOSPACE_GPUS=0 \
$T/run_srcearth_tests.sh --run-root $R/gpu1_double --amps srcGeospace/scripts/amps_gpu_select.sh \
    --exclude 'test/SA0[123]/' --np 1 --nt 32 -- -j 1
```
* Subsets: `--only '<regex>'` / `--exclude '<regex>'` of `run_srcearth_tests.sh`
  (`test/tools/sublist.py`).
* Compare two run roots: `srcGeospace/test/tools/compare_outputs.py REF NEW --out x.csv`,
  or `test/tools/compare_builds.py` for a report with plots (`test/tools/README.md`).
* One runner by hand: `srcGeospace/scripts/make_run_root.sh <root> <amps>`, then start the
  runner from `<root>`, e.g. `python3 srcEarth/test/C1/run_C1.py --mode gridless -np 4 -nt 16`.
* SA01-SA03 need the external data package given in their list entries.  Never add
  `--rebuild` to them: it deletes the repository `build/` and runs `make -j` there.
* Heavy tests (C19, C2 3d MESH; ~120 GB of memory) must not run concurrently.
