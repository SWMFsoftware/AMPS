# AMPS resonance-fluorescence models

`resonance_fluorescence` is a geometry-, application-, MPI-, and PIC-independent
core module for velocity-dependent cometary resonance-fluorescence models.  It
keeps published reference data and the equations that interpret those data in
one reusable location instead of embedding case tables in an application.

The public interface is `resonance_fluorescence.h`.  All public quantities use
SI units.  Positive radial velocity means motion away from the Sun.  Evaluation
is strict: an invalid distance, missing solar irradiance, or velocity outside a
published interval returns a non-`Ok` status; values are never clamped or
extrapolated.

## Implemented models

- `H_LYA_IKL2020` evaluates Equations 8--11 and 13 of
  Kowalska-Leszczynska et al. (2020), using the authors' nine Table-1
  coefficients.  The caller supplies the 1-AU solar Ly-alpha photon flux.  The
  module converts the published radiation-pressure ratio into a photon
  scattering rate by explicit momentum balance.
- `OH_AX00_SA1988_LOWELL_UNQUENCHED` linearly interpolates the reviewed Lowell
  Observatory OH A-X (0-0) energy-g-factor table.  The source quantity is an
  integrated energy rate.  The core deliberately does not invent a photon
  wavelength: applications that need photon counts must supply and document an
  effective wavelength.

The source files, retrieval snapshot, hashes, domains, units, and citations are
recorded in `reference/g_factors/manifest.json`.  `generated/reference_tables.inc`
is the immutable production representation; runtime evaluation performs no
file I/O.

## Regeneration and tests

From this directory:

```sh
python3 tools/generate_reference_tables.py --check
make test
```

To intentionally regenerate after reviewing a reference-data update:

```sh
python3 tools/generate_reference_tables.py
```

An upstream refresh is intentionally guarded against accidental replacement:

```sh
python3 tools/download_reference_data.py --replace-reviewed-data --workers 4
```

Then update the generated hash in the manifest and run all tests.  A changed
source hash is a scientific-data change and must not be accepted as a routine
formatting update.

The AMPS top-level build compiles this module into `libAMPS.a`.  Application
adapters select a named model and provide case conditions; no comet, observing
campaign, mesh, or species-index constants belong in this module.
