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
- Three atomic-oxygen models of the O I 130.4-nm triplet (130.217, 130.486,
  130.603 nm; upper level 3s 3S1).  Each returns the rate of triplet photons
  per O atom (all three lines together, as an instrument that does not resolve
  the triplet sees them) and the corresponding power.  They are additive
  components, not alternatives to one another's physics:
  - `O_I_1304_KVB2022_DIRECT_J2`: direct solar pumping at 130.217 nm, from the
    Killen, Vervack & Burger (2022, ApJS 263:37) Figure-14 data.  Those g
    values use Gladstone (1992) full-disk line shapes (two offset Gaussians per
    line) scaled to SORCE SOLSTICE on 2009-03-01, and are given at 0.352 au;
    the module converts them to 1 AU.  The caller supplies the 1-AU triplet
    photon flux; the table is scaled by its ratio to the triplet flux the
    Killen et al. tables imply, 8.412e9 photons cm^-2 s^-1 (SOLSTICE gives
    8.12e9 for the same day).  Domain -48..+48 km/s: the generated Gladstone
    profiles end near +-48.3 km/s.
  - `O_I_1304_SUMER2001_DIRECT_J2`: the same process with the line shape of the
    SOHO/SUMER quiet-Sun disk-centre atlas (Curdt et al. 2001), served by
    BASS2000 at 0.043-A sampling.  The continuum-subtracted 130.217-nm line is
    normalized to the Killen et al. share of the supplied triplet flux
    (0.3248).  The atlas includes the measured line wings and the underlying
    continuum.  Domain -100..+100 km/s.
  - `O_I_1304_SUMER2001_LYB_BOWEN_J2`: solar Lyman-beta pumping of the 3d 3D
    level at 102.576 nm, followed by the 1128.7-nm and 844.6-nm cascade into
    3s 3S1 (Bowen fluorescence).  The SUMER atlas shape is normalized to the
    supplied 1-AU flux in 102.4--102.8 nm; the cascade branching ratio (0.288)
    comes from NIST ASD transition probabilities.  Domain -100..+100 km/s.

  All three assume that the absorbing atoms are in the 2p4 3P2 ground level.
  The 3P1 and 3P0 levels decay radiatively in 3.1 h and 16 h, much faster than
  an O atom is photoionized and much faster than it is re-pumped (10^-5 s^-1),
  so freshly produced atoms reach 3P2 within ~10^5 km of their parent.  Use the
  sum of one direct model and the Bowen model for the total 130.4-nm emission.
  The two direct models bracket the solar line shape: they agree within 30%
  for |v| < 10 km/s, but at |v| = 30--40 km/s the Gaussian wings of the
  Gladstone profile fall 4--50 times below the SUMER wings.  There, the
  Bowen term (about 4e-7 s^-1 at 1 AU) is comparable to or larger than the
  direct term.  Electron-impact excitation, prompt 130.4-nm emission from
  photodissociative excitation of H2O, CO, or CO2, and optical depth are
  outside these models.

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

`--h-only`, `--oh-only`, or `--o-only` restrict the refresh to one model
family.  The O refresh retrieves the Killen et al. archive (only its four
expected regular files are accepted), the two BASS2000 atlas windows, and the
two NIST ASD line lists; on 2026-11-19 it reproduced every reviewed file
byte for byte.

Then update the generated hash in the manifest and run all tests.  A changed
source hash is a scientific-data change and must not be accepted as a routine
formatting update.

The AMPS top-level build compiles this module into `libAMPS.a`.  Application
adapters select a named model and provide case conditions; no comet, observing
campaign, mesh, or species-index constants belong in this module.
