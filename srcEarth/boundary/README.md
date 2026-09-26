# Boundary spectra — Roadmap Step 6

`spectrum.h` implements the five boundary spectrum families used by cutoff-to-flux
products: `POWER_LAW`, `POWER_LAW_CUTOFF`, `LIS_FORCE_FIELD`, `BAND`, and `TABLE`.
`spectrum.cpp` owns the single `gSpectrum` instance and writes `spectrum_input.dat`.

## Coordinate and uncertainty metadata

Every family accepts these optional keys in `#SPECTRUM`:

| Key | Values | Legacy default |
|---|---|---|
| `SPEC_ENERGY_BASIS` | `PER_PARTICLE`, `PER_NUCLEON` | `PER_PARTICLE` |
| `SPEC_INTENSITY_UNIT` | `PER_MEV`, `PER_MEV_PER_NUCLEON` | `PER_MEV` |
| `SPEC_MASS_NUMBER` | finite value greater than zero | `1` |
| `SPEC_RELATIVE_UNCERTAINTY` | finite value greater than or equal to zero | `0` |

For `PER_NUCLEON`, the public energy and channel bounds are MeV/nucleon, and intensity
is per MeV/nucleon. The trajectory adapters multiply the coordinate by mass number
before computing speed or rigidity. Spectrum evaluation and energy integration remain
in the declared coordinate. Declaring `PER_NUCLEON` with `PER_MEV`, or
`PER_PARTICLE` with `PER_MEV_PER_NUCLEON`, is an input error.

The uncertainty is a relative multiplicative interval. For boundary intensity `J` and
uncertainty `u`, the interval is `[max(0,1-u)J,(1+u)J]`. It is combined with the
independent access lower/upper interval only in `BoundaryProducts`, so numerical
non-resolution and source uncertainty remain distinguishable in the code.

## Time-dependent TABLE format

A legacy static table remains two positive columns:

```text
1      1000
10      100
100      10
```

A time-dependent table has one shared coordinate grid and one positive intensity per
grid point on every UTC row:

```text
ENERGY_MEV: 1 10 100
2026-01-01T00:00:00  1000 100 10
2026-01-01T00:10:00  2000 200 20
```

Selection is configured with:

| Key | Supported values |
|---|---|
| `SPEC_TIME_INTERPOLATION` | `LOG_INTENSITY` (`LOG_LINEAR` alias) |
| `SPEC_TIME_MAX_GAP_S` | `0` for unlimited, otherwise a finite nonnegative limit |
| `SPEC_TIME_GAP_POLICY` | `INTERPOLATE_FLAG`, `HOLD_NEAREST`, `FAIL` |
| `SPEC_TIME_OUT_OF_RANGE` | `CLAMP`, `ZERO`, `FAIL` |

Interpolation is linear in time and logarithmic in intensity. All rows must have the
same number of columns, UTC epochs must be unique, and intensities must be positive and
finite. `HOLD_NEAREST` resolves an exact midpoint tie to the earlier row. Every product
records the resulting status, gap flag, and interpolation fraction in `AUXDATA`; each
gridless trajectory spectrum zone records its own selection. Step 9 converts the SWMF
reference epoch plus authoritative PT clock into one absolute UTC used by the field,
coordinate transforms, ephemerides, and boundary spectrum. Fractional seconds through
nanosecond text precision are parsed rather than truncated. The lower-level
`SetEvaluationEpochUTCOffset()` API remains covered as an equivalent reference, but the
coupled adapter does not apply a second offset to the already absolute epoch.

`spectrum_input.dat` retains its historical two data columns and adds `AUXDATA` for
energy basis, mass number, intensity unit, uncertainty, temporal status, and gap flag.

## Validation

Run:

```bash
./test/UBoundaryProducts/run_test.sh
```

The suite compiles `spectrum.cpp` and compares all analytic families and table
interpolation with independent closed forms. It also exercises the production
key/value adapter, contradictory-unit rejection, a gap-crossing time table, ZERO and
FAIL policies, per-nucleon energy/Jacobian behavior, and uncertainty bounds.
