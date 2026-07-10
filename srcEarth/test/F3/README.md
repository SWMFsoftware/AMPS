# F3 — Dipole cutoff-filtered flux

F3 validates the first end-to-end magnetic-access flux calculation in the SEP
validation campaign.  It links the dipole cutoff benchmark to the density/flux
folding path by running gridless density/spectrum in a centered aligned dipole:

```text
FIELD_MODEL               DIPOLE
DS_TRANSMISSION_MODE      SCAN
DS_TRANSMISSION_SCAN_N    100
SPECTRUM_TYPE             POWER_LAW
SPEC_GAMMA                3.5
SPEC_EMIN/SPEC_EMAX       1 / 1000 MeV
```

The validation document describes F3 as a 9000 km shell latitude profile.  The
runner implements this with explicit `POINTS` rather than `SHELLS`, so it can use
the standard `gridless_points_density.dat`, `gridless_points_spectrum.dat`, and
`gridless_points_flux.dat` files.  The default point set is:

```text
alt = 9000 km
lon = 0, 90 deg
lat = -70, -60, -45, -30, 0, 30, 45, 60, 70 deg
```

## Reference solution

For each point, `run_F3.py` computes the vertical Størmer cutoff rigidity

```text
Rc(lambda,r) = 14.9 cos^4(lambda) / r_RE^2 GV
```

then converts that rigidity into a proton kinetic-energy cutoff, `Ecut`.  The
reference flux in each energy channel is

```text
F_ch = 4*pi*T_open*int_{max(E1,Ecut)}^{E2} J0*(E/E0)^(-gamma) dE
```

where `T_open` is the analytic straight-line open-sky fraction outside the inner
absorbing sphere.  This is an approximate external reference: the AMPS result is
computed from a full directional access map and contains penumbra structure,
while the Størmer reference collapses access to a vertical hard cutoff.

The runner also applies tighter internal checks that should be exact apart from
ASCII output precision:

- `J_local(E) = T(E) J_boundary(E)`
- density and flux reconstructed from the saved spectrum match the AMPS files
- centered-dipole longitude symmetry
- centered-dipole north/south symmetry
- total flux increases toward high absolute latitude

## Running

From the directory containing the `amps` executable:

```bash
python srcEarth/test/F3/run_F3.py -np 4 -nt 16
python srcEarth/test/F3/run_F3.py --lons 0,90,180,270
python srcEarth/test/F3/run_F3.py --analytic-flux-tol 0.75 --rc-tol 0.75
python srcEarth/test/F3/run_F3.py --dry-run --workdir test_output/F3_dryrun
python srcEarth/test/F3/run_F3.py --skip-run --workdir test_output/F3_gridless
```

## Files

Repository files:

```text
AMPS_PARAM_F3_gridless.in
reference_F3_dipole_cutoff.csv
run_F3.py
README.md
```

Run-directory outputs:

```text
AMPS_PARAM_F3.in
F3_amps.log
F3_summary.csv
F3_result.json
reference_F3_dipole_cutoff_used.csv
gridless_points_density.dat
gridless_points_spectrum.dat
gridless_points_flux.dat
```
