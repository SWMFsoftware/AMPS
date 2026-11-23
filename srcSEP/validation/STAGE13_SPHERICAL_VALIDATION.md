# Stage 13 spherical SWCME3D observational cases

## Build and select

Use a clean AMPS build directory for the shared proton field-line executable.
The compile-time `input/sep_stage13.conf` specifies the species, fixed-global-
step layout, and mover interface. `input/sep_stage13.in` uses the generic AMPS
named-case overlay, as in I3. It declares the physical boundary values and
the inputs for each comparison. One linked executable serves every case.

```sh
./Config.pl -install -compiler=gfortran,gcc_mpicc -spice-path=nospice
./Config.pl -application=sep_stage13 -amps-test=off
make -j
./amps -i input/sep_stage13.in --list-validation-cases
./amps -i input/sep_stage13.in -case gle-2012-05-17 --validate-runtime-input
```

Use a fresh evidence directory. Paths in each named case resolve against the
file that defines them. Use `-case` to select one named case or `--all` to
check the complete observational campaign. On hosts where MPI uses the
launcher, run the campaign driver as a single allocated rank so the legacy
OV01 linked process inherits that allocation. The numerical OV01 stage uses
MPI and the comparison plotting stage requires Matplotlib.

```sh
python3 -B srcSEP/validation/run_stage13.py --amps ./amps \
    -i input/sep_stage13.in --all --validate-runtime-input \
    --output-dir test_output/stage13-inputs
mpirun -np 1 python3 -B srcSEP/validation/run_stage13.py --amps ./amps \
    -i input/sep_stage13.in -case sep-2013-04-11 \
    --output-dir test_output/stage13-near-earth
```

The first command verifies the native selection and each declared artifact's
size and SHA-256. The scientific invocation records explicit SKIP and exits
77 while a required spherical comparison cannot run. An observation package
or model execution error exits 2. A scored failed criterion exits 1. The
April 2013 case runs the existing OV01 continuity comparison with the original
inputs, model, normalization and criteria. Its nested verdict applies to
OV01. The spherical comparison has its own status. JSON, JUnit, and native
command logs preserve the verdict and configuration in the evidence directory.

## Physical core and native admission

The core evaluates the spherical front and the moving open field lines at the
same epoch in HCI. It reuses canonical SWCME kinematics and the Stage 14
PFSS–Parker provider with mass per magnetic flux conserved on open tubes.
For the candidate wind it reports continuity, induction, momentum, and the heat
required by the isothermal closure. Its variational trajectory exports the
deformation tensor and Jacobian for independent material-volume tests. The
front activates a local fast shock only when the upstream inflow exceeds the
oblique fast-mode speed and the shared Rankine–Hugoniot solver returns an
admissible solution.

The shared manufactured test field permits testing at 1.05 solar radii, at
5 AU, and on either side of both physical boundaries without extrapolating
the wind. Those guards support geometry stencils. An open-field-line sample is
restricted to the provider's actual topology and magnetic sector. In the
detached SEP adapter, the sample behind a shock encounter cannot be requested
as a smooth transport state. The adapter exports a local downstream jump
primitive when admitted by the solver. Downstream species temperatures and
a downstream volume remain unavailable.

The earlier radial isothermal speed requires a sustaining force when applied
to the manufactured nonradial coronal tube. The [resolved thermal wind](../../src/models/sep_corona_swcme/docs/STAGE13_SPHERICAL_SHOCK_ONLY.md)
solves the field-aligned momentum equation across a continuous band of
axisymmetric tubes. Controlled tests verify the Cartesian derivatives, local
spherical shock queries and material map over that band. Event-specific wind
qualification and native field-line, source, escape and restart gates remain open.
The [native lifecycle audit](STAGE13_NATIVE_LIFECYCLE_AUDIT.md) identifies the
existing machinery for those gates. Observational scoring also requires the
case-specific model inputs and measurement operators below.

For a production replay, freeze the wind, proton scattering and source
with independent constraints. Use the full prescribed pitch-angle diffusion
law for the directional return check, including the model's 90-degree
behavior. Compare the Parker solver using the corresponding parallel diffusion
coefficient in the diffusion limit. Resolve particle escape at the physical
5 AU sphere and retain a provider guard outside it. The Ulysses deck declares
an 8 AU comparison control separately. A sensitivity comparison from 1 AU
must use the same source support in both domains. All AMPS host steps and
sampling intervals use the same approved fixed step within each run.

## Observation roles and prerequisites

The [Stage 13 data-source README](cases/stage13/README.md) records the case
sources and the tested archive routes for candidate wind realizations.

| Runtime name | Reference and comparison |
|---|---|
| `eclipse-2020-12-14` | [K-Cor polarized brightness and eclipse density analysis](cases/stage13/eclipse-2020-12-14/README.md). The available radial pB profile needs a reviewed synthetic line-of-sight operator and exposure/directional uncertainty. |
| `gle-2012-05-17` | [GOES data and two low-corona height constraints](cases/stage13/gle-2012-05-17/README.md). Resolve provenance of the particle-release clock and inference uncertainties before scoring early shock/release predictions. |
| `sep-2013-04-11` | [OV01 continuity](cases/stage13/sep-2013-04-11/README.md). Preserve the original comparison and obtain a separate dated shock and wind realization for the spherical model. |
| `sep-2001-09-24` | [SOHO/ERNE and published returning-particle constraints](cases/stage13/sep-2001-09-24/README.md). Sector responses and spacecraft magnetic-field orientation are needed for a directional comparison. |
| `sep-1978-09-23` | [Helios, IMP 8 and Voyager event](cases/stage13/sep-1978-09-23/README.md). The release archive has gaps near the SEP event peak and reservoir. Find calibrated channel-matched data for the selected upstream intervals. |
| `sep-1979-03-01` | [Voyager 2 event at 4.76 AU](cases/stage13/sep-1979-03-01/README.md). The checked Voyager flux products lack the selected event interval. |
| `ulysses-2003-oct-nov` | [COSPIN/HET proton spectra in the larger-domain control](cases/stage13/ulysses-2003-oct-nov/README.md). Combine the daily product with plasma, field, and transient-selection evidence for the selected intervals. |

The `srcSEP` regression includes the seven `IPW01`–`IPW07` numerical wind
cases, twelve controlled core cases `S13SC01`–`12`, the linked bridge case
`S13AD01`, and the seven `.in` selector and reference integrity checks.
The historical observational campaign is an explicit command because its
authentic inputs and interpretation have their own qualification.
