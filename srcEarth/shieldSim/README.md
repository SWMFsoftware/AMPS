# shieldSim

`shieldSim` is a compact Geant4 application for modeling energetic-particle transport through a planar shielding slab and estimating transmitted spectra and dose behind the shield.

The package version is split into normal Geant4-style source files while preserving the logic from the single-file prototype:

- CLI-selectable Geant4 physics list (`FTFP_BERT`, `FTFP_BERT_HP`, `Shielding`, `QGSP_BIC_HP`),
- source normalization using `J(E) dE`,
- `beam` and `isotropic` source modes controlled by CLI,
- local-coordinate shield rear-face scoring,
- dose accumulation in Geant4 internal units,
- detailed help text describing input/output units.

## Directory structure

```text
shieldSim/
├── CMakeLists.txt
├── README.md
├── shieldSim.cc
├── include/
│   ├── CLI.hh
│   ├── DetectorConstruction.hh
│   ├── EventAction.hh
│   ├── GCRSpectrum.hh
│   ├── OutputUtils.hh
│   ├── PrimaryGeneratorAction.hh
│   ├── RunAction.hh
│   ├── ShieldSimConfig.hh
│   ├── SpecBins.hh
│   └── SteppingAction.hh
├── src/
│   ├── CLI.cc
│   ├── DetectorConstruction.cc
│   ├── EventAction.cc
│   ├── GCRSpectrum.cc
│   ├── OutputUtils.cc
│   ├── PrimaryGeneratorAction.cc
│   ├── RunAction.cc
│   └── SteppingAction.cc
├── examples/
│   └── sep_spectrum.dat
└── macros/
    └── example_commands.sh
```

## Build

Source your Geant4 environment first, then build with CMake:

```bash
mkdir build
cd build
cmake ..
make -j
```

If CMake cannot find Geant4 automatically, provide `Geant4_DIR`, for example:

```bash
cmake .. -DGeant4_DIR=/path/to/geant4/lib/Geant4-11.*/
```

## Run examples

Single 2 mm Al shield with the built-in GCR-like spectrum and the default physics list:

```bash
./shieldSim --source-mode=beam --shield=G4_Al:2 --events=50000
```

Isotropic source over the upstream plane using high-precision neutron transport:

```bash
./shieldSim --physics-list=FTFP_BERT_HP --source-mode=isotropic \
            --shield=G4_Al:2 --events=50000
```

Tabulated SEP-like source spectrum:

```bash
./shieldSim --source-mode=isotropic --spectrum=examples/sep_spectrum.dat \
            --shield=G4_Al:2 --events=100000
```

Dose-vs-thickness sweep:

```bash
./shieldSim --sweep --source-mode=isotropic --sweep-material=G4_Al \
            --sweep-tmin=0.5 --sweep-tmax=30 --sweep-n=15 \
            --sweep-log --events=20000
```

## Physics-list selection

Use `--physics-list=<name>` to select the Geant4 reference physics list at runtime. Supported values are:

```text
FTFP_BERT
FTFP_BERT_HP
Shielding
QGSP_BIC_HP
```

The default is `FTFP_BERT`. For shielding calculations where low-energy secondary neutrons can affect dose behind the absorber, compare against one of the high-precision neutron options, especially `FTFP_BERT_HP`, `QGSP_BIC_HP`, or `Shielding`.

Example:

```bash
./shieldSim --physics-list=Shielding --source-mode=isotropic \
            --shield=G4_Al:5 --events=100000
```

## Spectrum input units

A tabulated spectrum file has three columns:

```text
E[MeV]   protonSpectrum   alphaSpectrum
```

`E` is total kinetic energy per particle. For alpha particles, `E` is total alpha kinetic energy, not MeV/nucleon.

In `beam` mode, spectrum columns are interpreted as differential pencil-beam source rates, for example:

```text
particles / s / MeV
```

In `isotropic` mode, spectrum columns are interpreted as differential intensities, for example:

```text
particles / cm2 / s / sr / MeV
```

The code applies the `pi` angular factor for isotropic plane-crossing flux.

## Output files

### `shieldSim_spectra.dat`

Tecplot-format differential spectra. Energy is MeV total kinetic energy per particle.

- `*_MC` columns: raw Monte Carlo spectra in `counts/(MeV primary)`.
- `*_Norm` columns: source-normalized spectra.

For beam mode, normalized units are `particles/s/MeV` if the input spectrum is `particles/s/MeV`.

For isotropic mode, normalized units are `particles/(cm2 s MeV)` if the input spectrum is `particles/(cm2 s sr MeV)`.

### `shieldSim_dose_sweep.dat`

Written only in sweep mode.

- `Dose_*_perPrimary` columns are `Gy/primary`.
- `DoseRate_*` columns are `Gy/s` when the source spectrum has physical units.

### `shieldSimOutput*.csv`

Geant4 analysis histograms for exit energies of protons, alphas, and neutrons.

## Notes and limitations

- The default physics list is `FTFP_BERT`, and the runtime option `--physics-list=` supports `FTFP_BERT`, `FTFP_BERT_HP`, `Shielding`, and `QGSP_BIC_HP`.
- The isotropic source uses a finite upstream plane, not an infinite half-space source. Very oblique trajectories can interact with the finite side boundaries.
- The default detector transverse size is 5 cm × 5 cm. Dose per primary depends on the scoring mass and therefore on this finite detector area.
