# shieldSim

`shieldSim` is a compact Geant4 application for modeling energetic-particle transport through a planar shielding slab and estimating transmitted spectra and dose behind the shield.

The package version is split into normal Geant4-style source files while preserving the logic from the single-file prototype:

- CLI-selectable Geant4 physics list (`FTFP_BERT`, `FTFP_BERT_HP`, `Shielding`, `QGSP_BIC_HP`),
- CLI-selectable shield materials using either Geant4/NIST names or built-in aliases for metals, hydrogenous polymers, composites, water, and regolith,
- CLI-selectable detector/absorber/target materials for downstream scoring slabs, including tissue proxies and electronics/semiconductor media,
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
│   ├── MaterialCatalog.hh
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
│   ├── MaterialCatalog.cc
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
./shieldSim --source-mode=beam --shield=Al:2 --events=50000
```

Isotropic source over the upstream plane using high-precision neutron transport:

```bash
./shieldSim --physics-list=FTFP_BERT_HP --source-mode=isotropic \
            --shield=Al:2 --events=50000
```

Tabulated SEP-like source spectrum:

```bash
./shieldSim --source-mode=isotropic --spectrum=examples/sep_spectrum.dat \
            --shield=Al:2 --events=100000
```

Dose-vs-thickness sweep over HDPE:

```bash
./shieldSim --sweep --source-mode=isotropic --sweep-material=HDPE \
            --sweep-tmin=0.5 --sweep-tmax=30 --sweep-n=15 \
            --sweep-log --events=20000
```

## Shield-material selection

Use `--shield=<material>:<thickness_mm>` for a single run, or `--sweep-material=<material>` in sweep mode.  The material name can be either a normal Geant4/NIST material such as `G4_Al`, `G4_WATER`, or `G4_Si`, or one of the built-in keys below.

Print the full material table, including aliases, composition notes, and reference notes, with:

```bash
./shieldSim --list-materials
```

The same table is also printed inside:

```bash
./shieldSim --help
```

Built-in shielding keys:

```text
Structural Metals
  Al               Aluminum (Al), alias to G4_Al
  Cu               Copper (Cu), alias to G4_Cu
  W                Tungsten (W), alias to G4_W
  Ta               Tantalum (Ta), alias to G4_Ta

Hydrogenous Polymers
  HDPE             Custom polyethylene approximation, CH2, density 0.95 g/cm3
  BPE              Custom borated polyethylene, 95 wt% HDPE + 5 wt% B
  Kevlar           Custom aramid/Kevlar approximation, C14H10N2O2
  Kapton           Polyimide/Kapton, alias to G4_KAPTON

Composites
  CFRP             Custom 60 wt% carbon fiber + 40 wt% epoxy approximation
  SiCComposite     Custom 70 wt% SiC + 30 wt% epoxy/plastic approximation

Water / Regolith
  Water            Water, alias to G4_WATER
  LunarRegolith    Approximate bulk lunar regolith
  MarsRegolith     Approximate bulk Martian regolith
```

Examples:

```bash
./shieldSim --shield=HDPE:10 --events=50000
./shieldSim --shield=BPE:10 --physics-list=Shielding --events=50000
./shieldSim --sweep --sweep-material=CFRP --sweep-tmin=1 --sweep-tmax=20 --sweep-n=10
./shieldSim --shield=LunarRegolith:50 --source-mode=isotropic --events=100000
```

The composite and regolith materials are approximate trade-study materials defined in `src/MaterialCatalog.cc`.  Replace their density/composition there if a project-specific measured material is required.


## Detector / absorber / target material selection

Use `--scoring=<material>:<thickness_mm>,...` to define downstream scoring slabs.  The alias `--target=<material>:<thickness_mm>,...` is equivalent and is often clearer when the slabs represent tissue, detector, or electronics media.

Print the full detector/target catalog with:

```bash
./shieldSim --list-target-materials
```

Built-in detector/target keys:

```text
Crewed Mission - Tissue Targets
  Skin             Skin/epidermis proxy; use e.g. Skin:0.07 for 0.07 mm
  EyeLens          Crystalline eye-lens tissue proxy
  BFO              Blood-forming-organ soft-tissue proxy; BFO:50 for 5 cm
  CNS              Brain/CNS tissue proxy
  SoftTissue       ICRU-style 4-element soft tissue / muscle proxy

Electronics - Silicon Family
  Si               Silicon, alias to G4_Si
  SiO2             Silicon dioxide / gate-oxide proxy
  SiC              Dense stoichiometric silicon carbide

Electronics - III-V and Compound Semiconductors
  GaAs             Gallium arsenide
  InGaAs           In0.53Ga0.47As proxy
  Ge               Germanium, alias to G4_Ge
  H2O              Water reference, alias to G4_WATER
```

Examples:

```bash
# Tissue and electronics scoring behind 2 mm aluminum
./shieldSim --shield=Al:2 --target=BFO:50,Si:1 --events=100000

# Shallow skin and eye-lens scoring
./shieldSim --shield=HDPE:5 --target=Skin:0.07,EyeLens:1 --source-mode=isotropic

# Electronics detector materials
./shieldSim --shield=Al:2 --target=Si:1,SiO2:0.01,GaAs:1,InGaAs:1,Ge:1
```

Organ/tissue names define material composition only.  `shieldSim` does not apply NASA limit checks, quality factors, LET/RBE weighting, organ weighting, NIEL conversion, displacement-damage conversion, or device response functions internally.  Those quantities require post-processing using the relevant radiation-protection or electronics-damage standard.

## Material-property references and adding new materials

The material-property references used by the built-in catalogs are documented in `src/MaterialCatalog.cc` and summarized by `./shieldSim --list-materials` and `./shieldSim --list-target-materials`.  In brief:

- `Al`, `Cu`, `W`, `Ta`, `Kapton`, and `Water` are aliases to Geant4/NIST `G4_*` materials.
- `HDPE`, `BPE`, `Kevlar`, `CFRP`, `SiCComposite`, `LunarRegolith`, and `MarsRegolith` are custom ShieldSim trade-study approximations.
- Composite and regolith definitions are intentionally approximate. Replace their density and mass fractions with measured project-specific values before using them for final quantitative analysis.
- Tissue/organ target definitions are simplified material proxies. They are not anatomical phantoms and do not include biological weighting factors.
- Electronics target definitions are transport media only. NIEL, displacement damage, LET, charge collection, and device response must be applied in post-processing.

To add a new material:

1. Edit `src/MaterialCatalog.cc` and add a `MaterialCatalogEntry` in `ShieldMaterialCatalog()` for a new shield/absorber material, or in `DetectorMaterialCatalog()` for a new detector/target scoring material.
2. If the material exists in Geant4/NIST, set `canonicalName` to the corresponding `G4_*` name and set `isCustom=false`.
3. If the material is custom, implement a `Build...()` function in `src/MaterialCatalog.cc`, use explicit Geant4 units such as `g/cm3`, register the builder in `BuildCustomMaterial()`, and document the density/composition/reference above the builder.
4. Add useful aliases. The parser ignores case and punctuation, so aliases can be user-friendly names.
5. Rebuild. The `--help`, `--list-materials`, and `--list-target-materials` output is generated from the catalogs automatically.

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
            --shield=Al:5 --events=100000
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
- The built-in material catalog is intended for shielding trade studies.  Use `G4_` NIST names or replace the custom definitions if an exact alloy, polymer formulation, composite layup, or regolith composition is required.
- The isotropic source uses a finite upstream plane, not an infinite half-space source. Very oblique trajectories can interact with the finite side boundaries.
- The default detector transverse size is 5 cm × 5 cm. Dose per primary depends on the scoring mass and therefore on this finite detector area.
