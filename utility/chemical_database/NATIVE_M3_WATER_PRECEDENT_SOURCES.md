# How other models handle water photofragment recoil

The optional dataset `water-photofragment-model-precedents-20261124` records two different modelling objectives. EXOSPHID follows fragments in a lunar exosphere. Shinnaka et al. modelled neutral hydrogen to interpret three Lyman-α images of comet 67P. Both supply explicit assumptions applicable to their respective observables. The ions, electron and neutral born in a *single* water dissociative photoionization event are the missing joint measurement for I22.

## EXOSPHID lunar particle model

At revision `62e670f9d0344ced1c1414a86b40fbba43b373ef`, the Julia water file defines the following inputs for its shortest wavelength bin, 0–662 Å. Each fraction is conditional on the model selecting this bin. The energy column is the total appearance threshold used in the program. The thesis credits Kronebusch and Berkowitz and explicitly simplifies wavelength intervals and channel weights. The shared database retains those numbers as model inputs, with no uncertainty assigned. [Pinned water parameters](https://github.com/charligolea/EXOSPHID/blob/62e670f9d0344ced1c1414a86b40fbba43b373ef/src/database/species/H2O.jl#L24-L105) · [TUM thesis, §§3.2.1 and 4.2](https://mediatum.ub.tum.de/doc/1784196/1784196.pdf)

| Water products | Threshold (eV) | Shortest-bin fraction |
|---|---:|---:|
| H₂O⁺ + e⁻ | 12.60 | 0.63 |
| H + OH⁺ + e⁻ | 18.11 | 0.31 |
| H₂ + O⁺ + e⁻ | 18.65 | 0.043 |
| H⁺ + OH + e⁻ | 18.72 | 0.017 |

The `DiPI_logic` implementation first creates an intermediate H₂O⁺ velocity from the incoming photon, parent velocity and a sampled electron direction. It then passes that intermediate velocity and the **same photon energy** to a separate dissociation calculation using the total appearance threshold. Both primitive functions add the incident photon momentum. The public event returns the two heavy-fragment velocities. The electron velocity used to obtain the intermediate recoil is discarded. Reconstructing that electron and combining it with the two heavy products would count the incident photon momentum twice under the primitive functions' own mass convention. The thesis describes isotropic emission. The pinned implementation folds sampled trial directions into the forward photon hemisphere. An application requiring electron-resolved transport would have to establish a single conserved photon budget and a measured or independently bounded energy-sharing distribution. [Dispatch](https://github.com/charligolea/EXOSPHID/blob/62e670f9d0344ced1c1414a86b40fbba43b373ef/src/photoreactions/photodestruction_logic.jl#L229-L249) · [Ionization](https://github.com/charligolea/EXOSPHID/blob/62e670f9d0344ced1c1414a86b40fbba43b373ef/src/photoreactions/SimplePhotoionisation.jl#L43-L108) · [Dissociation](https://github.com/charligolea/EXOSPHID/blob/62e670f9d0344ced1c1414a86b40fbba43b373ef/src/photoreactions/SimplePhotodissociation.jl#L79-L156)

## Comet 67P comparison

Shinnaka et al. Table 2 assigns the neutral H and OH produced at 1357–1860 Å speeds of 17.5 and 1.05 km/s with branch fraction 0.670. The aggregate ionization branch at wavelengths below 984 Å has fraction 0.059 and no product speeds in either column. The images constrain the extended neutral H coma. They do not resolve charged daughter energy sharing. The optional catalogue records one neutral row and the aggregate ionization row from the native PDF. [Shinnaka et al. 2017](https://doi.org/10.3847/1538-3881/153/2/76)

## Reproduce and query the source review

The recipe verifies SHA-256 for all six original files, reads the Julia declarations and the original PDF text layer, and writes a source packet accepted by the shared database's normal importer. Source IDs start with `water-model:`. Each entry has `physical_activation=false`. The catalogue publishes the facts and source identities under `metadata-only`. The online step requires curl and the replay requires Python and `pdftotext`. Both commands use new output directories.

```bash
work=$(mktemp -d)
python3 -B utility/chemical_database/native_m3_water_precedent.py \
  --fetch --acquisition "$work/downloads" --output-dir "$work/packet"
python3 -B utility/chemical_database/native_m3_water_precedent.py \
  --acquisition "$work/downloads" --output-dir "$work/replayed"
cmp "$work/packet/source_packet.json" "$work/replayed/source_packet.json"
```

The original release preparation retains the acquisition and manifest under `test_output/chemical_database/water-precedent-Pcq5D0ih/retrieved` and `packet/` in that same parent directory. A newly downloaded publisher PDF with a different cover will require a new edition review. Use the retained acquisition to reproduce the current release. Pass the manifest as an explicit source packet when building a new database release. In the shared regression suite, `AMPS_CHECK_WATER_PRECEDENT_RELEASE=1` selects the catalogue checks for a release containing this dataset.

```python
dataset = database.describe("water-photofragment-model-precedents-20261124")
facts = {entry["id"]: entry["quantity_metadata"] for entry in dataset["files"]}
print(facts["water-model:H2O.jl"]["branches"])
print(facts["water-model:photodestruction_logic.jl"]["returned_electron"])
print(facts["water-model:shinnaka2017.pdf"]["ionization"])
```
