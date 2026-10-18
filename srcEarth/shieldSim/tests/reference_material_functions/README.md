# Material/quantity reference-data audit

This directory records what physical data exist for each approved absorber and
quantity. It deliberately distinguishes **transport can score a primitive**
from **the exact material/quantity pair is qualified for production**. The
machine-readable record is
`material_quantity_data_status.v1.json`; `material_reference_data` authenticates
that record and every candidate table during normal offline tests.

## Prepared in the 2026-10-18 audit

- Exact NIST PSTAR and ASTAR tables for natural Ge and stoichiometric SiO2 were
  added to `tests/reference_tid/nist/`. These supplement the existing Al, Si,
  and water tables. Raw upstream responses, canonical CSVs, byte counts, and
  SHA-256 digests are pinned. The data support future proton/He-4 stopping and
  range comparisons for TID/LET qualification; adding tables does not silently
  widen the production validation domain.
- An ASTM E722-19 GaAs neutron response and SR-NIEL-7 v11.1 proton/He-4 GaAs
  responses are in `tests/reference_ddd/gaas/`. The Ga and As displacement
  thresholds are 21.5 eV. The GaAs 1-MeV-neutron reference is 70 MeV-mb and is
  kept separate from silicon's 95 MeV-mb reference.
- SR-NIEL proton/He-4 responses for natural Ge (40.5 eV displacement threshold)
  are in `tests/reference_ddd/germanium/`. No neutron table or Ge `n_eq`
  denominator was found, so the package is intentionally incomplete.
- SR-NIEL proton/He-4 responses for the catalog's Ga0.47In0.53As alloy are in
  `tests/reference_ddd/ingaas/`. The 15 eV In and 10 eV Ga/As thresholds are
  taken from Beaumel et al. (2020), and the result pages use the official
  calculator's predefined In_xGa_(1-x)As compound with x=0.53. The located
  publication does not provide a downloadable neutron response or standardized
  `n_eq` denominator.

Every package is offline and checksum-pinned. Candidate status means the data
still need a generalized target-specific loader, predeclared numerical gates,
transport tests, uncertainty review, and documentation before production use.

## Why no table was invented for every cell

TID and LET are transport scorers and in principle operate in every resolved
material. Independent qualification still needs an exact reference. The NIST
STAR material list contains exact Ge and SiO2 entries, but it does not contain
the model's exact Skin, EyeLens, BFO proxy, Brain/CNS, ICRU-4 soft tissue,
6H-SiC, GaAs, or Ga0.47In0.53As definitions. A nearby tissue, another SiC
polytype, or a different III-V alloy is not silently substituted.

DDD requires target- and incident-species-specific NIEL. Tissue-composition
slabs and liquid water are not crystal/device displacement-damage targets;
their environment effects need dosimetric quantities, not a fabricated NIEL
fold. SiO2, exact 6H-SiC, and complete Ge/InGaAs neutron response definitions
remain open. In particular, an accessible calculation used beta-SiC and a
single effective 25 eV displacement energy, whereas the material catalog says
6H-SiC. That value was not transferred across polytypes.

One-MeV-neutron-equivalent fluence is not a universal renaming of DDD. ASTM
E722 requires a response for the particular semiconductor. The standard and
the official calculator supply Si and GaAs responses; no approved equivalence
was found for the other ten targets. `H100/10` is different: it is a
shield-rear-face proton spectral ratio and is absorber-independent.

## Reproduction

NIST data can be deliberately refreshed with:

```bash
python3 scripts/download_nist_stage4.py --output tests/reference_tid/nist
```

The SR-NIEL result URLs are transient. Save result pages made with the exact
manifest parameters, then canonicalize them with:

```bash
python3 scripts/generate_gaas_damage_tables.py \
  --proton-html PROTON.html --alpha-html ALPHA.html \
  --neutron-html NEUTRON.html \
  --si-neutron-groups data/niel/astm_e722_19_si_neutron_groups.csv \
  --output tests/reference_ddd/gaas

python3 scripts/generate_germanium_damage_tables.py \
  --proton-html PROTON.html --alpha-html ALPHA.html \
  --output tests/reference_ddd/germanium

python3 scripts/generate_ingaas_damage_tables.py \
  --proton-html PROTON.html --alpha-html ALPHA.html \
  --output tests/reference_ddd/ingaas
```

Run the offline integrity gate with:

```bash
cd test_output/build
ctest --output-on-failure -R '^material_reference_data$'
cd ../..
```

## Sources

- [NIST STAR material list](https://physics.nist.gov/PhysRefData/Star/Text/table2.html)
  and [PSTAR/ASTAR methods](https://physics.nist.gov/PhysRefData/Star/Text/programs.html)
- [ASTM E722-19(2026) scope and significance](https://store.astm.org/e0722-19r26.html)
- [SR-NIEL charged-particle calculator](https://www.sr-niel.asi.it/index.php/sr-niel-web-calculators/niel-calculator-for-electrons-protons-and-ions/protons-ions-niel-calculator)
  and [ASTM neutron calculator](https://www.sr-niel.asi.it/index.php/sr-niel-web-calculators/niel-calculator-for-neutrons-using-astm-standards)
- [Beaumel et al., InGaAs displacement damage](https://doi.org/10.1109/TNS.2020.2984133)
- [Iwamoto et al., semiconductor DDD methods](https://doi.org/10.1371/journal.pone.0276364)
