# Pinned NIST Stage-4 reference tables

This directory contains the exact external reference bytes used by the Stage 4
tests and canonical CSV renderings made from those bytes. `raw/` is immutable
test evidence; `manifest.sha256.json` pins every raw and derived file. Normal
tests use these local files and make no network requests.

The tables were retrieved on 2026-10-17 from the official NIST Physical
Measurement Laboratory endpoints:

- PSTAR/ASTAR stopping-power and range tables:
  `https://physics.nist.gov/cgi-bin/Star/apdata.pl`
- XCOM photon cross-section/attenuation database:
  `https://physics.nist.gov/cgi-bin/Xcom/`
- X-Ray Mass Attenuation tables:
  `https://physics.nist.gov/PhysRefData/XrayMassCoef/`

The material selections are NIST material 013 (aluminum), 014 (silicon), and
276 (water). STAR energy is total kinetic energy in MeV for the named particle.
The XCOM CSVs contain the three exact Stage-4 photon nodes and both total-
attenuation definitions. The X-ray CSV additionally retains both mass
attenuation (`mu/rho`) and mass energy-absorption (`mu_en/rho`) columns,
including duplicate absorption-edge rows.

Stage 4 compares Geant4 with XCOM total `mu/rho` including coherent scattering.
The separate static X-ray tables preserve `mu_en/rho` provenance but are not
the photon acceptance oracle. The model does not substitute
`mu_en/rho` for finite-slab dose and does not claim that XCOM is an absorbed-
dose calculation. Photon transport is an attenuation verification; validated
TID remains limited to proton and helium-4 incidence.

To update the fixture deliberately, run from the repository root:

```bash
python3 scripts/download_nist_stage4.py
```

An update must be reviewed as reference-data provenance, because changed
checksums change the validation evidence even when source code is unchanged.
