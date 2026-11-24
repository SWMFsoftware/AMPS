# Helium ion line emission following electron impact

## Source-admission contract

The source-packet importer and the existing catalogue `describe` operation
publish two individually reported experimental observables and one theoretical
comparison as source metadata. The measurements use a neutral helium target
and a monoenergetic 200 eV incident electron beam. The observable is the
absolute cross section for production of a photon on the He II 2p → 1s line.
Neither paper supplies a 10 eV Maxwellian excited-ion fraction or a recoil
distribution. The independent oracles are the source authors' numerical
statements and the two corresponding entries in an independent assessment.
An absent, changed or incorrectly identified acquisition must stop preparation.
Both experiment records retain the reported absolute error in native units
and the method that established their scale. There is no stated confidence
level for these absolute errors. The shared source-access policy holds the
article bytes and full data tables pending permission for managed distribution.

| Publication | Reviewed physical input and limits |
|---|---|
| [Bloemen et al., 1981](https://doi.org/10.1088/0022-3700/14/4/019) | A modified absolute photoionization chamber gives a 30.4 nm cross section of (5.6 ± 1.4) × 10⁻¹⁹ cm² at 200 eV. Section 4 discusses contamination by adjacent radiation and the 25% overall error. The separate 200 keV neon-ion result in this paper has a different projectile. |
| [Forand, Becker and McConkey, 1985](https://doi.org/10.1088/0022-3700/18/7/020) | Relative calibration against He I emissions, with wavelength-dependent filter transmission and detector quantum efficiency, gives (6.2 ± 1.7) × 10⁻¹⁹ cm² at 200 eV for the 304 Å line. Table 2 contains 20 emission cross sections at incident energies 80–1000 eV normalized at 200 eV. Section 3.1.2 documents its uncertainty. Cascade feeding and emission polarization are uncorrected. The source's 65.4 eV calibration threshold is documented in Fig. 4. |
| [Zatsarinny and Bartschat, 2014](https://doi.org/10.1088/0953-4075/47/6/061001) | Separate angle- and energy-integrated B-spline R-matrix predictions for residual He⁺ 2p (Fig. 1) and 2s (Fig. 2). The authors' comparison supports the original Bloemen and Forand absolute scale. They report no 2s experimental benchmark. Agreement of two numerical methods at spot checks does not bound unmodelled electronic states or recoil. |

The NIST-hosted critical assessment independently lists both He II emission
cross sections at 200 eV and the two experimental normalization methods.
Its values appear in its own native scale of 10⁻¹⁸ cm². The assessment is a
cross-check of the original publications. [Table 8](https://www.nist.gov/system/files/documents/srd/jpcrd373.pdf)

## Reproduction and access

`native_m3_helium_emission.py` consumes the retained three publisher PDFs
and their HTTPS retrieval manifest. It checks the independently pinned file
identities, reads the original text layers and records short authored
descriptions of the two experimental results. The printed exponents and
uncertainties were checked against the original abstracts. Photon yields
must retain their beam energy and normalization when used for comparison.
A probability per ionization requires additional data on the total rate,
other residual states and cascades. Physical I23 remains SKIPPED.

```bash
case_dir=test_output/chemical_database/m3-helium-emission-kt4oVr
python3 -B utility/chemical_database/native_m3_helium_emission.py --acquisition "$case_dir"
python3 -B utility/chemical_database/native_m3_helium_emission.py --acquisition "$case_dir" --check
chemdb --lock "$CHEM_DB_LOCK" describe he-electron-ion-emission-20261124
```

The first command creates a new packet file. The second compares an offline
reduction against its bytes. The retained manifest contains the exact UTC,
URL, byte count and SHA-256 of each acquisition. An IOP download can receive
a fresh time-stamped cover page. A changed digest needs re-review against
the original scientific pages before a new edition can be used.

| Held scientific original | SHA-256 |
|---|---|
| Bloemen et al. PDF | `f368830ffbb532f1ca4dbe05876e25516aae691ae20823f033100cf1d8e769cd` |
| Forand et al. PDF | `4f3d9e3465cef11a7ec49fe532472222c6c5de36a3f90c941884f03273502c75` |
| Zatsarinny and Bartschat PDF | `c9c5a3efa5bc7e627585431e64457d5bdd9a76f5b96331c6481e68923d25acd6` |

The reviewed three-record packet has SHA-256
`f29732958978c4a26cf3f577b4ec462448a1fc2b7a67cbc435a000d66fb6ec45`.
Enable `AMPS_CHECK_M3_HELIUM_EMISSION=1` with the full shared suite and
the selected release lock to check both measurements, the curve domain,
the separate theoretical source and managed article-access restrictions.
