# Water and hydroxyl partial photo rates and the argon photoelectron energy benchmark

## Source and admission contract

Use the existing source-packet importer and the Python/C++ catalogue. Publish
two independently selectable tables with identifiable parents and permissions.
Each exposed number must replay from a checksum-verified original. Preserve
the published number strings, column relationships and row order. The oracle
for water is printed Table 1 on page 6. The oracle for argon is the publisher's
workbook, with the units independently checked against the horizontal axis of
Extended Data Fig. 4. All entries are catalogue inputs. No event sampler,
default lock, lunar test selection, or rate law changes in this work.

- Moore and Mapaye report 1 au unattenuated rates calculated with 1 nm
  TIMED/SEE v12 irradiance and PHIDRATES cross sections across 2004–2014.
  Read their Table 1 rows 6–16 for OH and H₂O, including all three charged
  water breakup branches. The independent variable is the numerical solar
  radio flux proxy in sfu, half the sum of daily F10.7 and its 81-day average.
  For each reaction the published fit gives the frequency in s⁻¹ as
  A0 × F10.7P^A1. The power of ten applies to both the A0 mantissa and its
  quoted uncertainty. Preserve the reported A1 uncertainty separately. The
  paper reports no joint parameter covariance in this table. These
  channel-specific estimates share cross-section lineage with D04. A new
  irradiance history supplies a useful comparison, while the fitted rates
  alone do not determine daughter energies, momenta or states.
  [Moore and Mapaye](https://doi.org/10.3847/PSJ/add6a2)
- Laurell et al. report two Ar⁺ fine-structure channels contributing to their
  narrow-band photoelectron signal. Project the 800 pairs of Gaussian center
  estimates in the Extended Data Fig. 4 workbook, in laboratory electron
  kinetic energy (eV). Each row is a Bayesian fit sample from the same laser
  experiment. The harmonic-19 XUV pump and bichromatic IR probe matter
  for interpreting that energy axis. Compare the mean center separation
  against the paper's 173 ± 9 meV. Preserve the published uncertainty. The
  centers provide neither channel yields nor broadband solar weighting.
  Their Eq. 2 uses squared amplitudes of 1/3 and 2/3 in its two-state model.
  Its Schmidt coefficients and the resulting purity have different physical
  meanings from a measured branching ratio. The archive for raw experiment
  data is separately identified by DOI 10.5878/m4fc-9a58.
  [Laurell et al.](https://doi.org/10.1038/s41566-024-01607-8)
  [Swedish National Data Service](https://doi.org/10.5878/m4fc-9a58)

The software checks cover exact source and rights identities, channel order,
product labels, printed coefficients, both API languages and the reported
splitting. Changed source bytes or an incomplete table are errors. The
scientific admission rules and predecessors for I21–I24 still apply.

## Replay and source identities

```bash
python3 -B utility/chemical_database/native_m3_photon_benchmarks.py \
  --acquisition test_output/chemical_database/m3-physical-r11-fS2iFdak \
  --output-dir /tmp/chemdb-photon-benchmarks-new
```

The retained Moore PDF SHA-256 is
`00c8d5da05b0c458f98d1d28d193189d05c5d71eb85df74695d4cada96853a0c`.
The argon workbook SHA-256 is
`9685df33108f39aeeff02edb45e02a33812a1327371fd4c434c7517edfd5a516`.
The displayed values are shared-string cells in the source XLSX. Using
the numeric string indices would produce incorrect energies. The mean
separation from its paired rows is 173.3706 meV. The paired sample standard
deviation is 6.3780 meV. Use the authors' published ±9 meV when citing
their experimental result.

The rounded Table 1 coefficients for intact H₂O⁺ production give
6.46 and 7.33 × 10⁻⁷ s⁻¹ at the stated example endpoints of 130 and
163 sfu. The authors' worked example reports 6.63–7.47 × 10⁻⁷ s⁻¹.
The packet uses the printed coefficients and retains their precision.
