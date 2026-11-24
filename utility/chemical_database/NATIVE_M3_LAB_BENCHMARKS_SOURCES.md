# OH surface loss and helium excitation source follow-up

This package continues the database-first searches for I22–I24.
The AMPS source revision is `e5f2978ac0576f42dc510e0da110dcbe0029d69a`.
The source roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`
and the canonical roadmap SHA-256 is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The AMPS roadmap adds the approved early implementation and resource
amendments. Both retain the same formal M3 qualification requirements.
The validation-data rules control how native data and missing values
are interpreted.

## Publication contract

Use the existing source-packet importer, immutable publisher and shared
read-only Catalogue APIs. The new release reuses all r9 source packets
and model definitions. Admit two exact original publications with pinned
source-specific permissions. Their projections contain the eight OH
conditions in Wang and Bruggeman Tables 2–3 and the three experimental
ratios in Sutton and Kay Table III. Use the published cell text as the
oracle, including blank entries, reported uncertainties and native
scale factors. Bind recipe bytes and receipt identities. Test every row
against independently read primary tables and compare C++ and Python
views. Register the relevant unresolved archive leads as metadata.
Application reaction selection, qualification criteria and native
validation archives are outside the publication contract.

## Admitted numerical sources

| Source | Scope | Native SHA-256 |
|---|---|---|
| [Wang and Bruggeman, Plasma Sources Sci. Technol. 32, 085016](https://doi.org/10.1088/1361-6595/aceff9) | Table 2 gives six dry/humid He conditions with OH diffusive flux in units of `10^14 cm^-2 s^-1`, thermal flux in `10^16 cm^-2 s^-1`, and dimensionless loss probability β. Table 3 gives two spatial resolutions at the same dry location and also reports density in `10^12 cm^-3`. The article grants CC BY 4.0 for its original content. | `d64a3b7d454b79efb389f812116c387bd62a89429a7a5d990d25518aeaf540ac` |
| [Crossref publisher record](https://api.crossref.org/works/10.1088/1361-6595/aceff9) | Confirms the DOI and the CC BY 4.0 licence of the version of record. Credit Jianan Wang and Peter J. Bruggeman, the article title, journal, and DOI. | `120609b87f3ec8c472fe7201d9124a01445c3ab3a43521a55b4b8d5a0ce4ba76` |
| [Sutton and Kay, NASA-TM-X-66137](https://ntrs.nasa.gov/citations/19730006011) | Table III, printed p. 22. Gas-phase neutral He bombarded by 200 eV electrons. Three measured ratios relative to the He⁺ 4s cross section, with reported uncertainty limits. This preprint's comparison columns cite earlier theory. | `ce5eb4aba0ce7512b81d87767b7748844ef1a352417799c6ceb6f29676369db9` |
| [Official NASA record](https://ntrs.nasa.gov/api/citations/19730006011) | Public distribution, government public use permitted. Credit J. F. Sutton and R. B. Kay, NASA Goddard Space Flight Center. | `748c46e84cd6d8e4f534f533495d5a089752a4e97adda1a2372a6d098120e335` |

The quartz source derives overall OH loss from diffusive and
thermal flux at an atmospheric-pressure helium plasma jet on fused
silica. Its three dry-position values range from 0.010 to 0.011.
The native humid 1.0 mm diffusive-flux and loss entries are absent.
The two Table 3 resolutions use the same experiment. The reported β
changes from 0.010 to 0.023 as the resolution is coarsened from
10 µm to 50 µm. The study constrains reactive loss under those
specified conditions. Incident sticking, retained OH, residence, and
energy accommodation in lunar vacuum require separate evidence.

The Sutton/Kay preprint differs from the 1972 NASA single-author
report already in r8. Its Table III gives `16.00 +10.0/-14.0`,
`1.10 +00.4/-00.6`, and `0.74 ±00.2` for 4p, 4d and 4f relative
to 4s. Keep the two reports separately identified. They concern the
same experimental programme. The ion ground state, total excited
yield and ion recoil are outside the three reported observables.
The imported ratios are the experimental column. Table I combines
earlier measurements and theory and supplies no new independent
Sutton/Kay measurements over its energy range.

## Remaining evidence

- The [Curtin CCC archive](https://atom.curtin.edu.au/ccc-db/index.html)
  supplies integral helium electron-impact data. Its current He I
  archive was examined in r6. Bellm et al. measure selected
  electron-angle and energy differential ratios for ionization with
  He⁺ left in `n=1–4`. Dogan et al. report a line-emission cross
  section from threshold to 300 eV. Dürr et al. give a 102 eV and
  1 keV coincidence benchmark. Public numerical arrays and permissions
  for a differential input covering the 10 eV Maxwellian case were
  not established. [Bellm et al.](https://doi.org/10.1103/PhysRevA.75.042704)
  [Dogan et al.](https://doi.org/10.1088/0953-4075/31/7/023)
  [Dürr et al.](https://doi.org/10.1103/PhysRevLett.96.243202)
- Lindle et al. supply He⁺ `n=2` photon partials, satellite ratios and
  photoelectron asymmetry in the 67.5–90 eV band. The Berkeley
  download attempts returned empty HTTP 202 responses. Positive
  retrieval and resolved rights would supply a finite-band anchor.
  Full He/Ne/Ar solar-weighted state probabilities are still open.
  [Lindle et al.](https://doi.org/10.1103/PhysRevA.31.714)
- The NASA water PEPICO record confirms source- and timescale-specific
  parent survival in H₂O and D₂O and is metadata-only. The RSC
  Suárez supplementary request again returned HTTP 403. Its source
  describes a preselected cation-state wave packet. Neither record
  supplies a joint electron/ion/neutral law for the three D04 charged
  water channels. [NASA water PEPICO](https://ntrs.nasa.gov/citations/19920033386)
  [Suárez et al.](https://doi.org/10.1039/C8CP03725F)
- I24 retains the D119990 rate-state discrepancy and angle-data gap.
  The Le et al. calculation covers 1.7–40 eV COM endothermic capture
  into excited H, outside the 4 keV/u ground-product channel. The
  Shingal and Fritsch papers identify further state-specific theory.
  This search did not recover a matching source-qualified numerical
  rate and differential table. [Le et al.](https://doi.org/10.1103/PhysRevA.68.012705)
  [Shingal et al.](https://doi.org/10.1088/0022-3700/19/3/014)
  [I24 archive record](https://db-amdis.org/collisiondb/datasets/119990/)

## Replay

The exact acquisition manifests and primary files are in
`test_output/chemical_database/m3-source-r10-CMOILD/`. The reduction
uses the PDF text layer and the visually checked printed tables. It
does not digitize a plotted curve, infer a missing measurement, assign
a confidence level, fit a model or convert native units.

```bash
python3 -B utility/chemical_database/native_m3_lab_benchmarks.py \
  --acquisition test_output/chemical_database/m3-source-r10-CMOILD \
  --output-dir /tmp/amps-m3-r10-new-packet
```

The final author-extracted views are shared reference inputs. Formal
I21–I24 remain SKIPPED under the existing gates.
