# Native atomic and thermochemical sources

These optional source packets provide catalogue/table access and reproducible
source preparation. Physical state profiles, reaction providers, and formal M3
qualification keep their separate admission gates. Source bytes, licences,
versions, and native missing fields remain explicit in each packet.

## NASA CEA thermochemical coefficients

The acquired official [NASA CEA repository](https://github.com/nasa/cea) is
pinned to commit `2f79a647e85737742fb08a127da00d1ea1d069da`.
The native [thermo.inp](https://raw.githubusercontent.com/nasa/cea/2f79a647e85737742fb08a127da00d1ea1d069da/data/thermo.inp)
has SHA-256 `fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363`.
Its global header date is `9/8/2021`; individual records retain their own
reference dates. The repository commit identifies the acquired version.

The repository [Apache-2.0 licence](https://raw.githubusercontent.com/nasa/cea/2f79a647e85737742fb08a127da00d1ea1d069da/LICENSE.txt)
and [NOTICE](https://raw.githubusercontent.com/nasa/cea/2f79a647e85737742fb08a127da00d1ea1d069da/NOTICE.txt)
explicitly identify the thermodynamic and transport data paths. Preserve both
documents and the original citations within each source record. This grant
applies to the published NASA product; direct acquisitions from its cited
providers retain their own permission policies.

The file contains **2,111 records**. The optional packet includes a native
nineteen-gas-species projection, table ID `nasa-selected-species`, and its
**56 interval rows** of nine coefficients, table ID
`nasa-selected-coefficients`. Original Fortran-D numeric spellings, order,
temperature bounds, phase codes, citations, and source line identities remain
unchanged. `native_cea.py` supplies `parse_thermo_bytes`,
`selected_gas_record`, `selected_species_csv`, and `coefficients_csv`;
its owned `Catalogue` returns defensive copies.

[NASA/TP-2002-211556](https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf),
Appendix C, defines the source columns and equations. The pinned
[database compiler](https://raw.githubusercontent.com/nasa/cea/2f79a647e85737742fb08a127da00d1ea1d069da/source/database_compile.f90)
also controls the current file layout: names use `A15`, and coefficient columns
33–48 are skipped with `16X`. Some current records put a printed zero in that
slot; the index preserves it without adding another polynomial term.

Temperature is in K, source molecular weight in u, and the assigned formation
enthalpy at 298.15 K and thermal increment in J/mol. NASA assigns
`H(298.15)=Hf(298.15)` and therefore
`H(0)=Hf(298.15)-[H(298.15)-H(0)]`. Assigned `H(0)` and `Hf(0)` have different
definitions. The electron-gas reference has its own thermal increment, which
must enter a balanced ionization-energy calculation. A zero thermal-increment
field is an unavailable-value marker documented by Appendix C. Keep it as
missing information.

Source molecular weights use the original elemental-weight convention. For
example, the file supplies `Na=22.9897700 u` and `Na+=22.9892214 u`; the
electron value also belongs to the source's older constants. These quantities
do not establish the fixed-isotope physical rest masses required by P02.
No uncertainty or covariance is supplied by these records or inferred from
their decimal precision. No fit evaluation or reaction/state mapping is
performed by the lexical API.

Eleven native condensed-phase intervals have nonincreasing temperature
limits, including `Br2(cr): 300.000–265.900 K`. They remain unchanged with
`native_domain_quality_status=FAIL`. The selected nineteen gas records have
increasing limits. The parser controls verify this separation; they do not
qualify the entire thermochemical database or an individual reaction model.

## ExoAtom sodium state lists

The official [ExoAtom service](https://www.exomol.com/exoatom/) and
[master file](https://www.exomol.com/exoatom/exoatom.all.json) provide a
separately published atomic-data product. The acquired master is version
`20260813`; its sodium NIST-derived datasets are version `20240901`.
The [ExoMol data licence](https://www.exomol.com/data/licence/) states
CC BY-SA 4.0. Preserve attribution, ShareAlike terms, native definitions, and
the NIST lineage described by the [ExoAtom publication](https://doi.org/10.1093/rasti/rzaf065).
The existing policy for direct NIST ASD acquisitions remains separate.

The packet retains the original neutral
[definition](https://www.exomol.com/exoatom/data/download/Na/NIST/Na__NIST.adef.json)
and [423-state list](https://www.exomol.com/exoatom/data/download/Na/NIST/Na__NIST.states),
and the ion
[definition](https://www.exomol.com/exoatom/data/download/Na/NIST/Na_p__NIST.adef.json)
and [186-state list](https://www.exomol.com/exoatom/data/download/Na/NIST/Na_p__NIST.states).
Energy and its supplied uncertainty are in cm⁻¹. Each spectrum has its own
ground-relative energy zero. The neutral ground row identifies `2p6.3s`,
`2S`, `J=0.5`; the ion ground row identifies `2s2.2p6`, `1S`, `J=0`.

Both definitions report `mass_in_Da=22.98976928` and `Ionisation=null`.
Retain those source facts: equal neutral/ion masses cannot supply a physical
Na⁺ rest mass, and the missing ionization field cannot supply an absolute
neutral/ion formation-energy offset. The `Na`/`Na_p` descriptors do not
explicitly resolve the approved `23Na` isotope mapping. Lifetimes are marked
unavailable. A relative-state catalogue alone provides no complete cascade,
photoionization energy budget, or daughter distribution. Follow the native
availability flags when reading optional columns.

## Remaining source and model decisions

The frozen D04 threshold product retains `2412.58 Å`, with its processed
conversion `5.13907097104346 eV`. A later same-day native-source check corrected
the earlier note: the solar-rate HTML supplies rate and mean excess energy,
not a separate `5.13907 eV` threshold. That truncated display must not be treated
as an independently measured or native source datum. Current
[NIST ASD output](https://physics.nist.gov/cgi-bin/ASD/ie.pl?at_num_out=on&biblio=on&e_out=0&el_name_out=on&format=1&level_out=on&seq_out=on&shells_out=on&spectra=Na&unc_out=on&units=1)
gives `5.13907696 eV` with a supplied `0.00000025 eV` uncertainty.
D04 cites [Baugh et al. (1998)](https://doi.org/10.1103/PhysRevA.58.1585);
the earlier [Ciocca et al. (1992)](https://doi.org/10.1103/PhysRevA.45.4720)
reports `41449.451(2) cm⁻¹`. No undocumented rounding interval, threshold
replacement, enlarged acceptance bound, or effective ion state is adopted
to reconcile these distinct numerical representations.

The later same-day research in repository document
`srcMoon/test/M3_PHYSICS_LITERATURE_20261119.md` establishes a source-derived
Na+ mass and proposes an explicit, energy-accounted source crosswalk. Those
are source/model admission tasks, not evidence that a separate ion-mass
measurement is required. The installed r2 snapshot is immutable and retains
its original documentation; this corrected source note belongs to a future
publication. No canonical threshold or installed file is edited.

A usable physical profile still needs an authenticated isotope/state and
formation-energy crosswalk with compatible source conventions and reviewed
error treatment. Molecular H₂O/OH products additionally need their actual
electronic/rovibrational identities and channel thresholds. The current NASA
thermochemical rows do not select those states. Direct ATcT home/API requests
in this acquisition returned HTTP403 challenge pages; they remain unavailable
acquisitions and are excluded from numerical source packets.

Run the source parser controls without rebuilding AMPS:

```sh
/usr/bin/python3 -B test/chemical_database/test_native_cea.py \
  --source "$NASA_THERMO" --output-dir test_output/NEW-native-cea-controls
```

Set `NASA_THERMO` explicitly to the acquired native file. Its pinned-source
run evaluates 32 software/source-integrity checks. Administrative fixtures
are labelled parser controls; they are not observations. Missing actual data
retains SKIPPED status, while malformed source/build infrastructure retains
ERROR. The native domain-quality failures and formal physical qualification
remain separately visible.
