# November 24 chemistry archive inputs

## Source admission and verification contract

The shared entry points are `build_database.py --source-packet` and
`source_packets.load_packets`. The published tables use the existing
Python, CLI and C++ Catalogue readers under an explicit lock. The new
immutable snapshot extends the r5 input family. It stores source material
without selecting a rate, species profile or particle interaction.

- **Zenodo 4555716.** Import the author's original ZIP and record. The
  deposit-specific CC BY 4.0 grant must bind the exact ZIP name, size and
  producer MD5 as well as the independently pinned native SHA-256. Expose
  two lossless CSV views for the `He0He1ct` branch. One indexes the 661
  center-of-mass energy files and their native header values. The other
  contains the original 768 angle, differential and weight tokens per
  file. The source differential includes the solid-angle Jacobian for
  the neutral angle defined in Appendix A. Retain all original files in
  the ZIP and identify the other branches without using them as CEX.
  [Author deposit](https://doi.org/10.5281/zenodo.4555716)
  [Appendix A](https://pmc.ncbi.nlm.nih.gov/articles/PMC8861972/)
- **IAEA D119990.** Add the current native HTML as corroborating evidence
  after exact comparison through the existing decoder. The numeric
  projection remains the frozen 12-knot table. The live container has
  its own SHA-256. Keep original citation, initial-state uncertainty
  and angle-support gaps explicit.
- **Other sources.** Retain metadata, source links and native identity
  hashes for the acquired NORAD, Curtin, Verner–Yakovlev, NIFS, NIST,
  Groningen and UPM records. NORAD's data grant is CC BY-NC 4.0 and
  does not cover general-purpose payload access. The other numerical
  inputs and article originals require a source-specific grant for this
  storage operation. Unacquired sources remain bibliography and gaps.
  [NORAD licence notice](https://norad.astronomy.osu.edu/)

Use the producer record and original ZIP as the independent reference
for archive identity, parser schema, exact lexical rows, and printed
quadrature integrals. Refuse missing or changed native bytes, an unrelated
grant, ambiguous ZIP members, malformed grids, invalid angle/weight
values, and disagreement with the printed integrals. Check that the
installed CLI and C++ SDK expose the identical index, and confirm that
held metadata do not grant payload access. These checks cover database
integrity and reader behavior. They confer no observational or linked
qualification and do not relax I21–I24 or their predecessors.

The source roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
The external canonical roadmap SHA-256 is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source implementation and resource amendments govern development
order. The common scientific and predecessor contracts control formal
qualification. This offline work requires no AMPS executable or D04–D06
processing.
