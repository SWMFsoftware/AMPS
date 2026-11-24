# Water fragment-ion yields at 21.2 eV

The NIST Chemistry WebBook references Kronebusch and Berkowitz's photodissociative-ionization study. Their own conference report D-1 gives three measurements from isolated water vapour at 584 Å (21.2 eV). The authors detected OH⁺ at 0.29, H⁺ at 0.05 and O⁺ at **0.0009** of total ionization. Keep that denominator. The reference counts ions and discriminates against ion kinetic energy above 1 eV. The authors regard signals for affected fragments as lower limits. An uncertainty per signal is unavailable. [NIST reference inventory](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=1DE1) [ASMS 1975 proceedings, paper D-1, printed p. 114](https://www.asms.org/docs/default-source/proceedings-archive/1975_asms_23rd_conference.pdf#page=117) [Journal paper](https://doi.org/10.1016/0020-7381(76)80088-5)

The authors link O⁺ to the molecular H₂ cofragment in this low-energy regime and report the opening of O⁺ + 2H at 23.13 eV. They did not detect the neutral in coincidence or measure a joint electron/ion/neutral distribution. The report supplies a finite-energy comparison for the yields of the three ions. It leaves the recoil, electron energy and daughter-state probabilities unresolved. Earlier shared Dibeler data at 584 Å use a different denominator (the parent ion signal).

The 712-page ASMS archive has SHA-256 `40d1c49380c06395bf22ce0a228d8ec1efe21fb110955db4f9ce5a1687352678`. Visual inspection of the scanned printed page 114 (PDF page 117) is essential for the O⁺ value. An online transcription dropped one zero. The database records the archive locator and a structured, source-attributed assessment of this single photon energy. General redistribution permission for the archived measurements remains unresolved. The journal metadata list a text-and-data-mining licence. Use `db.describe("kronebusch-water-photofragmentation")` and select file `kronebusch-water:asms1975`. Its quantity metadata preserve the three reported fractions, the denominator and the collection limit. The uncertainty field is null. Source `kronebusch-water:crossref` identifies the journal paper. The shared release exposes these records as metadata and withholds public access to the original PDF and numerical tables.

## Replay

~~~bash
python3 -B utility/chemical_database/native_m3_water_bibliography.py \
  --acquisition test_output/chemical_database/oh-ice-source-Mdgbrs/water-acquisition \
  --output-dir /tmp/amps-water-yields-bibliography
~~~

The 1976 journal record has SHA-256 `1b1428ddbb66426121c9a17bbb68619ead35d3be7439e54a4be320fdf68b7f5a`. The separate packet has SHA-256 `6cb3e5eec0d2bfd46b0db205507ef3886616c9e406dbf322f7b6dd1f6d278380`. Both are in `test_output/chemical_database/oh-ice-source-Mdgbrs/water-packet-final/`. An offline replay reproduced every packet file. Set `AMPS_CHECK_WATER_YIELDS_RELEASE=1` in the real-source suite to check the independently read values and held-source access.

The follow-up through DOI-bearing references also distinguished initial isotope and ionization mechanism. Powis and Reynolds measured D₂O, and Hult Roos et al. measured D₂O after valence and inner-valence photoionization. Mass and vibrational differences preclude transferring their quantitative recoil distributions to H₂O without an isotope model. Tan et al. derived dipole oscillator strengths with electron impact and do not supply joint recoil measurements for single-photon H₂O events. Bodi et al.'s H₂O threshold-electron coincidence measurements cover the onset of OH⁺, as recorded previously. [Powis and Reynolds](https://doi.org/10.1039/FT9918700921) [Hult Roos et al.](https://doi.org/10.1063/1.5055593) [Tan et al.](https://doi.org/10.1016/0301-0104%2878%2985080-0) [Bodi et al.](https://doi.org/10.1039/C4SC00791C)
