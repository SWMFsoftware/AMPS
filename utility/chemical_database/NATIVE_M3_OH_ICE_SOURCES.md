# OH uptake on ice in the IUPAC evaluation

The shared table `iupac2010-oh-ice:uptake` reproduces the two laboratory entries and one evaluator recommendation on printed page 9075, section V.A1.3, of Crowley et al. (2010). Two entries characterize the same experiment. The third uses that experiment as its basis. The experimental surface was water ice in a coated-wall flow tube, generally at 1 Torr. Surface conditioning changed the uptake. The evaluator describes deactivation after about 20 minutes of OH exposure. Its recommendation applies to clean ice. The removed gas-phase radical can react, and the experiment did not measure the energy of a reflected radical. Lunar sticking, retained adsorbed OH and kinetic-energy accommodation still require independent inputs. [IUPAC evaluation](https://doi.org/10.5194/acp-10-9059-2010) [Original experiment](https://doi.org/10.1021/jp952142z)

## Extraction and test contract

The source packet extends the shared importer. Consumers use the ordinary Python and C++ table APIs. All probabilities are dimensionless. The experimental temperature interval is 205–230 K. The evaluation extends its recommendation to 200–230 K. The fresh-surface value of approximately 0.1 has no reported absolute uncertainty. The experimental steady-state value is 0.03 ± 0.02. The evaluator's ±0.5 reliability is in log gamma. Its meaning is distinct from an additive probability error. The two uncertainty columns retain missing cells. Extraction uses the PDF text layer on PDF page 17, keyed to the single section, its column headings, formula symbols and the original experiment. There is no numerical fit or unit conversion. The source scope is the 2010 edition.

The independent oracle is the printed experimental and preferred-value tables. The shared regression checks the temperature limits, each coefficient, its uncertainty column, and the scope metadata through the normal catalogue. A C++ query must return the same cells as Python. An altered source, rights endorsement, or incompatible layout is an error. A missing external release has SKIPPED status. Successful extraction and access provide source-level evidence for U24. I22 remains subject to its existing predecessors and physics requirements.

The publisher applies CC BY 3.0 to the article. Attribution is to J. N. Crowley et al., *Evaluated kinetic and photochemical data for atmospheric chemistry. Volume V – heterogeneous reactions on solid substrates*, Atmospheric Chemistry and Physics 10 (2010), 9059–9223. The new CSV reformats its two adjacent tables. Preserve the citation to P. L. Cooper and J. P. D. Abbatt, J. Phys. Chem. 100 (1996), 2249–2254, for the measured entries. Only the evaluated article is redistributed under its licence. [Licence](https://creativecommons.org/licenses/by/3.0/) [Version-of-record metadata](https://api.crossref.org/works/10.5194/acp-10-9059-2010)

## Reproduce and query

Run from the AMPS checkout. Use a new output directory to preserve previous packets.

~~~bash
python3 -B utility/chemical_database/native_m3_oh_ice.py \
  --acquisition test_output/chemical_database/oh-ice-source-Mdgbrs/acquisition \
  --output-dir /tmp/amps-oh-ice-replay
~~~

For a fresh acquisition, pass `--download` and a new source directory. The script checks the exact bytes of the published PDF, its printed licence and its retrieval receipt. The original packet had SHA-256 `b5ae1c8cb8136f2756296bf91761566520e58db2995a012b1c77ff64453d5e57`. A new retrieval time changes the packet receipt. The projected CSV has SHA-256 `6361b1844e3edb7980bf0f7d5b9947f541191292a1954f5666cc36caa50d9232`. The original acquisition and a byte-identical offline replay are under `test_output/chemical_database/oh-ice-source-Mdgbrs/`. Pass the packet in the normal shared database build and set `AMPS_CHECK_OH_ICE_RELEASE=1` when testing the containing release.

The KIDA entry for OH includes desorption energies for several substrates. Its astrophysical ice and olivine entries measure or estimate other properties. The primary amorphous silicate experiment derives a lower bound of 143 meV and a conditional higher value that depends on interpretation of a mass-20 desorption signal. Transfer to lunar sticking requires an incident radical experiment with known substrate, energy, and recovery of the scattered fraction. PSI's PEPISCO archive and published coincidence work also did not yield a source-qualified joint distribution of the electron, ionic and neutral fragments for all three charged water branches in this search. The water precursor, its charge state and the photon preparation must match. [KIDA OH](https://kida.astrochem-tools.org/species/38/OH) [He and Vidali](https://doi.org/10.1088/0004-637X/788/1/50) [PEPISCO](https://pepisco.psi.ch/) [Valence ion photofragmentation](https://doi.org/10.1039/C8CP03725F)
