# Surface and relaxation native-source readers

These readers provide application-independent lexical access. They do not
select a material/state profile, enable a Moon reaction or surface law, infer
a missing uncertainty, or qualify an observational comparison. Their output
belongs to the existing shared chemical-database source-input test scope; no
new lunar U/I identifier is assigned.

## Source inventory and interpretation

* Zenodo `10.5281/zenodo.5044539`, Jones et al.2021: Apollo specimen15221,
  four native tab-separated figure tables and an explicit CC BY4.0 record.
  The two TPD files contain separate `Experimental TPD`, `Void TPD`, `Vacuum
  TPD` and `Sum` columns. Only the first is the measured signal; model
  components must remain distinguishable. Native temperature is K and TPD
  rates are ML/s. Coverage/energy and site-probability results are inferred
  under the authors' coupled desorption/transport model. Their native J/mol
  values do not establish a universal first-order residence-time law. The
  `Probability` measure is not inferred or renormalized. Signed TPD signals,
  repeated coverage and a tiny negative final coverage are preserved.
* Zenodo `10.5281/zenodo.14588634`, Li2025: six XLSX workbooks for soil
  adsorption and water migration. Every worksheet and represented cell is
  exported using the existing `native_workbooks.py` decoder. Materials include
  simulated lunar soil, loess and montmorillonite. Native units, formulas,
  calculated fields, labels and gaps control interpretation. The deposited
  files precede the final publication; figure numbering is not repaired.
* DOI `10.1038/s41557-021-00730-1`, Zhao et al.2021: three publisher source-data
  archives for OH(A)+H2 quenching. Figure2 includes reproduced experimental
  OH distributions and calculated distributions; Figure3 contains calculated
  H2 marginals; Figure4 contains wave-function densities. Quantum results use
  zero total nuclear angular momentum and declared collision energies.
  They are not Maxwellian rate coefficients, all-partial-wave cross sections,
  a complete outgoing joint law, or probabilities that can be renormalized
  after discarding the elastic/inelastic channel. Source PES symmetry labels
  must not be silently equated to physical Lambda-doublet populations.

The new packets retain explicit CC BY4.0 evidence and original source hashes.
LAMDA/BASECOL/KIDA source records acquired in the same search remain
metadata-only while managed-storage permission is unresolved. Ground-state
OH rotational relaxation and O(3P) fine-structure data do not supply
OH(A) electronic quenching or O(1D)/O(1S) relaxation.

## Reader APIs and exact commands

`new_surface_apollo15221.py` authenticates a source SHA256, enforces the
explicit native header for one of four named table kinds, and exports all
native string fields. Tab boundaries preserve explicit empty cells; changed
headers/ragged rows fail before output publication.

```sh
python3 -B utility/chemical_database/new_surface_apollo15221.py \
  --source-file NATIVE.txt --source-sha256 REVIEWED_SHA256 \
  --kind water-tpd --output-file NEW.csv
```

`new_surface_workbooks.py` calls the existing XLSX decoder. Both the native
file and the exact sibling decoder must match supplied hashes. The wrapper
and `native_workbooks.py` must be kept together when archiving the recipe.

```sh
python3 -B utility/chemical_database/new_surface_workbooks.py \
  --source-file NATIVE.xlsx --source-sha256 REVIEWED_SHA256 \
  --decoder-sha256 REVIEWED_DECODER_SHA256 --output-file NEW.cells.csv
```

`new_surface_relaxation.py` checks the complete fixed archive member set and
exports one long-cell CSV per source table. Native duplicate headers and `--`
tokens keep their positions. Figure4's original blank lines are explicit
records with column0/empty token; actual data retain columns1--6. This keeps
source row identities without erasing separators or inventing missing science
values. No extraction to source-controlled filesystem paths is performed.

```sh
python3 -B utility/chemical_database/new_surface_relaxation.py \
  --source-file SOURCE_DATA.zip --source-sha256 REVIEWED_SHA256 \
  --archive-number 5 --output-dir NEW
```

Each command refuses existing outputs. Return0 means source projection PASS;
77 means a missing required file;2 means corruption, parser or IO ERROR.
Authentication never substitutes for parsing validity. No command interprets
source numbers as validated physical coefficients.

## Verification contract

Purpose: preserve original source cells/tokens, missingness, units labels,
measured/model distinctions and immutable byte identities. Production entry
points are the three `project` functions used by packet preparation. Inputs
are explicitly authenticated source files; no random state, MPI, time/frame
transform or unit conversion is used. Independent oracles are raw source
TSV/ZIP/XML tokens and producer MD5 values, plus exact CSV round trips.

```sh
python3 -B -m unittest discover -s test/chemical_database \
  -p 'test_new_surface_*.py' -v
```

The20 source-only protocol tests exercise hash/decoder mismatches, explicit
missing prerequisites, corrupt XLSX, duplicate headings, missing tokens,
signed exponent strings, repeated coordinates, blank block separators,
ragged rows, complete archive membership, no overwrite and source-byte
preservation. Their miniature TSV/ZIP fixtures are administrative parser
fixtures, never observations. Actual-source QA is a separate explicit run
against downloaded/pinned bytes, including all23 new CSV tables.

Acceptance: exact string/row/cell identity and correct refusal/status; no
numerical tolerance. A changed source schema is an ERROR requiring review,
not permission to infer a column or weaken the parser. Source-only PASS does
not close missing sticking/accommodation, model-state/Q, quenching-domain,
uncertainty, convergence or lunar observational-validation gates.
