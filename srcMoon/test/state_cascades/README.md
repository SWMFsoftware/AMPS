# O I radiative-source acquisition for P05

`prepare_sources.py` downloads the official NIST Atomic Spectra Database O I
levels and lines, their native metadata/help, and the bibliographic records
for the selected transitions. It prepares **research evidence**. It changes
no C++ code, physical state, application input or canonical validation data.
The analysis-before-change contract is [CONTRACT.md](CONTRACT.md).

Run from the AMPS repository root, using a new output directory:

```sh
/usr/bin/python3 -B srcMoon/test/state_cascades/prepare_sources.py \
  --download \
  --output-dir test_output/srcMoon/oxygen-cascade-NEW \
  --wall-seconds 3600
```

Network access to `https://physics.nist.gov` is required only for acquisition.
There is no third-party Python dependency. To reproduce processing offline
from an existing complete acquisition packet, without touching its raw files:

```sh
/usr/bin/python3 -B srcMoon/test/state_cascades/prepare_sources.py \
  --source-packet test_output/srcMoon/oxygen-cascade-NEW \
  --output-dir test_output/srcMoon/oxygen-cascade-REPLAY \
  --wall-seconds 3600
```

The replay first checks every acquired file's SHA-256, MD5 and byte length,
and the exact request set. Both modes refuse to overwrite an existing output
directory. The one-hour ceiling includes all requests, parsing, QA and
reporting; a response exceeding the stated resource limit causes an error,
never an accepted truncated export. No application build or MPI run occurs.

The offline global-runner entry executes the maintained source/parser controls
and uses the retained research packet by default; it never makes network
requests. Override that packet explicitly when testing another acquisition:

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_oxygen_sources.py \
  --output-dir test_output/srcMoon/oxygen-source-controls-NEW
```

`--source-packet PATH` is optional for this wrapper. Missing acquired source
means SKIPPED/exit77. Executed criterion mismatches are FAIL/exit1, source or
infrastructure errors ERROR/exit2, and successful source controls PASS/exit0.
Its 15 controls include native row/component/missing-value preservation and
copied-manifest corruption/query refusal. No original raw file is modified.

## Requests and interpretation

The primary line query has no wavelength, intensity, strength or accuracy
cut. `line_out=0` retains missing A values, `output=0` requests the entire
result, and both allowed and forbidden transitions are selected. Vacuum
wavelengths in nm, level energies in cm⁻¹, and A values in s⁻¹ are requested
explicitly and checked against native headers. The complete level query is
also retained. CSV spreadsheet-text wrappers are decoded as text, never
executed; source significant figures, qualifiers and empty cells survive.

The selected states are identified by the actual configuration, term and J:
`2s2.2p4 / 1D / 2` and `2s2.2p4 / 1S / 0`. Each upper/lower level is mapped
to the level export. A second source query with maximum upper energy equal
to the exact native selected-state energy independently checks every listed
exit. This is an all-wavelength response check, not an alternate wavelength
selection. The untruncated primary export remains the source of record.

The 2026-11-18 acquisition returned ASD version **5.12**, **614 level rows**
and **910 transition rows**. Fifty-six transition rows have no A value; all
remain in the processed complete table. Eight components leave the selected
states and all eight have A values:

| Upper state | Lower state | Type | Native A, s⁻¹ | Native accuracy |
|---|---|---|---|---|
| 1D₂ | 3P₂ | M1 | 5.63e-03 | B+ |
| 1D₂ | 3P₂ | E2 | 2.11e-05 | C+ |
| 1D₂ | 3P₁ | M1 | 1.82e-03 | B+ |
| 1D₂ | 3P₁ | E2 | 3.39e-06 | C+ |
| 1D₂ | 3P₀ | E2 | 8.60e-07 | B+ |
| 1S₀ | 3P₂ | E2 | 2.42e-04 | C+ |
| 1S₀ | 3P₁ | M1 | 7.54e-02 | B+ |
| 1S₀ | 1D₂ | E2 | 1.26e+00 | B+ |

M1 and E2 components for the same state pair are distinct records and are
preserved. No one-branch assumption is imposed on 1S₀. The selected native
reference codes are T4539 and T5081:

- [Froese Fischer and Saha (1983), Phys. Rev. A 28, 3169](https://doi.org/10.1103/PhysRevA.28.3169).
- [Baluja and Zeippen (1988), J. Phys. B 21, 1455](https://doi.org/10.1088/0953-4075/21/9/007).

NIST defines B+ as an estimated accuracy within 7% and C+ within 18%. These
predominantly systematic critical-assessment estimates are **not** assigned
a Gaussian confidence level or independent-error covariance. Missing level
and wavelength uncertainties remain missing. The green line is 557.8887 nm
in this **vacuum** export; an air wavelength must not be substituted without
an explicitly sourced conversion.

The complete source also contains four 3P fine-structure descendant components:
3P₀→3P₂ E2 at 1.34e-10 s⁻¹, 3P₁→3P₂ M1 at 8.91e-05 s⁻¹ and E2 at
1.66e-11 s⁻¹, and 3P₀→3P₁ M1 at 1.75e-05 s⁻¹. Their stored graph has five
explicit levels and twelve components. They are retained instead of silently
collapsing every 3P_J product to one ground-state label. A node without a
stored exit receives no inferred physical lifetime.

## Completeness and physical status

Source-integrity QA can PASS while physical cascade qualification remains
**SKIPPED / NOT QUALIFIED**. NIST's all-lines option enumerates all *stored*
records. It is not a statement that every physically possible radiative exit
is represented. The acquired references assess M1/E2 transitions within 2p⁴;
this packet supplies no quantitative bound on omitted higher-multipole,
multi-photon or isotope/hyperfine-induced processes. The export also does
not specify isotope-resolved applicability, source product-state populations,
alignment/polarization, or collisional quenching and its omission domain.

Accordingly, `qa_report.json` gives only a **sum of listed available A values**
as a source-record diagnostic. It explicitly marks that sum as incomplete
for a total-decay claim. Total decay rate, lifetime and normalized branching
probabilities remain `null`; no reciprocal lifetime or unit branch is
generated. A pair with no stored transition is labelled as a source gap,
never assigned a zero rate. A future reviewed, explicitly bounded M1/E2-only
model would be a separate physical decision and is not activated here.

## Artifacts and status

`raw/` and `acquisition.json` preserve request URLs, original response bytes,
HTTP headers, retrieval intervals, byte lengths, SHA-256 and MD5. Processing
creates:

- `processed/oi_levels.csv` and `oi_transitions.csv`: complete native records.
- `processed/selected_exit_components.csv`: selected exact state identities,
  multipole components, A values, accuracy codes and native uncertainties.
- `processed/selected_lower_level_coverage.csv`: every catalogued lower-level
  pair, with missing-transition status preserved.
- `processed/stored_cascade_graph.json` and `stored_cascade_components.csv`:
  exact stored descendants, including the 3P fine-structure components.
- `qa_report.json`, `provenance.json`, `result.json` and `run_manifest.json`:
  interpretation, gaps, source/processor identities and actual outcomes.

The CLI returns 0 for completed source QA and 2 for an acquisition/parser/
hash/infrastructure error. A source PASS is never a radiative-cascade PASS.
All raw research packets remain immutable after acquisition; replay writes
to a new directory. No source record is changed to improve agreement with
AMPS, and no observational residual is used by this procedure.
