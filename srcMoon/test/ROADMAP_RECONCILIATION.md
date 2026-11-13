# Roadmap reconciliation through the 2026-11-12 M2 development work package

The required roadmaps are not byte-identical after this work package:

- current source-tree roadmap SHA-256:
  `2761b447b46539edb21395550cd1b5d9ff6add4839b0db5a7e6d734cf2a4f02a`
- canonical validation-data roadmap SHA-256:
  `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`

The earlier content difference was the dated Phase-2 implementation-status
paragraphs in the source-tree copy. They record implemented kernels, local
U09--U13 results, D01/D02/D03 readiness (including missing D03 package
metadata while recording that both selected days of raw RDR inputs are now
present), the now input-selected Na
accommodation law, and the remaining linked-test reasons why I14--I20 are not
yet qualified. The status now also records the implemented Paige C3--C9 and
grey-emission observation operators, finite-apparent-disc eclipse forcing,
observation-free full-column periodic spin-up, frozen I17 processor/selection
contract, and June-only production comparison/plot path. The pre-holdout
threshold, same-campaign I13 evidence, D03 package qualification, and linked
execution were unresolved in that earlier checkpoint; it established no I17
validation. A second dated status paragraph under R4.4 records the partial
uniform-global 40Ar source implementation and explicitly preserves U21 as
SKIPPED while D10, transient, and I28 prerequisites remain incomplete. No test
ID, acceptance threshold, dataset mapping, processing rule, or milestone
definition differs.

On 2026-10-14 the source R1.1 gained an explicitly user-approved operational
amendment: default I11 targets 480/240/120/60 km for the 200 GB budget, with
legacy FULL 40/20/10/5 km still selectable. I13 uses the selected finest grid.
This source amendment controls campaign selection and input derivation; the
shared numerical/statistical thresholds, test meanings, predecessor gates,
baseline review and observational qualification requirements remain binding.
The external roadmap remains unchanged (read-only for this work package).

On 2026-10-16 the source R1.1 gained a second explicitly user-approved
execution amendment: I11 completes one sample every 500 particle steps, stops
each requested mesh after at most 5000 particle steps, and writes no particle
restart files. This amendment controls execution length and checkpoint output
only. The independence, equilibrium, confidence-interval, and cross-grid
acceptance criteria remain those of the common roadmap. Ten possible samples
permit one nine-batch evaluation only when raw samples are independent;
reaching the cap without every gate is an evaluated `FAIL`, not qualification,
and cannot open I13 or observational gates. The external roadmap remains
unchanged and its hash above is current.

On 2026-11-01 the source R1.1 gained a non-qualifying physical-duration
diagnostic amendment after the retained bounded campaign established unequal
physical durations and no accepted means. It permits only the approved 480 km
member to continue beyond 5000 steps under unchanged statistical gates and
recorded resource guards, in order to propose a duration for investigator
review. This pilot does not alter the formal bounded I11 contract, cannot
report I11 PASS, and cannot apply its candidate automatically to finer grids.
Suppressing non-scoring visualization writers changes retained output policy,
with the intent of preserving production sampling. The 2026-11-11 audit found
that the common PIC cell reset was incorrectly conditional on those writer
flags; the November pilot therefore accumulated stale sampled windows. The
common reset is repaired and a production-path deterministic prerequisite is
added to I11. This correction does not change forces, sources, statistical
thresholds, the bounded 500/5000 contract, or the source-roadmap hash above.
An optional concurrent scheduler and shared initialization receipts change
execution cost only; their qualification still uses the independent existing
I11/I12/I13 scorers. The external roadmap remains unchanged.

The common roadmap text controls physics requirements, D01--D15 processing,
U/I meanings, acceptance criteria, and milestone completion. The added
source-tree paragraph controls only current implementation/status reporting.
The external-data `AGENTS.md` and native dataset labels continue to control
immutability, readiness, units, projections, frames, and archive metadata.

The November 12 user explicitly amended DEVELOPMENT order: M2 implementation,
local analytical tests, supporting production-path fixtures and data QA may
proceed before M1 closure. Formal qualification order is unchanged: I11-I13
must PASS before formal I14-I20, and observational scoring retains every
numerical/data/operator/hold-out gate. Supporting results have separate counts
and cannot install an I13 baseline or become formal M2 PASS. The source
amendment and matching AGENTS.md exception control this development decision.

The user separately approved continuous finite-inventory release, with exact
competing first-order depletion and explicitly weighted packets. The source
roadmap records this representation correction and its scope: thermal/PSD,
constant coefficients within the production time split, unchanged physical
parameters, prescribed baseline preserved, and no invented first-order law for
population-independent sputtering. The pre-amendment source hash was
`08d53e71daa8e27bde364f4d46a82323dabeaf5b518440d9072491df1d93f18a`;
the intermediate development-order-only hash was
`3b8b511aa34a52698546a9aa6d3de9b8b0054f7c4bc3fa4f84b8a27321e20a94`.
The current hash is given above. The external roadmap hash remains unchanged.
Earlier run manifests retain the contract hashes that actually governed them;
they are not retroactively rewritten.

The November 12 M2 supporting campaign now passes its scoped Na/Ar ledger,
terrain/mover, surface-distribution and U08-U13 checks; formal I14-I20 remain
SKIPPED. D03 selected-input readiness is now reproducibly READY after all 214
selected raw files match official manifests. New metadata/GCP/PRP products were
published; immutable raw files and frozen I17 selections were untouched. D02
remains PARTIAL for unresolved calculation epochs/operator requirements. These
are data/code status changes, not changes to either roadmap's scientific
acceptance contract. Historical source-roadmap status paragraphs retain their
dates; current evidence is in STATUS.md section 5.7. The external roadmap,
as distinct from newly published package metadata, remains unchanged.

On November 12 the supporting I11 event audit identified a Na production
chemistry dispatch bypass. The corrected dispatcher uses the already-declared
Moon lifetime/shadow law and is checked by actual linked sunlit/shadow survival
and exact event bookkeeping. This changes physical production behavior and
therefore requires nightly reference review; it does not add a rate coefficient
or qualify D04/I21. Both roadmap hashes above are unchanged. The pilot gains
an explicit cooperative wall-budget stop and supporting diagnostic artifacts,
while the formal 500/5000 contract, all thresholds, IDs, and I13 review remain
binding. External roadmaps/data remain untouched.
