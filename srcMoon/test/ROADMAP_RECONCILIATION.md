# Roadmap reconciliation for the 2026-10-13 M2 work package

The required roadmaps are not byte-identical after this work package:

- source-tree roadmap SHA-256:
  `9f86272737bac602d943c92d164f4497e3d10d38722fcd58053232d0281ce703`
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
execution remain unresolved gates; therefore no I17 comparison plot or
validation claim currently exists. A second dated status paragraph under R4.4 records the partial
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

The common roadmap text controls physics requirements, D01--D15 processing,
U/I meanings, acceptance criteria, and milestone completion. The added
source-tree paragraph controls only current implementation/status reporting.
The external-data `AGENTS.md` and native dataset labels continue to control
immutability, readiness, units, projections, frames, and archive metadata.
