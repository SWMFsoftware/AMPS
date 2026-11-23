# M3 physical-input follow-up — 2026-11-20

The database-first/primary-publication search is recorded in the shared
[physical-input assessment](../../utility/chemical_database/PHYSICAL_INPUT_STATUS_20261120.md).
The complete source ledgers, actual HTTP outcomes, hashes and per-lane reports
are in `test_output/chemical_database/m3-targeted-inputs-20261119T235459Z/`.
This document distinguishes new-data needs from implementation and model-choice
work; it does not alter the frozen roadmap or I21–I24 acceptance criteria.

## What the search changes

* The He/Ne/Ar atomic ground-ion thresholds already have retained isotope-
  specific primary evidence. Compatible ion physical masses are derivable.
  Complete source-bound profiles and source-effective energy crosswalks rather
  than requesting separately measured ion masses.
* Molecular dissociation/ionization/appearance energies have concrete primary
  anchors. NIST JRes121.005 supplies new native/machine-readable evaluation
  cells; the source-specific rights and original conventions remain attached.
  Existing r4 OH/OH+/H2 and CHIANTI O data supply internal-state offsets. A
  consistent eleven-channel formation-energy crosswalk is derivable once its
  explicit reference states are selected.
* H2+ is not a product of the frozen H2O/OH family. It is not a required input
  to this M3 network.
* The approved inherited-heavy-velocity reservoir does not need a full joint
  outgoing-electron distribution for signed bulk exchanges. Accuracy claims
  for recoil variance/heating/drag need additional evidence or bounds.
* CEX reference Q and rigorous all-angle error bounds are derivable. Small
  forward-scattering error and the missing D119990 experimental provenance are
  not supplied by that derivation.
* Named first-order OH and coverage-dependent water-release models exist.
  Their chosen material/domain/uncertainty and actual production integration
  remain work. A source-matched20K water sticking dataset is newly available.

## Remaining explicit choices before affected physical activation

| Choice | Concrete reference direction supported by the research | Scope/condition to retain |
|---|---|---|
| Molecular reference states | Use named H/O electronic states and lowest specified H2/OH/OH+/H2O+ reference levels, with an explicit water nuclear-spin treatment | This is a declared reduced population, not a measured nascent-state distribution; ambiguous source quantum labels are not guessed. |
| Missing printed molecular source threshold | Bind physical reference thresholds and source means through a separately named convention | Preserve D04 rates/quiet-active means; derived physical thresholds must not be falsely labelled literal PHIDRATES thresholds. |
| H2O04–06 dissociative ionization | Specify a conservative electron/ion/neutral partition with complete states and signed bath accounting | The already approved three-neutral microcanonical model is not automatically a charged/electron three-body law. Research alone does not select that partition. |
| U24 surface reference | Select one named material, coverage regime, mobility model, incoming conditions and accommodation/sticking closure | Compare20K sticking only under matching calibration/material conditions. Pair release energies/prefactors from the same final source version. Unsupported OH impacts remain explicit. |
| Collisionless versus nonzero bath | A prescribed zero-density reference gives exactly zero collisional quenching; otherwise use source-supported density/rate bounds | Preserve radiative cascades and all populated-state exits. Do not call zero-density reference evidence a measured lunar environment. |
| Approximation/error scope beyond Na | Compute source/state/speed/domain bounds and declare observable-specific permitted errors | Existing P02 numerical limits are Na-specific; source uncertainties, model error and numerical roundoff remain separate. |

These are scoped follow-up choices, not a request to reapprove P01–P09 or the
already approved numerical/statistical protocol. No choice in this table is
silently installed by the new database release. A complete measured outgoing
joint distribution is not introduced as a blanket gate for an approved reduced
reference model.

## Evidence that remains unresolved

D119990 still lacks its exact original experimental citation, initial
preparation, uncertainty and prompt-versus-eventual H1s interpretation. The
approved Na3s reference is distinct from those measurement facts. No small
source-specific angular error was established. The recovered Igenbergs source
law yields104PASS/2FAIL printed-coefficient comparisons; the two halfway
rounding discrepancies remain preserved, not widened away or repaired.

A complete matching incoming-OH sticking/accommodation law was not acquired.
Nonzero-bath electronic quenching requires applicable source/state/domain
bounds. Native download/managed-reuse holds for some primary material are
recorded separately from the existence and research use of its published facts.

## Development and formal status

Continue supported profile preparation, approved reduced-provider integration,
runner/scorer wiring, conservation and local/linked development verification.
Unresolved choices must be specified before their dependent channels activate.
Formal M3 additionally requires full frozen numerical/statistical matrices,
budgets and applicable same-campaign M1/M2 predecessors. Source/database tests
and a limited reference-model run cannot substitute for those gates.

No application physics, input default, scientific threshold, calibration/
hold-out membership or permanent U/I identifier changes in this search package.
