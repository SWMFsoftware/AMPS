# M3 photon registry and total-capture rate kernels — implementation contract

This implements the user-approved rate/registry phase only. Parent/root owns
runtime input/particle integration and global registration. Molecular/CEX daughter
kinematics, excited-state evolution, measured plasma forcing and full U17/I24
qualification remain unavailable. No reference or acceptance threshold is changed.

## Inventory before implementation

`srcMoon/species/ElectronImpact.*` is an independent native coefficient service;
there is no equivalent versioned D04 registry or lunar total-capture service.
`Moon.cpp::Chemistry::ProcessParticle` still provides guarded Na neutral loss.
The generic Huebner photochemistry arrays collapse oxygen electronic states and
sample fixed mean energies/isotropic directions. The shared charge-exchange API
supports hydrogen lifetime models, not the He/Na channels requested here. These
existing production modes and references will not be modified by this work.

D04 native package: `/data/vtenishe/moon_validation_data/photochemistry`.
README/provenance/QA and actual CSV schemas were inspected; all eleven processed
products and 33 native inputs match their archived hashes. D04 has four atomic,
seven H2O and four experimental-OH reactive channels, plus two absorption-total
catalogue rows. Rates have exactly quiet/active endpoints at optical depth zero
and 1 AU. Threshold metadata exists for the four atoms only. Missing molecular
initial/internal-state distributions and blank references are preserved as gaps.

CEX primary candidates are the retained original ORNL-6086/V1 A-32 He fit/table
and Morgan et al. (1985) Table 2 all-state Na capture values. The original PDF
bytes, formula, units, supported range and independent source evaluation must be
unambiguous before a kernel is enabled. A secondary webpage's rounded or broken
example is not the controlling law. The user approved named log-log interpolation
between positive populated Na knots; outside-range extension remains prohibited.

## Local supporting test contract

| Field | Declared contract |
|---|---|
| IDs/purpose | Supporting rate/registry subchecks within U15/U17 meanings; no new permanent ID; no claim of full product/linked qualification |
| Production entries | `Photochemistry::Registry` and `ChargeExchange::Registry` in `srcMoon/species`; PIC-independent production libraries used directly by a thin native probe |
| Configuration | Explicit frozen schema/source hashes, quiet or active solar case; CEX explicit channel/target-rest eV/u and named interpolation; no automatic production activation |
| Inputs | Complete D04 channel/rate/threshold/cross-section CSVs; original CEX reference PDFs and their verified hashes |
| Independent oracle | Archived D04 endpoint values and exact stoichiometry; direct primary ORNL polynomial/table evaluation; original Na source knots and independent log-linear interpolation |
| Units | Photon rates s^-1, distance AU, normalized visible photon fraction; cross section m2 from native cm2; CEX energy eV/u in target frame; prescribed cold-beam speed m/s and density m^-3 |
| State/frame/time | Preserve explicit O(3P), O(1D), O(1S) source labels; unknown molecular parent/product state stays unspecified. No new kinematic PDF, velocity-frame transform or epoch is inferred |
| Controls | Deterministic metadata/table/boundary fixtures; no fabricated observations or model-output reference |
| Procedure | Verify hashes/schema before transactional load; require complete unique families/endpoints, source identities, valid numeric domains and H/O/charge multiplicities; query production APIs; reject changed valid-looking rows, malformed/mixed metadata, unsupported cases and out-of-domain energies |
| Metrics | Exact channel/endpoint/product counts; rate/scaling/branch identities; primary cross-section residuals; untouched output/state on failure; source/hash refusal |
| Acceptance | Preserve existing U15 exact transcription/storage precision and declared branch contract; U17 direct-kernel relative error <=1% where applicable. Source accuracy is separate from numerical error. No lifetime/product gate is replaced by a kernel PASS |
| Expected PASS | Only the explicitly evaluated PIC-independent production rate/registry checks pass |
| Failure modes | Altered source/hash, omitted or duplicate branch, wrong atomic/state/multiplicity/units, total absorption sampled as destruction, silent solar mixing, invalid visible fraction, extrapolated capture, guessed kinetics |
| Required data/hashes | Fixed schema anchors the reviewed native CSV/PDF identities; absence is SKIPPED, contradictory bytes/schema/runtime ERROR, evaluated oracle mismatch FAIL |
| Predecessors | Formal U/I prerequisites remain parent/global-runner responsibilities; these supporting checks do not close M1/M2/M3 |
| Artifacts | Native compiler/probe command and logs, source/data hashes, result.json, comparison CSVs where evaluated; output under `test_output/srcMoon` or an explicit owned output directory |
| Status semantics | Local PASS is rate/registry verification only; full molecular/CEX transport remains OFF/NOT QUALIFIED; missing authoritative kinetics is never converted to PASS |

All pure-layer failing calls preserve output arguments and any previously loaded
registry. A visible-photon fraction is a supplied normalized gas operator result,
not a surface incidence cosine. CEX kernels accept source-frame energy directly;
they do not silently choose a projectile mass, a nonrelativistic velocity mapping,
a Maxwellian tail treatment, a solar-wind bath or product recoil.
