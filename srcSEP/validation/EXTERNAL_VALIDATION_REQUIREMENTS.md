# Inputs needed to complete external validation

Code and controlled mathematical/native tests can continue without new
observations. Their passing results do not complete the original cross-model,
SWMF or observational requirements. See [STATUS.md](../STATUS.md) for exact
build evidence and [development_requirements.json](development_requirements.json)
for the original requirement inventory.

## Inputs already located

| Input | Available evidence | Scope |
|---|---|---|
| Local MFLAMPA archives | `SWMF_data/SP/MFLAMPA/data/input/test_mflampa/MH_data_e20120123.tgz` and matching `.ref.gz` outputs; hashes in [reference_sources.json](reference_sources.json) | Real numerical inputs/reference outputs. The pinned upstream two-rank regression passes; AMPS equivalence remains unverified. |
| Pinned MFLAMPA source | Revision `38e9defc9c525bd3d49669bc2e7d8014987de488`, pinned dependencies and upstream parameter file | Reproducible comparator build/run. |
| Existing XM02/XM03 references | Published figure-derived curves and observational points with provenance | Reduced reconstructions, not numerical exports of the same initialized problem. |
| CLEAR event lists | Operational/calibrated event-summary CSVs with source locators and hashes | Candidate cohort/summary metadata; detector profiles and responses remain necessary. |
| Public AWSoM/MFLAMPA candidate | Metadata, SP input and output index for `Lulu_Zhao_071525_SH_1` | Possible comparison source. Publication equivalence is unverified; the recorded sample numeric-output request failed. |

The inventory distinguishes retained metadata from large result downloads. No
bulk observational/result archive was downloaded by this work.

## Equivalent-input AMPS versus MFLAMPA experiment

The local numerical comparator is sufficient to start preparing this
experiment. Before agreement can be scored, the two runs need a reviewed
mapping of:

- Line geometry/orientation, tube measure, observer location, frame and time
  history.
- Initial distribution and measure/units; source location, timing, species,
  spectrum and physical normalization.
- Spatial and momentum/energy boundaries and source/shock interfaces.
- Actual scattering/diffusion law, all radial/momentum dependences and regions
  where it applies. The pinned upstream nominal 0.3-AU setting scales as
  `r/AU * (p*c/1 GeV)^(1/3)` in its specified upstream region; it is not a
  constant-mean-free-path experiment.
- Plasma/field evolution, output grids, energy conventions and differential
  versus integral products with explicit unit conversions.
- Resolution, numerical uncertainty and comparison quantities declared before
  any fitted normalization or timing adjustment.

Preserve the archived MFLAMPA temperature-energy convention in keV and its
reader conversion. It differs from the live PT Kelvin interface. The local
archive has numeric snapshots through the upstream test's 600-second endpoint;
additional list/header time labels are not additional numeric states. A
current checkout does not establish the historical archive's producer mode.

## Full SWMF replay

A qualifying replay needs producer/run identity and configuration, ordered
epochs, line coordinates and stable material labels, vector B/U,
density/composition meaning, species pressure/temperature and units, and
actual grid/vertex motion. Focused transport also requires the Cartesian
velocity-gradient tensor or a provider able to evaluate it. A line's density
history does not determine that tensor.

The adopted standard BATSRUS wave contract is explicit: raw Wave1/Wave2 carry
total Alfvén-wave energy density; the immutable magnetic view is half, with
Wave1 along B and Wave2 opposite B. That contract does not supply L/R
polarization, spectral shape or composition. See
[the producer audit](physics/SWMF_FIELD_LINE_CONTRACT_AUDIT.md).

A controlled synthetic replay can test import software with independently
declared fields. Retain that label; it cannot replace actual coupled/event
validation.

## Observations and event qualification

The original OV/EV gates require time-tagged instrument flux, channel/energy
responses, differential/integral definitions, units, calibration/quality and
coverage flags, uncertainty, observation location/connectivity, and the event
inputs needed to predict the same quantities. Selection dossiers need
backgrounds, contaminating events and metric windows. Summary peaks alone
cannot determine onset, fluence, profiles, spectra or uncertainty coverage.

EV01's nine-event pilot and EV02's two holdouts do not satisfy the original
15–30-event/non-event calibration and sealed release campaign. A reviewed
cohort/split, objective/priors, response/processing contract, frozen metrics
and release/unseal provenance are needed before tuning or held-out scoring.
Already inspected/tuned events cannot become unseen by relabelling them.

Missing external evidence blocks its dependent validation claim while code
development and controlled tests continue. It must not be replaced by an
assumed source, background or synthetic observation.
