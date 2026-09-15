# OV01–OV05 observational validation implementation

## Evidence contract

All five cases are registered under `observational-validation` and execute the
binary supplied with `--amps`. Python writes a case-qualified native argument
manifest and invokes `amps --test OVxx`; the linked callback advances the
production `SEP::Transport::AdvanceParker` SDE and transactionally publishes a
CSV. Python then reads that saved output, copies the immutable observations,
scores the declared quantities, and writes CSV/JSON/JUnit plus PNG/EPS figures.
It never uses a reference curve as model output.

The inputs are fixed publication reconstructions. `--case-input` is rejected
for OV01–OV05 so a command line cannot silently change event, source, or
reference. Every plot title contains the short citation and exact figure/panel.
Every case directory includes a publication-input manifest and reference
provenance with the PDF SHA-256, extraction method, uncertainty meaning,
exclusions, reported parameters, reduced assumptions, and missing inputs.

`OV01` and `OV02` are release gates. `OV03`–`OV05` are diagnostic-only, as
specified by the campaign plan. Their agreement metrics remain in reports with
`gating=false`; this is not the same as omitting them. Reference coverage,
schema integrity, linked execution, and artifact publication always gate.

## Native observational model

`validation/cases/cross_model_validation_models.cpp` contains the common native
callback. A run declares one observer, a reviewed energy list, a source CSV,
Parker geometry, solar wind, radial/rigidity mean-free-path scaling, source
momentum index, timestep, cadence, duration, and keyed random seed. It supports:

- `time-profile`: first-passage weights are binned by arrival time for each
  instrument-channel center; and
- `integrated-spectrum`: all arrivals are accumulated at the requested energy
  centers over the declared event interval.

The source-radius mode is explicit. These observational reconstructions use
`inner-boundary`, because their papers identify low-coronal release but do not
publish an interplanetary shock trajectory on the selected line. This avoids
the physically incorrect alternative of extrapolating a constant-speed CME to
1 AU and then placing later compound-event injections at the observer.

All calculations use an absorbing 2.5-solar-radius inner boundary, first
passage at the observer, spherical-wind adiabatic cooling, an equatorial Parker
spiral, and the production spatial-diffusion provider. The implementation is a
controlled one-dimensional reconstruction, not a replay of unpublished PFSS,
MHD, shock, or detector-response state.

## OV01 — 2013 April 11 near-Earth release gate

Why: this is the campaign's best constrained near-Earth benchmark and spans
three epochs plus ACE, GOES, and SOHO energy coverage.

Reference: 80 Earth measurements from Liu et al. (2025), Figure 12(a–c), DOI
`10.3847/1538-4357/adc4e3`. OV01 deliberately reuses the single reviewed XM03
reference table and Figure 12(d) source history rather than maintaining a copy.

Input: 2.5 solar-radii injection, Earth at 1 AU, 363 km/s wind, 675 km/s event
speed, 15-minute connection delay, and
`lambda_parallel=0.3 AU (r/AU) (pc/GeV)^(1/3)`, with the paper's `p^-5` source.

Comparison: one amplitude across 4, 12, and 36 hours and every instrument; the
unpublished shock/collection area is the only normalized degree of freedom.
The formal score begins at 1 MeV because ACE/EPAM ion contamination is most
important below that energy.

## OV02 — 2020 May 29 PSP/STEREO-A radial release gate

Why: simultaneous observations at 0.33 and 0.96 AU probe radial evolution with
a 0.63 AU separation. Separate fits by spacecraft would destroy that test.

Reference: PSP/EPI-Hi 2.2 and 12.3 MeV observations from Cheng et al. Figure 3
and STEREO-A/LET 1.8–3.6 and 4.0–6.0 MeV traces from Figure 6, DOI
`10.3847/1538-4357/acac21`. Published simulation curves and the background-only
SEPT trace are excluded.

Input: the paper's 0.33/0.96 AU radii, 337 km/s event speed, 07:38–08:10 UT
low-coronal shock interval, 0.0465 AU one-GV radial mean free path, and momentum
index 6.16. The same source/coefficient set is used at both observers.

Comparison: one global logarithmic amplitude across both spacecraft and all
four channels. This retains radial intensity, dispersion, rise, and decay.

## OV03 — 2013 May 22 interacting-CME diagnostic

Why: the event tests whether a preceding CME seed episode changes the response
of the fast second CME, and exposes the limit of a one-line model at 141-degree
Earth/STEREO-A separation.

Reference: GOES-15 and STEREO-A/HET integral proton traces digitized from Ding
et al. (2014), Figure 1(b,d), DOI `10.1088/2041-8205/793/2/L35`.

Input: single-CME2 and twin-CME histories are separate immutable CSVs. They use
the reported 08:48/13:25 UT CME times, 1439 km/s CME2 speed, and 444 km/s solar
wind. Both hypotheses are compared against exactly the same spacecraft trace.

Comparison: unit-peak per integral channel preserves time shape without
claiming an absolute conversion from integral counts to a differential source.
RMSE/correlation are diagnostic, while coverage/execution gate.

## OV04 — 2014 January 6 connectivity diagnostic

Why: the behind-limb GLE is strongly connectivity- and anisotropy-sensitive.
Three fixed connection delays expose that uncertainty instead of hiding it in
a fitted best case.

Reference: PAMELA red measurement markers from Bruno et al. (2018), Figure 4
panel labelled 2014/01/06, DOI `10.3847/1538-4357/aacc26`. The blue fit is not
used as data. The digitized range is 90–1100 MeV; the printed fit parameters
(`gamma=2.14`, rollover 240.5 MeV) are preserved as provenance checks.

Input: one coronal release and early/nominal/late connection delays of 0, 30,
and 90 minutes. Every other physics/numerical parameter is unchanged.

Comparison: event-integrated first-passage spectra with one global amplitude
across all realizations. Metrics are diagnostic because a Parker line omits the
event's viewing-direction and detector-response physics.

## OV05 — September 2017 compound-event diagnostic

Why: the September 4, 6, and 10 eruptions, wide longitudinal separation, ICMEs,
and high-speed streams test compound-source persistence and expose missing
cross-field/transient physics.

Reference: STEREO-A LET/HET profiles at 4.25, 11, 31.5, and 80 MeV from Bruno
et al. (2019), Figure 2, DOI `10.1029/2018SW002085`. Dense trace overlap is
represented by larger digitization uncertainties in the reference CSV.

Input: three separately tagged pulses at the paper's flare peak times, the
paper's 450 km/s nominal Parker wind, 128-degree STEREO-A separation metadata,
and the 3163 km/s September 10 CME speed. Source placement stays at the inner
boundary because the evolving wide shock is not published as a 1-D trajectory.

Comparison: one global amplitude across four energies and the full 17-day
interval. Shape metrics are diagnostic; coverage and native execution gate.

## Commands and outputs

Run all five cases:

```sh
python3 test/run_tests.py --amps ../amps \
  --validation-case OV01 --validation-case OV02 \
  --validation-case OV03 --validation-case OV04 \
  --validation-case OV05 --output-dir test_output/OV01-OV05
```

Run source/provenance checks and, when `SEP_EXECUTABLE` is set, the linked
campaign:

```sh
make test-ov01-ov05-unit SEP_EXECUTABLE=../amps
```

The general isolated-process campaign includes these registry cases:

```sh
test/run_tests.py --amps ../amps --all --output-dir test_output/all
```

Each case directory under the output root contains resolved input, immutable
reference, native argument manifest, linked log/JSON/JUnit/model CSV, comparison
CSV, SHA-256 provenance, and PNG/EPS overlays. The aggregate report lists every
PASS, FAIL, SKIP, and ERROR and continues after individual failures.
