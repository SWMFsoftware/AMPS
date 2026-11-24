# I23 — Electron-impact ionization rate coefficients

The complete specification is [reference/acceptance.json](reference/acceptance.json).
Original acceptance: Coefficient relative error <=2%; ensemble lifetime error <=3%; explicit zero/invalid density/temperature handling. Existing U16 <=1e-12 criterion remains.

The [approved2026-11-18 addendum](../../M3_QUALIFICATION_ADDENDUM_20261118.md)
allocates relative0.003 numerical bias and fixes endpoint-survival lifetime
as primary. Existing direct/local coefficient gates are unchanged. The owner
records approval; missing physical electron-product/bath/source profiles and
actual native/predecessor evidence still block formal qualification.

This individual owner is currently a **formal capability/predecessor gate**.
It preserves source/input/executable and same-campaign predecessor evidence,
but full particle execution/scoring is unavailable while the recorded physical
kinetic policies/production event dispatch are missing. It reports SKIPPED,
never a scalar-rate PASS relabeled as chemistry qualification. The capability
register and [completion proposal](../../M3_COMPLETION_PROPOSAL.md) identify
the required implementation/data; new source recommendations are separate.

From the repository root (an existing dedicated campaign directory):

```sh
python3 -B srcMoon/test/integration-tests/I23/test.py \
  --amps "$PWD/amps" --input srcMoon/examples/m3_reaction_rates.in \
  --mpi-ranks 8 --campaign-mode regression \
  --output-dir test_output/srcMoon/M3-formal-NEW
```

The runner creates `<output>/I23/` with its real result/manifest,
stdout/stderr and an exact input snapshot when available. No `receipt.txt`
or comparison figure exists for a run that did not execute. Exit codes:
`0=PASS`, `1=FAIL`, `2=ERROR`, `77=SKIPPED`. Do not overwrite prior evidence.

The owner now retains `protocol_preparation.json`, binding the concrete
[physical model](../../M3_PHYSICAL_MODEL_SPECIFICATION.md) and
[qualification protocol](../../M3_QUALIFICATION_PROTOCOL.md). Matching
identities and a structural preparation PASS do not launch particle chemistry
or change this formal gate. Missing inputs remain SKIPPED; corrupt identities
are ERROR. Full source-family, product and predecessor requirements remain.

## Stationary helium pilot

Use a successful [atomic-family preparation](../../../../utility/chemical_database/state_profiles/ATOMIC_FAMILY_CONTRACT.md)
and its output `profile.json`.

```bash
/usr/bin/python3 -B srcMoon/test/integration-tests/I23/run_helium_pilot.py \
  --profile /path/to/profile.json --profile-sha256 PROFILE_SHA256 \
  --output-dir test_output/srcMoon/i23-helium-NEW
```

The runner compiles the production helium provider, authenticates D05 and
the physical state records, and compares the coefficient, density response,
fit slope, incident mean and thresholds at 10 eV. It writes
`comparison.csv`, the profile snapshot, logs, `run_manifest.json` and
`result.json`. The guarded composer is included in the verification.
Linked I23 still needs authenticated model-error evidence before a
helium event can commit. [Pilot contract](HELIUM_PILOT_CONTRACT.md)

Add `--kim-1992-pdf /path/to/jresv97n6p689_A1b.pdf` to run the optional
independent fit comparison. The runner verifies the exact PDF digest and
extracts its Table 1 neutral-helium row using `pdftotext`. SciPy integrates
the published fit over its finite range and records the result separately
in `kim_model_comparison.json`. It compares the predicted rate and incident
mean, with numerical quadrature errors distinct from physical data and
fit uncertainty. [Kim, 1992](https://doi.org/10.6028/jres.097.032)

Shared r11 records two independent absolute He II 30.4 nm emission
measurements at 200 eV and distinguishes them from the direct He⁺ 2s and
2p production calculations. Use the recorded photon observable, cascade
and polarization conditions when comparing those results. All-state
Maxwellian excitation and recoil budgets remain open.
[Helium evidence](../../../../utility/chemical_database/NATIVE_M3_HELIUM_EMISSION_SOURCES.md)

Turner et al.'s multi-code helium discharge benchmark assumes isotropic
scattering in the center-of-mass frame and divides the residual ionization
energy equally between the two electrons. Its purpose is to test whether
different particle codes solve the same prescribed problem. The published
WarpX binary-ionization routine also specifies the angular and energy-sharing
assumptions and includes ion recoil when solving the energy budget. Both
provide examples for event-conservation tests. The populated electronic state
of the residual He⁺ requires separate physical evidence. [Turner et al., section I](https://doi.org/10.1063/1.4775084)
[WarpX implementation](https://warpx.readthedocs.io/en/26.06/_static/doxyhtml/_split_and_scatter_func_8_h_source.html#l00504)

Peterson, Beaty and Opal measured the emitted electron's energy and angle
for incident electron energies of 100 eV to 2 keV. The incident-energy range
and single-electron observable limit its use for the controlled 10 eV
Maxwellian ensemble. Applying it there would require the collision-weighted
coverage of that energy interval and a separate residual-state assignment.
[Primary helium experiment](https://doi.org/10.1103/PhysRevA.5.712)

With `--kim-1992-pdf`, the pilot additionally evaluates the fraction of the
helium ionization rate at `Te = 10 eV` whose incident energy lies between 100 eV
and 2 keV, the limits of the Peterson et al. experiment. It uses the electron
flux times the Kim total single-ionization cross section and normalizes within
Kim's threshold-to-10-keV domain. The resulting fraction is about **0.5321%**.
The JSON records the contributions below, inside and above the window, plus an
independent whole-domain quadrature check. The eligible interval includes
unsampled laboratory energies. This result quantifies the incident-energy overlap. A residual He⁺ state
assignment still requires separate evidence.
[Kim total-cross-section fit](https://doi.org/10.6028/jres.097.032) ·
[Peterson, Beaty and Opal](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.5.712)

The r14 catalogue exposes the Opal secondary-energy marginal and the
parent-velocity heavy-ion reservoir as separate model definitions. The
approved reduced reservoir requires only a bounded heavy-particle observable.
It can account for one incident and two outgoing electrons collectively,
without reconstructing a fully differential outgoing-electron spectrum. At
10 eV its known Voronov coefficient and logarithmic slope determine a
**model** mean. Admitting the physical helium row still requires a
source/state threshold mapping, uncertainties on the temperature-dependent
shape or conditional moment (including the tail), an excited-state bound,
and a tolerable omitted-ion-impulse/mass-convention bound for the declared
observable. The linked exposure also needs the prescribed electron bath.
The 200 eV He II emission and discrete-energy electron distributions can test
limited parts of that construction. The Turner and WarpX recipes support
separate code comparisons for a chosen partition. A switch to one requires
an explicit new model choice and a state assignment. The full I23 numerical
and weighted event-budget gates remain open.
[Reduced heavy-product contract](../../M3_RESERVOIR_PRODUCTS_CONTRACT.md)
