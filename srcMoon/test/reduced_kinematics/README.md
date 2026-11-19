# Source-explicit reduced product kinematics

This component verifies approved P02–P04 numerical services in
`srcMoon/species/ReducedKinematics.{h,cpp}`. Read [CONTRACT.md](CONTRACT.md).
The kernels require explicit masses, formation/excitation values, source
identities, inertial frame and numerical/model-error bounds. There are no
embedded isotope masses, reaction energies, internal-state mappings, random
draws, PIC data or automatic physical-provider activation.

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_reduced_kinematics.py \
  --output-dir test_output/srcMoon/reduced-kinematics-NEW
```

Use a new directory. `--cxx` selects a C++17 compiler. This command compiles
only the maintained production kernel and thin adapter; it never configures
the application or edits `build/`. It writes source/compiler/executable
identity, `compile.log`, native stdout/stderr, named `result.json`, detailed
`comparison.csv`, ensemble CDF CSV and PNG/EPS figures. Python3 and a C++17
compiler are required; absent optional Matplotlib skips figures separately.
Exit codes:0 PASS,1 evaluated FAIL,2 infrastructure ERROR,77 missing compiler.

`ResolveState` uses supplied additive neutral-isotope/electron reference
masses and copies source provenance into an immutable state. `FormationQ`
counts the explicit incoming/outgoing formation and excitation occurrences.
Separate signed components are retained through a floating expansion and
copied into the receipt; a small excitation must survive large-formation
cancellation. Accepted near-unit vectors are normalized after checking their
explicit bound in[0,1), and zero vectors are always refused. Tiny open polar
uniforms retain their representable transverse recoil using a stable sine.
`EvaluateModelErrors` separately reports mass, force-acceleration and classical
kinetic errors against supplied limits; it selects no physical error budget.
The reference inventory identity must agree across an event. A later normal
AMPS build must use those same model masses/charges in the mover and ledger.

`GeneralTwoBodyFinalState` takes explicit total momentum, available COM
translation and relative direction. `PhotonTwoBody` calls that same solver
after including photon momentum and explicit state/threshold energy. It is
parent-rest inertial only: a nonzero parent momentum is refused. This is not
yet a moving-parent/SPICE-LSO event provider. Threshold and mean energy remain
separate in its receipt so cancellation cannot erase a supplied small excess.
No missing internal partition is inferred, and negative COM energy is refused.

`ThreeNeutralMicrocanonical` projects supplied `sqrt(m_i)*g_i` momenta into
the zero-total-momentum subspace and fixes the COM energy. Finite normals are
supplied by the caller; zero energy is deterministic, and an invalid positive-
energy draw is refused without resampling. All three input occurrences are
preserved, including equal H states. It supports the approved reduced neutral
phase-space model; it supplies no dissociative-ionization or measured breakup
distribution. Its heaviest-mass dependent coordinate avoids cancellation and
does not choose a physical partition parameter.

Independent controls use high-precision Decimal two-body formulas and a
Jacobi-coordinate construction, plus a fixed-seed40000-event Beta(3/2,3/2)
ensemble oracle. Native normal generation belongs to the adapter, never the
production service. All convenient kg/J fixtures are **analytical controls,
not isotope measurements or lunar product energies**. The17-statistic99%
simultaneous bound is frozen in the contract. Refusal cases check unchanged
outputs and source/state guards. Local PASS neither sets
`KinematicPolicy::qualified` nor establishes physical M3 qualification; these
approved reduced models remain **EXPERIMENTAL** pending data and linked gates.
