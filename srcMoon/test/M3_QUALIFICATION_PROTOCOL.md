# M3 concrete qualification protocol and preparation runner

The user approved the numerical/statistical choices on2026-11-18. The current
effective overlay is [the approved addendum](M3_QUALIFICATION_ADDENDUM_20261118.json)
and [qualification infrastructure](qualification/README.md). The ORIGINAL JSON
below remains frozen, including its preparation-only/pending fields, so old
attempts keep their original interpretation. Current I21–I24 owners additionally
bind/snapshot the approved addendum and remove only its two obsolete approval
blockers. Missing physical providers, source/state data, numerical/ledger proofs
and formal predecessors still require SKIPPED. An approved protocol is not a
native physical execution or milestone closure.

This specification extends the existing I21–I24 contracts without assigning a
new ID or changing a frozen threshold. Its machine-readable counterpart is
[`M3_QUALIFICATION_PROTOCOL.json`](M3_QUALIFICATION_PROTOCOL.json). The approved
initial direction is a common physical-time cohort. New statistical and
numerical choices are review proposals, identified separately from approved
choices and missing authoritative physics. The analysis-before-change contract
is [`M3_INDEPENDENT_SOFTWARE_CONTRACT.md`](M3_INDEPENDENT_SOFTWARE_CONTRACT.md).

```bash
/usr/bin/python3 -B srcMoon/test/m3_protocol.py \
  --output-dir test_output/srcMoon/m3-protocol-NEW
```

The runner verifies the identities and unchanged criteria of all four per-ID
contracts, the physical specification, source roadmap and exact consumed source
files. It produces a preparation result and manifest; it launches no particle
run. Current expected result is **SKIPPED**, structural checks PASS, because
physical product policies, numerical error budgets, stream/statistical details
and formal predecessors remain incomplete. Corruption or a falsely promoted
specification is ERROR, never a qualification PASS. An unavailable referenced
source is SKIPPED; altered bytes are ERROR. Exit codes are 0/1/2/77.

For an explicit planning calculation (not authorization to run):

```bash
/usr/bin/python3 -B srcMoon/test/m3_protocol.py \
  --output-dir test_output/srcMoon/m3-protocol-plan-NEW \
  --plan-comparisons 48
```

The declared independent unit is an original parent trajectory. At exposure
`lambda*T=log(2)` the true survival probability is 1/2. For a lifetime relative
error `d`, the accepted survival interval is
`[exp(-log(2)/(1-d)), exp(-log(2)/(1+d))]`. Its smaller distance from 1/2 and
the two-sided Hoeffding bound give a conservative cohort size divisible by
eight at proposed family planning risk `0.01/K`. This differs deliberately
from an exact-binomial minimum: it is a documented sufficient planning bound.
It adds no acceptance threshold. `K` must count the entire frozen comparison
family, including layouts, weights, time steps and seeds; 48 is an example,
not a selected physical campaign. Rare molecular branches and I22's rate
estimator require their own unresolved precision design. The U17 minimum
100,000 branch events remains; 100,000 finite-time parents do not guarantee
100,000 events. Siblings/split records are correlated, not extra trials.

The proposed campaign has homogeneous tests using actual production chemistry,
then actual-mover tests at three steps and the same physical horizon, both mover
orders, equal and unequal weights, fractional corrections and explicit one/eight
rank layouts. Relative common time is distinct from frozen environmental ET.
The proposed steps are `T/8`, `T/16`, `T/32` for one frozen representable
final horizon `T`; they are not fractions of the lifetime that would require
an incompatible variable-duration last common step. The actual represented
`lambda*T` must be used when freezing the final statistical design. The
half-survival calculation above is explicitly an analytical planning example.
All actual species/block steps must match for common-time admission. Mixed
local steps and asynchronous occupancy/LOS normalization are not qualified.

Rate tables do not supply the joint product phase space. Physical state/isotope
populations, masses, formation/internal energy and bath exchanges must follow
the separately reviewed physical specification. Rounded neutral/ion equality
and zero-Q mechanical fixtures cannot qualify real ionization. Exact nucleus
and charge identities remain distinct from a physical energy-defect convention.
The complete H2O/OH family and literal U24 dependency remain required; I23's
measured-event branch retains R4.1/D06. Missing data must not be repaired by
clipping, renormalization, silent interpolation or a chosen recoil law.

Formal execution requires same-campaign I13–I20 and local per-ID predecessors,
complete source/product policies, a reviewed numerical/statistical protocol,
and an implemented production executor/scorer. Development support cannot
replace those. Development packages share one 3,600-second deadline including
setup/build/plots/restoration/cleanup and 200 GB aggregate RAM. Longer physical
campaigns belong to an explicitly requested cluster tier. Results preserve
logs, input/build/source/data identities, parent/stream/MPI layout, event/budget
records, CSV/PNG/EPS comparisons and resource history without overwriting prior
evidence. This preparation runner creates no synthetic observation, calibration
fit, hold-out result or physical PASS.
