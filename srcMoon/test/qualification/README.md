# M3 approved qualification infrastructure

[CONTRACT.md](CONTRACT.md) records the implementation inventory and test
contract. [The approved addendum](../M3_QUALIFICATION_ADDENDUM_20261118.json)
preserves the original frozen protocol/criteria and records the user's
2026-11-18 choices. It changes no previous result or observational membership.

The production `species/QualificationRandom` is a stateless HMAC-SHA256
counter stream, using the maintained shared DataIO hash and the construction
defined in [RFC2104](https://www.rfc-editor.org/rfc/rfc2104.html). Its public
seed is a reproducibility key, not an authentication secret. Explicit cohort,
original-parent, genealogy, purpose and draw-index fields give the same history
under rank/reordering changes; physical daughters get distinct genealogy.
An exact52-bit midpoint gives strict open uniforms. This is a finite-grid
construction/statistical model; finite optical-depth tails and maximum draw
family remain explicit numerical/statistical evidence requirements. No mutable
generator, clock advance, PIC record/address or allocation-order seed is used.

`m3_qualification.py` implements independent endpoint-survival rate/lifetime
planning and scoring, exact conditional branch intervals, rare-branch precision
planning and numerical-bias checks. I22 uses total-removal `-log(S/N)/T`; branch
frequencies/partial rates and the event-time `D/exposure` estimator remain
separately named diagnostics. It never changes estimator after outputs.
Every complete source branch is retained, including source-zero branches.
Original parents have equal initial represented number within each binomial
cohort; daughters/splits/repeated samples are not extra trials. Unsupported
censoring or unequal initial weights are refused. Exact Clopper–Pearson/Binomial
planning requires SciPy; absence is SKIPPED, never a Gaussian substitute.

The readiness boundary is explicit:

```sh
/usr/bin/python3 -B srcMoon/test/m3_qualification_executor.py \
  --case-plan /path/to/frozen-cohort-plan.json \
  --campaign-root /path/to/same-campaign-predecessors \
  --output-dir test_output/srcMoon/m3-qualification-NEW
```

The plan must bind its approved protocol, complete source-channel list,
actual represented horizon, original cohort/seed and source-state profile.
The runner verifies actual artifacts/hashes and predecessor statuses, writes
result/manifest/input-plan/stdout/stderr and reports SKIPPED while complete
physical cohort providers/executors/numerical/ledger proofs are unavailable.
It invents no native CLI/receipt/figure and accepts no Boolean shortcut. Even
all-green infrastructure predecessor fixtures cannot create a native provider.
I21–I24 and full U24 remain unqualified when their required scope is missing.

Local supporting controls run independently of an AMPS rebuild:

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_qualification.py \
  --output-dir test_output/srcMoon/m3-qualification-controls-NEW
```

They compile the actual RNG/DataIO, compare independent Python HMAC vectors,
check strict-open endpoints, semantic key/shard identities and fixed statistical
controls, then verify the maintained independent scorer/planner on explicitly
analytical fixtures. These fixtures are not observations or source profiles.
CSV/PNG/EPS comparisons, compiler/source/executable identities and all outcomes
remain in the owned new output directory. Native application campaigns require
coordinated source freeze and one3600s/200GB package (larger cohorts need an
explicit cluster tier); these supporting controls do not configure AMPS.

The source-installed I21 adapter now connects this scorer to the actual normal
cohort owner. It still requires complete authoritative state/Q inputs, a fresh
source-matched compiled producer and frozen normal `.in` include identities.
`m3_qualification_executor.native_request` validates the concrete adapter schema;
`sodium_cohort/audit.py` reconciles every removal with per-rank committed event
records before independent ledger and first-hit audits. Three increasing common
step levels at both1/8MPI are mandatory. The whole actual attempt uses one
original resource clock and retains attempted-launch journals on failure.
The partial current Na profile remains SKIPPED; I22–I24 have no admitted complete
provider. Full matrix/milestone qualification is distinct from one cohort result.

## Frozen chemistry family

Use [run_m3_campaign.py](../run_m3_campaign.py) to run a whole declared set of
development cases from one command. The existing I21 and I24 owners retain
their native solvers and independent references. The new family runner pins
the selected inputs before launch, supervises all children under a single
budget, and indexes the independent checks in one CSV. The I21 auditor now
also compares the actual original-parent fate and reaction channel across
one and eight MPI ranks at each identical timestep. It reports the largest
rank-to-rank first-hit time difference in seconds alongside the already
approved independent timing criterion. [Campaign contract](CAMPAIGN_CONTRACT.md)

| Plan kind | Required input roles | What runs |
| --- | --- | --- |
| `i21-cohort` | `plan` | One existing frozen stationary, shadowed or moving sodium plan from `integration-tests/I21/prepare_cohort.py`. Both MPI layouts and all declared refinements are run. |
| `i24-sodium-ensemble` | `lock`, `mass_inventory` | The r14 source-bound 32,768-parent rate, survival and daughter event-plan comparison at one and eight ranks. |
| `i24-sodium-pic` | `lock`, `mass_inventory`, `producer` | The actual r14 Na/Na⁺/H SPICE producer with both movers and the fixed-horizon timestep and fixed-population weight matrix. |
| `i23-helium-ensemble` | `profile` | The approved controlled 10 eV electron-bath ⁴He ensemble (32,768 parents at 0, 1e6 and 1e7 m⁻³, 1 and 8 ranks) through the production provider and reviewed reservoir composer. |
| `i23-helium-pic` | `profile`, `producer` | The Na/He/He⁺ SPICE producer with both movers, the zero-density control and the fixed-horizon timestep and fixed-population weight matrix. |

Each pin is an object with `path` and `sha256`. Relative paths are resolved
against the plan's directory. For a typical family, run from the repository
root with paths to the three prepared I21 plans, the r14 database lock and
its NIST mass inventory, and the three-species producer receipt. See the
individual [I21 preparation](../integration-tests/I21/README.md) and
[I24 build](../source_qualification_build/README.md#build-the-i24-sodium-capture-case)
instructions. Use source-matched SPICE builds in disposable checkouts. A
build manifest binds the executable and library to the current source.

```bash
plan=/absolute/path/to/new-family.json
/usr/bin/python3 -B - "$plan" "$na_lit" "$na_dark" "$na_moving" "$db_lock" "$nist_mass" "$pic_producer" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

output, lit, dark, moving, lock, mass, producer = map(Path, sys.argv[1:])
def pin(path):
    path = path.resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

def case(name, kind, **inputs):
    return {"name": name, "kind": kind, "inputs": {key: pin(path) for key, path in inputs.items()}}

protocol = Path("srcMoon/test/M3_QUALIFICATION_ADDENDUM_20261118.json")
with output.open("x") as stream:
    json.dump({
        "schema": "srcMoon-M3-development-family-v1",
        "approved_addendum_sha256": pin(protocol)["sha256"],
        "cases": [
            case("na-illuminated", "i21-cohort", plan=lit),
            case("na-shadowed", "i21-cohort", plan=dark),
            case("na-moving", "i21-cohort", plan=moving),
            case("capture-ensemble", "i24-sodium-ensemble", lock=lock, mass_inventory=mass),
            case("capture-pic", "i24-sodium-pic", lock=lock, mass_inventory=mass, producer=producer),
        ],
    }, stream, indent=2)
    stream.write("\n")
PY
sha=$(sha256sum "$plan" | cut -d ' ' -f 1)
work=$(mktemp -d -p test_output/srcMoon m3-family-XXXXXXXX)
/usr/bin/python3 -B srcMoon/test/run_m3_campaign.py \
  --plan "$plan" --plan-sha256 "$sha" --output-dir "$work/results" --wall-seconds 1800
```

`result.json` reports counts for cases, independent checks and all four
formal IDs. `comparison.csv` lists case and check outcomes. `logs/` and
`resources/` retain per-child commands and measured resource histories.
`cases/` keeps the native output, committed ledgers, individual source
manifests and child comparison figures. `formal/` retains the I21–I24
readiness results. For existing predecessor receipts from the same
qualification campaign, pass their root via `--predecessor-root`. The
formal owners read those files in place. The normal Moon aggregate and
nightly retain their existing individual IDs. The global runner executes
the new administrative controls under U16's existing orchestration group.
Run them directly with `/usr/bin/python3 -B srcMoon/test/test_m3_campaign.py`.
