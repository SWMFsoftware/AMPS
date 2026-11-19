# U15 extension: frozen D04 photon-channel Monte Carlo selection

**Contract frozen before the first execution, 2026-11-17.** This is a supporting
extension of the existing U15 table/rate/sampler meaning, not a new test ID.
It activates no particle chemistry and does not qualify molecular daughter
kinetics, transported reaction products or linked I22.

## Purpose and production capability

Exercise the actual PIC-independent production
`Moon::Species::Photochemistry::Registry::SelectChannel` after the same strict
`LoadFrozenD04` used by the application. A thin native adapter supplies draws
and counts returned channel IDs; it contains no duplicate chemistry formula.
The existing eleven rate/registry checks retain deterministic endpoint,
branch-interval, transcription, scaling, refusal and source-identity coverage.
This extension adds the predeclared statistical control required by U15.

## Configuration, inputs and source identity

The complete ordered parent list is `He, Na, Ne, Ar, H2O, OH_exp`. For each
parent, quiet=0 then active=1 creates exactly12 cases. All15 reactive channels
are tested at both solar endpoints:30 predeclared positive-branch claims.
The two total-photoabsorption metadata rows are not reactive outcomes.
`OH_exp` explicitly selects the frozen experimental compilation; it does not
invent a resolved OH initial-state population.

All rates are the D04 1AU, optical-depth-zero endpoints in s^-1. There is no
spatial frame, observation epoch, sunlight geometry or fitted model parameter
in this sampler test. The oracle derives probabilities independently from the
native processed CSV rates as each rate divided by the complete parent/solar
family sum. It never infers probabilities from native output counts.

Required data are the four v1 production CSVs with their exact existing
Photochemistry registry hashes; the runner also checks the package README,
QA/provenance/raw-inventory identities, every used PHIDRATES raw file, every
PHIDRATES processed output and the committed processing-script hashes. The
data inventory D04 entry must explicitly say `Complete/data-ready`, the native
status currently supplied by `INVENTORY.md`, and package QA must be PASS. This
does not manufacture a separate machine-readable READY approval. Missing evidence
is SKIPPED; contradictory hashes/schema or ERROR QA is ERROR. Source-empty
metadata are retained and cannot be replaced with assumptions. CEX PDFs and
electron-impact data are not dependencies of this photon-only check.

## Stochastic controls frozen before output

- Exactly `N=100000` fresh pseudorandom draws for each of the12 cases.
- Master seed `731993`, `std::mt19937_64`.
- Case index `j=0..11` follows the ordered list above. Every case constructs
  its own generator with `std::seed_seq{731993u, uint32_t(j)}`. Streams are
  declared separately; no residual-based seed search, retry or seed change
  is permitted. PRNG streams are a reproducibility convention, not proof of
  mathematical independence between arbitrary finite pseudorandom sequences.
- One 64-bit engine output per draw. With `k=engine() >> 12`, use
  `u=(double(k)+0.5)*2^-52`. There are2^52 equiprobable midpoint values strictly
  inside(0,1); neither endpoint can round to1. Binary64 quantization is bounded
  by the grid spacing and is distinct from chemistry-branch errors.

## Independent oracle, metrics and acceptance

The familywise false-rejection budget is fixed at0.01 across all30 positive
channel claims. Each claim receives two-sided `alpha=0.01/30`; this Bonferroni
allocation is valid without assuming independence of the correlated
multinomial channel counts. Each marginal count is compared with the exact
binomial-distribution acceptance interval for its independently computed
source probability `p`:

```text
L = BinomialPPF(alpha/2, N, p)
U = BinomialPPF(1-alpha/2, N, p)
accept iff L <= count <= U
```

The interval is inclusive. Equivalently the rejection tails are
`P(X<L)=CDF(L-1)` and `P(X>U)=SF(U)`, each <=alpha/2, with discrete quantiles
evaluated by SciPy's binomial distribution rather than a Gaussian/chi-square
approximation. Degenerate p=1 branches must return N exactly. The runner
records both quantiles and their numerical CDF/SF tail checks **before**
starting the native probe, in `frozen_bounds.json`. No threshold is widened
after seeing samples. SciPy's documented CDF/PPF/SF definitions control the
implementation: [official binomial reference](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binom.html).

The degenerate p=1 and p=0 laws are evaluated analytically as X=N and X=0,
respectively. Quantile minimality is checked in addition to tail bounds so an
incorrect special-function result cannot quietly broaden an interval. A
five-fair-draw exact identity and a deliberately faulty quantile implementation
exercise this refusal before sampling. The preserved first implementation
attempt exposed the installed SciPy1.10.1 p=1 low-quantile issue; correction
enforces this contract's original exact-N condition, with unchanged seed,
sample size and alpha, and does not widen a statistical criterion.

All12 case identities, seed vectors, N, open-uniform extrema and channel sets
must match the frozen declaration; family counts must sum to N. Missing,
unknown, duplicated or malformed identities are infrastructure ERROR, never a
passed zero-count branch. A finite evaluated count outside its frozen interval
is FAIL. Every case and all30 channel criteria must pass for supporting PASS.

## Procedure, artifacts and status

1. Read and verify native metadata, source bytes and independent CSV oracle.
2. Write protocol/source hashes and all frozen statistical bounds.
3. Compile the thin adapter, Moon namespace compatibility translation unit,
   and maintained shared `DataIO.cpp`/`Photochemistry.cpp`, without Config.pl,
   generated AMPS trees, MPI or CEX sources. Hash the actual shared source and
   headers through the same dependency recipe used by compilation.
4. Execute the actual native selector once for all12 cases; preserve exact
   argv, stdout/stderr, compiler/executable/source identities and timings.
5. Score identities and the prewritten bounds; write comparison CSV, JSON,
   PNG and EPS. Plots display source frequencies and declared statistical
   intervals, not observations or physical daughter distributions.

The output directory must be new. Missing SciPy, NumPy, Matplotlib, compiler
or authoritative source evidence is SKIPPED (77); source corruption, malformed
metadata, compilation/runtime errors or missing output is ERROR (2);
evaluated statistical disagreement is FAIL (1); evaluated agreement is PASS
(0). No missing dependency is installed automatically. No original data or
acceptance file is edited. Expected PASS:12/12 families and30/30 declared
channel bounds, with exact source-family/channel identities retained.

Predecessor scope: the local table/selector implementation is tested with its
analytical/control companions. Parent U15 orchestration remains responsible
for U14/U25 and literal registry status. This check does not substitute for
production PIC dispatch, time-step/conservation tests, molecular state/angle
data, chemistry product transport, I22 or earlier formal milestone gates.
