# Reduced-model numerical verification and the product-floor certificate

This note documents the current shared `ReducedModels.h` arithmetic and its
independent verification. It does not change a source, model parameter,
observational gate, physical tolerance or application activation. Native guards
exercise the actual shared header; Decimal references are analytical fixtures,
not observations or frozen AMPS output.

## Why a zero COM-energy result needs a certificate

Two-body final momenta depend on
`K_COM=K_in+Q-P²/(2*(m1+m2))`. A tiny positive mass can disappear from the
floating total mass; cancellation can then round K_COM to zero while the light
product's true recoil velocity is representable. A rounded total-energy budget
alone can hide that error. The current service refuses a sum that entirely
loses a positive input mass, refuses an underflowed positive reduced mass, and
requires an exactness certificate for its special zero-COM case.

That certificate checks the mass/energy TwoSum tails, each squared-momentum
product's FMA residual, a canonical TwoSum tail for the squared-product sum, and
the FMA residual of `target_energy*(2*total_mass)-sum_squares`. It does not clamp
negative COM energy or add an epsilon. Exact dyadic zero-COM controls remain
accepted. If precision cannot establish the needed identities, the result is
typed `NumericalFailure` with the output unchanged.

## Machine-derived product range

For the tested binary floating environment, let p be
`numeric_limits<long double>::digits`, emin its `min_exponent`, and
`min_normal=2^(emin-1)`. With gradual subnormal support the least positive
subnormal is `2^(emin-p)`. The code uses

```
product_floor = scalbn(min_normal,p) = 2^F
F = emin-1+p
```

For a nonzero factor, a rounded product **strictly below** this floor is refused
as an exactness witness. Equality is permitted. Here is why the equality case
does not leave a half-subnormal-residual hole.

Write two nonzero floating magnitudes as normalized p-bit significands:
`a=A*2^(ra-p+1)`, `b=B*2^(rb-p+1)`, with integer
`2^(p-1)<=A,B<=2^p-1`. This representation also covers subnormal inputs by
normalizing them with trailing zero significand bits. Put `r=ra+rb`.
If `r<=F-2`, the largest possible exact product is

```
(2-2^(1-p))² * 2^(F-2)
= 2^F - 2^(F+1-p) + 2^(F-2p).
```

Its deficit from `2^F`, divided by half an ULP immediately below `2^F`, is
`4-2^(1-p)`, at least 3.5 for binary p>=2. Thus it cannot round up to the
floor in round-to-nearest arithmetic. A rounded product at or above the floor
therefore implies `r>=F-1`. The exact product's least possible nonzero bit has
exponent

```
r-2p+2 >= F-2p+1 = emin-p.
```

Every nonzero difference between this exact product and its rounded product is
therefore an integer multiple of the least subnormal, not a smaller residual
that rounds spuriously to zero. A single-rounding `std::fma(a,b,-rounded)` can
then certify the product. No prime-factorization assumption about p is needed.
The same reasoning applies to squared momentum and the generic target-energy
product. Explicit finite checks reject overflow. Nonzero factors whose product
rounds to zero are also rejected by the range guard.

The proof assumes ordinary binary round-to-nearest evaluation, correctly
single-rounded FMA and gradual underflow, as in the recorded strict C++17
build. It does not authorize fast-math/flush-to-zero or different floating
environment configurations as qualified implementations. Compiler and source
identities are retained with the actual test artifacts.

## Fixed defects and independent checks

- Opal's complementary-angle inverse loses relative accuracy for small support
  over width. Direct tangent is now used whenever that ratio is at most one;
  the complementary form is reserved for the large-ratio upper-quantile regime.
  Decimal atan/sine/cosine series independently verify the unchanged 1e-12
  numerical criterion. Native endpoint support is distinguished from the
  original decimal text's binary conversion; no physical support is widened.
- Thermal `nu*exp(-E/kBT)` can lose the exponential factor before multiplication
  by a large supplied prefactor. The actual rate uses the log-domain expression,
  while integrated exposure independently retains rate-times-duration in log
  space. Underflow/overflow of auxiliary exposure representation is explicit;
  no physical rate cutoff is selected.
- The inherited-velocity reservoir explicitly accounts for supplied charge,
  momentum, incoming bath energy and Q. A compensated four-term sum preserves
  small Q when large parent/heavy bulk kinetic terms cancel. Negative signed
  reservoir exchange remains explicit, not discarded.
- Ion binding/reference masses remain separately accessible. Radiative weak
  branches are not silently removed; an unrepresentable positive branch is a
  numerical refusal. Missing source values are never provided by these kernels.

The maintained standalone runner is
`python3 -B test/chemical_database/test_reduced_models.py --output-dir NEW`.
Final verification on the recorded source basis evaluates seven Python test
methods, 38 native domain/transaction/extreme controls, and 210 independent
Decimal metric comparisons. Earlier actual FAIL/ERROR attempts remain retained;
these controls establish mathematical software behavior, not physical M3
qualification or observational agreement.
