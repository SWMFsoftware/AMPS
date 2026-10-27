# Device math checks

Small, self-contained checks behind the `efmod` helper proposed for
`gpu/EngineConfig.h` (see `doc/BUILD.md`, section "Known build issue: fmod").

| file | what | how |
|---|---|---|
| `fmod_probe.cpp` | nvc++ 24.5 has no device `fmod`/`fmodf`; `remainder` is available | compile/link only (commands in the file header). **Never run the binary** (it initializes a GPU). |
| `efmod_exactness.cpp` | CPU proof that the `remainder()`-based `efmod` is bitwise identical to `std::fmod` | `g++ -O2 -std=c++17 -ffp-contract=off efmod_exactness.cpp -o efmod_exactness && ./efmod_exactness` |

Result of `efmod_exactness` (2026-10-21, g++ 9.4, 48.6 s):

```
efmod (remainder-based) vs std::fmod: 20000000 samples per type + edge cases: double mismatches=0 float mismatches=0
```

Samples per type: 2e7 engine-like arguments (|drift angle| in [0,1e6] rad, y=2*pi),
1.2e6 exact multiples/half-multiples of 2*pi and their floating-point neighbours,
2e7 random bit patterns (all exponents and signs, including NaN/inf/subnormals),
and all pairs of {±0, ±1, ±2*pi, ±inf, NaN, denorm_min, FLT/DBL_MIN, FLT/DBL_MAX}.
