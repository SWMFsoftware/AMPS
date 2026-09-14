# VAL01: manufactured Parker benchmark

`VAL01` compares the production `AdvanceParker` kernel with analytical
solutions that do not call the mover. For constant parallel diffusivity
`kappa`, constant field-aligned plasma velocity `U`, and constant velocity
divergence, the reference is

```text
E[s(t)-s(0)] = U t
Var[s(t)-s(0)] = 2 kappa t
p(t) = p(0) exp[-div(U)t/3].
```

The fixed keyed seed is `1501001`; 80,000 particles use `kappa=2.5e13 m2/s`,
`U=4e5 m/s`, `div(U)=1.5e-5 s-1`, and `t=40 s`. Acceptance requires the mean
within 4.5 analytical standard errors, variance relative error at most 0.02,
and momentum error within 128 machine epsilons. The test verifies one
manufactured constant-coefficient problem; it is not observational or coupled
evidence.
