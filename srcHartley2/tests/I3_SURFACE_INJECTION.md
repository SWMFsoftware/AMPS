# I3 surface-injection regression

The I3 parent gas is injected by the production AMPS routine
`PIC::Distribution::InjectMaxwellianDistribution`; the comet application does
not carry a second velocity sampler.  For a stationary surface at temperature
`T`, the crossing distribution is

```text
p(v_vector) proportional to v_n exp[-m |v|^2/(2 k T)],  v_n > 0.
```

Consequently, “`v` times a Maxwellian” is correct for the **speed marginal**:
if `q=m v^2/(2 k T)`, its CDF is `1-(1+q) exp(-q)`.  The complete boundary law
also contains the cosine angular distribution `p(mu)=2 mu`, not an isotropic
outward direction.  Each tangential component is a zero-mean Gaussian with
variance `k T/m`; the outward normal component has

```text
E[v_n]   = sqrt(pi)/(2 beta)
E[v_n^2] = 1/beta^2,
beta     = sqrt[m/(2 k T)].
```

`RadialMain --validate-surface-injection` samples the production routine for
H2O, CO, and CO2 at the agreed I3 source temperature of 175 K.  It verifies:

- every sample moves outward and the no-correction call returns unit weight;
- the analytical flux-energy and cosine-law CDFs with a KS distance;
- the first two normal-velocity moments;
- both tangential means and variances.

The random stream is fixed per MPI rank, making failures reproducible.  The
tolerances are eight asymptotic standard errors plus small minimum tolerances;
this avoids a flaky Monte Carlo gate without weakening sensitivity to an
ordinary volume-Maxwellian or an isotropic-hemisphere regression.

After compiling the shared I3 executable from `input/i3.conf`, run a selected
case from the combined runtime receipt:

```bash
srcHartley2/tests/run_i3_surface_injection_test.sh /absolute/path/to/amps
```

Optional environment controls are `I3_INJECTION_MPI_RANKS` (default 1),
`I3_INJECTION_PARTICLES` (default 200000), `I3_INJECTION_LOG`, `MPIEXEC`,
`I3_RUNTIME_INPUT` (default `input/i3.in`), and `I3_CASE` (default 1). The
runner always passes `-i ... -case ...`; it never rebuilds the executable.
