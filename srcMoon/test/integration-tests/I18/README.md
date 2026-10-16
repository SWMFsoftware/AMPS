# I18 — Species-specific sticking, residence time, and re-emission distributions

Purpose: exercise Na and Ar through the linked production impact/source
dispatch and compare the stochastic outcomes with independent analytical
probability laws. This is a verification test; it does not use lunar
observations as a substitute for a laboratory surface-interaction oracle.

I18 requires the dedicated Na+Ar qualification executable. The canonical
M0/M1 executable remains Na-only. During a normal aggregate run,
`run_all_test.py` derives a temporary Config.pl application that changes only
`SpeciesList=Na` to `Na,Ar` and enables the registered `RadiogenicArgon`
source block. It saves the executable and a hashed build receipt under the
campaign directory, removes the temporary input files, and then restores the
canonical Na-only build before running the general campaign. At runtime the
I18 runner derives `application.in` from `lola_surface.in` by changing only
`argon_source_enabled=false` to `true`.

Run this ID through the aggregate runner so its I13 predecessor is evaluated
in the same immutable campaign directory:

```bash
python3 srcMoon/test/run_all_test.py --campaign-mode qualification --only I11 --only I12 --only I13 --only I18 --mpi-ranks 8
```

For a deliberate build-reuse run, both executables must be named explicitly:

```bash
python3 srcMoon/test/run_all_test.py --skip-build \
  --amps ./amps --m2-amps /path/to/amps-moon-m2-na-ar \
  --campaign-mode qualification \
  --only I11 --only I12 --only I13 --only I18 --mpi-ranks 8
```

The runner always writes `stdout.log`, `stderr.log`, `application.in`,
`run_manifest.json`, and `result.json`. If I13 is not `PASS`, this test is
`SKIPPED` with `execution=not-run`; a U-series pass is not promoted to linked
evidence. The physical inputs are frozen in the reviewed `.in` receipt:
Tenishev et al. (2013), Eq. 14/Fig. 3 for Na and Kegerreis et al. (2017) for
the Ar modified-Arrhenius/Maxwell-flux model. The uniform continuous
`3.8e20 atoms s^-1` Ar source is the declared lower-bound Killen (2002) M2
test driver; it is not an I28 geography qualification.

Every MPI rank selects an actual LOLA facet and prefers a locally owned one;
I14 separately owns mover/facet MPI wiring. For Na, 4096 incident particles
per rank call
`ParticleSurfaceInteraction_SurfaceAccomodation`, the same callback used by
the production mover, with the facet normal and epoch-correct 6x6 SPICE state
transforms. For Ar, every rank makes 32768 calls to the registered production
`RadiogenicArgon` generator; only attempts whose randomly selected source
element is locally owned become samples. This preserves normal MPI ownership
weighting and permits valid decompositions in which a rank owns only far-field
blocks. The retained per-rank CSV tables contain SI
temperature and speeds, the outward direction cosine, and the production
status code. The Python scorer evaluates independent CDFs and never supplies
sample values to the model.

The adapter also installs a temporary known population in one ordinary
logical-surface element, calls the production thermal-desorption rate callback
for Na and Ar, and restores the original inventory. The scorer independently
evaluates the input-selected first-order Na Arrhenius law and modified-
Arrhenius Ar law. Each linked rate must agree within `1e-12` relative error.

I18 retains the approved 100 km LOLA surface mesh and all thermal/source/
surface parameters. Because it does not integrate a transport solution, its
derived receipt marks `mesh_profile=custom` and uses the bounded regression
volume targets only to distribute facet ownership over MPI ranks. A FULL
5--10 km transport volume would not change the sampled laws and is not used as
a proxy convergence test. The generated receipt states this distinction and
cannot satisfy I11-I13 or any M2 spatial-convergence gate.

Frozen acceptance gates are:

- the Na re-emission count lies inside the two-sided normal 95% binomial
  interval computed from each sample's input-table sticking probability;
- each one-sample Kolmogorov-Smirnov distance is at most `1.36/sqrt(N)` for
  Na accommodated energy, Na outward cosine, Ar Maxwell-flux energy, and Ar
  cosine-law direction;
- the production Ar source receipt is exactly `3.8e20 atoms s^-1` to a
  relative tolerance of `1e-15`;
- both production residence/desorption rates agree with the independently
  evaluated input laws within `1e-12` relative error;
- Na sample counts are exactly 4096 per MPI rank; the aggregate Ar row count
  must equal the linked accepted count and be at least 4096; all values are
  finite and physically bounded.

The Na energy variable is
`[v_f^2-(1-omega)v_i^2]/[omega 2 k T/m]`, whose reference is the unit
Maxwell energy CDF. The Ar variable is `m v^2/(2 k T)`, whose outward
Maxwell-flux CDF is `1-(1+q) exp(-q)`. Direction CDFs are `mu` for the Na
uniform hemisphere and `mu^2` for the Ar cosine law.

An executed comparison must additionally preserve the exact row-level
`comparison.csv` and render those unchanged rows as `comparison.png` and
`comparison.eps`. A skipped/non-executed test creates none of these three
files; empty or synthetic comparison graphics are prohibited.

`PASS` means all linked stochastic gates ran and passed. A missing Ar build,
source registration, terrain state, or I13 predecessor is `SKIPPED` when
identified before launch and `ERROR` when it contradicts the supplied build
or receipt. A statistically evaluated miss is `FAIL` and its samples remain
in the campaign directory.
