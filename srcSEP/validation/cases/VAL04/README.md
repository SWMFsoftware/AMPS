# VAL04: coupled-background replay

The implemented `VAL04-SWCME` case constructs the real header-only SWCME 1-D
model in `SHOCK_ONLY`/`SOURCE` mode, prepares 24 states at 300 s cadence, queries
SI background fields through `swcme::sep::Interface1D`, and publishes an
immutable model-owned srcSEP background snapshot at every epoch. The production
`AdvanceParker` kernel consumes SWCME velocity and divergence under the active
read phase. An independently accumulated characteristic checks arc length and
adiabatic momentum, while a production SWCME shock query must return an active
physical SEP source. Configuration fingerprint, seed, metrics, and tolerances
are retained in JSON and JUnit.

Zero diffusivity makes the replay deterministic and isolates coupling units,
signs, epochs, ownership, and source representation. It is a real model-to-core
integration, but does not exercise generated PIC particle storage or an MPI
driver.

`VAL04-SWMF` is a separate mandatory release gate. It requires a native coupled
run with real SWMF output, read-only import provenance, configuration, coupling
cadence, compiler, output schema, metrics, and checksum-verified inputs. The
source archive contains only a manifest template; absence of that external
record leaves coupled integration `INCOMPLETE` even when `VAL04-SWCME` passes.
