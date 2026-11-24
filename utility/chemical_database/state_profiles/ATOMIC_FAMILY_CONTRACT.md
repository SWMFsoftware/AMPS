# Ground-ion reference profiles for helium, neon and argon

I22 and I23 need isotope-specific state energies and masses before a product
provider can use their rates. The opt-in preparation takes the retained
primary-source ledger and a pinned existing state profile. It authenticates
the three PDFs and their exact retained source records, preserves the printed
wavenumbers and reported uncertainties, and produces additive `4He+:ground`,
`20Ne+:ground`, and `40Ar+:ground` reference states. The prior unresolved
entries remain in the profile with their original status.

The retained sources are CODATA 2018 for ⁴He, Eikema et al. (1994) for ²⁰Ne,
and Velchev et al. (1999) for ⁴⁰Ar. Exact SI conversion gives formation
energy `100 h c wavenumber` in joules. The existing P02 ground-neutral/free
electron zero gives an additive reference ion mass `m_neutral − m_e` and
physical ground ion mass `m_neutral − m_e + IE/c²` in kg. The native
`PhysicalStateProfile::Registry` probe independently reads the resulting
profile and checks formation, reference and physical masses against the
preparation calculation at the existing two-long-double-epsilon criterion.
No stochastic inputs, geometry, or epochs enter this source preparation.

The authenticated D04 quiet/active endpoints and printed thresholds are
retained beside the physical thresholds and their signed differences.
`source_effective_photon_candidates.json` prepares both endpoints for each
of the three isotopes. It records the printed threshold plus matched mean,
the remaining energy above the physical threshold before recoil, and each
native source record. Each row records a SKIPPED transport decision pending
source/state compatibility, ground-ion yield or a bounded reduced closure
and observable-error review. The worksheet is separate from the executable
state profile. A published database snapshot can consume the supported state
records through its versioned publication path.

An output records the exact source/recipe/compiler/executable hashes and
commands, native probe log, comparisons and per-check result. Missing source
or compiler is SKIPPED, changed bytes, malformed input or failed setup is
ERROR, and a native/oracle mismatch is FAIL. A PASS here verifies the
reference-state construction and loader only. I22/I23 retain their original
U-series, numerical, conservation and same-campaign predecessors.
