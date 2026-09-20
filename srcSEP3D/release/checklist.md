# Candidate release checklist

1. Freeze the input, source revision, compiler/MPI identity, executable SHA-256,
   campaign seed, mesh preflight, and provider fingerprints.
2. Run the complete standalone and production-build suite without weakening a
   tolerance or converting a missing prerequisite to PASS.
3. Run the selected V03 native profile and retain every native JSON artifact.
4. Run the V04 ladder in order: analytic, native, cross-model, observational,
   then live SWMF. A higher rung never substitutes for a missing lower rung.
5. Generate the V05 matrix. Release only if the requested profile is PASS.
6. Record all skips and blocks. Until R8 is implemented, the scientific-
   release profile must remain INCOMPLETE.
