# RH3D08

critical-Mach interpolation reproduces a versioned reference table over beta and obliquity, including the cold quasi-perpendicular limit; each declared `M_f`, `M_A`, or `M_An` convention compares the matching `M_chi` with `M_c^chi`, rejects cross-convention substitution, and handles the exact-`B_n=0` normal-Alfvén limit without NaN/Inf or false criticality.

Run from any directory with `python3 test.py`. The launcher delegates to the global registry, so individual and cumulative gates execute the identical implementation.
