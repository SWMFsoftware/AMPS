# CV10 — nonuniform turbulence advection

CV10 verifies conservative spectral wave transport for both characteristics
`U+V_A` and `U-V_A`. A cell-integrated sinusoidal wave-action profile is
advected periodically with `U=0` and `V_A=1 m/s`; every wavenumber bin carries
a known fixed fraction. Fixed-area, smoothly expanding-area, and mid-run
coarse-to-fine remap scenarios exercise common-face fluxes and moving geometry.

The independent reference integrates the shifted sinusoid over each final
cell. Resolutions 32, 64, and 128 provide automatic L2 refinement order. Gates
cover finest-grid error, wave-action conservation, non-negativity, both
propagation senses, and all spectral bins. The moving-grid reference includes
the expected conservative remap diffusion but no source, damping, reflection,
cascade, particle coupling, or shock injection.

Run `python3 test/run_tests.py --amps /path/to/amps --validation-case CV10
--output-dir test_output/CV10`.
