# Lunar neutral-surface interaction literature receipt

This note records the literature review used for the input-selectable M2
sticking and immediate non-sticking accommodation laws. It is a physics-input
receipt, not observational validation. U12 verifies the production kernels;
I18 must still exercise their linked stochastic distributions and budgets.

## Sodium

The selected source is V. Tenishev, M. Rubin, O. J. Tucker, M. R. Combi, and
M. Sarantos, “Kinetic modeling of sodium in the lunar exosphere,” *Icarus* 226
(2013) 1538–1549, doi:`10.1016/j.icarus.2013.08.021`. The reviewed local copy
is `/home/vtenishe/1-s2.0-S0019103513003692-main.pdf` (1,833,748 bytes,
SHA-256 `45edf66de2d6e9ed91c8bb7cae6d9ab9aac0038adb3592ba9725cb7a4d007fe8`).

- Section 3.3.1 and Eq. 14 define the outgoing speed of a non-sticking Na atom
  as `vf=sqrt(omega*vt^2+(1-omega)*vi^2)`, where `vi` is incident speed and
  `vt` is a Maxwellian thermal draw at the local surface temperature.
- The paper evaluates `omega=0.2`, cites Killen et al. (2012) and consistency
  with Mouawad et al. (2011), and distributes the outgoing direction
  isotropically.
- The cited lunar calculation by R. M. Killen et al., “The effect on the lunar
  exosphere of a coronal mass ejection passage,” *JGR Planets* 117 (2012),
  doi:`10.1029/2011JE004011`, independently documents the same Na choice above
  200 K: 20 percent thermal contribution and 80 percent rebound contribution.
  It also makes clear that this is a model selection, not a direct laboratory
  measurement of one universal regolith coefficient.
- Fig. 3 gives the temperature-dependent sticking coefficient and states that
  it is adapted from B. V. Yakshinskiy and T. E. Madey, “Temperature-dependent
  DIET of alkalis from SiO2 films: Comparison with a lunar sample,” *Surface
  Science* 593 (2005) 202–209, doi:`10.1016/j.susc.2005.06.062`.
- The numerical nodes in `srcMoon/examples/lola_surface.in` preserve the
  pre-existing AMPS digitization formerly embedded in `Moon.cpp`. The original
  digitization project/tool and graphical uncertainty are not present in the
  repository. Consequently the table is an explicit production parameter but
  is not claimed as independent validation evidence.

The example therefore selects `temperature-table` sticking,
`kinetic-energy` accommodation, and `omega=0.2`. Values outside the table's
100–494.99997 K support are clamped to the nearest endpoint, matching the
bounded interpolation contract tested by U12. Any alternative must be stated
in the input receipt with its own source string.

## Argon

The current residence-time source is J. A. Kegerreis et al., “Evidence for a
Localized Source of the Argon in the Lunar Exosphere,” *JGR Planets* (2017),
doi:`10.1002/2017JE005352`. The implemented modified-Arrhenius residence law is
`tau=h/(kT)*exp(Q/(RT))`, with the explicitly model-dependent `Q=28 kJ/mol`
selected in the example.

That paper supports residence/release modeling; it does not establish a
general measured immediate-impact sticking/accommodation curve for arbitrary
lunar regolith. The current M2 receipt therefore declares constant-unity Ar
sticking explicitly. Its accommodation model is also declared, but is
unreachable while sticking remains unity. This prevents the limiting choice
from being misreported as a laboratory-derived probability curve.

### Radiogenic source receipt

Killen, “Source and maintenance of the argon atmospheres of Mercury and the
Moon,” *Meteoritics & Planetary Science* 37 (2002) 1223–1231,
doi:`10.1111/j.1945-5100.2002.tb00891.x`, gives present lunar diffusion fluxes
of 990 and 1448 atoms cm^-2 s^-1 for its two diffusion-coefficient cases.
Integrated over a lunar sphere, Kegerreis et al. (2017) report the corresponding
range as `3.8e20`–`5.5e20 atoms s^-1`. The maintained input selects the lower
published bound exactly; it is not a midpoint, fit, or uncertainty reduction.
The NASA NTRS PDF retrieved 2026-10-13 had SHA-256
`1e4869bc418a07b2a01be6fdf2184f5354531d75988e0ea673e7ccd7cb4c6459`.

Kegerreis et al. implement a continuous source, choose random global surface
locations for their global case, and draw an accommodated local-temperature
Maxwell-Boltzmann flux distribution. Those are the implemented source
semantics. Their K-weighted, point, and top-hat source rates (`1.1e22`,
`1.9e21`, and `2.9e21 atoms s^-1`) are fit-dependent alternatives with
different inferred lifetimes; they are not copied into the uniform baseline.
The paper explicitly leaves the localized physical mechanism unresolved.

The local D10 Kaguya GRS package covers only -90 to -13 degrees latitude and
130 to 260 degrees east, has missing below-detection cells, and declares
`D10 DATA INCOMPLETE`. It therefore cannot define global source geography.
The parser accepts only `uniform-global` plus `map_file=none`; no missing area
is filled, mirrored, or assigned a mean K abundance. U21 preserves the exact
D10 metadata/proxy hashes but remains `SKIPPED` until global geography and its
linked stochastic/transient gates exist.

## Species scope

M2 currently promotes neutral Na and Ar surface inventories. He and Ne are
separate roadmap capabilities (U19/I26 and U20/I27 respectively), and ions
require neutralization/return branching data. This review does not assign Na
or Ar laws to those species. Missing species-specific evidence remains
`SKIPPED / NOT VALIDATED` rather than receiving an assumed coefficient.

The literature search also found a later, species-specific He constraint:
C. Grava et al., “LRO/LAMP observations of the lunar helium exosphere:
constraints on thermal accommodation and outgassing rate,” *MNRAS* 501 (2021)
4438–4451, doi:`10.1093/mnras/staa3884`, reports that its LAMP/model comparison
favours full He thermal accommodation (`alpha=1`). That result is recorded for
the later He work package, but it is not copied into the current Na/Ar M2
adapter: He source/wiring and its own verification/validation gates are not yet
operational. The reviewed Ne literature did not supply an equally explicit,
independently qualified sticking curve or accommodation coefficient, so no Ne
value is introduced here.

## Input selections

For each M2 neutral, the application receipt names independently:

- `*_sticking_model`: `constant` or `temperature-table`;
- `*_sticking_probability_by_temperature_k`;
- `*_accommodation_speed_model`: `linear-speed` or `kinetic-energy`;
- `*_accommodation_coefficient`, constrained to `[0,1]`;
- `*_surface_interaction_reference`.

The generic Exosphere callback retains `linear-speed` only when no application
runtime law is installed, preserving other AMPS applications. `srcMoon` M2
never relies on that implicit legacy selection.
