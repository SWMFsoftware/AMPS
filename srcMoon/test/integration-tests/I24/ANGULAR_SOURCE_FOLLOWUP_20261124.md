# I24 proton–sodium scattering source follow-up

The controlled I24 case selects
`H⁺ + Na(3s) → H(1s) + Na⁺(2p⁶)` at 4 keV/u in the target-rest
frame. Admission requires primary lineage for the integral cross section,
the two product states, and a differential law or a source-supported error
bound for the same reaction and energy.

## Archive and publication findings

- The retained IAEA CollisionDB D119990 record identifies an experimental
  JAEA/Kyoto ALADDIN compilation. Its numeric arrays contain 12 knots
  with projectile energy in eV/u and integral cross section in cm².
  The record names the outgoing states and leaves the initial Na state
  and original publication unspecified. The public page could not be
  retrieved through the browser during this follow-up. The previous
  successful acquisition supplies the byte-verified record. The NIFS
  CHART public search and OPEN-ADAS sodium catalogue did not produce a
  matching angular table. QST/JEAMDL identifies the Tabata metal-vapor
  fits as energy-dependent single capture with H as the product. The
  separate Rice absolute differential archive lists other targets.
  [IAEA D119990](https://db-amdis.org/collisiondb/datasets/119990/)
  [NIFS CHART](https://dbshino.nifs.ac.jp/nifsdb/chart/top)
  [OPEN-ADAS](https://open.adas.ac.uk/element/na)
  [QST/JEAMDL](https://www.qst.go.jp/site/jt60/6500.html)
  [Rice Collision Group](https://www.ruf.rice.edu/~atmol/charge_data.html)
- Knoop's Chapter 5, §5.2.3 identifies the outer-shell `H(n=1)`
  distribution from 7 keV/u upward. The published 25 keV/u trace in
  Fig. 5.10 is for the selected ground product channel and uses a
  milliradian laboratory angle inferred from recoil momentum. The
  4 keV/u panel shows the other channels. The curves are normalized to
  their own peaks. No public numerical bin table with angular acceptance,
  absolute normalization, uncertainties, and tail coverage was located.
  A plot boundary gives no upper bound on the true scattering angle.
  [Knoop, Chapter 5, Figs. 5.9–5.10](https://pure.rug.nl/ws/portalfiles/portal/2855778/c5.pdf)
  [Knoop, Chapter 4, Eq. 4.10](https://pure.rug.nl/ws/portalfiles/portal/2855777/c4.pdf)
- Courbin et al.'s published differential plots and the Daresbury EIKONXS
  examples use excited hydrogen products. They cannot supply ground
  capture at 4 keV/u. Knoop et al.'s 2004 inner-shell H(1s) measurements
  leave excited Na⁺. A common H(1s) label alone cannot identify the
  two-body state and energy balance.
  [Courbin et al.](https://oa.upm.es/61369/)
  [EIKONXS report](https://epubs.stfc.ac.uk/manifestation/5491/DL-SCI-TM69T.pdf)
  [Knoop et al.](https://doi.org/10.1103/PhysRevA.70.050702)

## Integral cross-section discrepancy

The preserved D119990 array gives `6.34e-15 cm²` at 4 keV/u.
Knoop's Fig. 5.6 places prompt outer-shell `H(n=1)` between
`1e-17` and `1e-16 cm²` at that energy. The D119990 value is of the
scale of the much stronger excited capture in the same figure. This
visual magnitude check challenges the proposed source crosswalk.
It does not establish the original D119990 process or justify changing
its label. No numerical data were digitized from the figures. Require
the original D119990 citation and state and measurement definitions, or
a reviewed replacement whose integral and angular measurements refer
to the same channel. [Knoop, Equation 5.1 and Fig. 5.6](https://pure.rug.nl/ws/portalfiles/portal/2855778/c5.pdf)

## Conditional isotropic comparison

The uniform fallback requested by the user is specified here as a
uniform **solid-angle** law in the incoming pair center-of-mass frame,
with independent `cos(theta) ~ Uniform(-1,1)` and
`phi ~ Uniform(0,2*pi)`. The outgoing H and Na⁺ carry equal and
opposite momenta there. Their momentum magnitude follows from the
resolved masses and reaction energy. The products then transform
back to the event frame. This law is suitable for an experimental
angular sensitivity case with a declared rate. It cannot resolve an
uncertain rate or qualify the physical recoil distribution.

For scale, assume carriers of 1 u and 23 u, a stationary target,
4 keV projectile kinetic energy, and zero reaction energy **only for
this calculation**. The target recoil is
`2*m1*m2/(m1+m2)^2 * E * (1 - cos(theta))`. It has mean
`319 eV` and maximum `639 eV` under isotropy. The forward limit
of this illustrative calculation gives zero recoil. This difference
is large enough that isotropy must be recorded as a sensitivity
model when interpreting pickup-ion energies.

## Provenance and status

AMPS revision `e5f2978ac0576f42dc510e0da110dcbe0029d69a`.
The source-roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
The validation-data roadmap SHA-256 is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The AMPS version authorizes early supporting work. The common
qualification contract governs I24.

The IAEA HTML SHA-256 is
`1015ec7a843e02ef7fc1ddf97f3589ac4dc4823f12bb28e9434313876dd90fd6`.
The 12-knot CSV SHA-256 is
`2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93`.
The local Chapter 5 PDF SHA-256 is
`7daca10c4ecf37bd79f0a5d353cfc9c6d6adc9c2b0277d74aa01d5381c5c1785`.
The Courbin PDF SHA-256 is
`030e49b566896e6839054a7fed5489e642032128f17c2a5b14670709666badf4`.
The rendered review pages and independent scalar calculation are under
`test_output/srcMoon/i24-angular-source-nCgP6h/`. The local plot
inspection used PDF pages 9 and 13 at their published labels
Fig. 5.6 and Fig. 5.10.

New angular data rows admitted **0**. The source review produced no
application build or linked test. The uniform case is a documented
experimental option. Formal I24 remains **SKIPPED**.

The shared chemical database publishes this source review and the
generic isotropic direction definition in r9. The original integral
table links to the source review. See the
[shared catalogue contract](../../../../utility/chemical_database/NATIVE_CEX_SCATTERING_SOURCES.md)
and the installed README section 15 for query examples.
