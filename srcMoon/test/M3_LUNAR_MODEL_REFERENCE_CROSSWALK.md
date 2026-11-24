# Reference solutions for lunar chemistry

## Existing AMPS output

I04 already compares a four-rank AMPS run with reviewed AMPS output.
Its frozen mean particle speed is 2211.773 m/s. Repeating the same case
checks the sampler, particle evolution and integration for changes. I21
records both native output and an independently computed event time. For
example, the eight-step, one-rank trajectory of parent 2 reacts at
18332.2712380 s in AMPS and 18332.2712397 s in the independent calculation.
These two references answer different questions about the same software.
[I04 reference](integration-tests/I04/reference/test_Moon.ref_np4)
[I21 event comparison](../../test_output/srcMoon/m3-followup-MfThQfny/trajectory-compact-aggregate/I21/trajectory/independent-audit/numerical_comparison.csv)

Comparing the archived one-rank and eight-rank AMPS outputs directly
gives identical parent fates and event times in **384 of 384 paired records**.
The pairs cover the same 128 parents at 8, 16 and 32 steps. They are matched
histories for a repeatability test. In a separate replay, the four-rank
sampler output again produced a zero-byte I04 diff against its approved
reference. Both comparisons reused existing simulation output.
[MPI comparison](../../test_output/chemical_database/lunar-source-comparison-pTDdeg8Y/mpi-repeatability.json)
[I04 comparator receipt](../../test_output/chemical_database/lunar-source-comparison-pTDdeg8Y/i04-recheck/test_Moon.reference.json)

An established AMPS application can also provide a **cross-model numerical
reference**. For water breakup, match the incident spectrum, channel, input
energy convention, compiled products and reference frame with the Europa
case. Compare reaction frequencies and the single-product speed distributions
in the parent's rest frame. Preserve the executable, input, seed and sampling
uncertainties with candidate results. Existing approval for the I04 reference
applies to the specified four-rank output. This review identifies the quantities
available for a new cross-model comparison. It freezes no new output.

## Which physical inputs other models use

| Process and model | What the model specifies | Use in AMPS |
| --- | --- | --- |
| Noble-gas photoionization in Peschel et al.'s lunar exosphere | In sunlight the neutral is removed. The paper's quiet-Sun rates for He, Ne and ⁴⁰Ar are 5.25 × 10⁻⁸, 5.90 × 10⁻⁷ and 3.06 × 10⁻⁷ s⁻¹. It studies the resulting neutral distributions over lunar terrain. | Useful reference for species survival and shadow gating. Its Ne rate is 3.21 times the D04 quiet value of 1.84 × 10⁻⁷ s⁻¹ and close to the active value of 5.87 × 10⁻⁷ s⁻¹. The source of that difference is unresolved. Preserve the published label. [Paper, Table 1 and Eq. 4](https://doi.org/10.3847/PSJ/ae4901) |
| Noble-gas daughters in EXOSPHID | The lunar code represents He, Ne and Ar each with one ionization channel, one threshold and unit branching to the singly charged ion plus an electron. The ion product has no electronic-state label. | Compare integrated destruction rates and this explicitly reduced channel model. For Ar the code uses 17.56 eV. NIST gives 15.7596119 eV from neutral ground state to ionic ground state. Substitution would change the energy available to the products by about 1.80 eV per event. [Pinned Ar code](https://github.com/charligolea/EXOSPHID/blob/62e670f9d0344ced1c1414a86b40fbba43b373ef/src/database/species/Ar.jl) [NIST Ar I](https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=argon&units=1) |
| Molecular ions in the lunar hybrid model of Poppe et al. | The model injects intact H₂O⁺ from sunlit neutral water to study solar-wind mass loading. Its tabulated conversion uses 7.0 × 10⁻⁷ s⁻¹ per H₂O. | Test ion production for this specific assumption. The output cannot determine the OH⁺, O⁺ and H⁺ breakup fractions or daughter recoil. [Authors' paper, §2.2 and Table 1](https://research.ssl.berkeley.edu/~poppe/papers/Poppe_et_al_2024b.pdf#page=3) |
| Charged water fragments in AMPS Europa and the published Callisto H corona | The Europa table uses 23.2, 39.8 and 30.5 eV of active-Sun mean excess energy for the OH⁺ + H, O⁺ + H₂ and H⁺ + OH channels, respectively. Each also produces an electron in the source definition. Callisto Table B1 uses the same three values. | These are energies for the entire channel. The Europa production code draws each transported product separately, and its shipped species list excludes electrons. Its output can test marginal heavy-product yields and speeds. The shared Huebner inputs and this sampling scheme leave the joint electron–ion–neutral distribution unresolved. [Europa channel table](../../src/models/photolytic_reactions/PhotolyticReactions_Huebner1992ASS.cpp) [Europa input](../../input/europa.input) [Shared transport](../../src/models/exosphere/Exosphere_Chemistry.cpp) [Callisto, Table B1](https://doi.org/10.1029/2022JE007294) |
| Water on lunar grains in Schörghofer's Planetary Code Collection | The thermal release flux depends on surface temperature and the number of adsorbed water molecules per square metre. The adopted monolayer coverage is 10¹⁹ m⁻² and the activated-grain binding-energy scale and width are 0.65 and 0.22 eV. | Numerically compare a finite water inventory below and above one monolayer at fixed temperature. This addresses bound molecular H₂O. Gas-phase OH hitting a grain needs its own sticking and accommodation law. [Versioned source](https://github.com/nschorgh/Planetary-Code-Collection/blob/ab8baaecb0c8f4414f651b7d31cad610339be0cb/Lunar/desorptionrates.f90#L35-L70) |
| Ion return in lunar Monte Carlo models | Smolka uses environment-dependent modeled recycling fractions, including about 30–39% for He and 43–54% for ⁴⁰Ar, from earlier calculated ion trajectories. | Compare the fraction of ions reaching the surface with a matched birth distribution and fields. Neutralization, implantation and release after impact act on the subset that returns. Applying a return fraction after explicitly tracing the same trajectories would count the transport loss twice. [Smolka, Table 3.1](https://mediatum.ub.tum.de/doc/1658109/document.pdf) [Original ion transport study](https://doi.org/10.1002/jgre.20148) |

## Electron impact and sodium charge exchange

Ambili and Choudhary's lunar ionosphere model gives a particularly clear
example of a different physical input. It uses an ion speed constrained by
Apollo SIDE and a characteristic distance to estimate the **bulk removal
rate** after interaction with the solar wind. Appendix A gives separate
source frequencies, including He, Ne and Ar photoionization at
4.51 × 10⁻⁸, 3.00 × 10⁻⁷ and 3.05 × 10⁻⁷ s⁻¹. Its ion continuity
equations do not propagate individual reaction partners and contain no Na
species. The SIDE speed can constrain population loss in that calculation.
It gives no birth velocity for Na⁺ in the AMPS charge-exchange mover.
[Model formulation and Appendix A](https://doi.org/10.1093/mnras/stab3734)

The archived AMPS Moon code under `/home/vtenishe/TEST_MOON/AMPS` contains
a legacy **electron-impact** estimate using 2 × 10⁻²⁰ m² for each of Na, Ne
and Ar and a speed of 400 km/s. The generated `build/pic/picGlobal.dfn`
disables the electron-impact ionization option. A historical run with that
configuration cannot supply an electron-impact rate or a proton–Na recoil
reference. For electron-induced ionization of He, Ren et al. measured
absolute triple-differential cross sections at 70.6 eV. Their benchmark gives
access to correlated kinematics at that energy. A 10 eV Maxwellian in I23
samples a range of incident collision energies. [Historical function](/home/vtenishe/TEST_MOON/AMPS/srcMoon/Moon.cpp:431)
[Archived configuration](/home/vtenishe/TEST_MOON/AMPS/build/pic/picGlobal.dfn:145)
[Ren et al.](https://doi.org/10.1103/PhysRevA.83.052711)

The located lunar pickup models prescribe total Na ion production from the
solar wind. Kaguya measurements include particles from both the surface and
the exosphere, with contributions from neighboring masses in the selected
Na⁺ band. These fluxes cannot isolate the differential cross section for
H⁺ + Na(3s) → H(1s) + ground-state Na⁺ at 4 keV/u. For the approved
isotropic sensitivity case an independent two-body calculation supplies the
reference recoil for each chosen angle. A reviewed AMPS histogram could then
detect numerical drift in the sampler. The physical angular distribution
for that channel is still open. [Kaguya analysis, §3](https://doi.org/10.1002/2013JE004529)
[I24 evidence](integration-tests/I24/SOURCE_ADMISSION_20261124.md)

The r14 shared catalogue now defines two selectable Na reference cases.
Poppe's case uses the three deposited median frequencies and starts pickup
ions at rest in the lunar frame. The I24 development case uses the single
Shingal A34 ground-product partial with modeled isotropic scattering in the
pair center of mass. Only the latter follows both individual collision
products. Comparing angular sensitivity requires the same production rate,
source positions and fields across runs. [Model IDs, source selectors and examples](../../utility/chemical_database/DATABASE_README.md#sodium-pickup-and-binary-charge-exchange-r14)
