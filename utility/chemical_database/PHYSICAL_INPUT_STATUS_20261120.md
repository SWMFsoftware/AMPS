# Targeted remaining physical inputs — 2026-11-20

This follows the user's request to search online databases first, then primary
publications for remaining M3 inputs. The source revision is
`ef73fd281f279e7538a7ead67620c02a629f0270` plus recorded working-source hashes.
The database/source-admission contract is
`PHYSICAL_INPUT_STATUS_TARGETED_CONTRACT_20261120.md`. Original records,
HTTP/UTC/hash/rights evidence and detailed reports are retained in
`test_output/chemical_database/m3-targeted-inputs-20261119T235459Z/`.

The key correction to the earlier gap list is that many required energies and
named reduced-model laws **are available**. They need consistent derivation,
explicit state/model selection and production integration. A full measured
joint electron/fragment distribution is not required by the already approved
bulk-reservoir, source-mean/isotropic or neutral microcanonical reference scope.
Their limited accuracy claims and actual admission gates remain binding.

## 1. New reusable source access

The separate snapshot is
`/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261120-r5`.
Select its own lock explicitly; original root/r2/r3/r4 locks and source bytes
remain unchanged. Its README section11 gives actual source/table API examples.
This release activates no chemistry or application profile.

* **Burgess2016 H/O thermochemistry**, [DOI10.6028/jres.121.005](https://doi.org/10.6028/jres.121.005):
  official PDF, native [PMC7339710 HTML](https://pmc.ncbi.nlm.nih.gov/articles/PMC7339710/),
  pinned NIST notice and two exact cell views. Table3 has208cells of OH
  formation-enthalpy material at **298.15K**, including historical methods;
  it is not an automatically selected0K state inventory. Table13 has375cells
  of ion-energy quantities with original units, references and uncertainties.
  Superscripts, footnotes and row/column spans survive. No value/uncertainty
  or isotope mapping is selected by the reader. The source-specific NIST
  technical-publication grant is separate from WebBook/ASD SRD permissions.
* **Laffon sticking source**, [Figshare13274840.v1](https://doi.org/10.6084/m9.figshare.13274840.v1),
  [primary article](https://doi.org/10.1038/s41550-020-01288-7): original
  Table1.txt and exact CC BY4.0 API record. Payload2136bytes, producer MD5
  `95fb477a081ae84bc1d7019c18dfe941`, SHA256
  `904badaf7780b14d85f3e0f8f03b3b2c58f1a6a56c4ea0950e8bcfe7d7600da2`.
  This is XPS-derived sticking **relative to gold assigned S=1**, with a
 20±1K surface and300K background gas. Water is the bare-substrate section;
  ice-coated sections concern CO2. It does not measure energy accommodation.
  The file's undeclared eight-bit encoding and irregular layout remain native
  bytes; no guessed UTF8/CSV view is published.

New publication/model-source records with unresolved managed-payload rights
are retained as metadata-only. Their original research facts can be cited;
managed full-paper redistribution and scientific usefulness are distinct.
Detailed packet/version/hash selection is in the publication plan/provenance.

## 2. Molecular energy and state evidence

Databases inspected include the existing ExoMol/CHIANTI products, ATcT and NIST.
ATcT direct requests returned403; these are acquisition limitations, not proof
that the evaluated energies do not exist. The literature supplies these
concrete anchors:

| Primitive | Primary evidence | Remaining mapping work |
|---|---|---|
| H2-16O lowest H+OH dissociation | Boyarkin2013, [10.1016/j.cplett.2013.03.007](https://doi.org/10.1016/j.cplett.2013.03.007),41145.92±0.12cm^-1 at95%confidence | Retain isotope and OH lambda-doublet/state convention. |
| OH binding | Sun2020, [10.1063/1674-0068/cjcp2003026](https://doi.org/10.1063/1674-0068/cjcp2003026),35580±15cm^-1 | Match the experimental state/asymptote; do not assign an unreported confidence law. |
| OH ionization | Wiedmann1992, [10.1063/1.463179](https://doi.org/10.1063/1.463179),104989±2cm^-1 | Match neutral/ion rotational references; NIST table13 preserves this source. |
| Water adiabatic ionization | Lauzin2015, [10.1080/00268976.2015.1049971](https://doi.org/10.1080/00268976.2015.1049971),101766.3(10)cm^-1 | Preserve nuclear-spin/rotational origin and uncertainty convention. |
| Water OH+ appearance threshold | Bodi2014, [10.1039/C4SC00791C](https://doi.org/10.1039/C4SC00791C), weighted18.1182±0.0009eV | Independent ion-cycle check; do not confuse the individual18.1183 datum or water-dimer channel with this result. |
| O electronic state in H2O three-neutral branch | Slanger/Black1982, [10.1063/1.444111](https://doi.org/10.1063/1.444111),2H(2S)+O(3P) at1216Å | This resolves the cited wavelength-specific lineage; it does not establish all-wavelength fine-structure populations. |

The retained primary-source report gives the complete eleven-channel symbolic
threshold crosswalk. For example, with source-defined dissociation energies
Dw(H2O),Do(OH),D2(H2), the H2+O(1D) threshold is
`Dw+Do-D2+E[O(1D)]`; OH+H+e is `Dw+IE(OH)`. This uses one separated-atom/
electron reference, not differences of rounded heavy rest energies.
Existing H2/O radiative/state sources supply internal offsets; physical ion
masses are derivable from compatible neutral masses and these energies.
Separate measured-ion-mass data are not a prerequisite.

**H2+ is not a product of the frozen eleven molecular channels.** Water05 is
O+ + H2 + e. It must not remain a required-input blocker for this network.
Existing isotope-specific4He/20Ne/40Ar threshold evidence is also available
in the preceding primary-source ledger; those profiles need implementation,
not another unsupported natural-average mass substitution. Fresh ExoAtom
He/Ne/Ar definitions have null ionization fields and do not replace those
primary threshold anchors.

Actual remaining choices are explicit product/reference levels and water spin
treatment, a named reference-threshold/source-mean convention where D04 prints
no molecular threshold, and a **separate dissociative-ionization partition**
for H2O04–06. The approved neutral-three-body P04 does not authorize applying
that law silently to an electron/charged three-body event. All original D04
rates, endpoint means and branches remain unchanged.

## 3. Electron impact and charge exchange

The IAEA/QST search was followed into Igenbergs etal2008,
[10.1016/j.adt.2008.06.001](https://doi.org/10.1016/j.adt.2008.06.001), whose
retained native preprint is dated2007. It is not the earlier bibliographic
shorthand Schweinzer1999. Equations2–5 and9 recover source-law definitions
for existing Na valence/subshell EI and proton SEC+ION composites. Source
`Elow` is a fit lower limit, not the ionization energy `I`; proton composite
energy in keV is not automatically the native mass-normalized eV/u coordinate.
The extra native composite A7 role is unresolved. Of106 printed-coefficient
rounding comparisons,104PASS and **2FAIL**, both halfway cases. Neither source
is repaired and the explanation is not assumed. These are retained source
crosswalk failures, distinct from new source-reader software checks.

The approved inherited-heavy-velocity **net reservoir** model can use the
existing scalar rate, physical state/Q, chosen heavy output and exact signed
number/charge/momentum/energy exchanges. It does not require a sampled joint
outgoing-electron law. For a stationary unpolarized target and isotropic,
zero-drift bath with no preferred collision axis, rotational symmetry gives
zero mean recoil vector. It does not determine variance/heating or moving-
target drag. Incoming moments can be derived mathematically from a rate's
Maxwellian transform, but differentiating a finite-domain approximate fit
does not produce measured moments or controlled uncertainty. The native Na
Voronov I=5.1eV is below the physical threshold, so an inferred microscopic
law cannot be silently truncated or renormalized. No new sampler is enabled.

For the approved ground-reference proton/Na partial, `Q=IE(H)-IE(Na)` and
conservative two-heavy-product kinematics are derivable from available state
energies. The all-angle bound is `|delta p|<=2*p_final`; an angular support
bound can tighten it to `2*p_final*sin(theta_max/2)`. The universal bound is
not evidence of small forward-scattering error. D119990's exact experimental
citation, initial preparation, uncertainty and prompt-versus-eventual-H1s
interpretation remain unresolved. A reference Na3s choice does not recover
those measurement facts. Total-capture and SEC+ION composites remain separate
from the selected H1s partial.

## 4. Surface and relaxation inputs

KIDA was followed into the **final** Minissale2022 review,
[10.1021/acsearthspacechem.1c00357](https://doi.org/10.1021/acsearthspacechem.1c00357).
It agrees with the paired OH first-order TST recommendation
`E/kB=5698K, nu=3.76e15s^-1` for the named ASW/olivine reference. The older
manuscript pair differs and is not mixed with the final version. These are
recommended/rederived kinetics, not new measured lunar impact probabilities.
KIDA zero-prefactor placeholders with instructions for derivation must not
become zero desorption rates.

Miyazaki2022, [10.3847/2041-8213/ac9d30](https://doi.org/10.3847/2041-8213/ac9d30),
gives1650±60K from64–76K diffusion-limited OH recombination on compact ASW.
This is not a desorption energy or a measurement of the absolute diffusion
prefactor, incoming-OH sticking or accommodation.

Schörghofer2023/2025,
[10.3847/PSJ/acf19b](https://doi.org/10.3847/PSJ/acf19b) and
[10.3847/PSJ/ade5b2](https://doi.org/10.3847/PSJ/ade5b2), supply a named
heterogeneous-site water-release law and identified author implementation.
Its coverage-dependent flux is available; implementation requires a selected
mobility/material/coverage regime and conservative population integration.
The source report derives a constant-temperature depletion oracle. This law
is not interchangeable with a single exponential residence time. The archive's
third-party code exceptions remain separate from reuse of an identified model.

An explicitly prescribed collisionless reference has exactly zero collisional
quenching. For nonzero baths, density/rate domain bounds and the populated
states can instead bound omission; no universal fixed error allowance is
invented. OH(X) rotational data do not determine OH(A) electronic quenching.
The approved ideal collector remains sufficient for its boundary reference
scope. A matching OH impact/accommodation law and extensions beyond measured
water sticking conditions still require evidence or a reviewed reduced model.

## 5. Consequences for implementation and qualification

Proceed with supported reference profiles, source/Q crosswalks, bulk reservoir
events, conservative reference CEX, selected surface inventory laws and actual
runner/scorer wiring. Freeze the remaining state/threshold/charged-breakup and
surface-model choices before activating the affected physical channels. The
new evidence does not authorize silent extrapolation, omitted channels,
instant relaxation, observational calibration or a changed acceptance limit.

Source/software verification, the104/2 coefficient audit, reduced-model
execution and formal I21–I24 have separate result scopes. Formal M3 still
requires original numerical/statistical/conservation gates and applicable
same-campaign M1/M2 predecessors. No M3 or observational PASS is claimed.
