# I24 first-channel source review

The approved priority is proton capture from ground-reference Na. The
transport gate requires a consistent source for the rate, populated states,
energy convention, uncertainty and differential information or a usable
angle bound. See the [I24 acceptance contract](reference/acceptance.json).

| Source | Supported use | Admission result |
| --- | --- | --- |
| [IAEA D119990](https://db-amdis.org/collisiondb/datasets/119990/) | Retained catalog and extraction contain 12 reported H(1s) partial-capture knots over 4–25 keV/u. The locally retained record remains the usable source. The current web lookup did not resolve the live data page. | SKIPPED for transport. Original experimental lineage, actual preparation and uncertainty, prompt versus eventual H(1s), and angle coverage are unresolved. |
| [Anderson et al. (1979)](https://doi.org/10.1016/0029-554X(79)90644-X) | The published abstract describes total charge transfer and charge fractions for proton impact on Na vapor at 1–25 keV. Use it for total-capture source comparison where its measurement definitions match. | The total-charge measurements do not establish the lineage or the state-specific interpretation of D119990. |
| [Shingal et al. (1986)](https://doi.org/10.1088/0022-3700/19/3/014) | Independent atomic-orbital calculation for Na(3s) capture into H(n,l) up to n=3, with a reported lab domain of 0.5–20 keV/u. Its partial-rate calculation is a cross-check lead for the overlap. | Numerical tables and a usable angular/error envelope were not acquired. |
| [Courbin-Gaussorgues et al. (1990)](https://oa.upm.es/61369/) | The author repository identifies state-to-state differential calculations with separately defined Na(3s) and Na(3p), over 0.5–5 keV. | Preferred same-target fallback. The archive search returned the abstract. The full source page timed out during retrieval, and no machine-readable differential data or model-error bounds were acquired. |
| [Schöffler et al. (2009)](https://atom.uni-frankfurt.de/publications/files/Schoeffler2009pra.pdf) | The primary article's Figure 2 reports ground-state single capture for proton–helium at 60, 100, 150 and 300 keV/u on **lab** angle axes. Its broader 60–630 keV/u scope includes other projectile species. | A controlled proton–He alternative at 60 keV/u merits review if the source supplies numeric bins, normalization and angular tail coverage. The figure's plotting range is not a bound on the physical distribution. Lab-to-COM mapping and relativistic truncation must be assessed for the chosen case. |

Keep the ORNL helium and Morgan sodium total-capture rate checks. Any
candidate differential source must be applicable to the *same initial and
final states and collision energy* as the selected partial rate. A documented
worst-case angular bound could be used with the already approved forward
model only when that bound fits reviewed heavy-product speed and energy
limits. No current source established such an error envelope in the requested
I24 regime. Carry the proposed alternative to case-scope review before
selecting an incident energy or using its numerical values.

The r5 native D119990 object matched SHA-256
`ee669388646f4dfc54720e9d53ded92c3a853a041e5b37d2ea9c39bfed61f590`.
Its extracted table matched SHA-256
`2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93`
from the [retained source decision](../../M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md).
The newer publications were inspected as source leads. This review produced
zero new numerical collision rows and made no runtime profile selection.
Formal I24 remains SKIPPED.
