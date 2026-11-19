# He+/isotope-exchanged silica conditional charge-survival operator

This component implements **only Tellez etal2012 Eq3 and its seven TableI
source-coordinate fits**, under approved P09. It exercises
`species/HeliumSilicaChargeSurvival.{h,cpp}`, using the actual production SHA256
service. It enables no Moon boundary callback, BCA simulation, implantation,
sticking, lunar retention or residence law.

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_helium_silica.py \
  --source-pdf test_output/srcMoon/m3-ion-surface-cex-data-sources-20261117T225600Z/tellez2012_silica_neutralization.pdf \
  --output-dir test_output/srcMoon/helium-silica-NEW
```

Use a new output directory. `--cxx` selects an installed C++17 compiler;
`--wall-seconds` bounds the entire local component, including source
processing/compilation/queries/plots. No application configuration or full
build is performed. SourcePDF, compiler, pdftotext or plotting dependency
absence is SKIPPED. Source corruption, build/runtime/receipt errors are ERROR;
an actually evaluated criterion mismatch is FAIL. CLI exit codes0/1/2/77 mean
PASS/FAIL/ERROR/SKIPPED. PASS means local published-law/transcription controls,
not observational agreement or full ion-surface qualification.

The complete preimplementation contract is [CONTRACT.md](CONTRACT.md).
`prepare_source.py` reproduces `table_i.csv` from the exact original native
PDF using both `pdftotext -raw` and `-layout`. It preserves missing r²/ratio
cells, source inequalities, fit intercepts and author-reported uncertainties.
The known PDF character mappings were checked against the original rendered
TableI; no OCR or plot digitization is used. The hash-verified sourceTableCSV
has SHAbffa29b2d3a608f05210b01da52bb70fc42e64b0f4399390c40764ce967d8576.

| Target/source row | Published fit support intersected with study | Vc ± reported value, m/s |
|---|---|---|
| Si CIN |600<=Ei<=7000eV |291000±6000 |
|16O CIN |2000<Ei<=7000eV |226000±29000 |
|18O CIN |2000<Ei<=7000eV |267000±28000 |
|O(mean) CIN |2000<Ei<=7000eV |247000±28000 |
|16O AN |600<=Ei<800eV |188000±3000 |
|18O AN |600<=Ei<800eV |191000±2000 |
|O(mean) AN |600<=Ei<800eV |190000±2000 |

The oxygen[800,2000]eV transition interval remains **unavailable**, including
both endpoints. The caption mentions600–800eV for estimating the AN line;
the table's strict `<0.8keV` inequality controls this bounded implementation.
O(mean) is the separate source fit, not an assumed isotope mixture. The
source has no Si AN row. Energies outside the600–7000eV study are unavailable;
no fit extrapolation, crossover selection or coefficient interpolation occurs.

`EvaluatePublishedCoordinate` evaluates `P+=exp(-Vc*x)` with an explicit
positive source inverse-velocity coordinate x[s/m]. The yield intercept a is
metadata from Eq4, not a factor in this probability. The parameter-envelope
endpoints propagate only published Vc±; the paper does not supply a confidence
level/covariance here, so they are not labelled1sigma or95percent intervals.
The reference/candidateCSV and PNG/EPS compare actual C++ output with an
independent high-precision Decimal evaluation at declared mathematical
coordinate queries. Those queries are not experimental observations or an
invented conversion from the selected incident energy to x.

**Important unavailable physical mapping:** the paper explicitly says velocity
coordinates depend on neutralization mechanism, referring to Brongersma2007.
General theory distinguishes nonlocal normal-velocity and local atom-relative
scaling, but the exact projection used to estimate each numerical TableI
coefficient is not established by the available source metadata. Consequently
`EvaluateFromVelocityVectors` explicitly returns Unavailable. No normal/full
speed projection, isotope/Si mass, elastic scattered energy or angular
acceptance is guessed. This laboratory coordinate law must not be installed
as an actual lunar particle/surface outcome provider without further evidence.

The source experiment used normal4He+ incidence and145degree analyzer
scattering on isotope-exchanged silica; fitted O yields were extrapolated to
100% isotope composition. Its P+ is an inferred charged fraction among all
backscattered projectiles in the selected class, including survival and
reionization. It is not the fraction of incident projectiles retained. A
positive-ion-only instrument count is not an unconditional implantation or
reflection probability. Native raw spectra/fit covariance and an independent
laboratory comparison are not available here; they are not synthesized.

Artifacts include `source_evidence/` nativeTableCSV/QA/commands/tool versions,
compiled source/executable identities, actual query stdout/stderr and return
codes, local result/manifest, unittest log, `comparison.csv`, `comparison.png`
and `comparison.eps`. The physical velocity adapter and complete lunar
ion-surface behavior remain SKIPPED/NOTQUALIFIED, separately from passing
local numerical controls. New global runner registration is owned by the
parent task; this component assigns no new U/I identifier.
