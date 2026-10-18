# Activity-enhancement coefficient file, version 1

This is the normative exchange format for the optional Stage-A3 C++ H2O
activity enhancement. It complements, and does not replace, the legacy
three-by-25 `ActivityFit.txt` format.

The file is UTF-8/ASCII text. Blank lines and text after `#` are ignored. All
other lines have `Key = Value` syntax. Metadata keys occur exactly once;
`Coefficient` is repeated. Unknown or duplicate fields, trailing tokens,
nonfinite numbers, missing/duplicate/unexpected modes, and incompatible
versions are fatal before pressure propagation begins.

```text
Format = liouville_cops.activity_enhancement
FormatVersion = 1
BasisID = liouville_cops.real_orthonormal_cs_v1
Species = H2O
BodyFrame = 67P/C-G_CK
AngleConvention = theta=acos(z/r);phi=atan2(y,x)
Ordering = increasing_l_then_m_minus_l_to_plus_l
Clipping = max(0,legacy_linear+enhancement)
MinimumDegree = 5
MaximumDegree = 5
Scale = 1.0
ScaleInterpretation = dimensionless_multiplier_of_normalized_basis
LegacyActivityIdentity = builtin:shap5
Regularization = not_applicable
Lambda = 0.0
Coefficient = 5 -5 0.0
Coefficient = 5 -4 0.0
# ...exactly one record for every order through (5,+5)
```

The basis equations and mode order are defined in `HARMONIC_BASIS.md`. The
applied field is

```text
DeltaA(x) = Scale * sum_(l=5..L) sum_(m=-l..l) b_lm Y_lm(x)
A(x) = max(0, A0_linear(x) + DeltaA(x)).
```

`Scale` is a positive, dimensionless numerical conditioning factor. It is not
a fitted physical observable. The coefficients and `DeltaA` are dimensionless
activity factors in the existing source multiplier. Clipping occurs only once,
after adding the residual to the unclipped legacy H2O field. CO2 and CO are not
accepted.

`MaximumDegree = 4` is the explicit zero-length control: it requires no
`Coefficient` records and produces exactly zero residual. Degrees 5, 6, 7, and
8 require respectively 11, 24, 39, and 56 records.

`LegacyActivityIdentity` prevents applying a residual to the wrong baseline.
It must equal `builtin:shap5` or `builtin:shap5_1` for an embedded map. For
`ActivityModel = file`, it must be `sha256:` followed by the lowercase SHA-256
of the resolved legacy activity file. `BodyFrame` must exactly match the
active SPICE body frame.

`Regularization` and `Lambda` describe how a fitter obtained the coefficients;
the forward model records but does not reinterpret or apply them. Lambda must
be finite and nonnegative. Stage A4 writes
`area_weighted_unclipped_residual_to_zero_relative_legacy_linear_power`, whose
exact definition is in `README.md`; any other producer must write its actual
regularization identifier rather than claiming a convention that was not used.
