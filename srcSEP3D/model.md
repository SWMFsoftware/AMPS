# Coronal-CME model specification moved

This application-local file is a **non-normative redirect**. The canonical,
generated model description shared by `srcSEP3D` and field-aligned `srcSEP` is
[`../src/models/sep_coronal_cme/model.md`](../src/models/sep_coronal_cme/model.md).

Do not add application-local copies of the physics, equations, configuration
contract, APIs, or validation requirements here. Edit the maintained modules
under `../src/models/sep_coronal_cme/model/`, update the structured requirement
registry when needed, and regenerate the canonical document with:

```sh
python3 src/models/sep_coronal_cme/tools/generate_model.py
python3 src/models/sep_coronal_cme/tools/generate_model.py --check
```

The dependency boundary is intentional: `srcSEP3D` links the shared
`sep_coronal_cme` model through its AMPS-facing adapter, whereas `srcSEP`
consumes an immutable, versioned field-line bundle through neutral
`sep_common` records and does not link or invoke the three-dimensional model.
