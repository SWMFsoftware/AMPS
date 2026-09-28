# DOCSCCM01

from `src/models/sep_coronal_cme/test`, validates `model/requirements.yaml`, enforces the registered numbered-section/module map, generates the review-disposition and trace tables, checks the existence of every declared section/configuration/API/test anchor, rejects duplicate registry IDs, requirements without a review-disposition link, duplicate or out-of-Section-17.2 canonical test definitions, unlinked or unknown test links, incomplete roadmap ranges, undeclared cross-stage reuse, and malformed source structure, regenerates `model.md`, and requires a byte-clean generated diff. It does not infer semantic relevance from token presence; that remains a scientific review responsibility. This model-neutral ID is owned by the shared library and is included in both applications' aggregate release evidence.

Run from any directory with `python3 test.py`. The launcher delegates to the global registry, so individual and cumulative gates execute the identical implementation.
