# Explicit physical-source packet import contract

2026-11-19, before implementation. Source basis: `8fbbd5e`.
The user authorized online acquisition of remaining physical inputs and their
integration into the existing shared chemical database. Existing raw sources,
source permissions, release locks and scientific activation gates remain intact.

## Inventory and production entry points

`source_inventory.inventory` currently discovers fixed retained research
packets; `amps_chemical_database.build.build_database` publishes immutable
catalogue/source/table objects; `install_database.py` installs only into a new
or empty directory. Python/C++ readers consume an explicit release/profile lock.
The current installation is catalogue release `amps-chemistry-20261118-r1`.
It must not be overwritten. New acquired records will be published into a new
snapshot below `/data/vtenishe/CHEMICAL_DATABASE/snapshots/` and selected by its
explicit lock. Existing default lock, numerical services and source selections
remain unchanged.

## Proposed interface and configuration

Add repeatable `--source-packet PATH` publication arguments. Each path identifies
one strict JSON `amps-chemical-source-packet-v1` manifest, with stable packet and
source IDs, provider/version, description and explicit per-file roles, original
HTTPS URL, retrieval UTC, HTTP status, source hash, relative native path,
scientific quantity metadata, limitations and a named reviewed permission policy.

The importer is offline. It verifies packet-contained paths and every declared
hash, rejects duplicate IDs, corrupt or absent declared successful payloads,
nonfinite/duplicate JSON and unsupported permission policies. Failure documents
may be retained as such, without being made active scientific payloads. No
record may claim physical activation or infer model capability from lexical
source/table access. Official licence/notice/metadata evidence is itself pinned.

Policy names select inspected provider-specific operations, not caller-provided
permission booleans. IAEA CC BY4.0, ExoMol CC BY-SA4.0 and the explicitly licensed
NASA CEA Apache2.0 data may be imported only under their corresponding verified
source domains/evidence. Zenodo file grants require the native record licence
and exact published file binding. Unknown or restricted providers remain
metadata-only. New policy support requires recorded primary evidence, not a
blanket licence inferred from a download button or an unrelated publication.

Native units, species/state/material labels, energy/reference conventions,
temperature/energy/angle support, method (measurement/theory), uncertainties and
missing fields remain documentary source metadata. Source numbers are never
clipped, fitted, smoothed, extrapolated, renormalized or silently converted.
A source thermochemical/absorption/integral-rate table does not supply a complete
state-resolved differential event law. Existing physical gates stay disabled.

## Independent references and tests

The oracle is the exact acquired native source and its independently pinned
publisher metadata/licence, not current model output. Test real packet imports,
source hash/URL/provenance retention, denied/unknown-source refusal, immutable old
lock, strict parser and path/symlink containment, duplicate IDs, attempted false
physical activation, corrupt/stale/missing evidence and cross-machine installed
Python/C++ access. Clearly labelled administrative fixtures test infrastructure;
they are not observations or physical qualification. Test IDs U01-U26/I01-I46
are not changed. Shared database protocol/API tests and installed examples are
used; no production-physics build is required when only acquisition/catalogue
code changes.

Artifacts: acquired bytes and acquisition logs, source packet manifests, primary
rights evidence, import/QA manifests, result.json/stdout/stderr, new snapshot lock
and deep integrity checks, old release/lock hashes, README/API examples and an
explicit gap matrix. PASS means the declared data/software criterion was met;
missing physical support stays SKIPPED/NOT QUALIFIED. Errors remain preserved.

## Governing roadmaps

Source SHA256: `712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`.
Canonical data SHA256: `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
Approved source development/resource amendments govern implementation order;
unchanged scientific criteria, native interpretation and formal predecessors
continue to govern qualification. No roadmap is silently reconciled.

## Explicit recipe dependencies in later packets

A derived source can declare `quantity_metadata.recipe_dependencies` as a
list of packet-relative `path`/`sha256` objects. Each is an actual decoder/helper
imported by the pinned main recipe. Admission authenticates and confines these
members exactly as it does the main recipe; duplicate/self paths, traversal,
missing or changed dependencies fail before publication. The builder archives
dependencies at their original relative locations beside the recipe, and the
release QA independently checks their member/packet/original hashes. Historical
packets without this optional field retain their behavior. This field grants
no scientific capability or third-party data permission.

Main recipes and dependencies are explicit Python software, parsed for syntax
without execution. Neither a recipe nor a dependency may alias a declared
native/derived/rights source by resolved path or identical byte hash, including
another selected packet. Renaming a held payload to `.py` does not clear its
managed-storage restriction. This check protects provenance archives as well
as the ordinary source-object store.

The targeted physical-input extension also adds two narrowly reviewed policies.
Figshare13274840.v1 requires its exact native article/version/CC-BY4 grant and
listed file25563521; a different record/file/licence fails. NIST JRes121.005
requires its technical-publication notice and only the reviewed article's
official PDF/native PMC rendition. Original source-specific SRD/other-provider
holds remain unchanged. Administrative controls test these boundaries without
representing their fixture metadata as downloaded measurements.
