"""Reproducible permission-aware publication of the AMPS chemical archive.

This builder consumes retained source bytes; it does not acquire, repair,
interpolate or reinterpret scientific numbers. Every source receives an
explicit disposition. Numerical storage is allowed only by the maintained
source inventory's evidence-bound policy. Unknown and restricted originals
remain outside the archive, with documentary identities recorded in JSON.
"""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import mimetypes
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from .api import (AVAILABILITIES, PERMISSIONS, SOURCE_COLUMNS, TABLE_COLUMNS,
                  Database, DatabaseError, _read_csv, member, sha256, strict_json)


def _json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False,
                               allow_nan=False) + "\n", encoding="utf-8")


def _csv(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(tuple(row.get(name, "") for name in columns) for row in rows)


def _load_inventory():
    """Find the maintained importer beside source tools or installed package."""
    local = Path(__file__).with_name("source_inventory.py")
    source = local if local.is_file() else Path(__file__).resolve().parents[2] / "source_inventory.py"
    if not source.is_file():
        raise DatabaseError("Unavailable", "maintained source_inventory.py is not installed")
    spec = importlib.util.spec_from_file_location("amps_chemical_source_inventory", str(source))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.inventory, source


def _load_acquisition():
    """Load the maintained metadata-only receipt reader, never a downloader."""
    local = Path(__file__).with_name("acquisition.py")
    source = local if local.is_file() else Path(__file__).resolve().parents[2] / "acquisition.py"
    if not source.is_file():
        raise DatabaseError("Unavailable", "maintained acquisition provenance reader is absent")
    spec = importlib.util.spec_from_file_location("amps_chemical_acquisition", str(source))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, source


def _load_iaea():
    local = Path(__file__).with_name("native_iaea.py")
    source = local if local.is_file() else Path(__file__).resolve().parents[2] / "native_iaea.py"
    if not source.is_file():
        raise DatabaseError("Unavailable", "maintained native IAEA reproduction recipe is absent")
    spec = importlib.util.spec_from_file_location("amps_native_iaea", str(source))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, source


def _iaea_processing(root, prefix, datasets, source_rows, table_rows, table_metadata):
    """Bind permitted derived CSVs to exact native responses and real recipes.

    D107581 is a composite fit/plot representation, with a Patch and a Line;
    D119990 is the distinct twelve-point measured partial. Reproduction checks
    preserve every renderer's native role and exact CSV bytes. They are source
    processing checks, never a fit, new observation or physical validation.
    """
    module, recipe = _load_iaea()
    files = {file["logical_name"]: file for dataset in datasets
             if dataset["id"] == "iaea-collisiondb-research" for file in dataset["files"]}
    source_by_id = {row["source_id"]: row for row in source_rows}
    result_members = {}
    for dataset_id, name in ((119990, "collisiondb_Na_H1s_119990.html"),
                             (107581, "collisiondb_Na3s_107581.html")):
        source = files.get(name)
        if not source or not source.get("imported_path"):
            continue
        directory = root / prefix / "provenance/reproductions" / str(dataset_id)
        records = module.reproduce(member(root, source["imported_path"]), directory, dataset_id)
        for record in records:
            original = files.get(record["output"])
            if not original or original.get("imported_sha256") != record["output_sha256"]:
                raise DatabaseError("CorruptSource", "native IAEA reconstruction does not match retained product: " + record["output"])
            table_id = original["table"]["id"]
            row = next(item for item in table_rows if item["table_id"] == table_id)
            row["source_ids"] = original["id"] + ";" + source["id"]
            row["notes"] += "; native renderer: {}; {}".format(
                record["native_renderer_role"], record["series_interpretation"])
            table_metadata[table_id]["processing"] = dict(record,
                source_id=source["id"], source_object_path=source_by_id[source["id"]]["path"],
                recipe_path=prefix + "/provenance/recipes/native_iaea.py")
        for path in directory.iterdir():
            result_members[str(path.relative_to(root))] = sha256(path)
    return result_members, recipe


def _copy_object(source, root, expected=None):
    """Copy exact bytes into content-addressed storage, checking disk space.

    Hash both sides. A declared expected native identity is checked before
    publication; a fresh observed digest alone is not called source approval.
    Reads/copies stream in bounded buffers and never follow old paths at run
    time. Original files are opened read-only and are never moved/deleted.
    """
    actual = sha256(source)
    if expected and actual != expected:
        raise DatabaseError("CorruptSource", "original identity changed: " + str(source))
    relative = "objects/sha256/{}/payload".format(actual)
    target = member(root, relative)
    if target.exists():
        if sha256(target) != actual:
            raise DatabaseError("CorruptSource", "existing content-addressed object is corrupt")
        return relative, actual
    size = Path(source).stat().st_size
    if shutil.disk_usage(root).free < size:
        raise DatabaseError("ResourceLimit", "insufficient free disk for {} ({} bytes). "
                            "Use a larger destination or a reviewed smaller source selection; "
                            "no source or table will be truncated.".format(source, size))
    target.parent.mkdir(parents=True)
    try:
        with Path(source).open("rb") as inp, target.open("xb") as output:
            shutil.copyfileobj(inp, output, 1024 * 1024)
    except OSError as error:
        raise DatabaseError("ResourceLimit", "source copy failed; preserved incomplete publication: " + str(error))
    if sha256(target) != actual:
        raise DatabaseError("CorruptSource", "copied source bytes do not match original")
    return relative, actual


def _source_revision(amps_root):
    result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(amps_root),
                            capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else "unavailable"


def _reaction_metadata(validation_root):
    """Index native channel identities without importing numerical rate data.

    This is a bibliographic/topological view of the existing reviewed channel
    catalogue. Product labels remain source strings, with repeated occurrences
    intact; no charge/state/isotope or energy crosswalk is inferred. In
    particular OH_exp identifies an experimental source variant, not an
    excited OH state. Chemical execution still requires explicit typed binding.
    """
    source = validation_root / "photochemistry/processed/phidrates/reaction_channels.csv"
    expected = "2d9fa1c8eb72d3bbba0b91ab4c9a912e6ebce2abc80b251096df2f41a4354061"
    if not source.is_file():
        return [], []
    if sha256(source) != expected:
        raise DatabaseError("CorruptSource", "frozen reaction identity catalogue changed")
    _, rows = _read_csv(source, ("species", "reaction_id", "reaction_type", "parent",
                                 "products", "phidrates_channel_label", "source_reference", "source_file"))
    species, reactions = {}, []
    for row in rows:
        variant, identifier, process, parent, products, label, reference, locator = row
        labels = products.split(" + ") if products != "TOTAL" else []
        for name in [parent] + labels:
            species.setdefault(name, {"id": "phidrates-label:" + name, "source_label": name,
                                      "isotope_composition": None, "resolved_state": None,
                                      "mass_definition": None, "qualification": "SOURCE_LABEL_ONLY"})
        reactions.append({"id": identifier, "source_variant": variant,
                          "process": process, "parent_source_label": parent,
                          "products_source_text": products,
                          "product_occurrences_source_labels": labels,
                          "reactive": process != "total_photoabsorption",
                          "native_channel_label": label, "reference": reference,
                          "original_source_locator": locator,
                          "source_catalogue_sha256": expected,
                          "rate_or_kinetic_data_inferred": False,
                          "qualification": "SOURCE_IDENTITY_METADATA_ONLY"})
    return list(species.values()), reactions


def _checked_permission(value):
    if value not in PERMISSIONS:
        raise DatabaseError("InvalidInput", "unknown source permission policy: " + str(value))
    return value


def _sdk(root, amps_root):
    """Install portable maintained source/headers and query tooling.

    A compiler-specific executable is deliberately not the only SDK artifact.
    Another machine can build the same standard-library C++ implementation.
    No symlink back to the developer checkout is introduced.
    """
    module = amps_root / "src/models/chemical_database"
    target = root / "sdk/src/models/chemical_database"
    if not module.is_dir():
        raise DatabaseError("Unavailable", "shared C++ implementation is not present")
    target.mkdir(parents=True)
    for path in sorted(module.iterdir()):
        if path.is_file() and path.suffix in (".h", ".cpp", ".md"):
            shutil.copy2(str(path), str(target / path.name))
    parser = amps_root / "src/general/input/amps_input.h"
    destination = root / "sdk/src/general/input/amps_input.h"
    destination.parent.mkdir(parents=True)
    shutil.copy2(str(parser), str(destination))
    examples = module / "examples"
    if examples.is_dir():
        shutil.copytree(str(examples), str(root / "sdk/examples"))
    package_source = Path(__file__).resolve().parent
    package_target = root / "python/amps_chemical_database"
    package_target.mkdir(parents=True)
    for path in package_source.glob("*.py"):
        shutil.copy2(str(path), str(package_target / path.name))
    _, importer = _load_inventory()
    shutil.copy2(str(importer), str(package_target / "source_inventory.py"))
    _, acquisition_source = _load_acquisition()
    shutil.copy2(str(acquisition_source), str(package_target / "acquisition.py"))
    # Source packet parsing is a maintained offline publication dependency,
    # never a downloader or a means to grant physical readiness at runtime.
    shutil.copy2(str(importer.with_name("source_packets.py")), str(package_target / "source_packets.py"))
    # New source-specific readers are also needed to reproduce the selected
    # lexical views after relocation. They do not run during ordinary queries
    # and cannot select/activate a physical process. Archive the same maintained
    # helpers that packet recipes import, including metadata-only amendments.
    source_helpers = (list(importer.parent.glob("native_*.py")) +
                      list(importer.parent.glob("new_*.py")) +
                      list(importer.parent.glob("amend_*.py")))
    for recipe in sorted(source_helpers):
        shutil.copy2(str(recipe), str(package_target / recipe.name))
    _, native_iaea_source = _load_iaea()
    shutil.copy2(str(native_iaea_source), str(package_target / "native_iaea.py"))
    rights = importer.with_name("RIGHTS.md")
    if rights.is_file():
        shutil.copy2(str(rights), str(root / "RIGHTS.md"))
    readme = importer.with_name("DATABASE_README.md")
    if not readme.is_file():
        raise DatabaseError("Unavailable", "detailed installation README template is absent")
    shutil.copy2(str(readme), str(root / "README.md"))
    # Ship the acquisition contract and source-specific caveats with the
    # portable installation. These are documentation, not a physical-profile
    # grant; users must not need the original checkout to discover data gaps.
    for note in [importer.with_name("SOURCE_PACKET_CONTRACT.md"), importer.with_name("PHYSICAL_MODELS_INTEGRATION_CONTRACT.md"), importer.with_name("PHYSICAL_INPUT_EXTENSION_CONTRACT.md")] + sorted(importer.parent.glob("NATIVE_*_SOURCES.md")) + sorted(importer.parent.glob("NEW_*_SOURCES.md")) + sorted(importer.parent.glob("PHYSICAL_INPUT_STATUS_*.md")):
        if note.is_file():
            shutil.copy2(str(note), str(root / note.name))
    bin_dir = root / "bin"
    bin_dir.mkdir()
    launcher = bin_dir / "chemdb"
    launcher.write_text("#!/usr/bin/env python3\n"
                        "# Portable installed entry; import only this installation.\n"
                        "from pathlib import Path\nimport sys\nsys.dont_write_bytecode=True\n"
                        "sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))\n"
                        "from amps_chemical_database.cli import main\n"
                        "raise SystemExit(main())\n", encoding="utf-8")
    launcher.chmod(0o755)


def build_database(amps_root, validation_root, output_root, release_id, source_packets=(),
                   include_model_definitions=False):
    """Publish a new immutable snapshot; refuse overwriting any existing root.

    Partial work survives under output_root on an ERROR for diagnosis. Only a
    verified release receives database-lock.json and publication_status=READY.
    Here READY means archive integrity, never scientific/model validation.
    """
    amps_root, validation_root = Path(amps_root).resolve(), Path(validation_root).resolve()
    root = Path(output_root).resolve()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", release_id):
        raise DatabaseError("InvalidInput", "release ID must be a portable identifier")
    if root.exists():
        raise DatabaseError("InvalidInput", "refusing to overwrite existing publication root: " + str(root))
    root.mkdir(parents=True)
    started = datetime.now(timezone.utc).isoformat()
    manifest = {"schema": "amps-chemical-publication-v1", "started_utc": started,
                "publication_status": "INCOMPLETE", "source_revision": _source_revision(amps_root),
                "release_id": release_id, "source_root": str(amps_root),
                "original_data_root": str(validation_root), "physical_activation": False,
                "scientific_validation": "NOT EVALUATED", "errors": []}
    _json(root / "publication.json", manifest)
    try:
        importer, importer_path = _load_inventory()
        acquisition, acquisition_path = _load_acquisition()
        acquisition_index = acquisition.build_acquisition_index(amps_root, validation_root)
        if acquisition_index.get("issues"):
            raise DatabaseError("CorruptSource", "retained acquisition metadata has unresolved parser/binding issues: " +
                                "; ".join(acquisition_index["issues"]))
        native_datasets = importer(amps_root, validation_root, source_packets=source_packets)
        # Authored software definitions are opt-in and have their own AMPS
        # software licence. They cannot clear or copy a held external source.
        # No download receipt is invented for locally authored JSON/CSV.
        models, authored_members = (), {}
        if include_model_definitions:
            from .model_definitions import materialize
            dataset, models = materialize(root, release_id, amps_root)
            native_datasets.append(dataset)
            authored_members = {str(path.relative_to(root)): sha256(path)
                for path in Path(dataset["expected_root"]).iterdir() if path.is_file()}
        manifest["include_model_definitions"] = bool(include_model_definitions)
        manifest["explicit_source_packets"] = [
            {"path": str(Path(path).resolve()), "sha256": sha256(path)} for path in source_packets]
        datasets, source_rows, table_rows, table_metadata = [], [], [], {}
        all_source_ids, all_table_ids = set(), set()
        for native_dataset in native_datasets:
            # Metadata remains documentary. Copy a defensive plain structure
            # so this builder cannot mutate the source inventory's policy.
            dataset = json.loads(json.dumps(native_dataset, default=str))
            dataset_files = []
            for original in dataset.get("files", []):
                source_id = original["id"]
                if source_id in all_source_ids or not source_id:
                    raise DatabaseError("InvalidInput", "duplicate inventory source ID: " + source_id)
                all_source_ids.add(source_id)
                source_path = Path(original["path"])
                storage = _checked_permission(original.get("storage_permission", "unknown"))
                evaluation = _checked_permission(original.get("evaluation_permission", "unknown"))
                present = source_path.is_file()
                observed = sha256(source_path) if present else None
                expected = original.get("sha256") or original.get("expected_sha256")
                if expected and present and observed != expected:
                    raise DatabaseError("CorruptSource", "inventory original hash mismatch: " + str(source_path))
                relative, digest = "", ""
                if present and storage == "allowed":
                    relative, digest = _copy_object(source_path, root, expected)
                    availability = "available"
                else:
                    availability = "metadata-only" if present else "unavailable"
                # Canonical metadata retains historical originals and observed
                # identity even when no numerical payload may be imported.
                entry = dict(original, original_path=str(source_path),
                             observed_sha256=observed, imported_path=relative or None,
                             imported_sha256=digest or None, availability=availability,
                             byte_count=source_path.stat().st_size if present else None)
                receipts = list(original.get("packet_acquisition_records", []))
                receipts += acquisition.lookup_acquisition(acquisition_index, source_path, observed) if observed else []
                entry["acquisition_records"] = receipts
                entry["acquisition_metadata_status"] = "native_receipts_retained" if receipts else "not_supplied_or_not_a_downloaded_object"
                # All native receipts survive, including mismatches/failures.
                # A unique exactly bound successful object's URL is useful in
                # the small runtime index; ambiguity never chooses a winner.
                urls = {receipt.get("final_url") or receipt.get("url") for receipt in receipts
                        if receipt.get("declared_sha256_matches_observed") is True}
                urls.discard(None)
                runtime_uri = next(iter(urls)) if len(urls) == 1 else original.get("source_uri", "")
                dataset_files.append(entry)
                source_rows.append({"source_id": source_id, "path": relative, "sha256": digest,
                                    "media_type": original.get("media_type", mimetypes.guess_type(str(source_path))[0] or "application/octet-stream"),
                                    "availability": availability, "storage_permission": storage,
                                    "evaluation_permission": evaluation,
                                    "capabilities": ";".join(original.get("capabilities", ["source-access"])),
                                    "source_uri": runtime_uri,
                                    "notes": "Original: {}; role: {}; imported without science-value edits: {}".format(
                                        original.get("logical_name", source_path.name), original.get("role", "reference"), bool(relative))})
                if original.get("table"):
                    table = original["table"]
                    table_id = table["id"]
                    if table_id in all_table_ids:
                        raise DatabaseError("InvalidInput", "duplicate table ID: " + table_id)
                    all_table_ids.add(table_id)
                    kind = table.get("format", "opaque")
                    columns, count = (), "unknown"
                    if kind == "csv" and present and availability == "available":
                        columns, rows = _read_csv(source_path, table.get("columns") or None)
                        if any(";" in column for column in columns):
                            raise DatabaseError("UnsupportedSchema", "CSV column cannot be encoded by runtime codec: " + table_id)
                        count = str(len(rows))
                        declared = table.get("row_count", "unknown")
                        if declared != "unknown" and int(declared) != len(rows):
                            raise DatabaseError("CorruptSource", "declared source CSV row count differs: " + table_id)
                    elif kind == "csv":
                        columns = tuple(table.get("columns", ()))
                        count = str(table.get("row_count", "unknown"))
                        if (count == "unknown" or not columns or len(set(columns)) != len(columns) or any(not column or
                                ";" in column or "\r" in column or "\n" in column for column in columns)):
                            # No row count is invented for an absent source.
                            # The same holds for held/restricted or incomplete
                            # responses: do not parse them as active tables or
                            # fabricate an empty valid CSV on their behalf.
                            kind, columns, count = "opaque", (), "unknown"
                    row = {"table_id": table_id, "path": relative, "sha256": digest,
                           "format": kind, "schema_id": table.get("schema_id", "native-source-v1"),
                           "source_ids": source_id, "availability": availability,
                           "evaluation_permission": evaluation, "capabilities": "table-access",
                           "columns": ";".join(columns), "row_count": count,
                           "notes": table.get("description", dataset.get("description", ""))}
                    table_rows.append(row)
                    table_metadata[table_id] = {key: value for key, value in table.items()
                                                if key not in ("id", "columns", "row_count", "format", "schema_id")}
                    table_metadata[table_id]["dataset_id"] = dataset["id"]
            dataset["files"] = dataset_files
            dataset["availability_summary"] = {name: sum(item["availability"] == name for item in dataset_files)
                                                for name in sorted(AVAILABILITIES)}
            datasets.append(dataset)
        if len({item["id"] for item in datasets}) != len(datasets):
            raise DatabaseError("InvalidInput", "duplicate dataset identifiers")
        prefix = "releases/" + release_id
        processing_members, native_iaea_source = _iaea_processing(
            root, prefix, datasets, source_rows, table_rows, table_metadata)
        catalogue_path = prefix + "/catalogue/catalogue.json"
        species, reactions = _reaction_metadata(validation_root)
        catalogue = {"schema": "amps-chemical-catalogue-v1", "datasets": datasets,
                     "species": species, "reactions": reactions, "table_metadata": table_metadata,
                     "scientific_validation": "NOT EVALUATED", "physical_activation": False,
                     "identity_policy": "native source assertions retained; no inferred state/isotope/reaction equivalence"}
        if include_model_definitions:
            catalogue["models"] = models
        _json(root / catalogue_path, catalogue)
        # The first publication is a comprehensive inventory/data-access
        # selection, not a chosen physical reaction network. No implicit
        # source/model ranking or complete-product capability is declared.
        profile_id = "catalogue-v1"
        profile_path = prefix + "/profiles/" + profile_id + ".json"
        profile = {"schema": "amps-chemical-profile-v1", "id": profile_id,
                   "scope": "source catalogue and permission-cleared native table access",
                   "datasets": [item["id"] for item in datasets],
                   "capabilities": ["catalogue", "source-access", "table-access"],
                   "topology_fingerprint": None, "physical_activation": False,
                   "scientific_validation": "NOT EVALUATED"}
        _json(root / profile_path, profile)
        members = {catalogue_path: sha256(root / catalogue_path),
                   profile_path: sha256(root / profile_path)}
        members.update(processing_members)
        members.update(authored_members)
        # Preserve the exact admitted packet declarations, not just their
        # historical filesystem paths. Source objects and recipes are stored
        # separately; this manifest records the original relative bindings,
        # source URLs, permissions and limitations for independent replay.
        # Check against the pre-import digest to detect a changed declaration
        # during publication instead of silently blessing its new contents.
        packet_recipes = []
        for index, packet in enumerate(manifest["explicit_source_packets"]):
            relative = prefix + "/provenance/source_packets/{:03d}.json".format(index)
            destination = root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(packet["path"], str(destination))
            if sha256(destination) != packet["sha256"]:
                raise DatabaseError("CorruptSource", "source packet changed during publication: " + packet["path"])
            packet["archived_path"] = relative
            members[relative] = packet["sha256"]
            # The declaration pins a recipe, but a historical absolute source
            # path is not a portable copy of that recipe. Retain its exact
            # reviewed bytes without executing it. Common maintained parser
            # dependencies are separately copied into provenance/recipes.
            # Source-packet admission has already checked path confinement and
            # the same SHA; member() and a second hash check close a race with
            # local edits between admission and publication.
            selected = strict_json(packet["path"])
            seen_recipe_paths = set()
            for source in selected["sources"]:
                # Preserve dependencies at their original packet-relative
                # locations as well as the main recipe. Admission has checked
                # every declared hash and confined each path to the packet.
                declarations = ([source["recipe"]] if source.get("recipe") else [])
                declarations += source.get("quantity_metadata", {}).get("recipe_dependencies", [])
                for recipe in declarations:
                    if recipe["path"] in seen_recipe_paths:
                        continue
                    seen_recipe_paths.add(recipe["path"])
                    original = member(Path(packet["path"]).parent, recipe["path"])
                    archived = prefix + "/provenance/packet_recipes/{:03d}/".format(index) + recipe["path"]
                    destination = root / archived
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(str(original), str(destination))
                    if sha256(destination) != recipe["sha256"]:
                        raise DatabaseError("CorruptSource", "source-packet recipe changed during publication")
                    members[archived] = recipe["sha256"]
                    packet_recipes.append(dict(packet_sha256=packet["sha256"],
                        original_path=recipe["path"], archived_path=archived, sha256=recipe["sha256"]))
        for source in source_rows:
            if source["path"]:
                members[source["path"]] = source["sha256"]
        provenance_path = prefix + "/provenance/import.json"
        # Preserve the exact publisher/policy implementation with the release,
        # not merely its current checkout location. These are AMPS software
        # sources, not a permission exception for a scientific provider's data.
        recipe_prefix = prefix + "/provenance/recipes"
        recipe_sources = list(Path(__file__).resolve().parent.glob("*.py")) + [
            importer_path, acquisition_path, native_iaea_source, importer_path.with_name("RIGHTS.md"),
            importer_path.with_name("DATABASE_README.md"), importer_path.with_name("source_packets.py"),
            importer_path.with_name("SOURCE_PACKET_CONTRACT.md")]
        recipe_sources += [importer_path.with_name("PHYSICAL_MODELS_INTEGRATION_CONTRACT.md")]
        recipe_sources += [importer_path.with_name("PHYSICAL_INPUT_EXTENSION_CONTRACT.md")]
        recipe_sources += list(importer_path.parent.glob("native_*.py")) + [importer_path.with_name("acquire_sources.py")]
        recipe_sources += list(importer_path.parent.glob("new_*.py")) + list(importer_path.parent.glob("amend_*.py"))
        recipe_sources += list(importer_path.parent.glob("NATIVE_*_SOURCES.md")) + list(importer_path.parent.glob("PHYSICAL_INPUT_STATUS_*.md"))
        recipe_sources += list(importer_path.parent.glob("NEW_*_SOURCES.md"))
        recipe_sources = list(dict.fromkeys(recipe_sources))
        for recipe in recipe_sources:
            destination = root / recipe_prefix / recipe.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(str(recipe), str(destination))
            members[str(destination.relative_to(root))] = sha256(destination)
        _json(root / provenance_path, {"schema": "amps-chemical-import-v1",
                                     "source_inventory_sha256": sha256(importer_path),
                                     "acquisition_reader_sha256": sha256(acquisition_path),
                                     "native_iaea_recipe_sha256": sha256(native_iaea_source),
                                     "builder_sha256": sha256(__file__),
                                     "source_revision": manifest["source_revision"],
                                     "explicit_source_packets": manifest["explicit_source_packets"],
                                     "packet_recipes": packet_recipes,
                                     "include_model_definitions": bool(include_model_definitions),
                                     "raw_originals_modified": False,
                                     "transformations": "byte-preserving payload import; source channel identity/occurrence metadata; independent IAEA native typed-array CSV reproduction; explicitly admitted native-token CSV projections with packet-pinned parent/recipe/output hashes; no scientific value conversion or repair",
                                     "missing_fields_inferred": False,
                                     "dataset_ids": profile["datasets"]})
        members[provenance_path] = sha256(root / provenance_path)
        release = {"schema": "amps-chemical-release-v1", "release_id": release_id,
                   "catalogue": {"path": catalogue_path, "sha256": members[catalogue_path]},
                   "members": dict(sorted(members.items())),
                   "source_revision": manifest["source_revision"],
                   "physical_activation": False, "scientific_validation": "NOT EVALUATED",
                   "release_scope": "native source inventory and explicitly permission-cleared payload only"}
        release_path = prefix + "/release.json"
        _json(root / release_path, release)
        release_digest = sha256(root / release_path)
        profile_digest = members[profile_path]
        export_prefix = "exports/{}/{}/csv-block-v1".format(release_digest, profile_digest)
        source_index_path = export_prefix + "/source_index.csv"
        table_index_path = export_prefix + "/tables_index.csv"
        _csv(root / source_index_path, SOURCE_COLUMNS, source_rows)
        _csv(root / table_index_path, TABLE_COLUMNS, table_rows)
        export_values = {"schema": "amps-chemical-runtime-v1", "codec": "csv-block-v1",
                         "release_id": release_id, "release_manifest_path": release_path,
                         "release_manifest_sha256": release_digest, "profile_id": profile_id,
                         "profile_path": profile_path, "profile_sha256": profile_digest,
                         "source_index_path": source_index_path,
                         "source_index_sha256": sha256(root / source_index_path),
                         "table_index_path": table_index_path,
                         "table_index_sha256": sha256(root / table_index_path),
                         "profile_capabilities": ";".join(profile["capabilities"]),
                         "topology_fingerprint": "none"}
        export_path = export_prefix + "/export.in"
        (root / export_path).write_text("#begin chemical_database: runtime\n" +
                                       "".join("  {} = {}\n".format(key, value) for key, value in export_values.items()) +
                                       "#end chemical_database\n", encoding="utf-8")
        lock = {"schema": "amps-chemical-database-lock-v1", "database_root": ".",
                "release_id": release_id, "release_manifest_sha256": release_digest,
                "profile_id": profile_id, "profile_sha256": profile_digest,
                "export_path": export_path, "export_manifest_sha256": sha256(root / export_path)}
        database = Database.open(root, release_id, release_digest, profile_id,
                                 profile_digest, export_path, lock["export_manifest_sha256"])
        verification = database.verify(deep=True)
        _json(root / prefix / "qa/integrity.json", verification)
        # QA is detached documentary evidence, not a member claiming its own
        # digest. The acyclic canonical release/export identity remains frozen.
        _sdk(root, amps_root)
        software_paths = [path for directory in (root / "python", root / "sdk", root / "bin")
                          for path in directory.rglob("*") if path.is_file()]
        software_paths += [root / "README.md", root / "RIGHTS.md"]
        _json(root / "software_manifest.json", {
            "schema": "amps-chemical-software-manifest-v1",
            "source_revision": manifest["source_revision"],
            "files": {str(path.relative_to(root)): sha256(path) for path in sorted(software_paths)},
            "scope": "installed reader/SDK/documentation source identities; not a scientific approval"})
        checksums = dict(members, **{release_path: release_digest,
                                    source_index_path: sha256(root / source_index_path),
                                    table_index_path: sha256(root / table_index_path),
                                    export_path: lock["export_manifest_sha256"]})
        (root / "SHA256SUMS").write_text("".join("{}  {}\n".format(digest, name)
                                                     for name, digest in sorted(checksums.items())), encoding="utf-8")
        manifest.update(publication_status="READY", finished_utc=datetime.now(timezone.utc).isoformat(),
                        archive_integrity=verification, dataset_count=len(datasets),
                        source_count=len(source_rows), table_count=len(table_rows),
                        permitted_payload_count=sum(bool(item["path"]) for item in source_rows),
                        metadata_only_count=sum(item["availability"] == "metadata-only" for item in source_rows),
                        lock=lock)
        _json(root / "publication.json", manifest)
        # The convenience lock is the final publication marker. An incomplete
        # SDK/README or failed verification never leaves a default runnable
        # lock that could be mistaken for a completed installation.
        _json(root / "database-lock.json", lock)
        return manifest
    except Exception as error:
        manifest["errors"].append({"code": getattr(error, "code", "InfrastructureError"), "message": str(error)})
        _json(root / "publication.json", manifest)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps-root", type=Path, required=True)
    parser.add_argument("--validation-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--source-packet", type=Path, action="append", default=[],
                        help="explicit authenticated offline source packet; repeat to supplement the historical catalogue")
    parser.add_argument("--include-model-definitions", action="store_true",
                        help="include authored application-independent model definitions; does not activate physical chemistry")
    args = parser.parse_args(argv)
    try:
        print(json.dumps(build_database(args.amps_root, args.validation_root,
                                        args.output_root, args.release_id,
                                        source_packets=args.source_packet,
                                        include_model_definitions=args.include_model_definitions), indent=2))
        return 0
    except (DatabaseError, OSError, ValueError, KeyError) as error:
        print(json.dumps({"status": "ERROR", "code": getattr(error, "code", "InfrastructureError"),
                          "message": str(error)}, indent=2), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
