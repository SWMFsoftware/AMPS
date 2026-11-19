"""Offline import of explicitly selected, authenticated physical-source packets.

Packets describe acquired source bytes, not executable physics.  Their policy
names select reviewed operations in source_inventory.py; arbitrary permission
booleans or a successful HTTP request never grant storage or model readiness.
Native quantities/units/state conventions remain strings/metadata.  This module
does not download, infer a missing state, convert units, evaluate a cross section
or modify an original.  See SOURCE_PACKET_CONTRACT.md for the pre-change contract.
"""
from __future__ import annotations

import copy
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path
from urllib.parse import urlsplit


_ID = re.compile(r"[A-Za-z0-9_.:-]+\Z")
_SHA = re.compile(r"[0-9a-f]{64}\Z")
# These are provider boundaries, not general licences for neighbouring sites.
# Each new reviewed policy also needs actual pinned licence/notice evidence.
_HOSTS = {
    "iaea-cc-by-4.0": {"db-amdis.org", "amdis.iaea.org"},
    "exomol-cc-by-sa-4.0": {"exomol.com", "www.exomol.com"},
    "nasa-cea-apache-2.0": {"raw.githubusercontent.com", "github.com", "nasa.github.io"},
    "zenodo-cc-by-4.0": {"zenodo.org"},
    "tuwien-regolith-cc-by-4.0": {"researchdata.tuwien.ac.at"},
    "chianti-use-attribution": {"sohoftp.nascom.nasa.gov", "db.chiantidatabase.org", "www.chiantidatabase.org"},
    "springer-cc-by-4.0": {"static-content.springer.com", "www.ebi.ac.uk"},
    "ntrs-19700032941-public-use": {"ntrs.nasa.gov"},
}


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path):
    def pairs(items):
        value = {}
        for key, item in items:
            if key in value:
                raise ValueError("duplicate source-packet JSON key: " + key)
            value[key] = item
        return value

    def invalid(token):
        raise ValueError("nonfinite source-packet JSON: " + token)

    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=pairs, parse_constant=invalid)


def _identifier(value, label):
    if not isinstance(value, str) or not _ID.fullmatch(value):
        raise ValueError("invalid source-packet " + label)
    return value


def _url(value):
    if not isinstance(value, str):
        raise ValueError("source URI must be an explicit HTTPS string")
    uri = urlsplit(value)
    if uri.scheme != "https" or not uri.hostname or uri.username or uri.password or uri.fragment:
        raise ValueError("invalid source-packet HTTPS URI: " + value)
    return uri


def _member(root, relative, expected):
    """Authenticate a real packet member and reject traversal/symlink escapes."""
    if not isinstance(relative, str) or not relative or "\\" in relative:
        raise ValueError("empty/nonportable source-packet member")
    part = Path(relative)
    if part.is_absolute() or ".." in part.parts:
        raise ValueError("source-packet member escapes its packet")
    path = (root / part).resolve()
    if root != path and root not in path.parents:
        raise ValueError("source-packet member symlink escapes its packet")
    if not isinstance(expected, str) or not _SHA.fullmatch(expected):
        raise ValueError("source-packet member requires lowercase SHA-256")
    if not path.is_file() or _sha(path) != expected:
        raise ValueError("missing/corrupt source-packet member: " + relative)
    return path


def _rights(root, evidence):
    if not isinstance(evidence, list):
        raise ValueError("source-packet rights evidence must be an explicit list")
    checked = []
    for item in evidence:
        if not isinstance(item, dict):
            raise ValueError("malformed source-packet rights evidence")
        _url(item.get("source_uri"))
        path = _member(root, item.get("path"), item.get("sha256"))
        checked.append(dict(item, original_path=str(path)))
    return checked


def load_packets(packet_paths, policies):
    """Produce ordinary inventory datasets without changing the default release.

    The caller supplies the maintained reviewed policy registry.  Duplicate
    source IDs, incorrect evidence, unknown policies or false activation fail
    before publication.  Unknown-policy *data* can use the explicit metadata-only
    policy, which retains its native identity but never copies numerical payload.
    """
    datasets, seen_packets, seen_sources = [], set(), set()
    for supplied in packet_paths:
        manifest = Path(supplied).resolve()
        root = manifest.parent
        value = _json(manifest)
        if not isinstance(value, dict) or value.get("schema") != "amps-chemical-source-packet-v1":
            raise ValueError("unsupported physical-source packet schema")
        if value.get("physical_activation", False) is not False:
            raise ValueError("source packet cannot activate physical chemistry")
        packet_id = _identifier(value.get("id"), "ID")
        if packet_id in seen_packets:
            raise ValueError("duplicate source-packet ID: " + packet_id)
        seen_packets.add(packet_id)
        for key in ("title", "provider", "version", "description"):
            if not isinstance(value.get(key), str) or not value[key]:
                raise ValueError("source packet requires explicit " + key)
        files = value.get("sources")
        if not isinstance(files, list) or not files:
            raise ValueError("source packet has no declared source records")
        dataset = {key: value[key] for key in ("id", "title", "provider", "version", "description")}
        dataset.update(expected_root=str(root), availability="present", files=[],
                       issues=list(value.get("limitations", [])), capabilities=["source-access"],
                       source_packet={"path": str(manifest), "sha256": _sha(manifest)},
                       rights=copy.deepcopy(policies["metadata-only"]))
        for row in files:
            if not isinstance(row, dict) or row.get("physical_activation", False) is not False:
                raise ValueError("source record cannot activate physical chemistry")
            source_id = _identifier(row.get("id"), "source ID")
            if source_id in seen_sources:
                raise ValueError("duplicate packet source ID: " + source_id)
            seen_sources.add(source_id)
            if not isinstance(row.get("logical_name"), str) or not row["logical_name"]:
                raise ValueError("source logical name is missing")
            uri = _url(row.get("source_uri"))
            policy_id = row.get("policy")
            if policy_id not in policies:
                raise ValueError("unreviewed source permission policy: " + str(policy_id))
            policy = copy.deepcopy(policies[policy_id])
            evidence = _rights(root, row.get("rights_evidence", []))
            if policy["storage_permission"] == "allowed":
                if not evidence or uri.hostname not in _HOSTS.get(policy_id, set()):
                    raise ValueError("allowed source requires matching provider and pinned rights evidence")
                if policy_id == "nasa-cea-apache-2.0" and not (
                    uri.path.startswith("/nasa/cea/") or
                    (uri.hostname == "nasa.github.io" and uri.path.startswith("/cea/"))):
                    raise ValueError("NASA CEA licence does not grant another repository's data")
                if policy_id == "ntrs-19700032941-public-use" and "19700032941" not in uri.path.split("/") and not any(
                        part.startswith("19700032941.") for part in uri.path.split("/")):
                    raise ValueError("NTRS public-use determination is specific to record 19700032941")
            role = row.get("role")
            if role not in ("source-native", "source-metadata", "rights-evidence", "derived-table"):
                raise ValueError("explicit source/metadata/derived role is required")
            path = _member(root, row.get("path"), row.get("sha256"))
            if role == "derived-table":
                if not row.get("derived_from") or not row.get("recipe"):
                    raise ValueError("derived source requires parent identities and a pinned recipe")
                recipe = row["recipe"]
                _member(root, recipe.get("path"), recipe.get("sha256"))
            else:
                if type(row.get("http_status")) is not int or row["http_status"] != 200:
                    # Retained HTTP errors belong in acquisition reports. A
                    # challenge/error page cannot become a successful source.
                    raise ValueError("native source did not return HTTP200: " + source_id)
                try:
                    timestamp = row["retrieved_utc"]
                    if not isinstance(timestamp, str):
                        raise ValueError("retrieval UTC must be a string")
                    instant = datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
                    if instant.utcoffset() is None:
                        raise ValueError("retrieval time lacks timezone")
                except (KeyError, TypeError, ValueError) as error:
                    raise ValueError("invalid source retrieval UTC") from error
            quantity = row.get("quantity_metadata")
            if not isinstance(quantity, dict):
                raise ValueError("explicit native quantity metadata is required, even when empty")
            entry = copy.deepcopy(row)
            entry.update(path=str(path), storage_permission=policy["storage_permission"],
                         evaluation_permission=policy["evaluation_permission"],
                         observed_sha256=row["sha256"], byte_count=path.stat().st_size,
                         capabilities=["source-access"], rights=policy,
                         rights_evidence=evidence, physical_activation=False,
                         scientific_status="SOURCE_INPUT_NOT_PHYSICAL_QUALIFICATION")
            if role != "derived-table":
                entry["packet_acquisition_records"] = [{
                    "url": row["source_uri"], "final_url": row.get("final_url", row["source_uri"]),
                    "retrieved_utc": row["retrieved_utc"], "http_status": 200,
                    "sha256": row["sha256"], "declared_sha256_matches_observed": True,
                    "manifest_path": str(manifest), "manifest_sha256": _sha(manifest),
                    "binding": "explicit-packet-relative-path-and-native-SHA256"}]
            if row.get("table"):
                table = row["table"]
                if not isinstance(table, dict):
                    raise ValueError("source table descriptor must be an object")
                _identifier(table.get("id"), "table ID")
                entry["table"] = copy.deepcopy(table)
            dataset["files"].append(entry)
        datasets.append(dataset)
    # Derived projections may retain source strings/rows, but cannot launder a
    # held parent's numerical data through an allowed output policy. All parent
    # identities must actually be in the explicitly selected packet family.
    entries = {entry["id"]: entry for dataset in datasets for entry in dataset["files"]}
    for entry in entries.values():
        if entry["role"] != "derived-table":
            continue
        parents = entry["derived_from"]
        if (not isinstance(parents, list) or not parents or
            any(not isinstance(parent, str) for parent in parents) or
            len(parents) != len(set(parents))):
            raise ValueError("derived table requires distinct original source IDs")
        for parent in parents:
            if parent not in entries or parent == entry["id"]:
                raise ValueError("derived table original source identity is missing/self-referential")
            if entries[parent]["role"] == "derived-table":
                raise ValueError("derived packet tables must cite native parents directly")
            if entry["storage_permission"] == "allowed" and (
                entries[parent]["storage_permission"] != "allowed" or
                entries[parent]["policy"] != entry["policy"]):
                raise ValueError("derived table cannot clear or change its native parent's source policy")
    return datasets
