"""Select one reviewed M2 contract without conflating three test meanings.

The .in may name a direct I15/I19/I20 document or a frozen contract-set whose
per-ID descriptors bind separate files. This is reference dispatch only: it
does not approve a child, certify data, or make a prerequisite pass.
"""
from __future__ import annotations

from datetime import datetime
import hashlib
import json
from pathlib import Path


class Unavailable(RuntimeError):
    """Missing/unreviewed scientific contract; callers map this to SKIPPED."""


def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):digest.update(block)
    return digest.hexdigest()


def read_json(path):
    """Reject duplicate keys and nonfinite numbers in scientific metadata."""
    def pairs(values):
        result={}
        for key,value in values:
            if key in result:raise ValueError("duplicate contract key: "+key)
            result[key]=value
        return result
    def invalid(value):raise ValueError("nonfinite contract number: "+value)
    return json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=invalid)


def resolve_contract(path,test_id,*,unavailable=Unavailable):
    """Return (selected path, immutable selection evidence).

    Relative child paths use the set document's directory. No default child
    path, implied review or path-only descriptor is accepted. A caller may
    supply its own unavailable exception to preserve existing status plumbing.
    Individual contracts still undergo their own full scientific preflight.
    """
    if test_id not in ("I15","I19","I20"):raise ValueError("unknown M2 formal contract owner")
    path=Path(path).resolve()
    if not path.is_file():raise unavailable("M2 formal contract unavailable: "+str(path))
    value=read_json(path)
    if not isinstance(value,dict):raise ValueError("M2 formal contract must be an object")
    evidence={"selected_test_id":test_id,"input_document":{"path":str(path),"sha256":sha256(path)}}
    if value.get("schema")!="srcMoon-M2-formal-contract-set-v1":
        # The selected test's own schema validator checks a direct document;
        # do not guess ownership from filenames or coerce a different schema.
        return path,evidence
    review=value.get("review",{})
    if not isinstance(review,dict):raise ValueError("M2 contract-set review must be an object")
    if value.get("status")!="APPROVED" or review.get("approved") is not True:
        raise unavailable("M2 formal contract set is not explicitly reviewed")
    if not review.get("reviewer") or not review.get("approved_utc"):
        raise ValueError("M2 contract-set review requires reviewer and approved_utc")
    stamp=datetime.fromisoformat(str(review["approved_utc"]).replace("Z","+00:00"))
    if stamp.tzinfo is None:raise ValueError("contract-set approval needs a time zone")
    contracts=value.get("contracts",{})
    if not isinstance(contracts,dict):raise ValueError("M2 contract-set contracts must be an object")
    record=contracts.get(test_id)
    if record is None:raise unavailable("M2 formal contract set has no "+test_id+" document")
    if not isinstance(record,dict) or set(record)!={"path","sha256"}:
        raise ValueError("M2 child descriptor requires exactly path and sha256")
    if not isinstance(record["path"],str) or not record["path"]:
        raise ValueError("M2 child path must be a nonempty string")
    target=Path(record["path"])
    target=(target if target.is_absolute() else path.parent/target).resolve()
    expected=record["sha256"]
    if not isinstance(expected,str) or len(expected)!=64 or any(c not in "0123456789abcdef" for c in expected):
        raise ValueError("invalid M2 child SHA-256")
    if not target.is_file():raise unavailable("M2 child contract unavailable: "+str(target))
    if target==path:raise ValueError("M2 contract-set self-reference")
    if sha256(target)!=expected:raise ValueError("M2 child contract SHA-256 mismatch")
    child=read_json(target)
    if not isinstance(child,dict):raise ValueError("M2 child contract must be an object")
    if child.get("schema")=="srcMoon-M2-formal-contract-set-v1":
        raise ValueError("nested M2 contract sets are not supported")
    evidence["selected_document"]={"path":str(target),"sha256":expected}
    return target,evidence
