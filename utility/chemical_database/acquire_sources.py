#!/usr/bin/env python3
"""Explicit HTTPS acquisition into a NEW research packet, never a runtime fetch.

The plan supplies exact URLs and portable output names. Native response bytes,
HTTP outcomes and timestamps are retained without scientific interpretation.
Failures remain failures, including server challenge/error bodies. Publication
is a separate offline, hash/rights-checked operation; downloading grants neither
managed-storage permission nor a complete physical reaction model.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen


def acquire(plan_path, output, timeout_seconds=30):
    """Retain bounded exact objects and a reproducible per-request receipt.

    ``max_bytes`` is a resource ceiling, not a scientific truncation rule. An
    incomplete body is explicitly ERROR/complete=False and is never admitted
    as a successful source. Existing packets/raw directories are not replaced.
    """
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    if output.exists() or "raw" in output.parts:
        raise ValueError("acquisition requires a new directory outside immutable raw trees")
    plan = json.loads(plan_path.read_text())
    if plan.get("schema") != "amps-chemical-acquisition-plan-v1":
        raise ValueError("unsupported acquisition plan")
    requests = plan.get("requests")
    if not isinstance(requests, list) or not requests:
        raise ValueError("explicit nonempty acquisition request list required")
    paths = set()
    for item in requests:
        uri = urlsplit(item["url"])
        path = Path(item["path"])
        limit = item.get("max_bytes", 50000000)
        if uri.scheme != "https" or not uri.hostname or uri.username or uri.password:
            raise ValueError("only explicit public HTTPS requests are supported")
        if path.is_absolute() or ".." in path.parts or "\\" in item["path"] or not item["path"]:
            raise ValueError("nonportable acquisition member path")
        if item["path"] in paths or type(limit) is not int or limit <= 0:
            raise ValueError("duplicate output path or invalid byte ceiling")
        paths.add(item["path"])
    output.mkdir(parents=True)
    shutil.copy2(plan_path, output / "acquisition_plan.json")
    records = []
    for item in requests:
        target = output / item["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        record = {"url": item["url"], "path": item["path"], "method": "GET",
                  "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                  "status": "ERROR", "complete": False, "physical_activation": False}
        digest, size = hashlib.sha256(), 0
        try:
            request = Request(item["url"], headers={"User-Agent": "AMPS-chemical-source-acquisition/1.0"})
            try:
                response = urlopen(request, timeout=timeout_seconds)
            except HTTPError as error:
                # Preserve the actual failed response; do not reinterpret a
                # 403/404/500 HTML body as the requested numerical table.
                response = error
            with response:
                record.update(http_status=response.code, final_url=response.geturl(),
                              headers={k: v for k, v in response.headers.items()
                                       if k.lower() in ("content-type", "content-length", "etag", "last-modified")})
                if urlsplit(record["final_url"]).scheme != "https":
                    raise ValueError("refuse a redirect from HTTPS to an insecure source")
                with target.open("xb") as stream:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        if size + len(chunk) > item.get("max_bytes", 50000000):
                            raise ValueError("source exceeded declared byte ceiling; reduce selection or increase explicit resource limit")
                        if shutil.disk_usage(output).free < len(chunk):
                            raise ValueError("insufficient storage; acquire into a larger filesystem")
                        stream.write(chunk)
                        digest.update(chunk)
                        size += len(chunk)
                record["complete"] = True
                expected = item.get("expected_sha256")
                if expected is not None and digest.hexdigest() != expected:
                    raise ValueError("download does not match independently supplied expected SHA-256")
                record["status"] = "PASS" if response.code == 200 else "ERROR"
        except (OSError, ValueError, URLError) as error:
            record["error"] = str(error)
        record.update(sha256=digest.hexdigest(), byte_count=size,
                      finished_utc=datetime.now(timezone.utc).isoformat())
        records.append(record)
        (output / "acquisition_manifest.json").write_text(json.dumps({
            "schema": "amps-chemical-explicit-acquisition-v1", "records": records,
            "source_plan_sha256": hashlib.sha256(plan_path.read_bytes()).hexdigest(),
            "status": "PASS" if all(r["status"] == "PASS" for r in records) else "ERROR",
            "scope": "HTTP/native-byte acquisition only; no physical qualification"}, indent=2) + "\n")
        print(item["path"] + ": " + record["status"], flush=True)
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=30)
    args = parser.parse_args(argv)
    if not math.isfinite(args.timeout_seconds) or args.timeout_seconds <= 0:
        parser.error("timeout must be positive")
    try:
        records = acquire(args.plan, args.output_dir, args.timeout_seconds)
        return 0 if all(r["status"] == "PASS" for r in records) else 2
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
