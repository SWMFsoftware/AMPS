"""Read-only discovery of native chemistry acquisition receipts.

Preimplementation contract (2026-11-18)
-------------------------------------
``build_acquisition_index(amps_root, validation_root)`` reads ONLY named
acquisition/retrieval metadata in the explicit retained chemistry packets.
``lookup_acquisition(index, source_path, observed_sha256)`` returns defensive
copies of ALL exact-path-matched records. Absence is an empty list, not a
guessed URL or timestamp. Duplicate/nonfinite JSON is documentary corruption
recorded in ``issues``; no entire mixed manifest or numerical array is copied.

Each returned record contains a small whitelist of native acquisition fields,
the manifest's absolute original path, SHA-256 and JSON pointer, and the exact
file-binding rule. A comparison with an observed source digest is diagnostic
only: it NEVER selects an authoritative expected hash, resolves permissions,
approves a source, repairs a file or establishes scientific readiness. HTTP
failure responses remain failure responses, even if their bytes match.

The independent checks use existing IAEA, NIST and ExoMol receipts and explicit
path/hash/absence invariants. Canonical files and source mtimes must remain
unchanged. These are provenance checks, not physical or observational tests.
There is no network, downloader, clock inference or write operation here.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
import re
from typing import Any


# Keep discovery bounded to sources already inventoried for this release.
# In particular, do not descend into the O I packet's intentionally corrupted
# source-controls/global-entry fixtures and mistake their fake receipts for
# acquisition evidence. These are original retained packet locations, not
# runtime paths that a published archive will subsequently follow.
_PACKETS = (
    'test_output/srcMoon/m3-atomic-online-databases-20261118T083400Z',
    'test_output/srcMoon/m3-molecular-databases-20261118T083551Z',
    'test_output/srcMoon/m3-ion-surface-cex-data-sources-20261117T225600Z',
    'test_output/srcMoon/m3-approved-products-20261118T163142Z/source-acquisition-v2',
    'test_output/srcMoon/m3-charge-exchange-research-20261117',
    'test_output/srcMoon/m3-kinematic-references-20261117',
)
_OXYGEN = 'test_output/srcMoon/m3-oxygen-cascade-sources-20261118T170845Z'
_DIGEST = re.compile(r'[0-9a-fA-F]{64}\Z')
_URL_FIELDS = ('url', 'final_url', 'original_url', 'requested_url',
               'request_url', 'source_url')
_FILE_FIELDS = ('file', 'path', 'filename', 'name', 'native_file', 'output')
_TIME_FIELDS = ('retrieved_utc', 'retrieval_utc', 'utc', 'started_utc',
                'finished_utc', 'completed_utc', 'created_utc',
                'start_utc', 'end_utc')
_FIELDS = _URL_FIELDS + _FILE_FIELDS + _TIME_FIELDS + (
    'http_status', 'status', 'content_type', 'method', 'sha256',
    'native_sha256', 'output_sha256')


def _strict_json(raw: bytes) -> Any:
    def unique(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate JSON field: ' + key)
            result[key] = value
        return result

    def invalid(token):
        raise ValueError('nonfinite JSON token: ' + token)

    return json.loads(raw.decode('utf-8'),
                      object_pairs_hook=unique, parse_constant=invalid)


def _is_metadata(path: Path) -> bool:
    """Only explicitly named receipts, never arbitrary scientific JSON."""
    name = path.name.lower()
    return name.endswith('.json') and (
        name.endswith('.metadata.json') or name == 'metadata.json' or
        name in ('acquisition.json', 'run_manifest.json',
                 'cross_sections_response.json', 'just_fotout_response.json') or
        re.fullmatch(r'brnout_[0-9]+\.response\.json', name) is not None or
        'retrieval' in name or ('download' in name and 'manifest' in name))


def _manifests(root: Path) -> list[Path]:
    if not root.is_dir():
        return []
    result = []
    for directory, children, names in os.walk(root, followlinks=False):
        children[:] = sorted(name for name in children if name not in
            {'.venv', '__pycache__', '.git', '.pytest_cache', 'node_modules'})
        result.extend(Path(directory) / name for name in sorted(names)
                      if _is_metadata(Path(name)))
    return result


def _pointer_part(text: str) -> str:
    return text.replace('~', '~0').replace('/', '~1')


def _records(value: Any, pointer: str = ''):
    """Yield dictionary locations without carrying arbitrary parent payload.

    A dictionary is a retrieval record only if it contains its own URL field.
    A parent's URL is never silently assigned to every child/source. Arrays,
    spectra, fit parameters, cookies, form bodies and header collections are
    traversed only to locate dictionaries; none are retained in the index.
    """
    if isinstance(value, dict):
        if any(isinstance(value.get(key), str) and value[key]
               for key in _URL_FIELDS):
            yield pointer, value
        for key, child in value.items():
            yield from _records(child, pointer + '/' + _pointer_part(key))
    elif isinstance(value, list):
        for number, child in enumerate(value):
            yield from _records(child, pointer + '/' + str(number))


def _whitelist(record: dict) -> dict:
    # Keeping native names avoids claiming that a package created_utc was the
    # retrieval time, or that a textual acquisition status was an HTTP code.
    selected = {key: record[key] for key in _FIELDS if key in record and
                isinstance(record[key], (str, int, float, bool, type(None)))}
    # Some receipts put only content type in their native headers. Copy that
    # one documentary scalar, never the whole headers/cookie/session object.
    if 'content_type' not in selected:
        headers = record.get('response_headers', record.get('headers'))
        pairs = headers.items() if isinstance(headers, dict) else headers
        if isinstance(pairs, (list, tuple)) or isinstance(headers, dict):
            values = [pair[1] for pair in pairs
                      if isinstance(pair, (list, tuple)) and len(pair) == 2 and
                      isinstance(pair[0], str) and pair[0].lower() == 'content-type' and
                      isinstance(pair[1], str)]
            if len(set(values)) == 1:
                selected['content_type'] = values[0]
    return selected


def _paths(record: dict, manifest: Path, packet: Path,
           amps: Path) -> list[tuple[str, Path, str]]:
    """Bind explicit native paths; do not recover paths from URLs or digests.

    Relative record paths may be manifest-relative or packet-relative, and
    the O I acquisition explicitly writes its ``file`` under ``raw/``.
    An AMPS-root-relative test_output path remains such a path. If different
    existing files satisfy those declarations, retain no ambiguous binding.
    Absolute historical paths are never retargeted by basename.
    """
    bindings = []
    for field in _FILE_FIELDS:
        text = record.get(field)
        if not isinstance(text, str) or not text or '\x00' in text:
            continue
        target = Path(text)
        if target.is_absolute():
            bindings.append((field, target.resolve(), 'explicit-absolute-path'))
            continue
        candidates = [(manifest.parent / target, 'manifest-relative-path'),
                      (packet / target, 'packet-relative-path')]
        if manifest.name == 'acquisition.json':
            candidates.append((manifest.parent / 'raw' / target,
                               'acquisition-raw-relative-path'))
        if target.parts and target.parts[0] == 'test_output':
            candidates.append((amps / target, 'amps-relative-path'))
        existing = {}
        for candidate, method in candidates:
            if candidate.is_file():
                existing.setdefault(candidate.resolve(), method)
        if len(existing) == 1:
            path, method = next(iter(existing.items()))
            bindings.append((field, path, method))
        elif len(existing) > 1:
            # The caller stores this explicit ambiguity as a manifest issue.
            bindings.append((field, None, 'ambiguous-existing-relative-path'))
        else:
            # A failed download can name an absent original. Preserve its
            # literal manifest-relative identity without inventing bytes.
            bindings.append((field, (manifest.parent / target).resolve(),
                             'literal-missing-manifest-relative-path'))
    return bindings


def _sidecar_paths(manifest: Path, canonical: Path) -> list[tuple[str, Path, str]]:
    """Use ONLY existing producer-defined canonical sidecar conventions.

    fetch_electron_impact.py and fetch_phidrates.py documentation use
    payload.ext.metadata.json. PHIDRATES also explicitly writes plot_NNN.png,
    brnout_NNN.raw and solar_rates/.../result.html beside their metadata.
    These rules follow the retained producer scripts, not URL extensions or
    matching numerical hashes. Other undocumented sidecars remain unbound.
    """
    try:
        relative = manifest.relative_to(canonical)
    except ValueError:
        return []
    if not relative.parts or relative.parts[0] != 'raw':
        return []
    candidates = []
    name = manifest.name
    if name.endswith('.metadata.json'):
        stem = name[:-len('.metadata.json')]
        candidates.append((manifest.with_name(stem), 'producer-extension-sidecar'))
        if re.fullmatch(r'plot_[0-9]+', stem):
            candidates.append((manifest.with_name(stem + '.png'), 'phidrates-plot-sidecar'))
        if re.fullmatch(r'brnout_[0-9]+', stem):
            candidates.append((manifest.with_name(stem + '.raw'), 'phidrates-brnout-sidecar'))
    elif name == 'metadata.json' and 'solar_rates' in relative.parts:
        candidates.append((manifest.with_name('result.html'), 'phidrates-solar-response-sidecar'))
    elif manifest.parent.name == 'network':
        # The source producer also saves the same response receipt under
        # network/. These exact names have explicit payload destinations in
        # fetch_phidrates.py; numerical response JSON is not parsed here.
        if name in ('cross_sections_response.json', 'just_fotout_response.json'):
            candidates.append((manifest.parent.parent / 'result.html',
                               'phidrates-network-response-sidecar'))
        elif re.fullmatch(r'brnout_[0-9]+\.response\.json', name):
            stem = name[:-len('.response.json')]
            candidates.append((manifest.parent.parent / 'data' / (stem + '.raw'),
                               'phidrates-network-brnout-sidecar'))
    existing = [(path.resolve(), method) for path, method in candidates if path.is_file()]
    if len(existing) != 1:
        return []
    path, method = existing[0]
    return [('sidecar', path, method)]


def build_acquisition_index(amps_root: Path, validation_root: Path) -> dict:
    """Build one deterministic index; missing known packets are documented.

    Returned paths are original documentary locators. Publication must retain
    them as provenance, never follow them during later runtime recovery.
    Only manifests are hashed here; the builder supplies the actual payload
    digest to ``lookup_acquisition`` independently.
    """
    amps = Path(amps_root).resolve()
    canonical = Path(validation_root).resolve() / 'photochemistry'
    roots = [amps / item for item in _PACKETS]
    roots.extend(amps / _OXYGEN / child for child in ('raw', 'acquired', 'final-replay'))
    roots.append(canonical)
    index = dict(schema='amps-chemical-acquisition-index-v1', records_by_path={},
                 scanned_manifests=[], issues=[])
    seen = set()
    for root in roots:
        if not root.is_dir():
            index['issues'].append('MISSING_ACQUISITION_PACKET: ' + str(root))
            continue
        for manifest in _manifests(root):
            manifest = manifest.resolve()
            if manifest in seen:
                continue
            seen.add(manifest)
            try:
                # Parse the exact byte snapshot whose digest is retained.
                # Two independent reads could bind parsed fields to an older
                # digest if another process replaced a receipt in between.
                raw = manifest.read_bytes()
                descriptor = dict(path=str(manifest), sha256=hashlib.sha256(raw).hexdigest())
                index['scanned_manifests'].append(descriptor)
                document = _strict_json(raw)
            except (OSError, UnicodeError, ValueError) as error:
                index['issues'].append('UNREADABLE_ACQUISITION_METADATA: ' +
                                       str(manifest) + ': ' + str(error))
                continue
            context = {key: document[key] for key in _TIME_FIELDS
                       if isinstance(document, dict) and key in document and
                       isinstance(document[key], str)}
            for pointer, record in _records(document):
                selected = _whitelist(record)
                bindings = _paths(record, manifest, root, amps)
                if not bindings and pointer == '':
                    bindings = _sidecar_paths(manifest, canonical)
                for field, target, method in bindings:
                    if target is None:
                        index['issues'].append('AMBIGUOUS_ACQUISITION_FILE: ' +
                            str(manifest) + '#' + pointer + ': ' + str(record[field]))
                        continue
                    receipt = dict(selected,
                        source_manifest=dict(descriptor, json_pointer=pointer),
                        source_path=str(target), file_binding_field=field,
                        file_binding_method=method)
                    if context:
                        receipt['manifest_time_context'] = dict(context)
                    group = index['records_by_path'].setdefault(str(target), [])
                    if receipt not in group:
                        group.append(receipt)
    index['scanned_manifests'].sort(key=lambda row: row['path'])
    index['issues'] = sorted(set(index['issues']))
    return index


def lookup_acquisition(index: dict, source_path: Path,
                       observed_sha256: str) -> list[dict]:
    """Return ALL matching receipts and a diagnostic declared-hash comparison.

    Match filename/path first. A matching digest never causes a receipt for a
    different file to be attached. Conversely, a mismatching receipt is kept
    (not silently replaced by a later match). No native URL/time/HTTP outcome
    is invented for missing provenance, and no caller can mutate the index by
    editing the returned dictionaries.
    """
    if not isinstance(index, dict) or index.get('schema') != 'amps-chemical-acquisition-index-v1':
        raise ValueError('unsupported acquisition index schema')
    if not isinstance(observed_sha256, str) or not _DIGEST.fullmatch(observed_sha256):
        raise ValueError('observed source SHA-256 must contain exactly 64 hex digits')
    result = copy.deepcopy(index['records_by_path'].get(str(Path(source_path).resolve()), []))
    for record in result:
        field = record['file_binding_field']
        declared = (record.get('output_sha256') if field == 'output' else
                    record.get('native_sha256') if field == 'native_file' else
                    record.get('sha256'))
        if declared is None:
            record['declared_sha256_matches_observed'] = None
        elif not isinstance(declared, str) or not _DIGEST.fullmatch(declared):
            record['declared_sha256_matches_observed'] = None
            record['declared_sha256_format'] = 'INVALID_NATIVE_DECLARATION'
        else:
            record['declared_sha256_matches_observed'] = declared.lower() == observed_sha256.lower()
    return result
