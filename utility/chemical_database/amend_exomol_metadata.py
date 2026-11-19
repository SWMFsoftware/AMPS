"""Reproduce a version-bound MYTHOS metadata correction in a NEW packet.

An earlier JSON projection compared ``states.max_energy`` with the largest
state energy. Matching-version native text instead identifies the quantity as
the higher energy with a complete set of transitions (cm^-1). These are
different quantities. This tool copies authenticated old members unchanged,
amends only their interpretation, and cites an independently acquired native
definition by source ID. It never edits an old packet/release or a science CSV.
"""
from __future__ import annotations

import argparse
import copy
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import shutil

try:
    from .source_packets import _json, _member
except ImportError:
    from source_packets import _json, _member


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def declaration(path, label):
    """Require one exact native label; never infer semantics from a JSON key."""
    found = []
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if '#' not in line:
            continue
        value, native_label = (part.strip() for part in line.split('#', 1))
        if native_label == label:
            found.append({'line': number, 'native_value': value, 'native_label': native_label})
    if len(found) != 1:
        raise ValueError('one exact native definition label required: ' + label)
    return found[0]


def amend(packet, expected_packet_sha256, evidence_packet, expected_evidence_sha256,
          evidence_source_id, output):
    """Authenticate both packets and copy old bytes before a metadata amendment.

    All paths remain packet-relative. Expected hashes are caller-supplied
    identities, not a fresh digest promoted to scientific approval. Native
    source admission/rights are still enforced by the ordinary publisher.
    """
    packet, evidence_packet, output = map(lambda p: Path(p).resolve(),
                                         (packet, evidence_packet, output))
    if output.exists():
        raise ValueError('amendment output must be a new directory')
    if sha(packet) != expected_packet_sha256 or sha(evidence_packet) != expected_evidence_sha256:
        raise ValueError('selected packet identity changed')
    old, proof = _json(packet), _json(evidence_packet)
    candidates = [row for row in proof['sources'] if row['id'] == evidence_source_id]
    if len(candidates) != 1:
        raise ValueError('exact native evidence source ID required')
    evidence = candidates[0]
    if evidence.get('role') not in ('source-native', 'source-metadata') or evidence.get('http_status') != 200:
        raise ValueError('evidence must be an actual successful native source')
    native = _member(evidence_packet.parent, evidence['path'], evidence['sha256'])
    version = declaration(native, 'Version number with format YYYYMMDD')
    completeness = declaration(native, 'Higher energy with complete set of transitions (in cm-1)')
    try:
        cutoff = Decimal(completeness['native_value'])
    except InvalidOperation as error:
        raise ValueError('invalid native completeness value') from error
    if not cutoff.is_finite():
        raise ValueError('nonfinite native completeness value')
    states = [row for row in old['sources'] if row['id'] == 'exomol-20261119:oh-mythos-states-bz2:lexical-csv']
    if len(states) != 1:
        raise ValueError('this amendment requires the original MYTHOS state view')
    checks = states[0]['quantity_metadata']['source_checks']
    if str(checks['source_dataset_version']) != version['native_value']:
        raise ValueError('native evidence and state dataset versions differ')
    if Decimal(checks['definition_max_energy_native_cm-1']) != cutoff:
        raise ValueError('native evidence and retained JSON quantity differ')
    # Validate/copy every selected member, including rights notices and recipe
    # bytes. No hardlinks: later scratch-file edits must not affect originals.
    members = {}
    for row in old['sources']:
        for member in [row] + list(row.get('rights_evidence', [])) + ([row['recipe']] if row.get('recipe') else []):
            path = _member(packet.parent, member['path'], member['sha256'])
            if member['path'] in members and members[member['path']][1] != member['sha256']:
                raise ValueError('inconsistent member hash in original packet')
            members[member['path']] = (path, member['sha256'])
    corrected = copy.deepcopy(old)
    corrected['version'] = str(old['version']) + '-metadata-correction-r4'
    corrected['limitations'] = [item for item in old['limitations']
        if not item.startswith('OH definition max_energy=')]
    corrected['limitations'].append(
        'Matching-version native .def identifies JSON max_energy as transition completeness, '
        'not the maximum state energy; states above it are retained and are not a defect.')
    selected = next(row for row in corrected['sources'] if row['id'] == states[0]['id'])
    amended = selected['quantity_metadata']['source_checks']
    amended.pop('maximum_energy_matches_definition', None)
    amended['transition_completeness_energy_native_cm-1'] = completeness['native_value']
    amended['native_definition_evidence_source_id'] = evidence_source_id
    corrected['metadata_amendment'] = {
        'scope': 'interpretation only; every native/derived science payload unchanged',
        'previous_packet_sha256': expected_packet_sha256,
        'native_evidence_packet_sha256': expected_evidence_sha256,
        'native_evidence_source_id': evidence_source_id,
        'native_evidence_sha256': evidence['sha256'],
        'version_declaration': version, 'completeness_declaration': completeness,
        'recipe_sha256': sha(__file__), 'physical_activation': False}
    output.mkdir(parents=True, exist_ok=False)
    for relative, (path, digest) in members.items():
        target = output/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(path), str(target))
        if sha(target) != digest:
            raise ValueError('source changed during amendment copy')
    target = output/'source_packet.json'
    target.write_text(json.dumps(corrected, indent=2) + '\n')
    return {'status': 'PASS', 'scope': 'metadata-only source correction',
            'packet': str(target), 'sha256': sha(target), 'unchanged_copied_members': len(members),
            'metadata_amendment': corrected['metadata_amendment'], 'physical_activation': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('packet', 'evidence-packet', 'output-dir'):
        parser.add_argument('--'+name, type=Path, required=True)
    for name in ('expected-packet-sha256', 'expected-evidence-sha256', 'evidence-source-id'):
        parser.add_argument('--'+name, required=True)
    args = parser.parse_args()
    result = amend(args.packet, args.expected_packet_sha256, args.evidence_packet,
                   args.expected_evidence_sha256, args.evidence_source_id, args.output_dir)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
