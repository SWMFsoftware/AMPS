#!/usr/bin/env python3
"""Exact native-source and reversible token access to Jones2020 Apollo TPD.

See NATIVE_APOLLO_SOURCES.md for the contract written before implementation.
These are real author-deposited measurements, inversions and model outputs;
their roles remain distinct. The tool does not implement adsorption physics,
repair source headers or assign missing physical units. It has no network or
AMPS dependencies. PASS means authenticated byte-preserving source access.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import stat


# These fixed identities were independently checked against the publisher-linked
# author deposit, including every producer MD5 and byte count. A later changed
# file/version requires source review; the current disk contents never retarget
# an authoritative expectation automatically.
SOURCE_RECORD = 'https://zenodo.org/records/3822365'
PUBLICATION_DOI = '10.1029/2019JE006147'
SOURCE_FILES = {
    "Figure1_left.txt": {
        "bytes": 22905,
        "md5": "4a1871d5275c2c5610c8cff9957b5a49",
        "sha256": "bd1a83743887e847ad62128b3dfcac5ade6ae5e8bf59f0890a7a281a79e5c134",
        "source_uri": "https://zenodo.org/api/records/3822365/files/Figure1_left.txt/content"
    },
    "Figure1_right.txt": {
        "bytes": 32723,
        "md5": "6c413e7d690f0857033095aadbb670dd",
        "sha256": "a55b62019217aa8f2d2916036625d093e83ea5f0fba65387df15ce68bab4eda5",
        "source_uri": "https://zenodo.org/api/records/3822365/files/Figure1_right.txt/content"
    },
    "figure2.txt": {
        "bytes": 36937,
        "md5": "10ca415b5e890652f6793591729e4c2b",
        "sha256": "037f289722e0e073a5f8dbfaa5fea42808ce388c8e4cffb717a464fb702ae415",
        "source_uri": "https://zenodo.org/api/records/3822365/files/figure2.txt/content"
    },
    "figure3.txt": {
        "bytes": 31528,
        "md5": "57d58e0a08dd5ee38ba6b7b6a2d2811d",
        "sha256": "b7d07e06961eeb3ee8c08fa86800e990fc26e0b65b775834adbb34d8b73e438d",
        "source_uri": "https://zenodo.org/api/records/3822365/files/figure3.txt/content"
    },
    "figure3inset.txt": {
        "bytes": 35355,
        "md5": "a1e9518950d8f1db4857d8228763d001",
        "sha256": "6ab7c237e0858fd9c86a40964cc650006b1524b1e082b6676c90964488e7dee3",
        "source_uri": "https://zenodo.org/api/records/3822365/files/figure3inset.txt/content"
    },
    "figure4.txt": {
        "bytes": 141735,
        "md5": "2349411c787f4977f64264c422d17a8b",
        "sha256": "5c77af11d573848e702695274d498cc26ec01c2e01fa7be6bf7603570fc57c82",
        "source_uri": "https://zenodo.org/api/records/3822365/files/figure4.txt/content"
    },
    "figure5bottom.txt": {
        "bytes": 91536,
        "md5": "23857ea9fc29f99b6da000a73044d457",
        "sha256": "53b9bcf61b2d8fe4a37634ee520211fee615254d49609e2229a7ff9de05f8ac5",
        "source_uri": "https://zenodo.org/api/records/3822365/files/figure5bottom.txt/content"
    },
    "figure5middle.txt": {
        "bytes": 73270,
        "md5": "9dc3c9a41ea814cab255ec909def7fe5",
        "sha256": "f4ba54c3323026a5910444a0eec2076bab7532372eebdbcefca5212821bd78b1",
        "source_uri": "https://zenodo.org/api/records/3822365/files/figure5middle.txt/content"
    },
    "figure6.txt": {
        "bytes": 113259,
        "md5": "1c7bcb62ddf80fbef84bed2213012c60",
        "sha256": "3d3789bdcd883927dc2cf4abbeeeb8ee2054c36831bfbe96066d4123e6961060",
        "source_uri": "https://zenodo.org/api/records/3822365/files/figure6.txt/content"
    },
    "figure7.txt": {
        "bytes": 30298,
        "md5": "de42e2b1834d7f5134b4428a4f427659",
        "sha256": "413ff90cd4fe8ac2d7caf8d461211c820b4393b189304f42fa6c7c80f4356d94",
        "source_uri": "https://zenodo.org/api/records/3822365/files/figure7.txt/content"
    },
    "figure7inset.txt": {
        "bytes": 45687,
        "md5": "65e9efc8bcb9f5e97a9233e9f6dfd86a",
        "sha256": "ce89053a6953e93fca798819bc1a26247ebe2de97fcfb54a645f35bf02d2e7c1",
        "source_uri": "https://zenodo.org/api/records/3822365/files/figure7inset.txt/content"
    },
    "figure8bottom.txt": {
        "bytes": 527,
        "md5": "1cb3916ff9a89939035dd4ed621af11f",
        "sha256": "162b6547b4238522de00d22b7b1c7af9feafe092dcfa3769be128c9217448358",
        "source_uri": "https://zenodo.org/api/records/3822365/files/figure8bottom.txt/content"
    },
    "figure8top.txt": {
        "bytes": 483,
        "md5": "116a04cd0690f7507c9a389519a07af6",
        "sha256": "9676284e498f49b1f501015158407fb8fc7299aaa746f79114335e6a854b6eff",
        "source_uri": "https://zenodo.org/api/records/3822365/files/figure8top.txt/content"
    },
    "readme.txt": {
        "bytes": 551,
        "md5": "cf2e73e1fbd85c64bfa9bd907e2ef181",
        "sha256": "62ce1864dde0dbfc4981e97f3c5192619f8fcf6940e97d5274b491aa33770011",
        "source_uri": "https://zenodo.org/api/records/3822365/files/readme.txt/content"
    }
}

CSV_COLUMNS = ('source_line_1based', 'source_field_1based', 'line_role',
               'native_token', 'native_line_ending')
SCHEMA = 'apollo-tpd-native-line-field-tokens-v1'
_ENDINGS = {'LF': '\n', 'CRLF': '\r\n', 'CR': '\r', 'NONE': ''}
_ENDING_NAMES = {value: key for key, value in _ENDINGS.items()}
_ROLES = {
    'Figure1_left.txt': 'measured-background-and-background-subtracted-TPD',
    'Figure1_right.txt': 'extracted-TPD-and-fitted-low-temperature-signal',
    'figure2.txt': 'measured-TPD-at-distinct-water-exposures',
    'figure3.txt': 'model-inverted-coverage-dependent-activation-energy',
    'figure3inset.txt': 'model-inverted-site-energy-distribution',
    'figure4.txt': 'mixed-measured-TPD-and-fitted-component-simulations',
    'figure5middle.txt': 'mixed-anneal-TPD-and-readsorption-recombination-simulations',
    'figure5bottom.txt': 'mixed-anneal-TPD-and-readsorption-recombination-simulations',
    'figure6.txt': 'mixed-measured-TPD-and-fitted-component-simulations',
    'figure7.txt': 'model-inverted-coverage-dependent-activation-energy',
    'figure7inset.txt': 'model-inverted-site-energy-distribution',
    'figure8top.txt': 'lunar-cycle-model-prediction-not-observation',
    'figure8bottom.txt': 'lunar-cycle-model-prediction-not-observation',
}


class NativeApolloError(ValueError):
    """Changed source identity or malformed lexical projection; no repair."""


class NativeApolloUnavailable(FileNotFoundError):
    """Only missing required authoritative sources, distinct from runtime I/O."""


def lexical_csv(raw: bytes) -> bytes:
    """Expose original tab fields with positions and native line endings.

    Source number strings are never parsed as float, so signed values, leading
    zeros, exponent spelling and blanks survive. The first native line remains
    data in the projection: we do NOT make its ambiguous names into CSV columns.
    Only CR/LF terminate native lines; Python's broader Unicode splitlines()
    convention would change valid token bytes such as vertical tabs.
    """
    if not isinstance(raw, bytes):
        raise NativeApolloError('native source must be bytes')
    text = raw.decode('utf-8')
    output = io.StringIO(newline='')
    writer = csv.writer(output, lineterminator='\n')
    writer.writerow(CSV_COLUMNS)
    line = 0
    for match in re.finditer(r'([^\r\n]*)(\r\n|\r|\n|$)', text):
        # The regex has an empty terminal match; this is not an extra source
        # line. A real trailing blank line has a native newline and is kept.
        if match.start() == len(text):
            continue
        line += 1
        native, ending = match.groups()
        role = 'native-header' if line == 1 else 'native-blank' if native == '' else 'native-data'
        for field, token in enumerate(native.split('\t'), 1):
            writer.writerow((line, field, role, token, _ENDING_NAMES[ending]))
    return output.getvalue().encode('utf-8')


def restore_native(projected: bytes) -> bytes:
    """Independently usable strict inverse of the documented lexical schema.

    This function does not authorize arbitrary bytes as scientific input: the
    source-specific reproduce() below additionally requires fixed identities.
    Changed indices, mixed line roles/endings, gaps, extra fields or a newline
    embedded in a source field are errors, not silently sorted/repaired rows.
    """
    if not isinstance(projected, bytes):
        raise NativeApolloError('native-token projection must be bytes')
    try:
        # strict=True refuses malformed quoting rather than treating an
        # incomplete transport object as a repaired lexical table.
        all_rows = list(csv.reader(io.StringIO(projected.decode('utf-8'), newline=''), strict=True))
    except csv.Error as error:
        raise NativeApolloError('malformed native-token CSV encoding') from error
    if not all_rows or all_rows[0] != list(CSV_COLUMNS):
        raise NativeApolloError('unexpected native-token CSV schema')
    output, tokens = [], []
    current_line, current_role, current_ending = 0, None, None

    def finish():
        if not tokens:
            return
        expected = 'native-header' if current_line == 1 else 'native-blank' if tokens == [''] else 'native-data'
        if current_role != expected:
            raise NativeApolloError('native line role contradicts its actual fields')
        if tokens == [''] and current_ending == 'NONE':
            raise NativeApolloError('spurious empty terminal source line')
        output.append('\t'.join(tokens) + _ENDINGS[current_ending])

    for values in all_rows[1:]:
        if len(values) != len(CSV_COLUMNS):
            raise NativeApolloError('ragged native-token CSV row')
        row = dict(zip(CSV_COLUMNS, values))
        if any(not re.fullmatch(r'[1-9][0-9]*', row[name])
               for name in ('source_line_1based', 'source_field_1based')):
            raise NativeApolloError('source indices must be canonical positive integers')
        line, field = int(row['source_line_1based']), int(row['source_field_1based'])
        ending, role, token = row['native_line_ending'], row['line_role'], row['native_token']
        if ending not in _ENDINGS or any(c in token for c in ('\r', '\n', '\t')):
            raise NativeApolloError('invalid line ending or delimiter inside native token')
        if line != current_line:
            if line != current_line + 1 or (current_line and current_ending == 'NONE'):
                raise NativeApolloError('source line order/gap or nonterminal NONE ending')
            finish()
            tokens = []
            current_line, current_role, current_ending = line, role, ending
        if field != len(tokens) + 1 or role != current_role or ending != current_ending:
            raise NativeApolloError('source field order or inconsistent repeated line metadata')
        tokens.append(token)
    finish()
    return ''.join(output).encode('utf-8')


def _authenticate(source_dir):
    """Read each regular original exactly once and verify the complete family."""
    source = Path(source_dir)
    if not source.is_dir():
        raise NativeApolloUnavailable('authoritative Apollo source directory unavailable: ' + str(source))
    missing = [name for name in SOURCE_FILES if not (source / name).exists()]
    if missing:
        raise NativeApolloUnavailable('authoritative Apollo source member unavailable: ' + ', '.join(missing))
    if {p.name for p in source.iterdir()} != set(SOURCE_FILES):
        raise NativeApolloError('unexpected native member; fixed author collection required')
    checked = {}
    for name, expected in SOURCE_FILES.items():
        path = source / name
        if not stat.S_ISREG(path.lstat().st_mode):
            raise NativeApolloError('native source must be a regular file, not a symlink: ' + name)
        raw = path.read_bytes()
        if (len(raw) != expected['bytes'] or hashlib.sha256(raw).hexdigest() != expected['sha256']
                or hashlib.md5(raw).hexdigest() != expected['md5']):
            raise NativeApolloError('changed native Apollo source: ' + name)
        checked[name] = raw
    return checked


def reproduce(source_dir: Path, output_dir: Path) -> list[dict]:
    """Validate all originals, then retain native bytes and make 13 token tables.

    Existing outputs are never overwritten. A write failure leaves incomplete
    evidence for diagnosis and must not be published as a complete release.
    No source file, mtime, database pointer or model switch is modified.
    """
    output = Path(output_dir)
    if output.exists() or output.is_symlink():
        raise FileExistsError('new output directory required: ' + str(output))
    originals = _authenticate(source_dir)
    prepared = []
    recipe_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for name in sorted(_ROLES):
        raw = originals[name]
        csv_bytes = lexical_csv(raw)
        if restore_native(csv_bytes) != raw:
            raise NativeApolloError('lexical projection is not byte-reversible: ' + name)
        rows = list(csv.reader(io.StringIO(csv_bytes.decode('utf-8'), newline='')))
        record = {
            'source_filename': name, 'source_sha256': SOURCE_FILES[name]['sha256'],
            'source_md5': SOURCE_FILES[name]['md5'], 'source_bytes': len(raw),
            'source_uri': SOURCE_FILES[name]['source_uri'],
            'output': 'tables/' + name + '.native_tokens.csv',
            'output_sha256': hashlib.sha256(csv_bytes).hexdigest(),
            'row_count': len(rows) - 1, 'columns': list(CSV_COLUMNS),
            'schema_id': SCHEMA, 'source_role': _ROLES[name],
            'native_header_tokens': raw.decode('utf-8').splitlines()[0].split('\t'),
            'units': {'source_line_1based': '1', 'source_field_1based': '1',
                      'line_role': 'lexical label', 'native_token': 'unconverted native string; heterogeneous units',
                      'native_line_ending': 'lexical label'},
            'native_units': {'temperature': 'K for native TPD temperature axes',
                             'desorption': 'ML/s for desorption series, per author readme',
                             'energy_coverage_distribution': None},
            'unit_gap': 'Per-figure energy/coverage/distribution units and native-model roles must be reviewed before numerical interpretation',
            'scientific_status': 'SOURCE_LEXICAL_ACCESS_NOT_PHYSICAL_QUALIFICATION',
            'physical_activation': False, 'recipe_sha256': recipe_sha,
        }
        prepared.append((record, csv_bytes))
    # Authenticate every source and check every reversible output before any
    # output directory is created. Readme is preserved even without a table.
    output.mkdir(parents=True, exist_ok=False)
    (output / 'native').mkdir()
    (output / 'tables').mkdir()
    for name, raw in originals.items():
        with (output / 'native' / name).open('xb') as stream:
            stream.write(raw)
    for record, payload in prepared:
        with (output / record['output']).open('xb') as stream:
            stream.write(payload)
    records = [record for record, _ in prepared]
    receipt = {'schema': 'apollo-tpd-source-reproduction-v1', 'source_record': SOURCE_RECORD,
               'publication_doi': PUBLICATION_DOI, 'native_files': SOURCE_FILES,
               'table_count': len(records), 'tables': records, 'physical_activation': False,
               'scope': 'authenticated source and lexical access only; no physical evaluator'}
    with (output / 'reproduction.json').open('x', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write('\n')
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        records = reproduce(args.source_dir, args.output_dir)
    except NativeApolloUnavailable as error:
        print(json.dumps({'status': 'SKIPPED', 'reason': str(error)}))
        return 77
    except (OSError, ValueError, UnicodeError) as error:
        print(json.dumps({'status': 'ERROR', 'reason': str(error)}))
        return 2
    print(json.dumps({'status': 'PASS', 'scope': 'source-lexical-access-only',
                      'table_count': len(records), 'rows': sum(r['row_count'] for r in records),
                      'physical_qualification': 'NOT EVALUATED'}, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
