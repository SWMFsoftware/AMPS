#!/usr/bin/env python3
"""Reproduce four frozen IAEA lexical CSV products from native HTML arrays.

Preimplementation contract, 2026-11-18
------------------------------------
This portable source-specific recipe reproduces ONLY the retained D119990
Scatter and D107581 Patch/Line arrays. It adapts the existing independent
source processors in sodium_partial_capture/prepare_source.py and the retained
atomic packet's extract_native_inventory.py without their unrelated ADAS or
Phelps operations. Raw HTML and each resulting CSV have frozen SHA-256 values.
Native Bokeh renderer identity, source order, axes, dtype, byte order and shape
are mandatory. Every native series is preserved in its original order.

``reproduce(source_html, output_dir, dataset_id)`` returns output provenance
records and refuses an existing output directory. No network, PIC, external
package, source-data write, curve digitization, interpolation, refitting,
renormalization, sorting, or missing-value repair occurs. The native Patch is
plotting geometry, not an extra measured branch. D107581 is an explicitly
fitted capture-plus-ionization composite, not a pure-CEX measurement. D119990
is a measured H(1s) product partial, not an all-state total or event law.

The independent oracle is exact agreement with all four already retained CSV
hashes, plus native source identity/shape invariants and source byte/mtime
preservation. Passing this recipe establishes source reproduction only; no
physical readiness, complete differential law or application validation follows.
"""
from __future__ import annotations

import argparse
import ast
import base64
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import struct


_HEADER = ['energy_eV_per_u', 'cross_section_cm2']
_SOURCES = {
    119990: dict(
        sha256='ee669388646f4dfc54720e9d53ded92c3a853a041e5b37d2ea9c39bfed61f590',
        url='https://db-amdis.org/collisiondb/datasets/119990/',
        roles=[('Scatter', 12)], selected=0,
        process='H+ + Na -> H(1s) + Na+(2p6); measured ground-product partial capture',
        initial_Na_state='not supplied by the native dataset label',
        interpretation='measured partial, not total capture or daughter kinetics',
        outputs=[('collisiondb_119990_embedded_arrays.csv', 0,
            '2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93')]),
    107581: dict(
        sha256='24cc0de43be8002fa3e29fa4ca61fb867fe0c2250e4c05c92d746d018ea8b34b',
        url='https://db-amdis.org/collisiondb/datasets/107581/',
        roles=[('Patch', 1000), ('Line', 500)], selected=1,
        process='HCX charge transfer + HIN ionization; explicitly composite fitted source',
        initial_Na_state='Na 3s, as stated by the native label',
        interpretation='native fit/plot representation, not independent measured pure-CEX data',
        outputs=[('collisiondb_107581_all_native_series_0.csv', 0,
            'd4f7779e19bf0ebe2f1e9c8032d8d64df0b58b4e49a5c816fafe1abde3e86022'),
            ('collisiondb_107581_all_native_series_1.csv', 1,
            '97a8bbb273a7410150e100f6289c33232146903674358033988f8ee79a125358'),
            ('collisiondb_107581_embedded_arrays.csv', 1,
            '97a8bbb273a7410150e100f6289c33232146903674358033988f8ee79a125358')]),
}


def _walk(value):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from _walk(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk(child)


def _json(text):
    def unique(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate native Bokeh JSON field: ' + key)
            result[key] = value
        return result

    def invalid(token):
        raise ValueError('nonfinite native Bokeh JSON token: ' + token)

    return json.loads(text, object_pairs_hook=unique, parse_constant=invalid)


def _decode(raw: bytes, dataset_id: int, specification: dict):
    """Decode native IEEE float64 arrays; never select a series by its values."""
    if hashlib.sha256(raw).hexdigest() != specification['sha256']:
        raise ValueError('native D%d HTML hash changed; a new response requires source review' % dataset_id)
    text = raw.decode('utf-8')
    matches = re.findall(r"const docs_json = ('(?:\\.|[^'])*');", text)
    if len(matches) != 1 or 'DataSet D%d' % dataset_id not in text:
        raise ValueError('missing or ambiguous native dataset/Bokeh identity')
    objects = list(_walk(_json(ast.literal_eval(matches[0]))))
    axes = {obj['attributes']['axis_label'] for obj in objects
            if obj.get('name') in ('LinearAxis', 'LogAxis')}
    if axes != {'E /eV.u-1', 'sigma /cm2'}:
        raise ValueError('native energy/cross-section axis units changed')
    frames = re.findall(r'<tr>\s*<td>Frame</td>\s*<td>([^<]+)</td>\s*</tr>', text)
    if frames != ['Target']:
        raise ValueError('native target-rest frame missing, changed or ambiguous')
    renderers = [obj for obj in objects if obj.get('name') == 'GlyphRenderer']
    sources = [obj for obj in objects if obj.get('name') == 'ColumnDataSource']
    if len(renderers) != len(specification['roles']) or len(sources) != len(renderers):
        raise ValueError('native renderer/source count changed')
    role_by_source = {}
    for renderer in renderers:
        attrs = renderer['attributes']
        source_id = attrs['data_source']['id']
        if source_id in role_by_source:
            raise ValueError('ambiguous native renderer for one data source')
        role_by_source[source_id] = attrs['glyph']['name']
    result = []
    seen_sources = set()
    for number, source in enumerate(sources):
        expected_role, expected_rows = specification['roles'][number]
        if source['id'] in seen_sources or role_by_source.get(source['id']) != expected_role:
            raise ValueError('native series order/role changed; no arbitrary renderer selection')
        seen_sources.add(source['id'])
        data = source['attributes']['data']
        if data.get('type') != 'map':
            raise ValueError('unsupported native source data map')
        columns = {}
        for key, array in data['entries']:
            if key not in ('x', 'y') or key in columns:
                raise ValueError('duplicate or unsupported native numeric axis')
            if (array.get('type') != 'ndarray' or array.get('dtype') != 'float64' or
                    array.get('shape') != [expected_rows] or array.get('order') != 'little' or
                    array.get('array', {}).get('type') != 'bytes'):
                raise ValueError('native float64 dtype/order/shape changed')
            payload = base64.b64decode(array['array']['data'], validate=True)
            if len(payload) != 8 * expected_rows:
                raise ValueError('native numeric byte count mismatch')
            columns[key] = struct.unpack('<%dd' % expected_rows, payload)
        if set(columns) != {'x', 'y'}:
            raise ValueError('missing native numeric axis')
        rows = list(zip(columns['x'], columns['y']))
        if any(not math.isfinite(value) for row in rows for value in row):
            raise ValueError('nonfinite native number; no fill or repair permitted')
        # A Patch deliberately doubles back in x. Do not sort it or apply a
        # monotone-cross-section interpretation to a closed plotting polygon.
        result.append(dict(role=expected_role, rows=rows))
    return result


def reproduce(source_html: Path, output_dir: Path, dataset_id: int) -> list[dict]:
    """Reproduce frozen CSV bytes into a NEW directory and return provenance.

    Validation and byte-hash checks complete before any output is created.
    A later disk/write failure leaves its partial output for diagnosis, never
    reuses it as success, and never changes the original HTML. The caller
    supplies a writable evidence/staging destination, not an immutable raw
    archive. There is no implicit network or source replacement fallback.
    """
    if isinstance(dataset_id, bool) or not isinstance(dataset_id, int) or dataset_id not in _SOURCES:
        raise ValueError('only frozen IAEA dataset IDs 107581 and 119990 are supported')
    source, output = Path(source_html).resolve(), Path(output_dir).resolve()
    if output.exists():
        raise FileExistsError('output directory must not exist: ' + str(output))
    specification = _SOURCES[dataset_id]
    raw = source.read_bytes()
    series = _decode(raw, dataset_id, specification)
    prepared = []
    script_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for filename, number, expected_hash in specification['outputs']:
        buffer = io.StringIO(newline='')
        writer = csv.writer(buffer, lineterminator='\r\n')
        writer.writerow(_HEADER)
        writer.writerows(series[number]['rows'])
        payload = buffer.getvalue().encode('utf-8')
        output_hash = hashlib.sha256(payload).hexdigest()
        if output_hash != expected_hash:
            raise ValueError('native CSV reproduction differs from frozen product: ' + filename)
        record = dict(dataset_id=dataset_id, source_url=specification['url'],
            native_file=str(source), native_sha256=specification['sha256'],
            output=filename, output_sha256=output_hash,
            row_count=len(series[number]['rows']), rows=len(series[number]['rows']),
            native_series_index=number, native_renderer_role=series[number]['role'],
            native_series_count=len(series),
            selected_native_series_index=specification['selected'],
            is_selected_series_alias=filename.endswith('_embedded_arrays.csv'),
            native_units=dict(energy='eV/u', cross_section='cm^2'),
            native_axis_labels=['E /eV.u-1', 'sigma /cm2'],
            native_frame='Target', frame_interpretation='target-rest projectile energy',
            native_dtype='float64', native_byte_order='little',
            process=specification['process'], native_initial_Na_state=specification['initial_Na_state'],
            source_interpretation=specification['interpretation'],
            series_interpretation=('plotting patch, not a measured reaction branch'
                if series[number]['role'] == 'Patch' else specification['interpretation']),
            extraction='native typed arrays; no digitization, interpolation, fit, sorting or renormalization',
            status='SOURCE_REPRODUCTION_ONLY_NOT_PHYSICS_QUALIFIED',
            script_sha256=script_sha)
        prepared.append((record, payload))
    output.mkdir(parents=True, exist_ok=False)
    for record, payload in prepared:
        with (output / record['output']).open('xb') as stream:
            stream.write(payload)
    records = [record for record, _ in prepared]
    (output / 'reproduction.json').write_text(json.dumps(records, indent=2, sort_keys=True) + '\n',
                                             encoding='utf-8')
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-html', type=Path, required=True)
    parser.add_argument('--dataset-id', type=int, required=True, choices=sorted(_SOURCES))
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        records = reproduce(args.source_html, args.output_dir, args.dataset_id)
    except FileNotFoundError as error:
        print(json.dumps(dict(status='SKIPPED', reason=str(error)))); return 77
    except (OSError, ValueError, KeyError, TypeError, SyntaxError, struct.error) as error:
        print(json.dumps(dict(status='ERROR', reason=str(error)))); return 2
    print(json.dumps(dict(status='PASS', scope='source-reproduction-only',
                         scientific_validation='NOT EVALUATED', outputs=records), indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
