#!/usr/bin/env python3
"""Acquire and audit native NIST ASD O I radiative-source records.

This prepares research evidence, never a physical cascade activation file.
In particular an unrestricted database query enumerates stored records; it
does not prove that every possible physical decay has an assessed A value.
No missing A value, lifetime, branch, uncertainty or isotope is synthesized.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import hashlib
import html
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = 'https://physics.nist.gov'
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LINE_FIELDS = [
    'obs_wl_vac(nm)', 'ritz_wl_vac(nm)', 'Aki(s^-1)', 'Acc',
    'Ei(cm-1)', 'Ek(cm-1)', 'conf_i', 'term_i', 'J_i', 'conf_k',
    'term_k', 'J_k', 'g_i', 'g_k', 'Type', 'tp_ref', 'line_ref', '',
]
LEVEL_FIELDS = ['Configuration', 'Term', 'J', 'g', 'Prefix',
                'Level (cm-1)', 'Suffix', 'Uncertainty (cm-1)', 'Reference']
TARGETS = [('O_I_1D_2', '2s2.2p4', '1D', '2'),
           ('O_I_1S_0', '2s2.2p4', '1S', '0')]

# These are source-interface choices, not observational cuts. In particular:
# output=0 is the ENTIRE result, line_out=0 preserves rows without A values,
# both transition types are requested, and no wavelength, strength, intensity
# or accuracy bound is supplied. page_size is irrelevant to entire output;
# it is nevertheless recorded exactly to make the HTTP request reproducible.
LINE_QUERY = dict(spectra='O I', limits_type='0', de='0', format='2',
    output='0', line_out='0', remove_js='on', page_size='1000',
    show_obs_wl='1', show_calc_wl='1', unc_out='1', order_out='0',
    en_unit='0', show_av='3', A_out='0', allowed_out='1', forbid_out='1',
    conf_out='on', term_out='on', enrg_out='on', J_out='on',
    g_out='on', bibrefs='1', unit='1')
LEVEL_QUERY = dict(spectrum='O I', units='0', format='2', output='0',
    page_size='1000', multiplet_ordered='0', conf_out='on', term_out='on',
    level_out='on', unc_out='1', j_out='on', g_out='on',
    lande_out='on', perc_out='on', biblio='on')


def utc_now():
    return datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z')


def digest(path, algorithm='sha256'):
    return hashlib.new(algorithm, Path(path).read_bytes()).hexdigest()


def strict_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate JSON field: ' + key)
            result[key] = value
        return result
    def nonfinite(value):
        raise ValueError('nonfinite JSON token: ' + value)
    value = json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                       parse_constant=nonfinite)
    if not isinstance(value, dict):
        raise ValueError('source manifest must be a JSON object')
    return value


class NativeText(HTMLParser):
    """Read source prose without executing embedded scripts/styles."""
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.parts = []
        self.suppressed = 0
    def handle_starttag(self, tag, attrs):
        if tag in ('script', 'style'):
            self.suppressed += 1
    def handle_endtag(self, tag):
        if tag in ('script', 'style') and self.suppressed:
            self.suppressed -= 1
        if not self.suppressed:
            self.parts.append(' ')
    def handle_data(self, value):
        if not self.suppressed:
            self.parts.append(value)


def plain_html(text):
    parser = NativeText()
    parser.feed(text)
    return ' '.join(''.join(parser.parts).split())


def source_url(endpoint, parameters):
    return BASE + endpoint + '?' + urllib.parse.urlencode(parameters)


def fixed_requests():
    result = [
        ('lines_form.html', BASE + '/PhysRefData/ASD/lines_form.html'),
        ('levels_form.html', BASE + '/PhysRefData/ASD/levels_form.html'),
        ('lines_help.html', BASE + '/PhysRefData/ASD/Html/lineshelp.html'),
        ('levels_help.html', BASE + '/PhysRefData/ASD/Html/levelshelp.html'),
        ('contents.html', BASE + '/PhysRefData/ASD/Html/contents.html'),
    ]
    for label, endpoint, query in (
        ('oi_lines_complete', '/cgi-bin/ASD/lines1.pl', LINE_QUERY),
        ('oi_levels_complete', '/cgi-bin/ASD/energy1.pl', LEVEL_QUERY),
    ):
        for extension, form in (('csv', '2'), ('html', '0')):
            parameters = dict(query, format=form)
            result.append((label + '.' + extension, source_url(endpoint, parameters)))
    for reference, database, number, kind in (
        ('T4539', 'tp', '4539', 'T'), ('T5081', 'tp', '5081', 'T'),
        ('L7288', 'el', '7288', 'E'),
    ):
        parameters = dict(db=database, db_id=number, comment_code='',
                          element='O', spectr_charge='1', type=kind)
        result.append(('reference_' + reference + '.html',
            source_url('/cgi-bin/ASBib1/get_ASBib_ref.cgi', parameters)))
    return result


def fetch(url, target, deadline):
    """Retain original bytes, headers, exact request and retrieval interval.

    The byte ceiling is a resource refusal, not a scientific truncation: an
    oversized response is retained up to the limit and the package errors.
    No partially read response is eligible for parsing or a source QA PASS.
    """
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError('whole source-package deadline exhausted before request')
    started = utc_now()
    request = urllib.request.Request(url, headers={
        'User-Agent': 'AMPS-srcMoon-P05-source-audit/1.0',
        'Accept': '*/*',
    })
    descriptor = dict(url=url, started_utc=started, method='GET',
                      request_headers=dict(request.header_items()))
    with urllib.request.urlopen(request, timeout=min(30, remaining)) as response:
        descriptor.update(status=response.status, final_url=response.geturl(),
                          response_headers=list(response.headers.items()))
        if response.status != 200:
            raise ValueError('official source returned non-200 status')
        data = response.read(64 * 1024 * 1024 + 1)
    target.write_bytes(data)
    descriptor.update(completed_utc=utc_now(), file=target.name,
                      bytes=len(data), sha256=digest(target), md5=digest(target, 'md5'))
    if len(data) > 64 * 1024 * 1024:
        raise ValueError('source exceeded 64MiB safety limit; response is not complete')
    if not data:
        raise ValueError('official source response is empty')
    return descriptor


def unwrap(value):
    """Remove only ASD's CSV spreadsheet-text wrapper ="...".

    This is not formula execution or spreadsheet processing. Internal state
    labels, qualifiers, blank values and native significant figures remain.
    """
    if value.startswith('="') and value.endswith('"'):
        return value[2:-1]
    return value


def native_csv(path, fields):
    with Path(path).open(newline='') as stream:
        reader = csv.reader(stream, strict=True)
        header = next(reader)
        if header != fields:
            raise ValueError('native CSV schema changed: ' + str(path))
        rows = []
        for index, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError('native CSV column count differs on line ' + str(index))
            if header[-1] == '' and row[-1] != '':
                raise ValueError('unexpected data in native trailing blank column')
            rows.append(dict((key, unwrap(value)) for key, value in zip(header, row) if key))
    return rows


def number(value):
    """Parse unqualified native numbers; blanks/qualifiers remain unavailable."""
    if not re.fullmatch(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', value):
        return None
    try:
        result = Decimal(value)
    except InvalidOperation:
        return None
    if not result.is_finite():
        return None
    return result


def metadata(path, kind):
    text = plain_html(Path(path).read_text())
    if 'NIST Atomic Spectra Database ' + kind + ' Data' not in text:
        raise ValueError('native HTML does not identify requested ASD data type')
    pattern = r'O\s+I\s*:?' + (r'\s*([0-9]+)\s+Lines of Data Found' if kind == 'Lines'
                              else r'\s*([0-9]+)\s+Levels Found')
    count = re.search(pattern, text)
    version = re.search(r'\(ver\.\s*([^\)]+)\)', text)
    if not count or not version:
        raise ValueError('native species/count/version metadata is unavailable')
    return dict(row_count=int(count[1]), version=version[1], native_prose=text)


def selected_levels(levels):
    selected = {}
    for identity, configuration, term, angular in TARGETS:
        matches = [row for row in levels if
            (row['Configuration'], row['Term'], row['J']) == (configuration, term, angular)]
        if len(matches) != 1 or number(matches[0]['Level (cm-1)']) is None:
            raise ValueError('selected upper-state identity is missing/ambiguous: ' + identity)
        if matches[0]['Prefix'] or matches[0]['Suffix']:
            raise ValueError('selected level has an unresolved native energy qualifier')
        selected[identity] = matches[0]
    return selected


def crosscheck_requests(levels):
    result = []
    for identity, row in selected_levels(levels).items():
        for extension, form in (('csv', '2'), ('html', '0')):
            # The upper-energy boundary is exactly the requested state's
            # native energy. This independent source query contains ALL
            # lower-state exits of that state, at every wavelength. It is
            # solely a database-response consistency check; the full export
            # remains the primary untruncated source and is always retained.
            parameters = dict(LINE_QUERY, format=form,
                              max_upp_enrg=row['Level (cm-1)'])
            result.append(('crosscheck_' + identity + '.' + extension,
                source_url('/cgi-bin/ASD/lines1.pl', parameters)))
    return result


def acquire(output, deadline):
    raw = output / 'raw'
    raw.mkdir()
    receipt = dict(schema='srcMoon-NIST-OI-acquisition-v1',
                   started_utc=utc_now(), requests=[], source_role='research-only')
    try:
        for name, url in fixed_requests():
            receipt['requests'].append(fetch(url, raw / name, deadline))
        levels = native_csv(raw / 'oi_levels_complete.csv', LEVEL_FIELDS)
        for name, url in crosscheck_requests(levels):
            receipt['requests'].append(fetch(url, raw / name, deadline))
        receipt['status'] = 'COMPLETE'
        receipt['completed_utc'] = utc_now()
    except Exception as error:
        receipt.update(status='ERROR', reason=str(error), completed_utc=utc_now())
        raise
    finally:
        (output / 'acquisition.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return output


def verify_source_packet(packet):
    receipt = strict_json(packet / 'acquisition.json')
    if receipt.get('schema') != 'srcMoon-NIST-OI-acquisition-v1' or receipt.get('status') != 'COMPLETE':
        raise ValueError('source packet has no complete acquisition receipt')
    seen = set()
    for descriptor in receipt['requests']:
        name = descriptor['file']
        if not isinstance(name, str) or name != Path(name).name or name in seen:
            raise ValueError('unsafe or duplicate raw source identity')
        seen.add(name)
        path = packet / 'raw' / name
        if (not path.is_file() or path.stat().st_size != descriptor['bytes'] or
                digest(path) != descriptor['sha256'] or digest(path, 'md5') != descriptor['md5']):
            raise ValueError('raw source hash/size mismatch: ' + name)
        if descriptor['status'] != 200:
            raise ValueError('source receipt has a non-200 request')
    levels = native_csv(packet / 'raw/oi_levels_complete.csv', LEVEL_FIELDS)
    expected = dict(fixed_requests() + crosscheck_requests(levels))
    actual = {d['file']: d['url'] for d in receipt['requests']}
    if actual != expected:
        raise ValueError('source request set differs from the documented unrestricted queries')
    return receipt


def accuracy_codes(lines_help):
    text = Path(lines_help).read_text()
    if 'id="OUTACC"' not in text:
        raise ValueError('NIST accuracy-code metadata anchor is unavailable')
    section = text.split('id="OUTACC"', 1)[1].split('</table>', 1)[0]
    result = {}
    for native in re.findall(r'<tr>(.*?)</tr>', section, re.S | re.I):
        cells = re.findall(r'<td[^>]*>(.*?)</td>', native, re.S | re.I)
        if len(cells) != 3:
            continue
        code, relation, percent = [plain_html(cell) for cell in cells]
        value = number(percent.rstrip('.').rstrip('%'))
        if relation not in ('≤', '>') or value is None:
            raise ValueError('NIST accuracy-code definition changed')
        result[code] = dict(relation=relation, percent=str(value),
                            fractional_bound=str(value / 100))
    if not result:
        raise ValueError('NIST accuracy-code mapping is empty')
    return result


def write_csv(path, rows, fields):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def select_exits(lines, upper):
    return [row for row in lines if (row['conf_k'], row['term_k'], row['J_k']) ==
            (upper['Configuration'], upper['Term'], upper['J'])]


def stored_cascade_graph(levels, lines):
    """Follow native stored exits, keeping every fine-structure occurrence.

    This is source-graph bookkeeping, not a complete physical rate matrix.
    IDs encode original level labels/energies without assigning an isotope,
    source population, extra transition or unlisted terminal-state lifetime.
    """
    def identity(level):
        return '|'.join(level[name] for name in ('Configuration', 'Term', 'J', 'Level (cm-1)'))
    pending = list(selected_levels(levels).values())
    nodes, edges = {}, []
    while pending:
        upper = pending.pop(0)
        upper_id = identity(upper)
        if upper_id in nodes:
            continue
        exits = select_exits(lines, upper)
        nodes[upper_id] = dict(source_level=upper,
            stored_exit_components=len(exits), physical_lifetime_s=None)
        for row in exits:
            matches = [lower for lower in levels if
                (lower['Configuration'], lower['Term'], lower['J'], lower['Level (cm-1)']) ==
                (row['conf_i'], row['term_i'], row['J_i'], row['Ei(cm-1)'])]
            if len(matches) != 1 or row['Ek(cm-1)'] != upper['Level (cm-1)']:
                raise ValueError('stored descendant component has an ambiguous native level identity')
            lower = matches[0]
            a = number(row['Aki(s^-1)'])
            if row['Aki(s^-1)'] and (a is None or a < 0):
                raise ValueError('stored descendant A has unsupported notation/sign')
            if (number(lower['Level (cm-1)']) is None or
                    number(lower['Level (cm-1)']) >= number(upper['Level (cm-1)'])):
                raise ValueError('stored spontaneous component is not a downward native transition')
            edges.append(dict(upper_native_level_id=upper_id,
                              lower_native_level_id=identity(lower), **row))
            pending.append(lower)
    return dict(schema='srcMoon-NIST-OI-stored-cascade-graph-v1',
        scope='all ASD-stored components reachable from selected upper states; not all physical exits',
        nodes=nodes, components=edges, physical_activation=False,
        complete_physical_cascade=False, physical_lifetime_s=None)


def audit(packet, output):
    receipt = verify_source_packet(packet)
    raw = packet / 'raw'
    levels = native_csv(raw / 'oi_levels_complete.csv', LEVEL_FIELDS)
    lines = native_csv(raw / 'oi_lines_complete.csv', LINE_FIELDS)
    level_meta = metadata(raw / 'oi_levels_complete.html', 'Levels')
    line_meta = metadata(raw / 'oi_lines_complete.html', 'Lines')
    if len(levels) != level_meta['row_count'] or len(lines) != line_meta['row_count']:
        raise ValueError('complete native HTML/CSV row counts differ')
    if line_meta['version'] != level_meta['version']:
        raise ValueError('line and level exports use different ASD versions')
    if 'vacuum' not in line_meta['native_prose']:
        raise ValueError('native line metadata does not confirm vacuum wavelengths')
    codes = accuracy_codes(raw / 'lines_help.html')
    processed = output / 'processed'
    processed.mkdir()
    write_csv(processed / 'oi_levels.csv', levels, LEVEL_FIELDS)
    write_csv(processed / 'oi_transitions.csv', lines, LINE_FIELDS[:-1])
    selected, lower_pairs, summaries = [], [], []
    for identity, upper in selected_levels(levels).items():
        exits = select_exits(lines, upper)
        if not exits:
            raise ValueError('no native exit records for selected state: ' + identity)
        keys = [(r['conf_i'], r['term_i'], r['J_i'], r['Type']) for r in exits]
        if len(set(keys)) != len(keys):
            raise ValueError('duplicate selected lower-state/multipole component')
        cross = native_csv(raw / ('crosscheck_' + identity + '.csv'), LINE_FIELDS)
        cross_meta = metadata(raw / ('crosscheck_' + identity + '.html'), 'Lines')
        if len(cross) != cross_meta['row_count'] or cross_meta['version'] != line_meta['version']:
            raise ValueError('independent upper-energy query metadata/count differs')
        def canonical(rows):
            return Counter(tuple(row[name] for name in LINE_FIELDS[:-1]) for row in rows)
        if canonical(select_exits(cross, upper)) != canonical(exits):
            raise ValueError('selected native exits differ between unrestricted and upper-energy queries')
        available_sum = Decimal(0)
        absent = 0
        for source in exits:
            if source['Ek(cm-1)'] != upper['Level (cm-1)']:
                raise ValueError('upper-level energy crosswalk differs from native level export')
            lower = [r for r in levels if
                (r['Configuration'], r['Term'], r['J'], r['Level (cm-1)']) ==
                (source['conf_i'], source['term_i'], source['J_i'], source['Ei(cm-1)'])]
            if len(lower) != 1:
                raise ValueError('selected lower-state mapping is missing or ambiguous')
            a = number(source['Aki(s^-1)'])
            if source['Aki(s^-1)'] and (a is None or a < 0):
                raise ValueError('selected transition probability has unsupported notation/sign')
            if a is None:
                absent += 1
            else:
                available_sum += a
            accuracy = codes.get(source['Acc'])
            # Preserve an unknown/qualified native accuracy code rather than
            # guessing its bound. B+/C+ definitions are source metadata, not
            # Gaussian 1sigma uncertainties or independent random errors.
            selected.append(dict(upper_state_id=identity, **source,
                lower_level_uncertainty_cm_1=lower[0]['Uncertainty (cm-1)'],
                upper_level_uncertainty_cm_1=upper['Uncertainty (cm-1)'],
                accuracy_relation=accuracy['relation'] if accuracy else '',
                accuracy_fractional_bound=accuracy['fractional_bound'] if accuracy else ''))
        # List every catalogued lower level, including pairs with no stored
        # transition. Missing pairs are NOT zero rates or fabricated branches.
        for lower in levels:
            energy = number(lower['Level (cm-1)'])
            if energy is None or energy >= number(upper['Level (cm-1)']):
                continue
            matches = [r for r in exits if
                (r['conf_i'], r['term_i'], r['J_i']) ==
                (lower['Configuration'], lower['Term'], lower['J'])]
            lower_pairs.append(dict(upper_state_id=identity,
                lower_configuration=lower['Configuration'], lower_term=lower['Term'],
                lower_J=lower['J'], lower_energy_cm_1=lower['Level (cm-1)'],
                listed_component_count=len(matches),
                source_status='LISTED_COMPONENTS' if matches else 'NO_STORED_TRANSITION_NOT_A_ZERO_RATE'))
        summaries.append(dict(upper_state_id=identity, source_level=upper,
            listed_exit_component_count=len(exits), listed_A_missing_count=absent,
            sum_of_listed_available_A_s_1=str(available_sum),
            sum_is_complete_decay_rate=False,
            complete_radiative_exit_evidence='UNAVAILABLE',
            total_decay_rate_s_1=None, lifetime_s=None,
            normalized_branching_probabilities=None))
    fields = ['upper_state_id'] + LINE_FIELDS[:-1] + [
        'lower_level_uncertainty_cm_1', 'upper_level_uncertainty_cm_1',
        'accuracy_relation', 'accuracy_fractional_bound']
    write_csv(processed / 'selected_exit_components.csv', selected, fields)
    write_csv(processed / 'selected_lower_level_coverage.csv', lower_pairs, list(lower_pairs[0]))
    graph = stored_cascade_graph(levels, lines)
    (processed / 'stored_cascade_graph.json').write_text(json.dumps(graph, indent=2) + '\n')
    write_csv(processed / 'stored_cascade_components.csv', graph['components'],
              ['upper_native_level_id', 'lower_native_level_id'] + LINE_FIELDS[:-1])
    references = {name: plain_html((raw / ('reference_' + name + '.html')).read_text())
                  for name in ('T4539', 'T5081', 'L7288')}
    gaps = [
        'ASD all-lines output enumerates all stored rows, not a certificate of every physical radiative exit.',
        'The selected transition references assess M1/E2 within 2p4; no acquired bound establishes omitted higher-multipole, multiphoton or isotope/hyperfine-induced contributions.',
        'O I export has no isotope-resolved applicability metadata; do not silently assign these entries to a selected oxygen isotope.',
        'Selected level uncertainties are blank; native line HTML states that explicit observed/Ritz wavelength uncertainties are unavailable for this selection.',
        'A accuracy codes are predominantly systematic critical-assessment estimates; covariance and statistical confidence are not provided.',
        'Source product-state populations, polarization/alignment, collisional quenching and the corresponding omission bounds are not supplied by this atomic export.',
    ]
    qa = dict(schema='srcMoon-NIST-OI-source-qa-v1', status='PASS',
        scope='source identity/schema/count/level/component/query-consistency checks only',
        database_version=line_meta['version'], native_line_rows=len(lines),
        native_level_rows=len(levels), selected_exit_components=len(selected),
        stored_cascade_levels=len(graph['nodes']),
        stored_cascade_components=len(graph['components']),
        complete_export_missing_A_rows=sum(not r['Aki(s^-1)'] for r in lines),
        selected_missing_A_rows=sum(not r['Aki(s^-1)'] for r in selected),
        untruncated_primary_queries=True, upper_energy_crosschecks_identical=True,
        vacuum_wavelength_unit='nm', level_energy_unit='cm^-1', A_unit='s^-1',
        accuracy_code_definitions=codes, selected_states=summaries,
        radiative_cascade_qualification='SKIPPED / NOT QUALIFIED',
        physical_activation=False, canonical_data_modified=False,
        missing_evidence=gaps)
    (output / 'qa_report.json').write_text(json.dumps(qa, indent=2, ensure_ascii=False) + '\n')
    source = dict(schema='srcMoon-NIST-OI-source-provenance-v1',
        role='RESEARCH_SOURCE_NOT_READY_FOR_PHYSICAL_CASCADE',
        source_packet=str(packet), acquisition_sha256=digest(packet / 'acquisition.json'),
        source_requests=receipt['requests'], database_version=line_meta['version'],
        source_doi='10.18434/T4W30F', bibliographic_native_text=references,
        processing_script=dict(path=str(Path(__file__).resolve()), sha256=digest(__file__)),
        transformations=['decode ASD spreadsheet-text CSV wrappers without executing them',
            'preserve source strings, missing cells, significant figures, accuracy and reference codes',
            'select exact native upper configuration/term/J identities before any AMPS residual',
            'map upper/lower identities to native levels without energy interpolation',
            'add native accuracy-code bounds as metadata without assigning confidence/covariance',
            'sum listed numeric A components only; do not call the sum a complete rate or derive lifetime/branches'],
        output_files={str(path.relative_to(output)): dict(sha256=digest(path), md5=digest(path, 'md5'))
                      for path in sorted(processed.iterdir())})
    (output / 'provenance.json').write_text(json.dumps(source, indent=2) + '\n')
    return qa


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True,
                        help='new research output directory; existing evidence is never overwritten')
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--download', action='store_true', help='retrieve exact official NIST requests')
    mode.add_argument('--source-packet', type=Path,
                      help='offline replay from an existing hash-verified acquisition packet')
    parser.add_argument('--wall-seconds', type=float, default=3600,
                        help='whole-package budget including retrieval/processing/reporting, maximum3600')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    result = dict(schema='srcMoon-NIST-OI-source-preparation-result-v1', status='ERROR',
                  scope='research source acquisition/QA only',
                  formal_cascade_status='SKIPPED / NOT QUALIFIED')
    script_hash = digest(__file__)
    try:
        if not 0 < args.wall_seconds <= 3600:
            raise ValueError('whole-package wall budget must be in(0,3600]seconds')
        deadline = start + args.wall_seconds
        packet = acquire(output, deadline) if args.download else args.source_packet.resolve()
        qa = audit(packet, output)
        if digest(__file__) != script_hash:
            raise ValueError('source-processing script changed during execution')
        if time.monotonic() > deadline:
            raise TimeoutError('whole source-package wall budget exceeded')
        result.update(status='PASS', source_qa_status=qa['status'],
                      native_line_rows=qa['native_line_rows'],
                      native_level_rows=qa['native_level_rows'],
                      selected_exit_components=qa['selected_exit_components'])
    except (OSError, ValueError, KeyError, csv.Error, StopIteration) as error:
        result.update(status='ERROR', reason=str(error))
    finally:
        result['elapsed_seconds'] = time.monotonic() - start
        (output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
        revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                  capture_output=True, text=True, timeout=5)
        (output / 'run_manifest.json').write_text(json.dumps(dict(
            command=sys.argv, source_revision=revision.stdout.strip(),
            processing_script_sha256=script_hash, python=sys.version,
            started_elapsed_clock='time.monotonic', canonical_data_modified=False,
            application_build_performed=False), indent=2) + '\n')
    print(json.dumps(result, indent=2))
    return 0 if result['status'] == 'PASS' else 2


if __name__ == '__main__':
    raise SystemExit(main())
