"""Read-only, rights-aware inventory of existing AMPS chemistry sources.

The inventory is the source-discovery layer, not a downloader/importer or an
approval registry. It never creates a missing physical datum, copies a source
file, chooses between conflicting sources, normalizes numbers or enables a
reaction. The builder must enforce each file's permissions before copying.
See RIGHTS.md for the source evidence and preimplementation contract.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import mimetypes
import os
from pathlib import Path
from typing import Any

ATOMIC = 'test_output/srcMoon/m3-atomic-online-databases-20261118T083400Z'
MOLECULAR = 'test_output/srcMoon/m3-molecular-databases-20261118T083551Z'
SURFACE = 'test_output/srcMoon/m3-ion-surface-cex-data-sources-20261117T225600Z'
MASSES = 'test_output/srcMoon/m3-approved-products-20261118T163142Z/source-acquisition-v2'
OXYGEN = 'test_output/srcMoon/m3-oxygen-cascade-sources-20261118T170845Z'
CEX = 'test_output/srcMoon/m3-charge-exchange-research-20261117'
KINEMATIC = 'test_output/srcMoon/m3-kinematic-references-20261117'


def _policy(storage: str, evaluation: str, redistribution: str,
            *evidence: dict[str, str], conditions: str = '') -> dict[str, Any]:
    """Keep operation permissions independent of scientific availability."""
    return dict(storage_permission=storage, evaluation_permission=evaluation,
                redistribution=redistribution, evidence=list(evidence),
                conditions=conditions)


IAEA = _policy('allowed', 'allowed', 'allowed-with-attribution',
    dict(url='https://db-amdis.org/collisiondb/licence/',
         note='CC BY4.0 for data/metadata unless explicitly stated otherwise.'),
    conditions='Retain original author/dataset citations, licence and modifications. Permission does not qualify states, domains or event kinetics.')
EXOMOL = _policy('allowed', 'allowed', 'allowed-with-attribution-sharealike',
    dict(url='https://www.exomol.com/data/licence/', note='All ExoMol data: CC BY-SA4.0.'),
    conditions='Preserve credit, licence link and changes; respect ShareAlike when sharing adaptations. No physical readiness follows.')
NIST_UNMARKED = _policy('allowed', 'allowed', 'allowed-with-source-credit',
    dict(url='https://www.nist.gov/copyrights-disclaimers',
         note='Use of NIST information permits copying of information not marked copyrighted.'),
    conditions='Only the inspected unmarked constants/isotope responses and their source pages. Not a public-domain assertion for every NIST product; exclude ASD-derived thresholds.')
NIST_ASD = _policy('unknown', 'unknown', 'not-cleared',
    dict(url='https://www.nist.gov/pml/atomic-spectra-database',
         note='ASD product homepage explicitly states copyright/all rights reserved.'),
    dict(url='https://www.nist.gov/open/license',
         note='SRD products are distinguished from the non-SRD reuse grant.'),
    conditions='Unmarked individual exports do not erase the ASD product notice. Keep metadata-only until managed-operation permission is established.')
NIST_OTHER = _policy('unknown', 'unknown', 'not-cleared',
    dict(url='https://physics.nist.gov/PhysRefData/Ionization/',
         note='Separate SRD107 product; its managed-operation permission was not established by this audit.'),
    dict(url='https://www.nist.gov/copyrights-disclaimers',
         note='Product-specific notices control; constants/isotope permission is not transferred.'),
    conditions='No conclusion about this product follows solely from another NIST product being public or copyrighted.')
ADAS = _policy('pending', 'pending', 'not-cleared',
    dict(url='https://open.adas.ac.uk/terms-and-conditions',
         note='Native ADF terms require written permission for insertion into a managed database.'),
    dict(local_terms=ATOMIC + '/adas_he_cex_2.dat', note='Native terms lines86-96; no permission obtained.'),
    conditions='No numerical payload in a managed database, including a private cache, pending required written permission. No provider contacted.')
UNKNOWN = _policy('unknown', 'unknown', 'not-cleared',
    conditions='No applicable managed-storage grant established; retain source metadata/links without copying numerical payload.')
PHIDRATES = _policy('unknown', 'unknown', 'not-cleared',
    dict(url='https://phidrates.space.swri.edu/', note='Scientific service; copyright notice; no managed-storage grant located.'))
VORONOV = _policy('unknown', 'allowed', 'not-cleared',
    dict(url='https://www.pa.uky.edu/~verner/info.html', note='Explicit research/code use and acknowledgement request; no managed-database licence located.'),
    dict(url='https://www.pa.uky.edu/~verner/col.html', note='Published native Voronov table and Fortran program.'))
LEIDEN = _policy('unknown', 'unknown', 'not-cleared',
    dict(url='https://home.strw.leidenuniv.nl/~ewine/photo/', note='Research numerical downloads/citation instructions; managed grant unresolved.'))
MAINZ = _policy('unknown', 'unknown', 'not-cleared',
    dict(url='https://www.uv-vis-spectral-atlas-mainz.org/uvvis/contact.html', note='Atlas imprint/citations; no source-specific managed-storage grant established.'))
UMIST = _policy('unknown', 'unknown', 'not-cleared',
    dict(url='https://umistdatabase.uk/downloads', note='RATE22 downloads; no explicit dataset grant located. New Mills notice is webpage software, not the data licence.'))
KIDA = _policy('unknown', 'allowed', 'not-cleared',
    dict(url='https://kida.astrochem-tools.org/help', note='Research use; complete reaction-set download and storing reaction lists are described.'),
    dict(url='https://kida.astrochem-tools.org/legal', note='Website IP reserved; no cleared shared managed-data licence established.'))
LXCAT = _policy('unknown', 'unknown', 'not-cleared',
    dict(url='https://nl.lxcat.net/instructions/redistribution.php', note='Contributor/third-party distribution restrictions; local managed-operation grant unresolved.'),
    dict(local_terms=ATOMIC + '/lxcat_export_plot.html', note='Captured terms; challenged numerical export was not bypassed.'))
TELLEZ = _policy('unknown', 'unknown', 'not-cleared',
    dict(url='https://pure.tue.nl/ws/portalfiles/portal/3503988/583864115491353.pdf', note='Publisher-PDF cover permits private research copying, restricts further distribution.'),
    dict(local_terms=SURFACE + '/tellez2012_silica_neutralization.txt', note='Cover terms line41; paper/table/release rights are separate.'))


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _files(root: Path) -> list[Path]:
    """List only native packet files, never virtual environments/bytecode.

    Existing research roots are explicit, not the whole test_output tree.
    Walk without following directory symlinks so accidental caches or a link
    cannot turn a source inventory into an unbounded unrelated tree scan.
    """
    if not root.is_dir():
        return []
    result = []
    for directory, subdirs, names in os.walk(root, followlinks=False):
        subdirs[:] = sorted(d for d in subdirs if d not in
            {'.venv', '__pycache__', '.git', '.pytest_cache', 'node_modules'})
        for name in sorted(names):
            path = Path(directory) / name
            if path.is_file() and path.suffix not in {'.pyc', '.pyo'}:
                result.append(path)
    return result


def _json(path: Path) -> Any:
    def unique(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate JSON field: ' + key)
            result[key] = value
        return result
    def invalid(token):
        raise ValueError('nonfinite JSON token: ' + token)
    return json.loads(path.read_text(), object_pairs_hook=unique, parse_constant=invalid)


def _metadata(root: Path, paths: list[Path]) -> tuple[dict, list[str]]:
    """Recover original URIs/expected hashes without treating them as review.

    Source acquisition manifests use several existing schemas. This lookup
    reads their ordinary path/url/hash descriptors and never executes an
    acquisition script. Ambiguous basenames are not silently resolved.
    """
    result: dict[str, dict] = {}
    issues = []
    by_name: dict[str, list[Path]] = {}
    for path in paths:
        by_name.setdefault(path.name, []).append(path)
    def walk(value):
        if isinstance(value, dict):
            yield value
            # Complete packet manifests also use a filename-keyed map.
            # Preserve the filename as its own descriptor before recursively
            # visiting ordinary nested dictionaries.
            if isinstance(value.get('files'), dict):
                for logical, descriptor in value['files'].items():
                    if isinstance(logical, str):
                        row = descriptor if isinstance(descriptor, dict) else {'sha256': descriptor}
                        yield dict(row, path=logical)
            for child in value.values():
                yield from walk(child)
        elif isinstance(value, list):
            for child in value:
                yield from walk(child)
    candidates = [p for p in paths if p.suffix == '.json' and
                  any(word in p.name for word in ('manifest', 'inventory', 'provenance', 'retrieval', 'retention', 'acquisition'))]
    def authority(path):
        # This ordering is declared independently of observed file values.
        # Initial/superseded extracts remain historical evidence, rather than
        # being allowed to change the expectation of a current snapshot.
        if path.name.startswith(('initial_', 'initial-series')):
            return (100, str(path))
        rank = {'packet_manifest.json': 0, 'retention_audit.json': 0,
                'raw_inventory.json': 0, 'acquisition.json': 0,
                'provenance.json': 1, 'native_inventory.json': 2}
        return (rank.get(path.name, 10), str(path))
    candidates.sort(key=authority)
    for path in candidates:
        try:
            value = _json(path)
        except (OSError, ValueError) as error:
            issues.append('UNREADABLE_SOURCE_METADATA: ' + str(path) + ': ' + str(error))
            continue
        for record in walk(value):
            uri = next((record.get(k) for k in ('source_url', 'url', 'final_url')
                        if isinstance(record.get(k), str) and record[k].startswith(('https://', 'http://'))), None)
            for key in ('path', 'file', 'filename', 'native_file', 'output'):
                name = record.get(key)
                if not isinstance(name, str):
                    continue
                target = Path(name)
                choices = [target] if target.is_absolute() and target in paths else []
                if not target.is_absolute() and root / target in paths:
                    choices = [root / target]
                if not choices and not target.is_absolute():
                    # Acquisition paths are relative to that receipt's own
                    # directory; NIST's raw responses are its raw/children.
                    choices = [p for p in (path.parent / target, path.parent / 'raw' / target)
                               if p in paths]
                if not choices and len(by_name.get(target.name, [])) == 1:
                    choices = by_name[target.name]
                for choice in choices:
                    entry = result.setdefault(str(choice), {})
                    if uri:
                        entry.setdefault('source_uri', uri)
                    # A derived CSV and its native HTML often occur in the
                    # SAME descriptor. Bind each path to its own hash; using
                    # native_sha256 for output would falsely report corruption.
                    if key == 'output':
                        expected = record.get('output_sha256') or record.get('sha256')
                    elif key == 'native_file':
                        expected = record.get('native_sha256') or record.get('sha256')
                    else:
                        expected = record.get('sha256')
                    if isinstance(expected, str) and len(expected) == 64:
                        if entry.get('expected_sha256') not in (None, expected):
                            issues.append('HISTORICAL_METADATA_HASH_DIFFERS_FROM_AUTHORITATIVE_SNAPSHOT: '
                                          + str(choice) + ' in ' + str(path))
                        entry.setdefault('expected_sha256', expected)
                        entry.setdefault('expected_sha256_from', str(path))
    return result, issues


def _role(path: Path) -> str:
    name = path.name.lower()
    if any(word in name for word in ('license', 'licence', 'terms')):
        return 'license'
    if (path.suffix in {'.json', '.py', '.log'} or
        any(word in name for word in ('manifest', 'provenance', 'qa_report', 'sha256sums'))):
        return 'provenance'
    if path.suffix.lower() == '.pdf' or name.endswith('.md'):
        return 'reference'
    return 'processed' if 'processed' in path.parts or path.suffix.lower() == '.csv' else 'raw'


UNITS = {
    'energy_eV_per_u': 'eV/u', 'cross_section_cm2': 'cm^2',
    'wavelength_nm': 'nm', 'photon_energy_eV': 'eV',
    'cross_section_m2': 'm^2', 'threshold_eV': 'eV',
    'rate_1au_s-1': 's^-1', 'excess_energy_eV': 'eV',
    'fit_A_cm3_s': 'cm^3/s', 'fit_A_m3_s': 'm^3/s',
    'temperature_min_eV': 'eV', 'temperature_max_keV': 'keV',
}


def _table(path: Path, file_id: str) -> dict | None:
    """Describe lexical source tables, preserving native bytes and headers.

    The IAEA CSV has its exact maintained extractor header. ExoMol photo
    files are headerless two-column native ASCII; column labels/units are
    the documented ExoMol source schema, not inserted text. Other CSVs are
    listed only when a genuine textual header exists. Headerless material
    spectra with unresolved axes receive no invented column names/table law.
    """
    if path.suffix.lower() == '.photo':
        rows = 0
        with path.open() as stream:
            for line in stream:
                if not line.strip():
                    continue
                if len(line.split()) != 2:
                    return None
                rows += 1
        return dict(id=file_id + ':table', name=path.name, format='native-ascii',
                    schema_id='exomol-photo-two-column-native-v1',
                    columns=[], row_count=rows,
                    units={'wavelength': 'nm', 'cross_section': 'cm^2/molecule'},
                    header_in_source=False,
                    column_labels_from='ExoMol native photo definition/format metadata')
    if path.name == 'codata2022.txt':
        return dict(id=file_id + ':table', name=path.name, format='native-ascii',
                    schema_id='nist-codata-2022-native-fixed-width-v1', columns=[],
                    row_count='unknown', units={},
                    interpretation='Native per-record quantities/units; no new fixed-width parser or derived SI table supplied here.')
    if path.suffix.lower() != '.csv':
        return None
    with path.open(newline='') as stream:
        reader = csv.reader(stream, strict=True)
        try:
            header = next(reader)
        except StopIteration:
            return None
        if not header:
            return None
        # A numeric first record is data, not a header. Never consume that
        # record merely because its filename has a .csv extension.
        try:
            for cell in header:
                float(cell)
        except ValueError:
            pass
        else:
            return None
        count = sum(1 for _record in reader)
    schema = 'source-native-headered-csv-v1'
    if header == ['energy_eV_per_u', 'cross_section_cm2']:
        schema = 'iaea-collisiondb-native-energy-cross-section-v1'
    return dict(id=file_id + ':table', name=path.name, format='csv',
                schema_id=schema, columns=header, row_count=count,
                units={name: UNITS[name] for name in header if name in UNITS},
                header_in_source=True)


def inventory(amps_root: Path, validation_root: Path) -> list[dict]:
    """Return all known datasets; missing packets remain explicit records.

    Paths are caller-owned and explicit. Dataset/file IDs depend only on
    stable source IDs and native logical names, never machine absolute paths,
    traversal order, source values or a newly observed fit/model residual.
    """
    amps, validation = Path(amps_root).resolve(), Path(validation_root).resolve()
    datasets: dict[str, dict] = {}
    file_ids = set()

    def dataset(identifier, title, provider, version, description, policy,
                expected_root, issues=(), source_uri=''):
        evidence = copy.deepcopy(policy)
        for record in evidence['evidence']:
            if 'local_terms' in record:
                record['local_terms'] = str(amps / record['local_terms'])
        value = dict(id=identifier, title=title, provider=provider, version=version,
            description=description, rights=evidence, capabilities=[],
            issues=list(issues), files=[], expected_root=str(expected_root),
            source_uri=source_uri,
            availability='present' if expected_root.is_dir() else 'missing')
        if not expected_root.is_dir():
            value['issues'].append('MISSING_PACKET: ' + str(expected_root))
        datasets[identifier] = value
        return value

    def add(ds, path, logical_name, meta=None, policy=None):
        path = path.resolve()
        permission = policy or ds['rights']
        logical = logical_name.replace(os.sep, '/')
        identifier = ds['id'] + ':' + hashlib.sha256(logical.encode()).hexdigest()
        if identifier in file_ids:
            raise ValueError('duplicate stable file identity: ' + ds['id'] + '/' + logical)
        file_ids.add(identifier)
        record = dict(id=identifier, path=str(path), role=_role(path),
            logical_name=logical, media_type=mimetypes.guess_type(path.name)[0] or 'application/octet-stream',
            storage_permission=permission['storage_permission'],
            evaluation_permission=permission['evaluation_permission'],
            source_uri=(meta or {}).get('source_uri', ds['source_uri']),
            capabilities=['source-access'], observed_sha256=_sha(path),
            byte_count=path.stat().st_size)
        if meta and 'expected_sha256' in meta:
            record.update(expected_sha256=meta['expected_sha256'],
                          expected_sha256_from=meta['expected_sha256_from'])
            if record['observed_sha256'] != record['expected_sha256']:
                ds['issues'].append('SOURCE_HASH_MISMATCH: ' + logical)
        try:
            table = _table(path, identifier)
        except (OSError, UnicodeError, ValueError, csv.Error) as error:
            ds['issues'].append('UNINTERPRETED_TABLE: ' + logical + ': ' + str(error))
            table = None
        if table:
            record['table'] = table
            record['capabilities'].append('table-access')
        ds['files'].append(record)

    canonical = validation / 'photochemistry'
    photo = dataset('phidrates-d04', 'PHIDRATES frozen D04 source package', 'PHIDRATES/SwRI',
        'provider version not supplied; canonical2026-10-07 snapshot',
        'Four atomic, seven H2O and four OH_exp reactive channels; two absorption totals. Quiet/active1AU endpoints remain distinct.',
        PHIDRATES, canonical, ['RATE_SOURCE_READY_IS_NOT_DAUGHTER_KINETIC_READY',
        'HE_AR_LONG_WAVELENGTH_THRESHOLD_STATE_INTERPRETATION_UNRESOLVED'], 'https://phidrates.space.swri.edu/')
    electron = dataset('voronov-d05', 'Voronov atomic electron-impact fits', 'Voronov/UKY',
        'Version2, March24,1997', '406 native atomic/ionic fits; separate four-neutral selection. Maxwellian coefficients, not beam cross sections or daughter distributions.',
        VORONOV, canonical, ['NO_MOLECULAR_ELECTRON_FITS', 'MANAGED_STORAGE_PERMISSION_UNRESOLVED'],
        'https://www.pa.uky.edu/~verner/col.html')
    common = dataset('canonical-chemistry-mixed', 'Shared canonical package metadata and combined index',
        'mixed PHIDRATES/Voronov', 'canonical2026-10-07 snapshot',
        'Mixed numerical summaries/index/provenance inherit unresolved source rights; no blanket own-code exemption.', UNKNOWN, canonical)
    paths = _files(canonical)
    metadata, problems = _metadata(canonical, paths)
    common['issues'].extend(problems)
    for path in paths:
        relative = path.relative_to(canonical).as_posix()
        if ('phidrates' in path.parts or path.name in ('fetch_phidrates.py', 'prepare_phidrates.py')):
            owner = photo
        elif ('electron_impact' in path.parts or path.name in ('fetch_electron_impact.py', 'prepare_electron_impact.py')):
            owner = electron
        else:
            owner = common
        add(owner, path, relative, metadata.get(str(path)))

    atomic_root = amps / ATOMIC
    atomic_groups = {
        'iaea': dataset('iaea-collisiondb-research', 'IAEA CollisionDB retained partial/composite sources', 'IAEA CollisionDB',
            'retained2026-11-18 responses', 'Native data and exact embedded arrays; explicit partial/composite/plotting roles retained.', IAEA, atomic_root,
            ['D119990_INITIAL_NA_STATE_BIBLIOGRAPHY_UNCERTAINTY_NOT_SUPPLIED',
             'D107581_COMPOSITE_CAPTURE_PLUS_IMPACT_IONIZATION_NOT_PURE_CEX',
             'D107581_DISPLAY_GRID_IS_NOT_AUTHOR_VALIDITY_DOMAIN',
             'PATCH_ARRAYS_ARE_PLOTTING_GEOMETRY_NOT_A_BRANCH'], 'https://db-amdis.org/collisiondb/'),
        'adas': dataset('open-adas-research', 'OPEN-ADAS ADF01/ADF07 research references', 'OPEN-ADAS',
            'native SCCS/source versions retained per file', 'Atomic Maxwellian and state-conditioned integral CEX alternatives; no joint angle law.', ADAS, atomic_root,
            ['WRITTEN_MANAGED_DATABASE_PERMISSION_REQUIRED', 'NO_NUMERICAL_CACHE_BYPASS'], 'https://open.adas.ac.uk/'),
        'nist': dataset('nist-electron-research', 'NIST electron total/differential research responses', 'NIST SRD107',
            'source snapshot2026-11-18', 'He energy-differential marginals and H2O BEB totals are not full product laws.', NIST_OTHER, atomic_root,
            ['PRODUCT_SPECIFIC_MANAGED_PERMISSION_NOT_ESTABLISHED', 'SDCS_IS_NOT_A_JOINT_ANGLE_ENERGY_DISTRIBUTION'], 'https://physics.nist.gov/PhysRefData/Ionization/'),
        'lxcat': dataset('lxcat-research', 'LXCat selection/terms and unavailable export evidence', 'LXCat',
            'source snapshot2026-11-18', 'Process identities and metadata; challenged export is not a downloaded numerical table.', LXCAT, atomic_root,
            ['NUMERICAL_EXPORT_CHALLENGED_NOT_BYPASSED'], 'https://nl.lxcat.net/'),
        'phelps': dataset('phelps-published-research', 'Original Phelps compilation and extracted Na tables', 'Phelps/JILA',
            'native publication retained', 'Na3s and3p beam total cross sections; not a Maxwellian/product sampler. Paper rights are separate from LXCat.', UNKNOWN, atomic_root,
            ['PAPER_AND_MANAGED_TABLE_RIGHTS_UNRESOLVED', 'NO_ANGLE_OR_UNMEASURED_TAIL_LAW'], 'https://nl.lxcat.net/notes/index.php?download=phelps5'),
        'mixed': dataset('atomic-research-mixed', 'Mixed atomic research manifests/processors/reports', 'mixed research providers',
            '2026-11-18', 'Mixed files may embed restricted ADAS numbers; no provider-wide numerical clearance is inherited.', UNKNOWN, atomic_root),
    }
    paths = _files(atomic_root)
    metadata, problems = _metadata(atomic_root, paths)
    atomic_groups['mixed']['issues'].extend(problems)
    for path in paths:
        name = path.name.lower()
        kind = ('adas' if name.startswith(('adas', 'openadas', 'prior_gayet')) else
                'iaea' if name.startswith('collisiondb_') or name == 'iaea_collisiondb.html' else
                'nist' if name.startswith('nist_') else
                'lxcat' if name.startswith('lxcat_') else
                'phelps' if name.startswith('phelps_') else 'mixed')
        add(atomic_groups[kind], path, path.relative_to(atomic_root).as_posix(), metadata.get(str(path)))

    molecular_root = amps / MOLECULAR
    molecular_groups = {}
    for key, title, provider, policy, uri, issues in (
        ('leiden', 'Leiden photon tables and aggregate partial arrays', 'Leiden photo', LEIDEN,
         'https://home.strw.leidenuniv.nl/~ewine/photo/', ['NEGATIVE_HDF5_SAMPLES_UNRESOLVED', 'PARTIALS_DO_NOT_EXACTLY_CLOSE', '2016_PARTIAL_AND2017_TEXT_NOT_SPLICED', 'NORMALIZED_SOLAR_IS_NOT1AU']),
        ('mainz', 'MPI-Mainz H2O absorption research data', 'MPI-Mainz UV/VIS Atlas', MAINZ,
         'https://www.uv-vis-spectral-atlas-mainz.org/uvvis/', ['ABSORPTION_IS_NOT_CHANNEL_YIELD_OR_SOLAR_RATE', 'NATIVE_AIR_VACUUM_AND_TEMPERATURE_METADATA_CONTROL']),
        ('exomol', 'ExoMol OH MYTHOS/H2O DTU photo tables', 'ExoMol', EXOMOL,
         'https://www.exomol.com/data/data-types/photo/', ['OH_ACTUAL19173_VS_METADATA19172_ROWS_UNRESOLVED', 'TOTAL_ABSORPTION_IS_NOT_DAUGHTER_KINETICS', 'ISOTOPOLOGUE_AND_SOURCE_TEMPERATURE_SCOPE_RETAINED']),
        ('umist', 'UMIST RATE22 retained network', 'UMIST', UMIST,
         'https://umistdatabase.uk/downloads', ['ISRF_NETWORK_RATES_ARE_NOT_SOLAR1AU', 'REDUCED_NETWORK_PRODUCTS_ARE_NOT_RESOLVED_BRANCHES']),
        ('kida', 'KIDA2024 retained network/help', 'KIDA', KIDA,
         'https://kida.astrochem-tools.org/networks', ['DRAINE_ISRF_IS_NOT_SOLAR1AU', 'RESEARCH_USE_STATEMENT_IS_NOT_CLEARED_SHARED_MANAGED_LICENCE']),
        ('mixed', 'Mixed molecular research manifests/processors/reports', 'mixed molecular providers', UNKNOWN,
         '', ['MIXED_SOURCE_NUMERICAL_SUMMARIES_INHERIT_UNRESOLVED_RIGHTS']),
    ):
        molecular_groups[key] = dataset('molecular-' + key + '-research', title, provider,
            'native source versions retained; OH MYTHOS20241206 where applicable', title,
            policy, molecular_root, issues, uri)
    paths = _files(molecular_root)
    metadata, problems = _metadata(molecular_root, paths)
    molecular_groups['mixed']['issues'].extend(problems)
    for path in paths:
        name = path.name.lower()
        kind = next((key for key in ('leiden', 'mainz', 'exomol', 'umist', 'kida') if name.startswith(key + '_')), 'mixed')
        add(molecular_groups[kind], path, path.relative_to(molecular_root).as_posix(), metadata.get(str(path)))

    mass_root = amps / MASSES
    constants = dataset('nist-codata-2022-unmarked', 'NIST CODATA2022 native constants', 'NIST/CODATA',
        '2022 adjustment; native ASCII retrieved2026-11-18', 'Unmarked public-information source response; source units/uncertainties retained.',
        NIST_UNMARKED, mass_root, source_uri='https://physics.nist.gov/cuu/Constants/Table/allascii.txt')
    isotopes = dataset('nist-isotopes-unmarked', 'NIST isotope relative masses/compositions', 'NIST',
        'native compilation/AME2012 source retained', 'Explicit isotope source, not a lunar abundance population or mixture-mass substitute.',
        NIST_UNMARKED, mass_root, source_uri='https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all')
    thresholds = dataset('nist-asd-ground-research', 'NIST ASD H/Na ground ionization responses', 'NIST ASD',
        'source snapshot2026-11-18', 'Source thresholds lack isotope-shift/uncertainty crosswalk; product copyright retained.',
        NIST_ASD, mass_root, source_uri='https://physics.nist.gov/PhysRefData/ASD/')
    mixed_mass = dataset('nist-reference-state-mixed', 'Derived mixed reference-state inventory and provenance', 'mixed NIST source classes',
        'srcMoon-isotope-reference-mass-evidence-v1', 'Contains both public-information constants/isotopes and ASD thresholds; cannot inherit blanket permission.', UNKNOWN, mass_root,
        ['MIXED_ASD_THRESHOLDS_KEEP_DERIVED_INVENTORY_ON_HOLD'])
    paths = _files(mass_root)
    metadata, problems = _metadata(mass_root, paths)
    mixed_mass['issues'].extend(problems)
    for path in paths:
        owner = constants if path.name == 'codata2022.txt' else isotopes if path.name == 'isotopes.html' else thresholds if path.name in ('hydrogen-ground.html', 'sodium-ground.html') else mixed_mass
        add(owner, path, path.relative_to(mass_root).as_posix(), metadata.get(str(path)))

    oxygen_root = amps / OXYGEN
    oxygen = dataset('nist-asd-oi-5.12-research', 'NIST ASD5.12 O I levels/transitions', 'NIST ASD',
        '5.12', '614 native levels/910 transitions; stored O1D/O1S descendants are not a complete physical lifetime certificate.',
        NIST_ASD, oxygen_root, ['MISSING_A_VALUES_PRESERVED', 'COMPLETE_PHYSICAL_RADIATIVE_EXIT_BOUND_UNAVAILABLE',
        'NO_ISOTOPE_POPULATION_ALIGNMENT_OR_QUENCHING_DEFAULT'], 'https://physics.nist.gov/PhysRefData/ASD/')
    paths = []
    for component in ('raw', 'acquired', 'final-replay'):
        paths.extend(_files(oxygen_root / component))
    paths.extend(p for p in (oxygen_root / 'REPORT.md', oxygen_root / 'completion_report.json') if p.is_file())
    metadata, problems = _metadata(oxygen_root, paths)
    oxygen['issues'].extend(problems)
    for path in sorted(paths):
        add(oxygen, path, path.relative_to(oxygen_root).as_posix(), metadata.get(str(path)))

    surface_root = amps / SURFACE
    surface_groups = {
        'tellez': dataset('tellez-silica-2012', 'Tellez2012 He+/silica conditional charge survival', 'Tellez etal/AIP/TUe',
            'APL101151602 (2012)', 'Published scalar-coordinate fits; not reflection, implantation or lunar retention.', TELLEZ, surface_root,
            ['FITTED_VELOCITY_PROJECTION_NOT_ESTABLISHED', 'OXYGEN800_TO2000EV_NO_FIT_GAP', 'NO_ION_RETENTION_LAW'], 'https://doi.org/10.1063/1.4758699'),
        'adas': dataset('open-adas-surface-packet-reference', 'OPEN-ADAS integral CEX in material research packet', 'OPEN-ADAS',
            'native ADF01 retained', 'High-energy state-integral reference, distinct from low-energy KVI source.', ADAS, surface_root,
            ['WRITTEN_MANAGED_DATABASE_PERMISSION_REQUIRED'], 'https://open.adas.ac.uk/'),
        'enstatite': dataset('zenodo-enstatite-15978325', 'Author-deposited enstatite material research data', 'Cymes etal/Zenodo',
            'DOI10.5281/zenodo.15978325', 'Three acquired source files; simulated depth/composition and experimental spectra are separate observables.', UNKNOWN, surface_root,
            ['MOST_DEPOSIT_FILES_NOT_ACQUIRED', 'SOURCE_AXES_AND_MEASUREMENT_SIMULATION_ROLES_NOT_FULLY_QUALIFIED'], 'https://zenodo.org/records/15978325'),
        'apollo': dataset('zenodo-apollo17-eels-10794539', 'Apollo17 EELS author-deposit metadata', 'Cymes etal/Zenodo',
            'DOI10.5281/zenodo.10794539', 'Metadata acquired;682.5MB numerical archive not acquired.', UNKNOWN, surface_root,
            ['NUMERICAL_ARCHIVE_NOT_ACQUIRED'], 'https://zenodo.org/records/10794539'),
        'other': dataset('ion-material-other-research', 'Other retained ion-material/differential references', 'mixed authors/repositories',
            'native versions retained', 'CaSiO3 request-only spectra, Rice finite-angle benchmarks and material leads; preserve all gaps.', UNKNOWN, surface_root,
            ['CASIO3_NUMERICAL_SPECTRA_REQUEST_ONLY_NOT_ACQUIRED', 'RICE0.002_VS0.02DEG_LABEL_GAP', 'NO_COMPLETE_LUNAR_ION_MATERIAL_KERNEL']),
    }
    # The native author-deposit record supplies an explicit CC BY4.0 grant.
    # Apply it only to that deposit's metadata/selected files, not neighbouring
    # publisher PDFs or unrelated material sources in the same directory.
    for key, filename in [('enstatite', 'cymes2026_enstatite_data_record.json'),
                          ('apollo', 'cymes2024_apollo17_eels_record.json')]:
        record = surface_root / filename
        if record.is_file():
            try:
                if _json(record).get('metadata', {}).get('license', {}).get('id') == 'cc-by-4.0':
                    surface_groups[key]['rights'] = _policy('allowed', 'allowed', 'allowed-with-attribution',
                        dict(local_terms=str(record), note='Native Zenodo metadata licence.id=cc-by-4.0.'),
                        dict(url=surface_groups[key]['source_uri'], note='Author deposit; preserve attribution/source roles.'))
            except (OSError, ValueError) as error:
                surface_groups[key]['issues'].append('UNREADABLE_NATIVE_LICENCE: ' + str(error))
    paths = _files(surface_root)
    metadata, problems = _metadata(surface_root, paths)
    surface_groups['other']['issues'].extend(problems)
    for path in paths:
        name = path.name.lower()
        kind = ('tellez' if name.startswith('tellez') else
                'adas' if name.startswith('openadas') else
                'enstatite' if name.startswith(('enstatite_selected_', 'cymes2026_enstatite_data_record.')) else
                'apollo' if name == 'cymes2024_apollo17_eels_record.json' else 'other')
        # The two numeric .csv spectra have no authoritative column labels in
        # the acquired files. _table refuses to consume their first data row.
        add(surface_groups[kind], path, path.relative_to(surface_root).as_posix(), metadata.get(str(path)))

    for relative, identifier, title, provider, uri, issues in (
        (CEX, 'published-total-cex-research', 'ORNL/Morgan source publications and acquisition evidence',
         'Barnett/ORNL and Morgan/NIST/AIP', 'https://www.osti.gov/servlets/purl/6570226',
         ['SOURCE_SPECIFIC_PAPER_MANAGED_RIGHTS_UNRESOLVED', 'TOTAL_RATES_NOT_DAUGHTER_KINETICS', 'FAILED_RETRIEVAL_BODIES_ARE_NOT_PUBLICATIONS']),
        (KINEMATIC, 'kinematic-publication-research', 'Molecular/atomic differential and state-publication references',
         'mixed source authors', '', ['PAPER_LICENCE_DOES_NOT_AUTOMATICALLY_CLEAR_LINKED_NUMERICAL_TABLES', 'FULL_SOURCE_JOINT_DISTRIBUTIONS_NOT_ACQUIRED']),
    ):
        base = amps / relative
        ds = dataset(identifier, title, provider, 'native publication versions retained', title, UNKNOWN, base, issues, uri)
        paths = _files(base)
        metadata, problems = _metadata(base, paths)
        ds['issues'].extend(problems)
        for path in paths:
            add(ds, path, path.relative_to(base).as_posix(), metadata.get(str(path)))

    for ds in datasets.values():
        ds['files'].sort(key=lambda row: row['logical_name'])
        ds['capabilities'] = sorted({cap for row in ds['files'] for cap in row['capabilities']})
        if ds['availability'] == 'present' and not ds['files']:
            ds['availability'] = 'missing'
            ds['issues'].append('MISSING_SOURCE_FILES: ' + ds['expected_root'])
        ds['issues'] = list(dict.fromkeys(ds['issues']))
    return [datasets[key] for key in sorted(datasets)]
