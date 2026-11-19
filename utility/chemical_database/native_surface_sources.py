#!/usr/bin/env python3
"""Reproduce the 18 native TU Wien lunar sputter tables without changing values.

Preimplementation contract, 2026-11-19
------------------------------------
The sole supported source is the author-deposited, CC-BY4.0 figure-data archive
at record hjn8d-hyh84 (version index2, DOI10.48436/hjn8d-hyh84), linked to
Broetzner et al.2025, DOI10.1038/s43247-025-02546-0. Exact producer size/MD5
and acquired SHA256 authenticate it before any output is created. All 18
members, original names, CSV headers, number strings, order and newline bytes
are preserved. Four files are experimental (17 rows); fourteen are model
outputs. No source value, uncertainty, interpolation or fitted curve is made.

The independent oracle is the original ZIP member bytes and producer checksum,
not current AMPS output. These tables measure/report sputtered target material,
not incident-ion reflection, charge survival, implantation or residence. The
publication's 68501/68051 sample-identifier inconsistency remains unresolved.
Temperature and uncertainty confidence levels absent from the source remain
absent. Successful reproduction qualifies source integrity only, never a
boundary law, independent model validation or a complete ion-surface kernel.

``reproduce(source_zip, output_dir)`` is offline, standard-library only and
refuses an existing output directory. Missing source is unavailable; changed
source, schema, dimensions or invalid numerical data is an error. A later
write failure retains partial evidence; it never modifies the original ZIP.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import io
import json
from pathlib import Path
import zipfile


SOURCE_SHA256 = 'de22ab6a26ef363349e11626e8bf98c35df5244ec58f528ae4a1b3f2630a4f1a'
SOURCE_MD5 = '3d72470af26187925596faf0cfbad47b'
SOURCE_BYTES = 22673
SOURCE_URI = ('https://researchdata.tuwien.ac.at/api/records/hjn8d-hyh84/files/'
              'As_published_in_Comm_Earth_Env.zip/content')
_MODEL = ['Angle_(deg)', 'Yield_(amu/ion)', 'Yield_(atoms/ion)']
_EXPERIMENT = ['Angle_(deg)', 'Yield_(amu/ion)', 'Error_(amu/ion)',
               'Yield_(atoms/ion)', 'Error_(atoms/ion)']
# These are source member identities and native dimensions, not selected
# residuals or an invented smooth grid. In particular,299 and300 are distinct.
_ROWS = {
    'Fig1a_Experiment.txt': 6, 'Fig1a_SDTrimSP_.txt': 300,
    'Fig1a_SDTrimSP_default.txt': 300, 'Fig1a_SDTrimSP_Szabo.txt': 300,
    'Fig1a_SpuBase.txt': 89, 'Fig1a_SRIM.txt': 15,
    'Fig1b_Experiment.txt': 6, 'Fig1b_SDTrimSP_.txt': 299,
    'Fig1b_SDTrimSP_default.txt': 299, 'Fig1b_SDTrimSP_Szabo.txt': 299,
    'Fig1b_SpuBase.txt': 89, 'Fig1b_SRIM.txt': 15,
    'Fig2a_Experiment_Rough.txt': 2, 'Fig2a_SDTrimSP__Regolith.txt': 8,
    'Fig2a_SDTrimSP__Rough.txt': 6,
    'Fig2b_Experiment_Rough.txt': 3, 'Fig2b_SDTrimSP__Regolith.txt': 8,
    'Fig2b_SDTrimSP__Rough.txt': 6,
}


def reproduce(source_zip: Path, output_dir: Path) -> list[dict]:
    """Authenticate and copy every native CSV-formatted .txt member verbatim.

    Parsing checks native headers, dimensions, finite nonnegative values and
    ordered source angles only. Parsed numbers are NEVER serialized back into
    the products. Both amu/ion and atoms/ion are already author-supplied columns;
    their bulk-stoichiometry conversion is not recalculated by this recipe.
    """
    source, output = Path(source_zip).resolve(), Path(output_dir).resolve()
    if output.exists():
        raise FileExistsError('output directory already exists: ' + str(output))
    if source.stat().st_size != SOURCE_BYTES:
        raise ValueError('native archive byte length changed')
    raw = source.read_bytes()
    if (hashlib.sha256(raw).hexdigest() != SOURCE_SHA256 or
            hashlib.md5(raw).hexdigest() != SOURCE_MD5):
        raise ValueError('native archive digest changed; source review required')
    prepared = []
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)) or set(names) != set(_ROWS):
            raise ValueError('native archive member identity/count changed')
        if archive.testzip() is not None:
            raise ValueError('native ZIP member CRC failure')
        for name in names:
            payload = archive.read(name)
            table = list(csv.reader(io.StringIO(payload.decode('utf-8'))))
            experimental = 'Experiment' in name
            header = _EXPERIMENT if experimental else _MODEL
            if not table or table[0] != header or len(table)-1 != _ROWS[name]:
                raise ValueError('native table header/row count changed: ' + name)
            previous = None
            for values in table[1:]:
                if len(values) != len(header):
                    raise ValueError('ragged native table: ' + name)
                try:
                    numbers = [Decimal(value) for value in values]
                except InvalidOperation as error:
                    raise ValueError('invalid native numeric token: ' + name) from error
                if any(not value.is_finite() or value < 0 for value in numbers):
                    raise ValueError('nonfinite/negative native source value: ' + name)
                if numbers[0] > 90 or (previous is not None and numbers[0] <= previous):
                    raise ValueError('native incidence-angle order/domain changed: ' + name)
                previous = numbers[0]
            record = {
                'output': name, 'native_archive_member': name,
                'source_uri': SOURCE_URI, 'source_sha256': SOURCE_SHA256,
                'source_md5': SOURCE_MD5, 'source_bytes': SOURCE_BYTES,
                'output_sha256': hashlib.sha256(payload).hexdigest(),
                'row_count': _ROWS[name], 'columns': header,
                'source_role': 'measured-sputter-yield' if experimental else 'model-output-sputter-yield',
                'projectile_label': 'H' if name.startswith(('Fig1a_', 'Fig2a_')) else 'He',
                'sample_form': 'flat-film' if name.startswith('Fig1') else
                    'rough-pellet' if experimental or name.endswith('__Rough.txt') else 'modeled-porous-regolith',
                'units': {'Angle_(deg)': 'degree', 'Yield_(amu/ion)': 'amu/ion',
                          'Yield_(atoms/ion)': 'atoms/ion',
                          **({'Error_(amu/ion)': 'amu/ion',
                              'Error_(atoms/ion)': 'atoms/ion'} if experimental else {})},
                'native_atoms_per_ion_conversion': 'Author supplied; bulk-stoichiometry assumption in publication; not recalculated',
                'uncertainty_confidence_level': None,
                'temperature_K': None,
                'material_id_status': 'Publication abstract/results68501 versus methods68051; unresolved',
                'incoming_ion_retention_fraction': None,
                'outgoing_charge_distribution': None,
                'status': 'SOURCE_REPRODUCTION_ONLY_NOT_PHYSICS_QUALIFIED',
                'recipe_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            }
            prepared.append((record, payload))
    # The complete archive/schema is checked before creating an output. Native
    # member names are fixed above, so no archive path can escape this directory.
    output.mkdir(parents=True, exist_ok=False)
    for record, payload in prepared:
        with (output/record['output']).open('xb') as stream:
            stream.write(payload)
    records = [record for record, _ in prepared]
    (output/'reproduction.json').write_text(json.dumps(records, indent=2, sort_keys=True)+'\n')
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-zip', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        records = reproduce(args.source_zip, args.output_dir)
    except FileNotFoundError as error:
        print(json.dumps({'status': 'SKIPPED', 'reason': str(error)})); return 77
    except (OSError, ValueError, UnicodeError, zipfile.BadZipFile) as error:
        print(json.dumps({'status': 'ERROR', 'reason': str(error)})); return 2
    print(json.dumps({'status': 'PASS', 'scope': 'source-reproduction-only',
                      'physical_qualification': 'NOT EVALUATED', 'outputs': records}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
