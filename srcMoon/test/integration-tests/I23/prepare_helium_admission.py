#!/usr/bin/env python3
"""Derive the reviewed I23 helium admission evidence from pinned sources.

The 2026-11-25 approval record (reference/helium_admission_approval_20261125.json)
states the user's decisions. This committed script turns them into the two
files that the production ``HeliumElectronImpact::Compose`` authenticates by
SHA-256: a threshold/state crosswalk and a model-error/observable record. It
introduces no physical number of its own. Every value is either read from an
authenticated source (D05, NIST state profile, the retained Kim 1992 fit
comparison, the retained CCC diagnostic) or is an approved tolerance, or is
an explicitly labelled evaluation-rounding allowance.

Units: SI (J, kg, m/s) with eV shown alongside. Frame: LSO. Time: SPICE ET.
Exit codes: 0 PASS, 2 ERROR, 77 SKIPPED (missing source).
"""
import argparse
from decimal import Decimal as D, localcontext
import csv
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
TEST = HERE.parents[1]
ROOT = TEST.parents[1]
sys.path.insert(0, str(TEST / 'U16_electron_impact'))
from voronov_verification import MissingPrerequisite, independent_rows, preflight, reference_log_slope

APPROVAL = HERE / 'reference/helium_admission_approval_20261125.json'
CASE = HERE / 'reference/helium_case.json'
EXIT = {'PASS': 0, 'ERROR': 2, 'SKIPPED': 77}
ENERGY_REFERENCE = 'P02:separated-neutral-ground-isotope-atoms-and-free-electron-reference-zero'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pinned(path, digest, label):
    path = Path(path)
    if not path.is_file():
        raise MissingPrerequisite(label + ' unavailable: ' + str(path))
    if sha(path) != digest:
        raise ValueError(label + ' SHA-256 mismatch: ' + str(path))
    return path


def binary64_spacing(value):
    """Unit in the last place of a normal binary64 value (53-bit significand)."""
    return long_double_spacing(value) * 2 ** 11


def long_double_spacing(value):
    """Unit in the last place of an x87 extended (64-bit significand) value."""
    value = abs(D(value))
    exponent = 0
    two = D(2)
    while value >= two:
        value /= two
        exponent += 1
    while value < 1:
        value *= two
        exponent -= 1
    return two ** (exponent - 63)


def derive(case, approval, profile_path, data_root):
    profile = json.loads(Path(profile_path).read_text())
    states = {}
    with Path(profile['states']['path']).open(newline='') as stream:
        for row in csv.DictReader(stream):
            states[row['state_id']] = row
    constants = {}
    with Path(profile['constants']['path']).open(newline='') as stream:
        for row in csv.DictReader(stream):
            constants[row['name']] = D(row['value'])
    for required in ('4He:ground', '4He+:ground', 'electron:free'):
        if required not in states or states[required]['status'] != 'available':
            raise MissingPrerequisite('state profile lacks available ' + required)
    ev = constants['elementary charge']
    me = D(states['electron:free']['physical_rest_mass_kg'])
    # Composer kinetic (reference) masses: neutral isotope minus z electrons.
    m_he = D(states['4He:ground']['physical_rest_mass_kg'])
    m_ion_kinetic = m_he - me
    package, hashes = preflight(Path(data_root))
    fit = independent_rows(package)['He']
    with localcontext() as context:
        context.prec = 70
        temperature = D(case['bath']['thermal_energy_eV'])
        slope = reference_log_slope(fit, temperature)
        fit_threshold_j = fit['energy'] * ev
        physical_threshold_j = D(states['4He+:ground']['formation_J'])
        difference = abs(fit_threshold_j - physical_threshold_j)
        # Rounding allowance only, not a physical tolerance. Production stores
        # the D05 decimal threshold (24.6 eV) as binary64, so one binary64
        # spacing of that value [eV] times e covers its representation error;
        # 64 extended-precision ulps of the larger threshold [J] cover the
        # composer's long-double evaluation of |Q + I_fit|.
        allowance = (binary64_spacing(fit['energy']) * ev
                     + 64 * long_double_spacing(max(fit_threshold_j, physical_threshold_j)))
        kim = json.loads(pinned(case['kim_comparison']['path'], case['kim_comparison']['sha256'],
                                'Kim 1992 comparison').read_text())
        slope_error = D(repr(kim['analytic_slope_equivalent_difference']))
        ccc = json.loads(pinned(case['ccc_diagnostic']['path'], case['ccc_diagnostic']['sha256'],
                                'CCC diagnostic').read_text())['curtin_he']
        excited = ccc['finite_table_integrals'][0]['single_ionization_fraction_above_HeII_n2_onset']
        ccc_mean = D(repr(ccc['finite_table_integrals'][0]['mean_incident_energy_eV']))

        def recoil(log_slope):
            e_in = temperature * ev * (D('1.5') + log_slope)
            e_out = e_in - physical_threshold_j
            impulse = (2 * me * e_in).sqrt() + (4 * me * e_out).sqrt()
            dv = impulse / m_ion_kinetic
            return dict(mean_incident_energy_J=e_in, mean_incident_energy_eV=e_in / ev,
                        outgoing_aggregate_energy_J=e_out, rms_impulse_bound_kg_m_s=impulse,
                        rms_speed_bound_m_s=dv, kinetic_bound_stationary_J=impulse * impulse / (2 * m_ion_kinetic),
                        kinetic_bound_stationary_eV=impulse * impulse / (2 * m_ion_kinetic) / ev)
        central, upper = recoil(slope), recoil(slope + slope_error)
        lower_mean = temperature * ev * (D('1.5') + slope - slope_error)
    rec = approval['admission']['recoil']
    v_limit = D(rec['maximum_heavy_velocity_error_m_s'])
    e_limit_j = D(rec['maximum_heavy_kinetic_error_eV']) * ev
    if upper['rms_speed_bound_m_s'] > v_limit or upper['kinetic_bound_stationary_J'] > e_limit_j:
        raise ValueError('approved recoil tolerances are below the composer upper-end bound; no relaxation')
    if lower_mean < physical_threshold_j:
        raise ValueError('admitted slope error reaches below the physical helium threshold')
    sources = {str(Path(profile_path).resolve()): sha(profile_path), **hashes,
               str(APPROVAL): sha(APPROVAL),
               case['kim_comparison']['path']: case['kim_comparison']['sha256'],
               case['ccc_diagnostic']['path']: case['ccc_diagnostic']['sha256']}
    common = dict(approval={'path': str(APPROVAL.relative_to(ROOT)), 'sha256': sha(APPROVAL),
                            'date': approval['approval_date']},
                  scope=approval['admission']['observable_scope'], frame='LSO',
                  time_scale=case['time_scale'], energy_reference=ENERGY_REFERENCE,
                  generator=str(Path(__file__).resolve().relative_to(ROOT)),
                  formal_I23_status='SKIPPED / NOT QUALIFIED')
    s = lambda value: format(value, 'f') if isinstance(value, D) else value
    crosswalk = dict(common, schema='srcMoon-I23-helium-state-crosswalk-v1',
        reaction='e(bath) + 4He(1s2 1S0) -> 4He+(1s 2S1/2) + 2 e(bath)',
        rate_fit={'source': 'D05 Voronov cfit.dat', 'row_line': fit['line'],
                  'threshold_eV': s(fit['energy']), 'threshold_J': s(fit_threshold_j)},
        physical={'state': '4He+:ground', 'source_id': states['4He+:ground']['formation_source_id'],
                  'source_record': states['4He+:ground']['formation_source_record'],
                  'formation_J': s(physical_threshold_j), 'formation_eV': s(physical_threshold_j / ev)},
        product_state_assignment=approval['admission']['threshold_state']['product_state_assignment'],
        excited_ion_yield={'status': 'unknown', 'ccc_fraction_above_HeII_n2_onset': excited,
                           'ccc_onset_eV': ccc['spectroscopic_excited_product_onset_eV'],
                           'interpretation': 'energy eligibility in an unpublished tabulated model, not a yield; not used in any computation'},
        source_threshold_difference_J=s(difference), source_threshold_difference_eV=s(difference / ev),
        evaluation_rounding_allowance_J=s(allowance),
        admitted_threshold_crosswalk_maximum_difference_J=s(difference + allowance),
        sources=sources)
    model_error = dict(common, schema='srcMoon-I23-helium-model-error-v1',
        observable=approval['admission']['observable_scope'],
        conditional_mean={'rule': approval['admission']['conditional_mean']['rule'],
                          'voronov_log_slope_at_10eV': s(slope),
                          'kim_minus_voronov_slope_equivalent': s(slope_error),
                          'ccc_mean_incident_energy_eV_recorded_only': s(ccc_mean),
                          'admitted_maximum_absolute_log_slope_error': s(slope_error),
                          'lower_end_mean_incident_energy_eV': s(lower_mean / ev),
                          'interpretation': approval['admission']['conditional_mean']['interpretation']},
        recoil={'rule': rec['rule'], 'central': {k: s(v) for k, v in central.items()},
                'upper_end': {k: s(v) for k, v in upper.items()},
                'admitted_maximum_heavy_velocity_error_m_s': s(v_limit),
                'admitted_maximum_heavy_kinetic_error_J': s(e_limit_j)},
        mechanical=approval['admission']['mechanical_limits_He'],
        sources=sources)
    return crosswalk, model_error


def run(profile, data_root, output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    result = {'status': 'ERROR'}
    try:
        case = json.loads(CASE.read_text())
        approval = json.loads(APPROVAL.read_text())
        if approval.get('approval_status') != 'APPROVED':
            raise MissingPrerequisite('I23 helium admission is not approved')
        pinned(profile, case['profile_sha256'], 'helium state profile')
        crosswalk, model_error = derive(case, approval, profile, data_root)
        for name, payload in (('helium_state_crosswalk.json', crosswalk), ('helium_model_error.json', model_error)):
            (output / name).write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
        result = {'status': 'PASS', 'files': {name: sha(output / name) for name in
                  ('helium_state_crosswalk.json', 'helium_model_error.json')}}
    except MissingPrerequisite as error:
        result = {'status': 'SKIPPED', 'reason': str(error)}
    except (OSError, ValueError, KeyError, TypeError, ArithmeticError) as error:
        result = {'status': 'ERROR', 'reason': str(error)}
    (output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', type=Path, required=True)
    parser.add_argument('--data-root', type=Path, default=Path('/data/vtenishe/moon_validation_data'))
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    value = run(args.profile, args.data_root, args.output_dir)
    print(json.dumps(value, indent=2))
    sys.exit(EXIT[value['status']])
