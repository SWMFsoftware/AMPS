#!/usr/bin/python3
"""Prepare the 2013-04-11 real-event inputs of the srcSEP nightly tests.

This is the only Python of the nightly case and it runs only to *set up the
event*.  The nightly Table recipe never calls it: it reads the checked-in
runtime deck ``event.in``.  Developers regenerate the checked-in products
from the pinned raw archives with

    prepare_event.py derive --raw RAW_DIR      # raw archives -> derived.json,
                                               #   observations.csv
    prepare_event.py deck                      # derived.json -> event.in
    prepare_event.py check --raw RAW_DIR       # both, byte-compared

``derive`` needs NumPy and cdflib (CDF reading).  ``deck`` uses only the
standard library.  Every raw file is verified against ``raw_manifest.json``
before it is read.  No observation is interpolated, smoothed, extrapolated
or gap-filled; fill values and out-of-range values are dropped and counted.

Units: SI unless the name says otherwise.  Times are UTC.  The model clock
(AMPS PIC time, SWCME launch-relative time) starts at the physical epoch of
the first shock knot, ``model_t0_utc`` in derived.json.
"""

import argparse
import csv
import hashlib
import io
import json
import math
import re
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
AU_M = 1.495978707e11            # IAU 2012 exact [m]
RSUN_M = 6.957e8                 # IAU 2015 nominal solar radius [m] (= SWCME)
C_M_S = 299792458.0              # [m/s]
KB = 1.380649e-23                # [J/K]
# Proton mass exactly as compiled into the AMPS species table by
# input/test/sep_reduced_front_protons.input (CODATA 2018).  The deck's
# spectrum_particle_mass_kg must name the compiled species (main.cpp check).
MP_KG = 1.67262192369e-27        # [kg]
MALPHA_KG = 6.6446573450e-27     # alpha mass, CODATA 2022 [kg]
EV_J = 1.602176634e-19           # [J]
UTC = timezone.utc


def utc(text):
    if not text.endswith('Z'):
        raise ValueError(f'expected an explicit UTC timestamp ending in Z: {text}')
    return datetime.fromisoformat(text[:-1]).replace(tzinfo=UTC)


def iso(t):
    return t.astimezone(UTC).isoformat(timespec='milliseconds').replace('+00:00', 'Z')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_json(name):
    return json.loads((HERE / name).read_text())


# ---------------------------------------------------------------------------
# Raw-archive integrity
# ---------------------------------------------------------------------------

def verified(raw, manifest, key):
    """Return the path of a pinned raw file after checking its SHA-256."""
    entry = manifest['files'][key]
    path = Path(raw) / key
    if not path.is_file():
        raise ValueError(f'missing raw file {path} (source: {entry["source"]})')
    digest = sha256(path)
    if digest != entry['sha256']:
        raise ValueError(f'checksum mismatch for {path}: {digest} != {entry["sha256"]}')
    return path


# ---------------------------------------------------------------------------
# CDF helpers (derive only)
# ---------------------------------------------------------------------------

def cdf_variable(cdf, name):
    """Values with FILLVAL / VALIDMIN / VALIDMAX applied as NaN (never zero)."""
    import numpy as np
    attributes = cdf.varattsget(name)
    values = np.asarray(cdf.varget(name), dtype=float)
    bad = np.zeros(values.shape, dtype=bool)
    if 'FILLVAL' in attributes:
        bad |= values == float(attributes['FILLVAL'])
    if 'VALIDMIN' in attributes:
        bad |= values < float(np.min(attributes['VALIDMIN']))
    if 'VALIDMAX' in attributes:
        bad |= values > float(np.max(attributes['VALIDMAX']))
    values[bad] = np.nan
    return values, attributes


def cdf_epochs(cdf, name='Epoch'):
    import cdflib
    import numpy as np
    stamps = cdflib.cdfepoch.to_datetime(cdf.varget(name))
    return [datetime.utcfromtimestamp(int(x) / 1e9).replace(tzinfo=UTC)
            for x in np.asarray(stamps, dtype='datetime64[ns]').astype('int64')]


def window_mask(epochs, start, stop):
    return [start <= t < stop for t in epochs]


def statistics(values, mask):
    import numpy as np
    selected = np.asarray(values)[np.asarray(mask)]
    good = selected[np.isfinite(selected)]
    if good.size == 0:
        raise ValueError('no valid samples in the averaging window')
    return {'mean': float(np.mean(good)), 'median': float(np.median(good)),
            'std': float(np.std(good)), 'min': float(np.min(good)),
            'max': float(np.max(good)), 'valid': int(good.size),
            'in_window': int(selected.size)}


# ---------------------------------------------------------------------------
# Horizons text records
# ---------------------------------------------------------------------------

def horizons_vectors(path):
    """Parse VEC_TABLE=3 records: list of (UTC datetime, RG [m], LT [s])."""
    text = Path(path).read_text()
    body = text.split('$$SOE', 1)[1].split('$$EOE', 1)[0]
    records = []
    stamp = re.compile(r'A\.D\. (\d{4}-[A-Za-z]{3}-\d{2} \d{2}:\d{2}:\d{2}\.\d+)')
    current = None
    for line in body.splitlines():
        match = stamp.search(line)
        if match:
            current = datetime.strptime(match.group(1), '%Y-%b-%d %H:%M:%S.%f').replace(tzinfo=UTC)
            continue
        if 'RG=' in line and current is not None:
            fields = dict(re.findall(r'(LT|RG|RR)\s*=\s*([-+0-9.Ee]+)', line))
            records.append((current, float(fields['RG']) * 1e3, float(fields['LT'])))
            current = None
    if not records:
        raise ValueError(f'no Horizons vector records in {path}')
    return records


def horizons_subobserver(path, column):
    """(UTC, value) rows of ObsSub-LON (column 2) or ObsSub-LAT (column 3).

    Target frame IAU_SUN, east-longitude positive: its prime meridian is the
    Carrington meridian, so ObsSub-LON is the sub-Earth Carrington longitude.
    """
    text = Path(path).read_text()
    if 'IAU_SUN' not in text or 'East-longitude positive' not in text:
        raise ValueError('sub-observer table is not in the IAU_SUN east-longitude frame')
    body = text.split('$$SOE', 1)[1].split('$$EOE', 1)[0]
    rows = []
    for line in body.splitlines():
        parts = line.split()
        if len(parts) >= 4:
            t = datetime.strptime(' '.join(parts[:2]), '%Y-%b-%d %H:%M').replace(tzinfo=UTC)
            rows.append((t, float(parts[column])))
    return rows


def horizons_subobserver_latitude(path):
    return horizons_subobserver(path, 3)


def ephemeris_at(records, t):
    """Linear interpolation of a smooth ephemeris between bracketing records."""
    for (t0, a0, *_), (t1, a1, *_) in zip(records, records[1:]):
        if t0 <= t <= t1:
            w = (t - t0).total_seconds() / (t1 - t0).total_seconds()
            return a0 + w * (a1 - a0)
    raise ValueError(f'{t} is outside the ephemeris table')


# ---------------------------------------------------------------------------
# derive
# ---------------------------------------------------------------------------

def derive_ambient(raw, manifest, settings):
    """1-AU ambient values from OMNI (V, n_p, T_p, |B|, polarity) and Wind/SWE
    (alpha/proton ratio, T_alpha from quality-10 fits, T_e)."""
    import cdflib
    import numpy as np
    start, stop = (utc(x) for x in settings['clocks']['ambient_window_utc'])

    omni = cdflib.CDF(str(verified(raw, manifest, 'omni/omni_hro_1min_20130401_v01.cdf')))
    t = cdf_epochs(omni)
    mask = window_mask(t, start, stop)
    result = {'window_utc': [iso(start), iso(stop)]}
    for name, key, unit in (('flow_speed', 'V_sw', 'km/s'), ('proton_density', 'n_p', 'cm^-3'),
                            ('T', 'T_p', 'K'), ('F', 'B_magnitude', 'nT'),
                            ('BX_GSE', 'Bx_gse', 'nT'), ('BY_GSE', 'By_gse', 'nT')):
        values, attributes = cdf_variable(omni, name)
        result[key] = dict(statistics(values, mask), unit=unit, source=f'OMNI_HRO_1MIN {name}',
                           catdesc=attributes.get('CATDESC', ''))
    for name in ('PLS', 'IMF'):
        values, _ = cdf_variable(omni, name)
        ids = sorted({int(v) for v in np.asarray(values)[np.asarray(mask)] if np.isfinite(v)})
        result[f'omni_{name.lower()}_source_ids'] = ids

    # Sector polarity from the field itself: B_r = -B_x(GSE) at L1.  An
    # outward (away) sector has mean B_r > 0, i.e. parker.radial_polarity=+1.
    bx, _ = cdf_variable(omni, 'BX_GSE')
    by, _ = cdf_variable(omni, 'BY_GSE')
    sel = np.asarray(mask) & np.isfinite(bx) & np.isfinite(by)
    result['polarity'] = {
        'mean_radial_field_nT': float(-np.mean(bx[sel])),
        'away_fraction': float(np.mean((bx[sel] < 0) & (by[sel] > 0))),
        'toward_fraction': float(np.mean((bx[sel] > 0) & (by[sel] < 0))),
        'rule': 'parker.radial_polarity = sign(mean B_r), B_r = -Bx_GSE',
        'radial_polarity': 1 if -np.mean(bx[sel]) > 0 else -1}

    h1 = cdflib.CDF(str(verified(raw, manifest, 'wind_swe/wi_h1_swe_20130411_v01.cdf')))
    t1 = cdf_epochs(h1)
    m1 = np.asarray(window_mask(t1, start, stop))
    flag, _ = cdf_variable(h1, 'fit_flag')
    np_fit, _ = cdf_variable(h1, 'Proton_Np_nonlin')
    na_fit, _ = cdf_variable(h1, 'Alpha_Na_nonlin')
    wa_fit, _ = cdf_variable(h1, 'Alpha_W_nonlin')
    good = m1 & (flag == 10)
    counts = {int(k): int(v) for k, v in zip(*np.unique(flag[m1][np.isfinite(flag[m1])], return_counts=True))}
    ratio = na_fit / np_fit
    result['alpha_to_proton'] = dict(statistics(ratio, good), unit='1',
        source='WI_H1_SWE Alpha_Na_nonlin/Proton_Np_nonlin, fit_flag==10', fit_flag_counts=counts)
    # W = sqrt(2 k T / M) (SWE most-probable thermal speed, CDAWeb WI_K0_SWE notes).
    talpha = MALPHA_KG * (wa_fit * 1e3) ** 2 / (2.0 * KB)
    result['T_alpha'] = dict(statistics(talpha, good), unit='K',
        source='WI_H1_SWE Alpha_W_nonlin, fit_flag==10, T = M W^2 / (2 k)')

    h5 = cdflib.CDF(str(verified(raw, manifest, 'wind_swe/wi_h5_swe_20130411_v01.cdf')))
    t5 = cdf_epochs(h5)
    te, _ = cdf_variable(h5, 'T_elec')
    result['T_e'] = dict(statistics(te, window_mask(t5, start, stop)), unit='K',
                         source='WI_H5_SWE T_elec (total electron temperature)')

    # SWCME multi_species reads ambient.density_1au as the electron density
    # and closes n_p = n_e / (1 + 2 n_alpha/n_p) (swcme_solarwind.hpp).  The
    # OMNI value is the proton density, so charge neutrality gives n_e.
    f = result['alpha_to_proton']['mean']
    result['deck'] = {
        'wind_speed_km_s': result['V_sw']['mean'],
        'electron_density_1au_cm3': result['n_p']['mean'] * (1.0 + 2.0 * f),
        'magnetic_field_1au_nT': result['B_magnitude']['mean'],
        'proton_temperature_K': result['T_p']['mean'],
        'electron_temperature_K': result['T_e']['mean'],
        'alpha_temperature_K': result['T_alpha']['mean'],
        'alpha_to_proton_ratio': f,
        'radial_polarity': result['polarity']['radial_polarity']}
    return result


def lario_table2_knots(raw, manifest):
    """Shock heights above the L1 footpoint (Lario et al. 2014, Table 2).

    The four L1-connected rows (near-Earth instruments) are transcribed and
    then verified against the pinned text of the preprint so a transcription
    error cannot pass silently.  Their times have the light travel time
    already added (table note a), i.e. they are Earth-frame times.
    """
    text = verified(raw, manifest, 'docs/lario2014_preprint.txt').read_text()
    rows = [
        ('Protons/SOHO/ERNE', '07:58', 9, 3.68, 1.27, r'Protons/SOHO/ERNE\s+07:58.9 min \(VDA\) 3\.68.1\.27'),
        ('Near-relativistic Electrons/WIND/3DP', '07:45', 5, 2.57, 0.90,
         r'Near-relativistic Electrons/WIND/3DP\s+07:45.5 min \(VDA\) 2\.57.0\.90'),
        ('175-315 keV Electrons/ACE/EPAM/DE', '07:47', 2, 2.73, 0.65,
         r'175-315 keV Electrons/ACE/EPAM/DE\s+07:47.2 min \(TSA\)\s+2\.73.0\.65'),
        ('0.25-0.70 MeV Electrons/SOHO/EPHIN', '07:35', 2, 1.52, 0.63,
         r'0\.25-0\.70 MeV Electrons/SOHO/EPHIN\s+07:35.2 min \(TSA\)\s+1\.52.0\.63'),
    ]
    knots = []
    for label, hhmm, sigma_min, height, sigma_h, pattern in rows:
        if not re.search(pattern, text):
            raise ValueError(f'Lario Table 2 row not found in the pinned preprint text: {label}')
        knots.append({'source': f'Lario et al. (2014) ApJ 797:8, Table 2, {label}',
                      'listed_utc_light_time_added': f'2013-04-11T{hhmm}:00Z',
                      'time_uncertainty_min': sigma_min,
                      'height_above_surface_rsun': height, 'height_uncertainty_rsun': sigma_h})
    return sorted(knots, key=lambda k: k['listed_utc_light_time_added'])


def derive_knots(raw, manifest):
    soho = horizons_vectors(verified(raw, manifest, 'horizons/vectors_soho.txt'))
    knots = []
    for k in lario_table2_knots(raw, manifest):
        listed = utc(k['listed_utc_light_time_added'])
        # Light travel time Sun -> SOHO at the listed epoch (Horizons LT).
        light_s = ephemeris_at([(t, lt) for t, _, lt in soho], listed)
        physical = listed - timedelta(seconds=light_s)
        knots.append(dict(k, light_time_subtracted_s=light_s, physical_utc=iso(physical),
                          radius_rsun=1.0 + k['height_above_surface_rsun'],
                          _physical=physical))
    # In-situ shock at Wind: CfA Interplanetary Shock Database (2013 list),
    # fast-forward shock at day-of-year 103.925955 = 2013-04-13 22:13:22.5 UT.
    cfa = verified(raw, manifest, 'catalogs/wi_2013.html').read_text(errors='replace')
    row = re.sub(r'<[^>]*>', ' ', cfa)
    row = re.sub(r'\s+', ' ', row)
    if not re.search(r'2013 04 13 103\.925955 2213 \S+ \S+ \S+ FF', row):
        raise ValueError('CfA Wind shock 2013-04-13 DOY 103.925955 (FF) not found in pinned list')
    icme = re.sub(r'<[^>]*>', ' ', verified(raw, manifest, 'catalogs/icmetable2.htm').read_text(errors='replace'))
    icme = re.sub(r'\s+', ' ', icme)
    if not re.search(r'2013/04/13 2254 .{0,160}2013/04/11 0724', icme):
        raise ValueError('Richardson-Cane association of 2013/04/13 2254 with the 2013/04/11 0724 CME not found')
    shock_t = datetime(2013, 1, 1, tzinfo=UTC) + timedelta(days=103.925955 - 1.0)
    wind = horizons_vectors(verified(raw, manifest, 'horizons/wind_at_shock.txt'))
    if len(wind) != 1 or abs((wind[0][0] - shock_t).total_seconds()) > 0.5:
        raise ValueError('Wind ephemeris is not evaluated at the CfA shock epoch')
    knots.append({'source': 'CfA Interplanetary Shock Database, Wind 2013 list, FF shock DOY 103.925955; '
                            'Richardson & Cane ICME list associates the 2013/04/13 22:54 disturbance '
                            'with the 2013/04/11 07:24 CME',
                  'physical_utc': iso(shock_t), 'time_uncertainty_min': None,
                  'radius_rsun': wind[0][1] / RSUN_M,
                  'wind_heliocentric_distance_m': wind[0][1], '_physical': shock_t})
    # The model clock starts on a whole UTC second (the reduced-front event
    # label run.reference_epoch has 1 s resolution) at the first knot rounded
    # to the nearest second; all knot epochs are rounded to whole seconds.
    # The reconstructed epochs are uncertain by 2-9 min (Table 2), so the
    # <= 0.5 s rounding is four orders of magnitude below the data precision.
    first = knots[0]['_physical']
    t0 = (first + timedelta(microseconds=500000)).replace(microsecond=0)
    for k in knots:
        k['physical_minus_model_clock_s'] = round((k['_physical'] - t0).total_seconds(), 6)
        k['model_time_s'] = int(round((k.pop('_physical') - t0).total_seconds()))
    if knots[0]['model_time_s'] != 0:
        raise ValueError('first knot does not round onto the model clock origin')
    knots[0]['model_t0'] = t0
    radii = [k['radius_rsun'] for k in knots]
    times = [k['model_time_s'] for k in knots]
    if any(b <= a for a, b in zip(times, times[1:])) or any(b < a for a, b in zip(radii, radii[1:])):
        raise ValueError('shock knots must have increasing time and nondecreasing radius')
    return knots


def derive_geometry(raw, manifest, settings, t0, t_end):
    out = {}
    for name in ('soho', 'ace', 'earth'):
        rec = horizons_vectors(verified(raw, manifest, f'horizons/vectors_{name}.txt'))
        inside = [r for t, r, _ in rec if t0 - timedelta(hours=1) <= t <= t_end + timedelta(hours=1)]
        out[f'{name}_heliocentric_distance_m'] = {
            'mean': sum(inside) / len(inside), 'min': min(inside), 'max': max(inside),
            'samples': len(inside), 'rule': 'mean of hourly DE441 vectors over the run interval'}
    subobs = verified(raw, manifest, 'horizons/sun_from_earth_subobs.txt')
    b0 = ephemeris_at(horizons_subobserver_latitude(subobs), t0)
    out['earth_heliographic_latitude_deg'] = b0
    # Sub-Earth Carrington longitude at the model origin (monotone, no wrap in
    # the 6-h table over this interval, so linear interpolation is exact to
    # the table's smoothness).
    out['earth_carrington_longitude_deg'] = ephemeris_at(horizons_subobserver(subobs, 2), t0)
    # The pinned ACE ephemeris starts at 06:00 UT, so the L1 sampling radius
    # of the ambient values is taken at the model origin; the hourly table
    # change is recorded to show it is negligible against the 1% scale of r.
    ace = horizons_vectors(verified(raw, manifest, 'horizons/vectors_ace.txt'))
    out['ace_heliocentric_distance_t0_m'] = ephemeris_at([(t, r) for t, r, _ in ace], t0)
    out['ace_heliocentric_distance_hourly_change_m'] = max(
        abs(b[1] - a[1]) for a, b in zip(ace, ace[1:]))
    # SWCME/srcSEP: parker.sin_theta is the sine of the heliographic colatitude.
    out['sin_colatitude'] = math.cos(math.radians(b0))
    return out


def erne_windows(raw, manifest, settings, t0):
    """ERNE HED proton intensities averaged over model-grid windows.

    Each requested window [centre - w/2, centre + w/2) is moved to the model
    output grid (multiples of time_step * output_cadence_steps from t0) by
    choosing the grid start whose window centre is nearest the requested
    centre.  A 1-min sample is used if its whole minute lies inside the
    window, its intensity is finite and inside VALIDMIN/VALIDMAX, and, for
    PH3-PH10, status bit b12 (HED HIT 3x, mask 0x0008 with b0 the most
    significant bit) is clear (export description section 7.2 note 2).  In
    b12 minutes the archive replaces those intensities by time-interpolated
    values, which are not measurements.
    """
    import cdflib
    import numpy as np
    numerics = settings['numerics']
    cadence = numerics['time_step_s'] * numerics['output_cadence_steps']
    reference = utc(settings['clocks']['comparison_reference_utc'])
    minutes = []
    for day in ('20130411', '20130412', '20130413'):
        cdf = cdflib.CDF(str(verified(raw, manifest, f'erne/soho_erne-hed_l2-1min_{day}_v01.cdf')))
        if cdf.varattsget('PH')['UNITS'] != '(cm^2 s sr MeV)^-1':
            raise ValueError('unexpected ERNE PH units')
        ph, _ = cdf_variable(cdf, 'PH')
        est = np.asarray(cdf.varget('est')).astype(int)
        for t, row, word in zip(cdf_epochs(cdf), ph, est):
            minutes.append((t, row, int(word)))
    windows = []
    for spec in settings['windows']:
        centre = reference + timedelta(hours=spec['hours_after_reference'])
        half = spec['width_s'] / 2.0
        k = round(((centre - t0).total_seconds() - half) / cadence)
        start_s = k * cadence
        stop_s = start_s + spec['width_s']
        if abs(stop_s - start_s - spec['width_s']) > 1e-9 or spec['width_s'] % cadence:
            raise ValueError('window width must be a whole number of model outputs')
        start, stop = t0 + timedelta(seconds=start_s), t0 + timedelta(seconds=stop_s)
        channels = []
        for index, (lo, hi) in enumerate(settings['erne_revised_f50_channels_mev']):
            values, rejected = [], {'status_b12': 0, 'invalid': 0}
            for t, row, word in minutes:
                if not (start <= t and t + timedelta(minutes=1) <= stop):
                    continue
                if index >= 2 and word & 0x0008:
                    rejected['status_b12'] += 1
                    continue
                if not np.isfinite(row[index]):
                    rejected['invalid'] += 1
                    continue
                values.append(float(row[index]))
            if not values:
                raise ValueError(f'no usable ERNE minutes for {spec["id"]} channel {lo}-{hi} MeV')
            mean = sum(values) / len(values)
            sd = math.sqrt(sum((v - mean) ** 2 for v in values) / max(1, len(values) - 1))
            channels.append({'channel': f'PH{index + 1}', 'energy_low_mev': lo, 'energy_high_mev': hi,
                             'mean_intensity_per_cm2_s_sr_mev': mean, 'sample_std': sd,
                             'minutes_used': len(values), 'rejected': rejected})
        windows.append({'id': spec['id'], 'role': spec['role'], 'requested_centre_utc': iso(centre),
                        'model_start_s': start_s, 'model_stop_s': stop_s,
                        'start_utc': iso(start), 'stop_utc': iso(stop), 'channels': channels})
    return windows


def wso_radial_coefficients(raw, manifest):
    """Parse the WSO radial-model PFSS coefficients (Harmonic.rad/CR2135).

    Definition (WSO 'Notes on PFSS Extrapolation', X. Sun, eq. 10): the
    photospheric field is B_r(R0,th,ph) = sum P_l^m(cos th)(g cos m ph +
    h sin m ph), independent of the source-surface radius, with WSO's
    Schmidt semi-normalized P_l^m WITHOUT the Condon-Shortley phase (eq. 24).
    Returns {(l, m): (g, h)} in the file's units (not stated in the file).
    """
    text = verified(raw, manifest, 'wso/Harmonic.rad.CR2135.txt').read_text()
    if 'CT2135:360' not in text:
        raise ValueError('WSO file is not the CR2135 radial-model expansion')
    g_part, h_part = text.split("h's", 1) if "h's" in text else (None, None)
    if g_part is None:
        raise ValueError('WSO file lacks the h coefficient block')
    # The h block starts at the last "0 0.000" line before the "h's" label.
    lines = [ln for ln in text.splitlines() if ln.strip()]
    rows = []
    for ln in lines:
        nums = re.findall(r'-?\d+\.\d+|-?\d+', ln.split("g's")[0].split("h's")[0])
        rows.append([float(x) for x in nums])
    if len(rows) != 20:
        raise ValueError(f'expected 2 x 10 coefficient rows, found {len(rows)}')
    coefficients = {}
    for block in (0, 1):
        for l in range(10):
            row = rows[block * 10 + l]
            if int(row[0]) != l or len(row) != l + 2:
                raise ValueError(f'malformed WSO coefficient row for l={l}')
            for m, value in enumerate(row[1:]):
                g, h = coefficients.get((l, m), (0.0, 0.0))
                coefficients[(l, m)] = (value, h) if block == 0 else (g, value)
    return coefficients


def wso_to_provider_rows(coefficients, scale_t_per_unit):
    """Convert WSO Schmidt (no CS phase) to the provider's basis.

    Provider (sep_coronal_cme/pfss_harmonics.{h,cpp}): real orthonormal with
    the Condon-Shortley phase, N_lm P_l^m cos(m ph), N_lm =
    sqrt((2l+1)(2-d_m0)(l-m)!/(4 pi (l+m)!)).  WSO eq. 24 gives
    P_wso = (-1)^m sqrt((2-d_m0)(l-m)!/(l+m)!) P_CS, hence
        b_lm = (-1)^m sqrt(4 pi / (2l+1)) g_lm * scale     (same for h).
    The l=0 monopole is dropped: the provider rejects degree < 1 (div B = 0)
    and, on a Parker-exterior line whose B_r = G(footpoint)/r^2, its effect
    is a pure amplitude factor that the single |B| calibration absorbs.
    """
    rows = []
    for (l, m), (g, h) in sorted(coefficients.items()):
        if l == 0:
            continue
        factor = ((-1) ** m) * math.sqrt(4.0 * math.pi / (2 * l + 1)) * scale_t_per_unit
        rows.append((l, m, g * factor, 0.0 if m == 0 else h * factor))
    return rows


def derive_goes(raw, manifest):
    """GOES-13 >10 MeV integral-flux event record (CCMC SEP campaign, OpSEP).

    OpSEP calculate_fluence multiplies Int J dt by 4 pi sr, so the published
    fluence [cm^-2] is 4 pi Int J dt (docs/opsep/operational_sep_quantities.py).
    """
    record = json.loads(verified(raw, manifest, 'ccmc/goes-13.json').read_text())
    submission = record['sep_observation_submission']
    if submission['observatory']['short_name'] != 'GOES-13' or \
            submission['observatory']['flux_type'] != 'integral':
        raise ValueError('unexpected CCMC observatory record')
    for item in submission['observations']:
        channel = item['energy_channel']
        if channel['min'] == 10.0 and channel['max'] == -1 and channel['units'] == 'MeV':
            if item['peak_intensity']['units'] != 'pfu' or item['fluences'][0]['units'] != 'cm^-2':
                raise ValueError('unexpected CCMC units')
            length = item['event_lengths'][0]
            return {'source': 'CCMC SEP Scoreboard campaign record (ccmc/goes-13.json), OpSEP measurement',
                    'peak_pfu': item['peak_intensity']['intensity'],
                    'peak_utc': item['peak_intensity']['time'],
                    'event_start_utc': length['start_time'], 'event_end_utc': length['end_time'],
                    'threshold_pfu': length['threshold'],
                    'fluence_cm2': item['fluences'][0]['fluence'],
                    'fluence_definition': '4 pi Int J dt over the event interval (OpSEP calculate_fluence)'}
    raise ValueError('no >10 MeV GOES-13 record in the CCMC file')


def derive(args):
    settings = load_json('event_settings.json')
    manifest = load_json('raw_manifest.json')
    knots = derive_knots(args.raw, manifest)
    t0 = knots[0].pop('model_t0')
    numerics = settings['numerics']
    end_needed = utc(settings['clocks']['run_end_not_before_utc'])
    steps = math.ceil((end_needed - t0).total_seconds() / numerics['time_step_s'])
    cadence_steps = numerics['output_cadence_steps']
    steps = cadence_steps * math.ceil(steps / cadence_steps)
    t_end = t0 + timedelta(seconds=steps * numerics['time_step_s'])
    derived = {
        'schema': 'srcsep-event-2013-04-11-derived-v1',
        'generator': 'srcSEP/test/event_2013_04_11/prepare_event.py derive',
        'settings_sha256': sha256(HERE / 'event_settings.json'),
        'raw_manifest_sha256': sha256(HERE / 'raw_manifest.json'),
        'model_t0_utc': iso(t0),
        'model_t0_note': 'first shock knot (light-time corrected) rounded to the whole second; AMPS and reduced-front time zero',
        'total_iterations': steps,
        'run_end_utc': iso(t_end),
        'ambient': derive_ambient(args.raw, manifest, settings),
        'shock_knots': knots,
        'geometry': derive_geometry(args.raw, manifest, settings, t0, t_end),
        'erne_windows': erne_windows(args.raw, manifest, settings, t0),
        'goes13_gt10mev': derive_goes(args.raw, manifest),
    }
    (HERE / 'derived.json').write_text(json.dumps(derived, indent=1, sort_keys=True) + '\n')
    write_observations(derived, settings)
    print('wrote derived.json and observations.csv')


# ---------------------------------------------------------------------------
# Observation table (comparison operator inputs)
# ---------------------------------------------------------------------------

def observer_id(instrument, lo, hi):
    """Stable deck observer ID for one instrument channel (no '.' allowed)."""
    tag = re.sub(r'[^a-z0-9]+', '-', instrument.lower()).strip('-')
    fmt = lambda x: f'{x:.6g}'.replace('.', 'p')
    return f'{tag}-{fmt(lo)}-{fmt(hi)}mev'


def write_observations(derived, settings):
    """observations.csv: one row per model-observation comparison.

    quantity=differential: intensity [(cm^2 s sr MeV)^-1] averaged over
    [model_start_s, model_stop_s); quantity=integral_peak: max of the window
    means of the >10 MeV integral intensity [pfu] over that interval;
    quantity=integral_fluence: 4 pi Int J dt [cm^-2] (OpSEP definition).
    """
    t0 = utc(derived['model_t0_utc'])
    numerics = settings['numerics']
    cadence = numerics['time_step_s'] * numerics['output_cadence_steps']
    rows = []
    for window in derived['erne_windows']:
        for ch in window['channels']:
            rows.append({'id': f'{window["id"]}-{ch["channel"]}', 'role': window['role'],
                         'instrument': 'SOHO/ERNE-HED', 'observer': observer_id('erne-hed', ch['energy_low_mev'], ch['energy_high_mev']),
                         'quantity': 'differential', 'energy_low_mev': ch['energy_low_mev'],
                         'energy_high_mev': ch['energy_high_mev'], 'model_start_s': window['model_start_s'],
                         'model_stop_s': window['model_stop_s'],
                         'value': ch['mean_intensity_per_cm2_s_sr_mev'], 'unit': '(cm^2 s sr MeV)^-1',
                         'note': f'{ch["minutes_used"]} usable 1-min samples'})
    liu = HERE.parent.parent / 'validation/cases/XM03/reference/liu_figure12_earth_observations.csv'
    reference = utc(settings['clocks']['comparison_reference_utc'])
    rule = settings['liu_holdout']
    with liu.open(newline='') as stream:
        for n, row in enumerate(csv.DictReader(stream)):
            hours = float(row['elapsed_hours'])
            lo, hi = float(row['energy_low_mev']), float(row['energy_high_mev'])
            instant = (reference + timedelta(hours=hours) - t0).total_seconds()
            k = math.floor(instant / cadence)
            role = 'holdout' if hours in rule['elapsed_hours'] and lo >= rule['minimum_energy_mev'] else 'excluded'
            rows.append({'id': f'liu-{int(hours):02d}h-{n:02d}', 'role': role,
                         'instrument': row['instrument'], 'observer': observer_id(row['instrument'], lo, hi),
                         'quantity': 'differential', 'energy_low_mev': lo, 'energy_high_mev': hi,
                         'model_start_s': k * cadence, 'model_stop_s': (k + 1) * cadence,
                         'value': float(row['differential_intensity_pfu_per_mev']),
                         'unit': '(cm^2 s sr MeV)^-1',
                         'note': 'Liu et al. (2025) Fig. 12 (vector-digitized, XM03); model output window containing the snapshot'})
    goes = derived['goes13_gt10mev']
    start, stop = utc(goes['event_start_utc']), utc(goes['event_end_utc'])
    a = math.ceil((start - t0).total_seconds() / cadence) * cadence
    b = math.floor((stop - t0).total_seconds() / cadence) * cadence
    for quantity, value, unit in (('integral_peak', goes['peak_pfu'], 'pfu'),
                                  ('integral_fluence', goes['fluence_cm2'], 'cm^-2')):
        rows.append({'id': f'goes13-gt10mev-{quantity.split("_")[1]}', 'role': 'holdout',
                     'instrument': 'GOES-13 (CCMC/OpSEP)', 'observer': observer_id('goes-13-integral', 10.0, 200.0),
                     'quantity': quantity, 'energy_low_mev': 10.0, 'energy_high_mev': 200.0,
                     'model_start_s': a, 'model_stop_s': b, 'value': value, 'unit': unit,
                     'note': 'model integrates to E_max=200 MeV (source cutoff); windows fully inside the OpSEP event interval'})
    fields = ['id', 'role', 'instrument', 'observer', 'quantity', 'energy_low_mev', 'energy_high_mev',
              'model_start_s', 'model_stop_s', 'value', 'unit', 'note']
    buffer = io.StringIO()
    table = csv.DictWriter(buffer, fieldnames=fields, lineterminator='\n')
    table.writeheader()
    for row in rows:
        table.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in row.items()})
    (HERE / 'observations.csv').write_text(buffer.getvalue())


# ---------------------------------------------------------------------------
# deck
# ---------------------------------------------------------------------------

def observers_for_deck(settings):
    """Unique observer channels with the heliocentric radius of their host."""
    with (HERE / 'observations.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    derived = load_json('derived.json')
    radius = {'SOHO/ERNE-HED': 'soho', 'SOHO/ERNE': 'soho', 'ACE/EPAM': 'ace',
              'GOES-13/EPEAD': 'earth', 'GOES-13 (CCMC/OpSEP)': 'earth'}
    seen = {}
    for row in rows:
        if row['role'] == 'excluded':
            continue
        host = radius[row['instrument']]
        r = derived['geometry'][f'{host}_heliocentric_distance_m']['mean']
        seen.setdefault(row['observer'], (r, float(row['energy_low_mev']), float(row['energy_high_mev'])))
    return seen


# ---------------------------------------------------------------------------
# Reduced-front event (shared core, schema shock-front-ambient-v1.3)
# ---------------------------------------------------------------------------
#
# The physics is evaluated by the C++ core (src/models/sep_corona_swcme).
# The field line is anchored at Earth at t0 (mesh.field_line.anchor_*): srcSEP
# traces the core field from Earth sunward through the Parker exterior and the
# open PFSS shell to the 1.1 Rs base at run time, so Python never evaluates
# the magnetogram.  Python only sizes the line (its length must reach past
# the outer sphere) and reports the exterior winding; for that it repeats the
# core's isothermal Parker (1958) wind with the core's constants.

GM_SUN = 1.32712440018e20   # sep_coronal_cme Constants::kSolarGravitationalParameterM3PerS2
OMEGA_SUN = 2.86533e-6      # provider rotation rate [rad/s] (sidereal Carrington)
# PFSS source surface [Rs]: the WSO standard (Hoeksema 1984; WSO radial-model
# coefficients are independent of Rs) and the footpoint method of Lario et al.
# (2014, Sect. 2.1, after Park et al. 2013).  User decision 2026-11-25.
RSS_RS = 2.5
WSO_T_PER_UNIT = 1.0e-6     # WSO microtesla; the single |B| calibration fixes the amplitude


def parker_y(r, rc):
    """(U/c)^2 on the transonic branch: y - ln y = 4 ln(r/rc) + 4 rc/r - 3."""
    rhs = 4.0 * math.log(r / rc) + 4.0 * rc / r - 3.0
    if abs(r - rc) <= 1e-13 * rc:
        return 1.0
    lo, hi = (1e-300, 1.0) if r < rc else (1.0, max(4.0, rhs + math.log(max(2.0, rhs)) + 2.0))
    while r > rc and hi - math.log(hi) < rhs:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f = mid - math.log(mid) - rhs
        # y - ln y decreases on (0,1) and increases on (1,inf).
        if (f > 0) == (r < rc):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def parker_speed(r, c):
    return c * math.sqrt(parker_y(r, GM_SUN / (2.0 * c * c)))


def fit_sound_speed(r_ref, u_ref):
    """Isothermal sound speed with U(r_ref) = u_ref (unique: U increases with c)."""
    lo = math.log(math.sqrt(GM_SUN / (4.0 * r_ref)))
    hi = lo + math.log(100.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if parker_speed(r_ref, math.exp(mid)) < u_ref:
            lo = mid
        else:
            hi = mid
    return math.exp(0.5 * (lo + hi))


def gauss_log_integral(f, r0, r1, panels=4000):
    """Int_r0^r1 f(r) dr with 4-point Gauss-Legendre panels uniform in ln r."""
    nodes = (-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526)
    weights = (0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538)
    a, b = math.log(r0), math.log(r1)
    h = (b - a) / panels
    total = 0.0
    for i in range(panels):
        mid = a + (i + 0.5) * h
        for x, wgt in zip(nodes, weights):
            r = math.exp(mid + 0.5 * h * x)
            total += wgt * f(r) * r
    return 0.5 * h * total


def winding(r, rss, c):
    """Parker-exterior azimuth lag of the source footpoint: phi_foot = phi + W(r)."""
    if r <= rss:
        return 0.0
    return gauss_log_integral(lambda x: OMEGA_SUN * (1.0 - (rss / x) ** 2) / parker_speed(x, c), rss, r)


def spiral_arc_length(r0, r1, rss, c, sin_theta):
    """Arc length of the Parker-exterior field line between radii r0 < r1."""
    return gauss_log_integral(lambda x: math.sqrt(1.0 + (
        OMEGA_SUN * sin_theta * (x - rss * rss / x) / parker_speed(x, c)) ** 2), r0, r1)


def model_frame_harmonics(raw, carrington_longitude_deg):
    """WSO radial-model coefficients in the model frame, provider convention.

    Model frame: inertial HCI with +z the solar rotation axis and +x the
    projection of the Sun->Earth direction at t0, so at t=0 (where the core's
    body frame equals HCI) body azimuth = Carrington longitude - L_E(t0).
    g cos m(phi+c) + h sin m(phi+c) = g' cos m phi + h' sin m phi with
    g' = g cos mc + h sin mc, h' = h cos mc - g sin mc.
    """
    manifest = load_json('raw_manifest.json')
    wso = wso_radial_coefficients(raw, manifest)
    c = math.radians(carrington_longitude_deg)
    rotated = {}
    for (l, m), (g, h) in wso.items():
        rotated[(l, m)] = (g * math.cos(m * c) + h * math.sin(m * c),
                           h * math.cos(m * c) - g * math.sin(m * c))
    return wso, wso_to_provider_rows(rotated, WSO_T_PER_UNIT)


def write_asset(rows, path):
    lines = ['# WSO radial-model PFSS coefficients, Carrington rotation 2135 (raw/wso/Harmonic.rad.CR2135.txt),',
             '# rotated to the model HCI frame (x toward Earth at t0) and converted to the',
             '# provider real-orthonormal Condon-Shortley basis; l=0 monopole excluded (div B = 0).',
             f'# Amplitude unit {WSO_T_PER_UNIT!r} T per WSO unit; the event calibration fixes |B| at L1.',
             'degree,order,cosine_T,sine_T']
    lines += [f'{l},{m},{gc!r},{hs!r}' for l, m, gc, hs in rows]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(lines) + '\n')


def deck(args):
    settings = load_json('event_settings.json')
    derived = load_json('derived.json')
    a = derived['ambient']['deck']
    g = derived['geometry']
    n = settings['numerics']
    s = settings['source']
    out = Path(args.output_dir) if getattr(args, 'output_dir', None) else HERE
    rss = RSS_RS * RSUN_M
    if settings['decisions']['source_surface_radius_rs'] != RSS_RS:
        raise ValueError('event_settings.json source surface differs from RSS_RS')
    # The user-decision record and the deck inputs must state the same values.
    if settings['decisions']['source_energy_mev'] != [s['minimum_energy_mev'], s['maximum_energy_mev']]:
        raise ValueError('event_settings.json: decisions.source_energy_mev differs from the source band')
    if settings['decisions']['release_distance_shock_radius_fraction'] != s['release_distance_shock_radius_fraction']:
        raise ValueError('event_settings.json: decisions and source release fractions differ')
    r_inner = settings['decisions']['field_line_inner_radius_rs'] * RSUN_M
    r_outer = n['outer_radius_au'] * AU_M
    r_cover = n['coverage_radius_au'] * AU_M
    knots = derived['shock_knots']
    b0 = math.radians(g['earth_heliographic_latitude_deg'])
    theta = 0.5 * math.pi - b0                       # colatitude of Earth at t0
    sin_theta = math.sin(theta)
    r_ref = g['ace_heliocentric_distance_t0_m']
    r_earth = g['earth_heliocentric_distance_m']['mean']
    u_ref = a['wind_speed_km_s'] * 1e3
    c_w = fit_sound_speed(r_ref, u_ref)
    # Exterior winding of the line through Earth at t0 (azimuth 0) down to the
    # source surface; the PFSS part below Rss is traced by srcSEP at run time.
    phi_foot = winding(r_earth, rss, c_w)
    # Line length: exterior spiral from Rss to line_end_radius plus the radial
    # depth of the PFSS shell. A non-radial PFSS path is longer, which only
    # shortens the reach past the outer sphere (margin line_end - outer); the
    # run reports the actual end radius.
    length = (rss - r_inner) + spiral_arc_length(rss, n['line_end_radius_au'] * AU_M, rss, c_w, sin_theta)
    # Geometric spacing a q^k with chord ~ eps*r: a = eps*r_base and q = 1+eps
    # give chord_k = a + eps*s_k; the point count closes the length with q~1+eps.
    eps = n['relative_segment_length']
    first_arc = eps * r_inner
    point_count = 1 + round(math.log(1.0 + length * eps / first_arc) / math.log(1.0 + eps))
    end_s = derived['total_iterations'] * n['time_step_s']
    wind_knot = knots[-1]
    _, rows = model_frame_harmonics(args.raw, g['earth_carrington_longitude_deg'])
    asset = out / 'assets' / 'wso_cr2135_model_frame.harmonics.csv'
    write_asset(rows, asset)
    unit = (math.cos(b0), 0.0, math.sin(b0))
    vec = lambda r: ','.join(repr(r * x) for x in unit)
    ev = []
    e = ev.append
    e('! Reduced shock-front event for the 2013-04-11 real-event nightly test.')
    e('! GENERATED by srcSEP/test/event_2013_04_11/prepare_event.py deck; do not edit.')
    e(f'! Model time 0 = {derived["model_t0_utc"]}; HCI x axis toward Earth at t0.')
    e('#begin physics: shock_front')
    e('  schema = shock-front-ambient-v1.3')
    e('  profile = heliocentric-sphere-pfss-parker-pchip-v1')
    label = utc(derived['model_t0_utc'])
    if label.microsecond:
        raise ValueError('model clock origin is not a whole UTC second')
    e(f'  run.reference_epoch = {label.strftime("%Y-%m-%dT%H:%M:%SZ")}')
    e('  run.time_scale = UTC-resolved-to-elapsed-SI-seconds')
    e('  run.coordinate_frame = HCI')
    e('  run.start_s = 0')
    e('  ! Core coverage runs to the last observed knot (Wind shock); the')
    e(f'  ! particle run ends earlier, at {end_s!r} s.')
    e(f'  run.end_s = {wind_knot["model_time_s"]!r}')
    e(f'  run.background_dt_s = {n["front_cache_cadence_s"]!r}')
    e('  run.particle_mode = disabled')
    e(f'  domain.solar_radius_m = {RSUN_M!r}')
    e(f'  domain.first_valid_plasma_radius_m = {RSUN_M!r}')
    e(f'  domain.coverage_radius_m = {r_cover!r}')
    e('  ! PFSS source surface at 2.5 Rs (WSO standard).  Below it the line runs')
    e('  ! along open PFSS field lines; every front point (>= 2.52 Rs) is in the')
    e('  ! Parker exterior.')
    e(f'  ambient.source_surface_radius_m = {rss!r}')
    e('  ! L1 (ACE) sampling radius of the OMNI ambient values.')
    e(f'  ambient.reference_radius_m = {r_ref!r}')
    e(f'  ambient.electron_density_at_reference_m3 = {a["electron_density_1au_cm3"] * 1e6!r}')
    e('  ! Hydrostatic closed-loop base density: used only on closed PFSS tubes,')
    e('  ! which no material vertex or front point may enter (such a query is')
    e('  ! rejected); 1e8 cm^-3 is the maintained value.')
    e('  ambient.closed_base_electron_density_m3 = 1e14')
    e(f'  ambient.rotation_rate_rad_per_s = {OMEGA_SUN!r}')
    e('  ambient.minimum_magnetic_field_t = 1e-15')
    e(f'  ambient.trace_step_m = {n["pfss_trace_step_m"]!r}')
    e('  ambient.gradient_relative_step = 1e-5')
    e('  ambient.trace_maximum_steps = 10000')
    e(f'  ambient.wind_table_points = {n["wind_table_points"]}')
    e('  ! Open PFSS tubes carry the same isothermal wind: field-aligned in the')
    e('  ! corotating frame, rho u/|B| constant, isothermal Bernoulli along the tube.')
    e('  ambient.open_wind_density_model = isothermal-flux-tube')
    e('  ! Isothermal Parker transonic wind fitted to the OMNI speed at r_ref;')
    e('  ! it sets U and rho only (plasma.* temperatures stay the measured values).')
    e(f'  ambient.wind_speed_at_reference_m_s = {u_ref!r}')
    e('  ! One magnetogram amplitude: |B| = OMNI F at L1 (Sun-Earth line) at t0.')
    e(f'  ambient.magnetic_calibration_position_hci_m = {vec(r_ref)}')
    e('  ambient.magnetic_calibration_epoch_s = 0')
    e(f'  ambient.magnetic_calibration_magnitude_t = {a["magnetic_field_1au_nT"] * 1e-9!r}')
    e('  plasma.gamma_adiabatic = 1.6666666666666667')
    e(f'  plasma.electron_temperature_k = {a["electron_temperature_K"]!r}')
    e(f'  plasma.proton_temperature_k = {a["proton_temperature_K"]!r}')
    e(f'  plasma.alpha_temperature_k = {a["alpha_temperature_K"]!r}')
    e(f'  plasma.alpha_to_proton_number_ratio = {a["alpha_to_proton_ratio"]!r}')
    e('  plasma.include_electron_mass = true')
    e('  geometry.model = heliocentric-sphere')
    e(f'  geometry.physical_inner_radius_m = {r_inner!r}')
    e('  ! Shock radius versus model time: Lario et al. (2014) Table 2 heights above')
    e('  ! the L1 footpoint (light-time corrected) and the Wind in-situ shock.')
    e('  history.model = data-driven-pchip')
    e('  history.knot_times_s = ' + ','.join(repr(float(k['model_time_s'])) for k in knots))
    e('  history.knot_radii_m = ' + ','.join(repr(k.get('wind_heliocentric_distance_m', k['radius_rsun'] * RSUN_M)) for k in knots))
    e('  history.validity_policy = report-only')
    e('  diagnostics.criticality_model = none')
    e(f'  numerics.polar_cells = {n["front_polar_cells"]}')
    e(f'  numerics.azimuth_cells = {n["front_azimuth_cells"]}')
    e('  numerics.surface_topology = triangular-sphere-v1')
    e('  numerics.weak_mach_tolerance = 1e-6')
    e('  numerics.rh_residual_tolerance = 1e-8')
    e('  ! Endpoint = the Wind shock knot; the chart pole is the Earth direction.')
    e(f'  endpoint.apex_radius_m = {wind_knot["wind_heliocentric_distance_m"]!r}')
    e(f'  endpoint.observer_position_hci_m = {vec(wind_knot["wind_heliocentric_distance_m"])}')
    e('  endpoint.require_fast_shock_at_observer = false')
    e(f'  assets.harmonics_file = assets/{asset.name}')
    e(f'  assets.harmonics_sha256 = {sha256(asset)}')
    e('#end physics')
    (out / 'event.event').write_text('\n'.join(ev) + '\n')

    nt = settings['nightly']
    nightly_eps = nt['relative_segment_length']
    nightly_first = nightly_eps * r_inner
    nightly_points = 1 + round(math.log(1.0 + length * nightly_eps / nightly_first) / math.log(1.0 + nightly_eps))
    mev = 1.0e6 * EV_J

    def runtime_deck(nightly):
        """Lines of event.in (long tier, nightly=False) or event.nightly.in (nightly tier)."""
        steps = nt['time_steps'] if nightly else derived['total_iterations']
        samples = nt['samples_per_step'] if nightly else n['samples_per_step']
        grid_eps, grid_points, grid_first = ((nightly_eps, nightly_points, nightly_first) if nightly
                                             else (eps, point_count, first_arc))
        lines = []
        w = lines.append
        if nightly:
            w('! srcSEP schema-6 runtime deck for the 2013-04-11 real-event test, NIGHTLY tier:')
            w('! the same event, physics and command line as event.in; only the grid, the source')
            w('! samples per step, the number of steps and the inner regression observers differ.')
        else:
            w('! srcSEP schema-6 runtime deck for the 2013-04-11 real-event test, LONG tier (full event).')
        w('! GENERATED by srcSEP/test/event_2013_04_11/prepare_event.py deck; do not edit.')
        w(f'! Model time 0 = {derived["model_t0_utc"]}.  Run with --reduced-shock-event event.event:')
        w('! the shared reduced core then owns B, U, rho, T and the shock everywhere.')
        w('#begin run')
        w(f'  time_step_s = {n["time_step_s"]!r}')
        w(f'  maximum_time_steps = {steps}')
        w('#end run')
        w('')
        w('#begin application: sep')
        w('  schema_version = 6')
        w('  ! The mover is selected on the command line (--particle-mover).')
        w('  transport = parker')
        w('  #begin turbulence')
        w('    ! Transport coefficients come from the prescribed lambda law on the')
        w('    ! command line; the turbulence state is prescribed and not evolved.')
        for key, value in (('source', 'prescribed'), ('representation', 'integrated'), ('coupling', 'disabled'),
                           ('advection', 'false'), ('reflection', 'false'), ('cascade', 'false'),
                           ('shock_injection', 'false'), ('inner_boundary', 'fixed-reservoir'),
                           ('inner_boundary_value_si', '0'), ('outer_boundary', 'fixed-reservoir'),
                           ('outer_boundary_value_si', '0'), ('reflection_coefficient', '0.6'),
                           ('cascade_coefficient', '0.8'), ('perpendicular_correlation_length_m', '1e7'),
                           ('electron_heating_fraction', '0.3'), ('spectral_k_min_per_m', '1e-10'),
                           ('spectral_k_max_per_m', '1e-2'), ('spectral_bins', '128'), ('cfl_safety', '0.8'),
                           ('operator_accuracy_safety', '0.25'), ('maximum_source_fraction', '0.25'),
                           ('maximum_cascade_fraction', '0.2'), ('minimum_substep_s', '1e-12'),
                           ('maximum_substeps', '1000000'), ('conservation_relative_tolerance', '1e-12'),
                           ('shock_source_efficiency', '0.02'), ('shock_source_plus_fraction', '0.5'),
                           ('diagnostic_cadence_steps', '0')):
            w(f'    {key} = {value}')
        w('  #end turbulence')
        w('#end application')
        w('')
        w('#begin mesh')
        w(f'  global_cell_size_m = {r_outer / 4.0!r}')
        w(f'  minimum_cell_size_m = {r_inner / 4.0!r}')
        w('  maximum_level = 12')
        w('  #begin field_line')
        w('    origin_x_m = 0')
        w('    origin_y_m = 0')
        w('    origin_z_m = 0')
        w('    ! The line through Earth at t0: srcSEP traces the core field from this')
        w('    ! point sunward to inner_radius_m and starts the line there.')
        anchor = [r_earth * x for x in unit]
        w(f'    anchor_x_m = {anchor[0]!r}')
        w(f'    anchor_y_m = {anchor[1]!r}')
        w(f'    anchor_z_m = {anchor[2]!r}')
        w(f'    ! Arc length: PFSS shell depth + exterior spiral to {n["line_end_radius_au"]!r} AU (beyond the outer sphere).')
        w(f'    length_m = {length!r}')
        w(f'    ! Geometric spacing, chord ~ {grid_eps!r} r (first arc {grid_eps!r} r_base, growth ~1+{grid_eps!r}).')
        w(f'    point_count = {grid_points}')
        w(f'    first_segment_length_m = {grid_first!r}')
        w('  #end field_line')
        w('  #begin domain')
        w(f'    inner_radius_m = {r_inner!r}')
        w(f'    outer_radius_m = {r_outer!r}')
        w('  #end domain')
        w('  #begin refinement: solar')
        w('    enabled = true')
        w(f'    surface_cell_size_m = {r_inner / 4.0!r}')
        w(f'    transition_outer_radius_m = {10.0 * RSUN_M!r}')
        w('    profile = smoothstep')
        w('    exponent = 1')
        w('  #end refinement')
        w('  #begin refinement: tube')
        w('    enabled = true')
        w(f'    reference_radius_m = {AU_M!r}')
        w(f'    radius_at_reference_m = {0.04 * AU_M!r}')
        w('    radius_mode = constant-angular-width')
        w(f'    center_cell_size_m = {0.02 * AU_M!r}')
        w('    transverse_profile = smoothstep')
        w('    transverse_exponent = 1')
        w('  #end refinement')
        w('#end mesh')
        w('')
        w('#begin physics: background')
        w('  provider = swcme')
        w('#end physics')
        w('#begin physics: shock')
        w('  authority = swcme')
        w('#end physics')
        w('')
        w('! Required by schema 6 and cross-checked against the line geometry (wind')
        w('! speed, base radius, colatitude).  With --reduced-shock-event the shared')
        w('! core replaces this state at every vertex; values are the measured ones.')
        w('#begin physics: swcme')
        w('  preset = fast')
        w(f'  ambient.wind_speed = {a["wind_speed_km_s"]!r} km/s')
        w(f'  ambient.density_1au = {a["electron_density_1au_cm3"]!r} cm^-3')
        w(f'  ambient.magnetic_field_1au = {a["magnetic_field_1au_nT"]!r} nT')
        w(f'  ambient.proton_temperature = {a["proton_temperature_K"]!r} K')
        w('  ambient.adiabatic_index = 1.6666666666666667')
        w(f'  ambient.alpha_to_proton_ratio = {a["alpha_to_proton_ratio"]!r}')
        w(f'  ambient.electron_temperature = {a["electron_temperature_K"]!r} K')
        w(f'  ambient.alpha_temperature = {a["alpha_temperature_K"]!r} K')
        w('  ambient.thermodynamic_closure = multi_species')
        w(f'  parker.radial_polarity = {a["radial_polarity"]}')
        # Colatitude of the anchor direction: srcSEP compares it with the radial
        # projection of the anchor on the inner sphere.
        w(f'  parker.sin_theta = {math.cos(b0)!r}')
        w(f'  parker.source_radius = {r_inner / RSUN_M!r} Rs')
        w('  cme.kinematics = data_driven')
        w('  cme.data_times = ' + ', '.join(f'{float(k["model_time_s"])!r} s' for k in knots))
        w('  cme.data_radii = ' + ', '.join(f'{k["radius_rsun"]!r} Rs' for k in knots))
        w(f'  cme.launch_radius = {knots[0]["radius_rsun"]!r} Rs')
        first_speed = (knots[1]['radius_rsun'] - knots[0]['radius_rsun']) * RSUN_M / (
            knots[1]['model_time_s'] - knots[0]['model_time_s']) / 1e3
        w(f'  cme.launch_speed = {first_speed!r} km/s')
        w('  cme.drag_coefficient = 0 1/km')
        w('  cme.extrapolation = outside_time')
        w('  shock.region_mode = shock_only')
        w('  shock.acceleration_mode = source')
        w('  shock.relative_source_weight_per_area = 1')
        for key, value in (('geometry.sheath_thickness_1au', '0.12 AU'), ('geometry.ejecta_thickness_1au', '0.22 AU'),
                           ('smoothing.shock_width_1au', '0.010 AU'), ('smoothing.leading_edge_width_1au', '0.020 AU'),
                           ('smoothing.trailing_edge_width_1au', '0.030 AU'), ('sheath.ramp_power', '2'),
                           ('sheath.leading_edge_speed_factor', '1.12'), ('ejecta.density_factor', '0.50'),
                           ('ejecta.speed_factor', '0.80')):
            w(f'  {key} = {value}')
        w('  event.launch_epoch = 0 s')
        w('  event.valid_from = 0 s')
        w(f'  event.valid_until = {end_s!r} s')
        w('#end physics')
        w('')
        w('#begin source: shock')
        w('  enabled = true')
        w(f'  samples_per_step = {samples}')
        w(f'  minimum_energy_j = {s["minimum_energy_mev"] * mev!r}')
        w(f'  maximum_energy_j = {s["maximum_energy_mev"] * mev!r}')
        w('  ! Net release fraction eta: nominal scale; the physical amplitude is the')
        w('  ! single global factor fitted on the ERNE 4 h window (score_observations.pl).')
        w(f'  injection_efficiency = {s["nominal_injection_efficiency"]!r}')
        w('  ! L_ref = f R_sh(t): the upstream release surface is the sphere of radius')
        w('  ! (1+f) R_sh (see README.md, "Release distance").')
        w(f'  release_distance_shock_radius_fraction = {s["release_distance_shock_radius_fraction"]!r}')
        w(f'  spectrum_reference_energy_j = {mev!r}')
        w(f'  spectrum_particle_mass_kg = {MP_KG!r}')
        w('  spectrum_charge_number = 1')
        w('  spectrum_normalization = relative-only')
        w('  spectrum_reference_intensity_per_m2_s_sr_j = 1')
        w('#end source')
        w('')
        w('#begin species')
        w(f'  macroparticle_weight = {n["macroparticle_weight"]!r}')
        w('#end species')
        w('')
        for oid, (r, lo, hi) in sorted(observers_for_deck(settings).items()):
            w(f'#begin observer: {oid}')
            w('  kind = field-line-radius')
            w(f'  heliocentric_radius_m = {r!r}')
            w(f'  minimum_energy_j = {lo * mev!r}')
            w(f'  maximum_energy_j = {hi * mev!r}')
            w('  energy_bins = 1')
            w('  energy_spacing = logarithmic')
            w('  pitch_angle_bins = 24')
            w('#end observer')
        w('! One broadband diagnostic observer at SOHO for the regression spectra.')
        w('#begin observer: soho-spectrum')
        w('  kind = field-line-radius')
        w(f'  heliocentric_radius_m = {g["soho_heliocentric_distance_m"]["mean"]!r}')
        w(f'  minimum_energy_j = {s["minimum_energy_mev"] * mev!r}')
        w(f'  maximum_energy_j = {s["maximum_energy_mev"] * mev!r}')
        w('  energy_bins = 40')
        w('  energy_spacing = logarithmic')
        w('  pitch_angle_bins = 24')
        w('#end observer')
        w('')
        if nightly:
            w('! Inner regression observers (nightly tier only): within the first hour the')
            w('! source particles fill the line out to about 1 AU, so the 1 AU observers above')
            w('! record almost nothing; these give the DiffNum products their signal.')
            w('! They are not compared with observations.')
            for r_au in nt['inner_observer_radii_au']:
                w(f'#begin observer: nightly-r{repr(r_au).replace(".", "p")}au')
                w('  kind = field-line-radius')
                w(f'  heliocentric_radius_m = {r_au * AU_M!r}')
                w(f'  minimum_energy_j = {s["minimum_energy_mev"] * mev!r}')
                w(f'  maximum_energy_j = {s["maximum_energy_mev"] * mev!r}')
                w(f'  energy_bins = {nt["inner_observer_energy_bins"]}')
                w('  energy_spacing = logarithmic')
                w('  pitch_angle_bins = 24')
                w('#end observer')
            w('')
        w('#begin output')
        w('  initialization_mesh_tecplot_file = event-mesh.dat')
        w('  initialization_field_line_tecplot_file = event-field-line.dat')
        w('  initialization_data_tecplot_file = event-native-data.dat')
        w('#end output')
        return lines

    (out / 'event.in').write_text('\n'.join(runtime_deck(False)) + '\n')
    (out / 'event.nightly.in').write_text('\n'.join(runtime_deck(True)) + '\n')

    # Command-line options shared by every mover and both tiers. srcSEP takes the
    # transport law and these numerics only on the command line; the Table recipes
    # read this one line ($$(cat event.args)), so the values have one source here.
    t = settings['transport']
    if settings['decisions']['source_sampling'] != 'log-uniform importance sampling':
        raise ValueError('event_settings.json: unsupported source_sampling decision')
    args = ['--reduced-shock-particles', 'on',
            '--shock-energy-sampling', 'log-uniform',
            '--coefficient-source', 'prescribed',
            '--mean-free-path-provider', 'radial-rigidity-power-law',
            '--mfp-reference-m', repr(t['lambda_reference_au'] * AU_M),
            '--mfp-reference-radius-m', repr(t['lambda_reference_radius_au'] * AU_M),
            '--mfp-reference-rigidity-v', repr(t['lambda_reference_rigidity_v']),
            '--mfp-radial-exponent', repr(t['lambda_radial_exponent']),
            '--mfp-rigidity-exponent', repr(t['lambda_rigidity_exponent']),
            '--spatial-diffusion-provider', 'from-mfp',
            '--pitch-angle-diffusion-provider', 'isotropic-from-mfp',
            '--no-particle-coupling',
            '--turbulence-source', 'prescribed']
    if n['shock_boundary_root_table']:
        args += ['--shock-boundary-root-table', 'on']
    args += ['--shock-boundary-tolerance-m', repr(n['shock_boundary_tolerance_m'])]
    (out / 'event.args').write_text(' '.join(args) + '\n')
    summary = {'wind_isothermal_sound_speed_m_s': c_w, 'wind_critical_radius_rsun': GM_SUN / (2 * c_w * c_w) / RSUN_M,
               'source_surface_radius_rsun': RSS_RS,
               'source_surface_foot_azimuth_deg': math.degrees(phi_foot),
               'source_surface_foot_carrington_longitude_deg': g['earth_carrington_longitude_deg'] + math.degrees(phi_foot),
               'line_base_radius_rsun': r_inner / RSUN_M, 'line_length_m': length,
               'line_point_count': point_count, 'line_first_segment_m': first_arc,
               'nightly_line_point_count': nightly_points, 'nightly_line_first_segment_m': nightly_first,
               'note': 'the 1.1 Rs footpoint is traced through the PFSS shell by srcSEP (REDUCED_ANCHOR line of the run log)'}
    (out / 'deck_summary.json').write_text(json.dumps(summary, indent=1, sort_keys=True) + '\n')
    print(f'wrote {out / "event.event"}, {out / "event.in"}, {out / "event.nightly.in"}, {out / "event.args"}, {asset}')


def check(args):
    """Regenerate every product in a temporary copy and byte-compare."""
    import shutil
    import tempfile
    names = ('derived.json', 'observations.csv', 'event.in', 'event.nightly.in', 'event.args', 'event.event', 'deck_summary.json',
             'assets/wso_cr2135_model_frame.harmonics.csv')
    saved = {name: (HERE / name).read_bytes() for name in names if (HERE / name).exists()}
    try:
        derive(args)
        deck(argparse.Namespace(output_dir=None, raw=args.raw))
        changed = [name for name in names if (HERE / name).read_bytes() != saved.get(name)]
    finally:
        for name, data in saved.items():
            (HERE / name).write_bytes(data)
    if changed:
        raise ValueError('regenerated products differ from the checked-in files: ' + ', '.join(changed))
    print('PASS: every checked-in input is reproducible: ' + ', '.join(names))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='action', required=True)
    d = sub.add_parser('derive', help='raw archives -> derived.json, observations.csv')
    d.add_argument('--raw', required=True, type=Path)
    k = sub.add_parser('deck', help='derived.json + WSO map -> event.event, event.in, asset (standard library only)')
    k.add_argument('--raw', required=True, type=Path)
    k.add_argument('--output-dir', default=None)
    c = sub.add_parser('check', help='regenerate and byte-compare all products')
    c.add_argument('--raw', required=True, type=Path)
    args = parser.parse_args()
    try:
        {'derive': derive, 'deck': deck, 'check': check}[args.action](args)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f'prepare_event: {error}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
