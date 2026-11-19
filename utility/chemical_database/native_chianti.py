"""Lossless CHIANTI level/radiative source access, independent of AMPS models.

Format authority: CHIANTI technical reports 1 (current, version-8+ elvlc) and 10
(wgfa), https://www.chiantidatabase.org/tech_reports/. This module preserves
the native tokens/row identities and never supplies a missing measurement,
physical state/mass, uncertainty, emissivity model, or application selection.

Levels use cm^-1, with observed and theoretical columns kept separate. A native
-1 energy is unavailable; it is never silently replaced by the other column.
Signed wavelengths in Angstrom retain the CHIANTI provenance convention:
negative means at least one endpoint is not experimentally known. Zero is a
special two-photon/autoionization-style entry, not a single emitted photon.
gf is dimensionless; the decay-rate column is s^-1. No uncertainty is present
in these fields, so precision is not promoted to an error bar. The numerical
column may include nonradiative/special rates; consumers must check its kind.

Permissions, native hashes and original theory credits are authenticated by
the source-packet/import layer. A lexical Catalogue is not physical readiness.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import hashlib
import io
import json
import re

_INTEGER = re.compile(r"[+-]?[0-9]+\Z")
_NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eEdD][+-]?[0-9]+)?\Z")


def _integer(text, label):
    token = text.strip()
    if not _INTEGER.fullmatch(token):
        raise ValueError("invalid native CHIANTI " + label + ": " + repr(token))
    return int(token)


def _number(text, label):
    token = text.strip()
    if not _NUMBER.fullmatch(token):
        raise ValueError("invalid native CHIANTI " + label + ": " + repr(token))
    try:
        value = Decimal(token.replace("D", "E").replace("d", "e"))
    except InvalidOperation as error:
        raise ValueError("invalid native CHIANTI " + label) from error
    if not value.is_finite():
        raise ValueError("nonfinite native CHIANTI " + label)
    return token, value


def _block(native, label):
    """Decode exactly one native data block; retain unrestricted footer text.

    CHIANTI defines everything after the first -1 as comments, including a
    possible final -1. We retain those bytes/text and do not reinterpret a
    numeric-looking footer line as scientific data. A missing first sentinel is
    a truncated/unsupported input, not a silently accepted EOF.
    """
    try:
        rows = native.splitlines(keepends=True)
        lines = [row.decode("ascii").rstrip("\r\n") for row in rows]
    except UnicodeError as error:
        raise ValueError("native CHIANTI " + label + " is not ASCII") from error
    boundary = next((i for i, line in enumerate(lines) if line.strip() == "-1"), None)
    if boundary is None:
        raise ValueError("missing native CHIANTI " + label + " data sentinel")
    if boundary == 0:
        raise ValueError("empty native CHIANTI " + label + " data block")
    return rows, lines, boundary


@dataclass(frozen=True)
class Catalogue:
    """Owned immutable snapshot, with defensive nested copies for all queries."""
    _snapshot: str

    def as_dict(self):
        return json.loads(self._snapshot)

    @property
    def levels(self):
        return tuple(self.as_dict()["levels"])

    @property
    def transitions(self):
        return tuple(self.as_dict()["transitions"])


def parse_ion_bytes(levels, radiative, *, maximum_bytes=None):
    """Parse current CHIANTI elvlc/wgfa bytes transactionally and without fits.

    Neither file is modified. Any malformed layout, duplicate identity, invalid
    reference or nonfinite number refuses the entire returned object. Optional
    maximum_bytes limits the combined input bytes; it never truncates data.
    Native energy-order/missing-energy defects are explicit source-quality
    issues, not repaired values. Indices may be out of energy order, as allowed
    by the CHIANTI format; numeric lower/upper *indices* are not energy ordering.
    """
    if not isinstance(levels, bytes) or not isinstance(radiative, bytes):
        raise ValueError("native CHIANTI inputs must be immutable bytes")
    if maximum_bytes is not None:
        if isinstance(maximum_bytes, bool) or not isinstance(maximum_bytes, int) or maximum_bytes <= 0:
            raise ValueError("maximum_bytes must be an explicit positive integer")
        if len(levels)+len(radiative) > maximum_bytes:
            raise ValueError("native CHIANTI bytes exceed maximum_bytes; no truncation performed")
    level_raw, level_lines, level_end = _block(levels, "elvlc")
    rate_raw, rate_lines, rate_end = _block(radiative, "wgfa")
    entries, by_index, issues = [], {}, []
    for i, text in enumerate(level_lines[:level_end]):
        if len(text) != 87:
            raise ValueError("unsupported elvlc fixed-column width at line " + str(i+1))
        index = _integer(text[:7], "level index")
        if index <= 0 or index in by_index:
            raise ValueError("invalid/duplicate CHIANTI level index")
        multiplicity = _integer(text[42:47], "spin multiplicity")
        j_token, j = _number(text[52:57], "J")
        observed, observed_value = _number(text[57:72], "observed energy")
        theory, theory_value = _number(text[72:87], "theoretical energy")
        if multiplicity <= 0 or j < 0 or 2*j != (2*j).to_integral_value():
            raise ValueError("invalid native multiplicity or non-half-integer J")
        if any(value < 0 and value != -1 for value in (observed_value, theory_value)):
            raise ValueError("unsupported negative energy; only native -1 means unavailable")
        row = {"level_index": index, "configuration_native": text[7:37].strip(),
               "label_native": text[37:42].strip(), "spin_multiplicity_native": text[42:47].strip(),
               "L_native": text[47:52].strip(), "J_native": j_token,
               "observed_energy_cm_inverse_native": observed,
               "observed_energy_status": "UNAVAILABLE_NATIVE_MINUS_ONE" if observed_value == -1 else "REPORTED",
               "theoretical_energy_cm_inverse_native": theory,
               "theoretical_energy_status": "UNAVAILABLE_NATIVE_MINUS_ONE" if theory_value == -1 else "REPORTED",
               "uncertainty_status": "NOT_SUPPLIED_IN_FILE", "native_line": i+1,
               "native_row": text, "native_row_sha256": hashlib.sha256(level_raw[i]).hexdigest()}
        if not row["configuration_native"] or not row["L_native"]:
            raise ValueError("empty native CHIANTI configuration or L field")
        if observed_value == theory_value == -1:
            issues.append({"code": "BOTH_LEVEL_ENERGIES_UNAVAILABLE", "level_index": index,
                           "native_line": i+1, "disposition": "preserve gap; no energy assigned"})
        entries.append(row)
        by_index[index] = row
    links, pairs = [], set()
    for i, text in enumerate(rate_lines[:rate_end]):
        if len(text) < 55:
            raise ValueError("truncated wgfa fixed-column row at line " + str(i+1))
        lower, upper = _integer(text[:5], "lower index"), _integer(text[5:10], "upper index")
        if lower == upper or lower not in by_index or upper not in by_index or (lower, upper) in pairs:
            raise ValueError("invalid/duplicate/unresolved CHIANTI transition endpoints")
        pairs.add((lower, upper))
        wavelength, wave = _number(text[10:25], "wavelength")
        gf_token, gf = _number(text[25:40], "gf")
        rate_token, rate = _number(text[40:55], "decay rate")
        if gf < 0 or rate < 0:
            raise ValueError("negative native gf/decay rate")
        row = {"lower_index": lower, "upper_index": upper,
               "signed_wavelength_angstrom_native": wavelength,
               "wavelength_status": ("OBSERVED_ENDPOINTS" if wave > 0 else
                                     "THEORETICAL_OR_UNOBSERVED_ENDPOINT" if wave < 0 else
                                     "NO_SINGLE_PHOTON_WAVELENGTH_NATIVE_ZERO"),
               "gf_native": gf_token, "decay_rate_s_inverse_native": rate_token,
               "rate_kind": "SINGLE_PHOTON_RADIATIVE" if wave else "UNRESOLVED_TWO_PHOTON_OR_AUTOIONIZATION",
               "uncertainty_status": "NOT_SUPPLIED_IN_FILE", "label_native": text[55:].strip(),
               "native_line": i+1, "native_row": text,
               "native_row_sha256": hashlib.sha256(rate_raw[i]).hexdigest()}
        # Compare only two genuinely reported observed energies. No theoretical
        # substitution or repaired wavelength/ordering is generated by this QA.
        lo, hi = by_index[lower], by_index[upper]
        if lo["observed_energy_status"] == hi["observed_energy_status"] == "REPORTED":
            if Decimal(lo["observed_energy_cm_inverse_native"]) >= Decimal(hi["observed_energy_cm_inverse_native"]):
                issues.append({"code": "NONINCREASING_OBSERVED_TRANSITION_ENERGY", "native_line": i+1,
                               "lower_index": lower, "upper_index": upper,
                               "disposition": "preserve native transition; no wavelength/level repair"})
        links.append(row)
    if not links:
        raise ValueError("empty native CHIANTI radiative block")
    value = {"schema": "amps-chianti-native-lexical-index-v1", "format": "CHIANTI8+-elvlc-and-wgfa",
             "level_source_sha256": hashlib.sha256(levels).hexdigest(),
             "radiative_source_sha256": hashlib.sha256(radiative).hexdigest(),
             "level_comments": level_lines[level_end+1:], "radiative_comments": rate_lines[rate_end+1:],
             "levels": entries, "transitions": links, "native_quality_issues": issues,
             "native_domain_quality_status": "FAIL" if issues else "PASS",
             "scientific_status": "SOURCE_INPUT_NOT_PHYSICAL_QUALIFICATION"}
    return Catalogue(json.dumps(value, sort_keys=True))


def _selection(catalogue, indices):
    if not isinstance(catalogue, Catalogue) or not isinstance(indices, (list, tuple)) or not indices:
        raise ValueError("explicit Catalogue and nonempty level-index sequence required")
    if any(type(i) is not int or i <= 0 for i in indices) or len(set(indices)) != len(indices):
        raise ValueError("explicit distinct positive integer level indices required")
    by_index = {row["level_index"]: row for row in catalogue.levels}
    if any(i not in by_index for i in indices):
        raise LookupError("requested CHIANTI level is unavailable")
    return by_index


def levels_csv(catalogue, level_indices):
    """Copy lexical level quantities in explicit caller order, with no fallback."""
    by_index = _selection(catalogue, level_indices)
    columns = ("level_index", "configuration_native", "label_native", "spin_multiplicity_native", "L_native", "J_native",
               "observed_energy_cm_inverse_native", "observed_energy_status", "theoretical_energy_cm_inverse_native",
               "theoretical_energy_status", "uncertainty_status", "native_line", "native_row_sha256")
    stream = io.StringIO(newline="")
    writer = csv.writer(stream); writer.writerow(columns)
    for index in level_indices:
        writer.writerow(by_index[index][name] for name in columns)
    return stream.getvalue()


def transitions_csv(catalogue, level_indices):
    """Project a downward-closed level set, preserving every native weak link.

    If a selected upper state has a known outgoing link to an omitted lower
    level, refusing prevents a partial source model being mistaken for a
    complete cascade. No rate/branch fraction is normalized by this operation.
    Source transition order is preserved, independent of requested level order.
    """
    _selection(catalogue, level_indices)
    selected = set(level_indices)
    if any(row["upper_index"] in selected and row["lower_index"] not in selected for row in catalogue.transitions):
        raise ValueError("selected CHIANTI levels omit a known downward endpoint; no implicit channel deletion")
    columns = ("lower_index", "upper_index", "signed_wavelength_angstrom_native", "wavelength_status", "gf_native",
               "decay_rate_s_inverse_native", "rate_kind", "uncertainty_status", "label_native", "native_line", "native_row_sha256")
    stream = io.StringIO(newline="")
    writer = csv.writer(stream); writer.writerow(columns)
    for row in catalogue.transitions:
        if row["upper_index"] in selected and row["lower_index"] in selected:
            writer.writerow(row[name] for name in columns)
    return stream.getvalue()
