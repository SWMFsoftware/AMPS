"""Lossless lexical access to the official NASA CEA thermochemical source.

This module indexes source coefficients; it does not evaluate a fit, provide
an isotope/rest-mass profile, infer a ground/excited state, or approve a
reaction threshold. Source authentication and operation-specific permissions
belong to source_packets.py. The caller supplies native bytes here.

Field authority: NASA/TP-2002-211556 Appendix C and the acquired official
nasa/cea database_compile.f90 FORMAT 99002/99003/99005/99006. The current
compiler reads A15 names and skips coefficient columns33--48 with16X;
nonblank text in that skipped slot must be preserved, not mistaken for a
new polynomial coefficient. Raw line/block identities retain that detail.

All quantities remain source strings: temperature K, source molecular weight
u, assigned heat of formation at298.15K J/mol, and H(298.15)-H(0) J/mol.
A zero final thermal-increment field is an explicit unavailable marker per
Appendix C. Assigned H(0) generally differs from heat of formation Hf(0).
Electron-reference conventions and source elemental-weight conventions must
be resolved before any physical state/reaction use. No uncertainty is inferred
from decimal precision. Native domain defects are retained as quality issues;
lexical access does not turn an invalid domain into an evaluable fit.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import hashlib
import io
import json


def _number(text):
    token = text.strip()
    try:
        value = Decimal(token.replace("D", "E").replace("d", "e"))
    except InvalidOperation as error:
        raise ValueError("invalid native CEA number: " + repr(token)) from error
    if not value.is_finite():
        raise ValueError("nonfinite native CEA number: " + repr(token))
    return token


@dataclass(frozen=True)
class Catalogue:
    """Immutable owned snapshot; every mapping returned is a defensive copy.

    Keeping the canonical JSON string immutable prevents a caller changing a
    returned coefficient/list from retargeting a later query. It also permits
    the original caller's input bytes/local dictionary to be released safely.
    """
    _snapshot: str

    def as_dict(self):
        return json.loads(self._snapshot)

    @property
    def records(self):
        return tuple(self.as_dict()["records"])

    @property
    def source_sha256(self):
        return self.as_dict()["source_sha256"]


def parse_thermo_bytes(native, *, maximum_bytes=None):
    """Index a complete native file, preserving every original numeric token.

    Malformed/truncated structure or nonfinite numbers raise ValueError before
    a Catalogue is returned. Optional maximum_bytes is an explicit resource
    admission control, never a physical cutoff/truncation. Out-of-order source
    temperature limits remain visible quality issues and are never repaired.
    Native reactant-only zero-interval records and duplicate names are retained.
    """
    if not isinstance(native, bytes):
        raise ValueError("native CEA input must be immutable bytes")
    if maximum_bytes is not None:
        if isinstance(maximum_bytes, bool) or not isinstance(maximum_bytes, int) or maximum_bytes <= 0:
            raise ValueError("maximum_bytes must be an explicit positive integer")
        if len(native) > maximum_bytes:
            raise ValueError("native CEA input exceeds explicit maximum_bytes; no truncation performed")
    raw = native.splitlines(keepends=True)
    try:
        lines = [row.decode("ascii").rstrip("\r\n") for row in raw]
    except UnicodeError as error:
        raise ValueError("native CEA source is not ASCII") from error
    positions = [i for i, row in enumerate(lines) if row.strip().lower() == "thermo"]
    if len(positions) != 1:
        raise ValueError("expected exactly one native thermo sentinel")
    first = positions[0]
    if any(row.strip() and not row.startswith("!") for row in lines[:first]):
        raise ValueError("unexplained source data before thermo sentinel")
    if first + 1 >= len(lines):
        raise ValueError("missing native global schedule/date")
    schedule = lines[first + 1].ljust(40)
    temperatures = [_number(schedule[offset:offset+10]) for offset in range(0, 40, 10)]
    i, section = first + 2, "products"
    records, issues = [], []
    ended_products = ended_reactants = False
    while i < len(lines):
        text = lines[i]
        if text.strip() == "END PRODUCTS":
            if ended_products or section != "products":
                raise ValueError("duplicate/out-of-order END PRODUCTS")
            ended_products, section, i = True, "reactants", i + 1
            continue
        if text.strip() == "END REACTANTS":
            if not ended_products:
                raise ValueError("END REACTANTS precedes END PRODUCTS")
            ended_reactants, i = True, i + 1
            break
        if not text.strip() or text.startswith("!"):
            i += 1
            continue
        start = i
        if i + 1 >= len(lines):
            raise ValueError("truncated native record metadata")
        metadata = lines[i+1].ljust(80)
        try:
            intervals, phase = int(metadata[:2]), int(metadata[50:52])
        except ValueError as error:
            raise ValueError("invalid native interval/phase at line " + str(i+2)) from error
        if not 0 <= intervals <= 99:
            raise ValueError("native two-column interval count is invalid")
        name = text[:15].strip()  # Actual pinned compiler's A15, not an inferred chemical identity.
        if not name:
            raise ValueError("empty native species name")
        record = {"name": name, "section": section, "phase": phase,
            "native_start_line": start+1, "native_header": text,
            "reference": text[15:].strip(), "reference_date_code": metadata[3:9].strip(),
            "composition": [], "interval_count": intervals,
            "molecular_weight_native_u": _number(metadata[52:65]),
            "hf298_native_J_mol": _number(metadata[65:80]), "intervals": []}
        for offset in range(10, 50, 8):
            symbol, amount = metadata[offset:offset+2].strip(), metadata[offset+2:offset+8]
            if symbol:
                record["composition"].append({"symbol": symbol, "count_native": _number(amount)})
            elif amount.strip() and Decimal(_number(amount)) != 0:
                raise ValueError("unnamed nonzero native composition entry")
        i += 2
        if not intervals:
            if i >= len(lines):
                raise ValueError("truncated reactant-only reference row")
            record["reactant_reference_temperature_K"] = _number(lines[i].ljust(80)[:11])
            record["reactant_native_reference_row"] = lines[i]
            i += 1
        for _ in range(intervals):
            if i+2 >= len(lines):
                raise ValueError("truncated native coefficient interval")
            row, left, right = [line.ljust(80) for line in lines[i:i+3]]
            if row[22:23] != "7":
                raise ValueError("unsupported native Cp coefficient count")
            lo, hi = _number(row[:11]), _number(row[11:22])
            exponents = [_number(row[k:k+5]) for k in range(23, 63, 5)]
            if [Decimal(v) for v in exponents] != [Decimal(v) for v in (-2,-1,0,1,2,3,4,0)]:
                raise ValueError("unsupported native CEA polynomial exponents")
            coefficients = [_number(left[k:k+16]) for k in range(0, 80, 16)]
            coefficients += [_number(right[k:k+16]) for k in (0,16,48,64)]
            delta = _number(row[65:80])
            ordered = Decimal(lo) < Decimal(hi)
            if not ordered:
                issues.append({"code": "NONINCREASING_NATIVE_TEMPERATURE_LIMITS",
                    "name": name, "native_line": i+1, "minimum_K": lo, "maximum_K": hi,
                    "disposition": "preserve source values; no evaluator/physical qualification supplied"})
            record["intervals"].append({"native_start_line": i+1,
                "temperature_min_K": lo, "temperature_max_K": hi,
                "temperature_limits_increasing": ordered, "Cp_exponents_native": exponents,
                "coefficients_a1_to_a7_b1_b2_native": coefficients,
                "native_compiler_ignored_columns_33_to_48": right[32:48],
                "h298_minus_h0_native_J_mol": delta,
                "h298_minus_h0_status": "UNAVAILABLE_NATIVE_ZERO_MARKER" if Decimal(delta)==0 else "REPORTED"})
            i += 3
        record["native_end_line"] = i
        record["native_block_sha256"] = hashlib.sha256(b"".join(raw[start:i])).hexdigest()
        records.append(record)
    if not ended_products or not ended_reactants or not records:
        raise ValueError("missing native end sentinel or species records")
    if any(row.strip() and not row.startswith("!") for row in lines[i:]):
        raise ValueError("unparsed trailing native records")
    data = {"schema": "amps-nasa-cea-native-lexical-index-v1",
        "source_sha256": hashlib.sha256(native).hexdigest(), "native_line_count": len(lines),
        "native_global_schedule_and_date": lines[first+1],
        "native_global_temperatures_K": temperatures,
        "native_global_date": schedule[40:].strip(), "records": records,
        "native_domain_quality_status": "FAIL" if issues else "PASS", "native_quality_issues": issues,
        "scientific_status": "SOURCE_INPUT_NOT_PHYSICAL_QUALIFICATION"}
    return Catalogue(json.dumps(data, sort_keys=True))


def selected_gas_record(catalogue, name):
    """Select exactly one original products/gas record; never merge duplicates."""
    if not isinstance(catalogue, Catalogue) or not isinstance(name, str) or not name:
        raise ValueError("explicit Catalogue and nonempty native species name required")
    matches = [row for row in catalogue.records if row["section"]=="products" and row["phase"]==0 and row["name"]==name]
    if len(matches) != 1:
        raise LookupError("missing/ambiguous native products gas species: " + name)
    return matches[0]


def selected_species_csv(catalogue, names):
    """Return native strings in caller-declared order, with source row identity.

    This projection copies no new physics: no H(0) conversion, charge inference,
    isotope remapping, numerical fit, or inferred uncertainty is included.
    The caller must preserve its recipe, source SHA, licence and NOTICE.
    """
    if not isinstance(names, (list, tuple)) or not names or len(set(names)) != len(names):
        raise ValueError("explicit distinct source species names required")
    stream = io.StringIO(newline="")
    writer = csv.writer(stream)
    writer.writerow(("name", "native_start_line", "native_end_line", "phase", "reference_date_code",
        "molecular_weight_native_u", "hf298_native_J_mol", "h298_minus_h0_native_J_mol", "source_reference"))
    for name in names:
        row = selected_gas_record(catalogue, name)
        increments = {item["h298_minus_h0_native_J_mol"] for item in row["intervals"]}
        if len(increments) != 1:
            raise ValueError("native thermal-increment field differs/missing across selected intervals")
        writer.writerow((row["name"], row["native_start_line"], row["native_end_line"], row["phase"],
            row["reference_date_code"], row["molecular_weight_native_u"], row["hf298_native_J_mol"],
            next(iter(increments)), row["reference"]))
    return stream.getvalue()


def coefficients_csv(catalogue, names):
    """Project actual nine-coefficient arrays with native domains and row IDs.

    Fortran-D strings are retained verbatim, so generic decimal-only evaluators
    must refuse them rather than silently choosing a number conversion. The
    polynomial's source equations determine the coefficient units; this helper
    does not calculate Cp, enthalpy, entropy, a reaction energy or a state mass.
    A source-domain FAIL remains explicit in each row and is never repaired.
    """
    if not isinstance(names, (list, tuple)) or not names or len(set(names)) != len(names):
        raise ValueError("explicit distinct source species names required")
    stream = io.StringIO(newline="")
    writer = csv.writer(stream)
    writer.writerow(("name", "native_species_start_line", "native_interval_start_line", "interval_index",
        "temperature_min_K", "temperature_max_K", "native_domain_quality",
        "a1_native", "a2_native", "a3_native", "a4_native", "a5_native", "a6_native", "a7_native", "b1_native", "b2_native",
        "h298_minus_h0_native_J_mol", "h298_minus_h0_status"))
    for name in names:
        record = selected_gas_record(catalogue, name)
        for index, interval in enumerate(record["intervals"]):
            writer.writerow((record["name"], record["native_start_line"], interval["native_start_line"], index,
                interval["temperature_min_K"], interval["temperature_max_K"],
                "PASS" if interval["temperature_limits_increasing"] else "FAIL",
                *interval["coefficients_a1_to_a7_b1_b2_native"],
                interval["h298_minus_h0_native_J_mol"], interval["h298_minus_h0_status"]))
    return stream.getvalue()
