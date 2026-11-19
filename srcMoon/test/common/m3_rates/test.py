#!/usr/bin/env python3
"""Production-backed supporting D04/CEX checks; never full product V&V.

Compile only the PIC-independent native sources. Do not configure AMPS or touch
its generated tree. Source-corruption fixtures live exclusively in the owned
output directory; canonical data and original research PDFs remain read-only.
Every invocation/response and source identity is retained, including failures.
"""
import argparse
from collections import Counter
import csv
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import time
import unittest

import primary_reference as reference

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0,str(REPO/'srcMoon/test'))
from chemical_database_build import shared_sources, shared_dependencies
CONTEXT = None
SENTINEL = 17.125


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


class RecordedResult(unittest.TextTestResult):
    """Retain actual method/fixture outcomes, including missing-data groups.

    unittest does not count a skipped/errored setUpClass as a started method.
    Counting records explicitly avoids negative PASS totals or pretending that
    missing-source methods ran. A fixture outcome is identified as such.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.records = []

    def addSuccess(self, test):
        super().addSuccess(test)
        self.records.append({"test": test.id(), "status": "PASS"})

    def addFailure(self, test, error):
        super().addFailure(test, error)
        self.records.append({"test": test.id(), "status": "FAIL"})

    def addError(self, test, error):
        super().addError(test, error)
        self.records.append({"test": test.id(), "status": "ERROR"})

    def addSkip(self, test, reason):
        super().addSkip(test, reason)
        self.records.append({"test": test.id(), "status": "SKIPPED", "reason": reason})


class Context:
    def __init__(self, args):
        self.args = args
        self.out = args.output_dir.resolve()
        self.out.mkdir(parents=True, exist_ok=False)
        self.photo = args.data_root / "photochemistry/processed/phidrates"
        self.capture = args.cex_reference_dir
        self.calls = []
        self.comparisons = []
        self.probe = self.out / "probe"

    def compile(self):
        sources = [HERE / "probe.cpp", REPO / "srcMoon/species/Photochemistry.cpp",
                   REPO / "srcMoon/species/ChargeExchange.cpp"]
        sources += shared_sources(REPO,('photochemistry','charge-exchange'))
        command = [self.args.cxx, "-std=c++17", "-O2", "-g", "-Wall", "-Wextra",
                   "-Werror", "-pedantic", "-I", str(REPO / "srcMoon/species")]
        command += [str(p) for p in sources] + ["-o", str(self.probe)]
        with (self.out / "compile.log").open("w") as log:
            done = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                  check=False, timeout=180)
        metadata = {
            "command": command, "returncode": done.returncode,
            "compiler": subprocess.check_output([self.args.cxx, "--version"], text=True).splitlines()[0],
            "source_sha256": {str(p): sha(p) for p in sources + [
                REPO / "srcMoon/species/Photochemistry.h",
                REPO / "srcMoon/species/ChargeExchange.h", Path(__file__),
                HERE / "primary_reference.py", HERE / "CONTRACT.md"]+
                shared_dependencies(REPO,('photochemistry','charge-exchange'))},
        }
        if done.returncode == 0:
            metadata["executable_sha256"] = sha(self.probe)
        (self.out / "build.json").write_text(json.dumps(metadata, indent=2) + "\n")
        if done.returncode:
            raise RuntimeError("native rate-probe compilation failed; see compile.log")

    def batch(self, mode, lines):
        directory = self.photo if mode == "photo" else self.capture
        command = [str(self.probe), mode, str(directory)]
        done = subprocess.run(command, input="\n".join(lines) + "\n", text=True,
                              capture_output=True, timeout=120, check=False)
        index = len(self.calls)
        record = {"command": command, "stdin": lines, "returncode": done.returncode,
                  "stdout": done.stdout, "stderr": done.stderr}
        self.calls.append(record)
        (self.out / ("call-%03d.json" % index)).write_text(json.dumps(record, indent=2) + "\n")
        if done.returncode:
            raise RuntimeError("production probe failed; retained call-%03d.json" % index)
        rows = [json.loads(line) for line in done.stdout.splitlines()]
        if len(rows) != len(lines):
            raise ValueError("native adapter response count differs from requested cases")
        return rows

    def mutation(self, name, content):
        directory = self.out / "source-corruption-fixtures"
        directory.mkdir(exist_ok=True)
        path = directory / name
        with path.open("xb") as stream:
            stream.write(content)
        return path


class NativeCase(unittest.TestCase):
    def close(self, actual, expected, tolerance=1e-12):
        expected = float(expected)
        self.assertTrue(math.isfinite(actual))
        self.assertLessEqual(abs(actual / expected - 1), tolerance)

    def refused(self, row):
        self.assertFalse(row["ok"])
        self.assertEqual(row["value"], SENTINEL)
        self.assertTrue(row["error"])


class HashIdentity(NativeCase):
    def test_streaming_sha256_matches_published_and_independent_vectors(self):
        # Published SHA vectors and an independent hashlib evaluation test
        # empty input, padding boundaries, multiple blocks and a long stream.
        for index, raw in enumerate((b"", b"abc", b"a" * 55, b"a" * 56,
                                     bytes(range(256)) * 1024, b"a" * 1000000)):
            path = CONTEXT.mutation("hash-%d.bin" % index, raw)
            done = subprocess.run([str(CONTEXT.probe), "sha", str(path)],
                                  text=True, capture_output=True, check=False)
            record = {"command": [str(CONTEXT.probe), "sha", str(path)],
                      "returncode": done.returncode, "stdout": done.stdout,
                      "stderr": done.stderr}
            index = len(CONTEXT.calls)
            CONTEXT.calls.append(record)
            (CONTEXT.out / ("call-%03d.json" % index)).write_text(json.dumps(record, indent=2) + "\n")
            self.assertEqual(done.returncode, 0, done.stderr)
            value = json.loads(done.stdout)
            self.assertTrue(value["ok"])
            self.assertEqual(value["digest"], hashlib.sha256(raw).hexdigest())


class PhotonRegistry(NativeCase):
    @classmethod
    def setUpClass(cls):
        required = ("reaction_channels.csv", "photo_rates_1au.csv",
                    "thresholds.csv", "cross_sections.csv")
        missing = [str(CONTEXT.photo / n) for n in required if not (CONTEXT.photo / n).is_file()]
        if missing:
            raise unittest.SkipTest("authoritative D04 products unavailable: " + repr(missing))
        cls.catalogue = read_csv(CONTEXT.photo / required[0])
        cls.rates = read_csv(CONTEXT.photo / required[1])
        cls.thresholds = read_csv(CONTEXT.photo / required[2])
        cls.cross = read_csv(CONTEXT.photo / required[3])

    def test_complete_catalogue_endpoint_and_source_identity(self):
        value = CONTEXT.batch("photo", ["inventory"])[0]
        self.assertEqual((value["catalogue"], value["reactive"]), (17, 15))
        actual = {row["id"]: row for row in value["channels"]}
        self.assertEqual(set(actual), {row["reaction_id"] for row in self.catalogue})
        for row in self.catalogue:
            native = actual[row["reaction_id"]]
            self.assertEqual(native["key"], row["species"])
            self.assertEqual(native["parent"], row["parent"])
            self.assertEqual(native["source_products"], row["products"])
            self.assertEqual(native["source_label"], row["phidrates_channel_label"])
            self.assertEqual(native["source_reference"], row["source_reference"])
            self.assertEqual(native["source_file"], row["source_file"])
            self.assertEqual(native["parent_state"], "not_supplied_by_source")
            self.assertFalse(native["daughter_kinetics"])
            self.assertEqual(native["experimental_oh"], row["species"] == "OH_exp")
            cross = [r for r in self.cross if r["reaction_id"] == row["reaction_id"]]
            self.assertEqual(native["sample_count"], len(cross))
            self.assertEqual(native["wavelength_nm"],
                             [min(float(r["wavelength_nm"]) for r in cross),
                              max(float(r["wavelength_nm"]) for r in cross)])
        for row in self.rates:
            native = actual[row["reaction_id"]]
            index = int(float(row["solar_activity"]))
            self.assertEqual(native["rates"][index], float(row["rate_1au_s-1"]))
            self.assertEqual(native["mean_energy_eV"][index], float(row["excess_energy_eV"]))
        for row in self.thresholds:
            native = actual[row["reaction_id"]]
            self.assertTrue(native["has_threshold"])
            self.assertEqual(native["threshold_eV"], float(row["threshold_eV"]))
        for row in actual.values():
            if row["key"] in ("H2O", "OH_exp"):
                self.assertFalse(row["has_threshold"])
        # Preserve the real unexplained OH tail; loading a catalogue must not
        # silently clip/reinterpret it as a newly qualified threshold model.
        self.assertEqual(actual["PHIDRATES:OH_exp:02:o_1s_h"]["wavelength_nm"][1], 6132)
        (CONTEXT.out / "photon_inventory.json").write_text(json.dumps(value, indent=2) + "\n")

    def test_products_states_multiplicities_and_exact_nucleus_charge_identities(self):
        actual = CONTEXT.batch("photo", ["inventory"])[0]["channels"]
        for row in actual:
            if not row["reactive"]:
                self.assertFalse(row["products"])
                continue
            expected = Counter()
            total = Counter()
            charge = 0
            for token in row["source_products"].split(" + "):
                state_match = re.fullmatch(r"(.+)\(([^)]+)\)", token)
                chemical, state = ((state_match[1], state_match[2]) if state_match
                                   else (token, "not_supplied_by_source"))
                q = -1 if chemical == "e" else int(chemical.endswith("+"))
                expected[(chemical, state, q)] += 1
            native = Counter((p["species"], p["state"], p["charge"])
                             for p in row["products"] for unused in range(p["multiplicity"]))
            self.assertEqual(native, expected)
            for product in row["products"]:
                for nucleus, count in product["nuclei"].items():
                    total[nucleus] += count * product["multiplicity"]
                charge += product["charge"] * product["multiplicity"]
            parent = Counter()
            for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", row["parent"]):
                parent[element] += int(count or 1)
            self.assertEqual(total, parent)
            self.assertEqual(charge, 0)

    def test_all_family_rates_branches_selector_and_inverse_square_visibility(self):
        calls, expected = [], []
        for key, count in (("He", 1), ("Na", 1), ("Ne", 1), ("Ar", 1), ("H2O", 7), ("OH_exp", 4)):
            for solar in (0, 1):
                rows = [r for r in self.rates if r["species"] == key and float(r["solar_activity"]) == solar]
                self.assertEqual(len(rows), count)
                total = sum((Decimal(r["rate_1au_s-1"]) for r in rows), Decimal(0))
                calls += ["rate %s %d 1 1" % (key, solar),
                          "rate %s %d 2 0.25" % (key, solar),
                          "rate %s %d 1 0" % (key, solar)]
                expected += [total, total / 16, Decimal(0)]
        for row, truth in zip(CONTEXT.batch("photo", calls), expected):
            self.assertTrue(row["ok"], row)
            if truth:
                self.close(row["value"], truth)
            else:
                self.assertEqual(row["value"], 0)
        for key in ("H2O", "OH_exp"):
            for solar in (0, 1):
                rows = [r for r in self.rates if r["species"] == key and float(r["solar_activity"]) == solar]
                total = math.fsum(float(r["rate_1au_s-1"]) for r in rows)
                branches = dict(CONTEXT.batch("photo", ["branches %s %d" % (key, solar)])[0]["branches"])
                self.assertEqual(set(branches), {r["reaction_id"] for r in rows})
                for r in rows:
                    self.close(branches[r["reaction_id"]], float(r["rate_1au_s-1"]) / total)
                self.assertAlmostEqual(math.fsum(branches.values()), 1, places=15)
                # Endpoint and interval-interior controls use the native source
                # row order independently of the C++ selector's internals.
                # This exact binary64 predecessor of 1 supports Python 3.8,
                # which predates math.nextafter; it still exercises u<1 at
                # the final representable inverse-CDF endpoint.
                uniforms = [0.0, 1.0 - 2.0 ** -53]
                identities = [rows[0]["reaction_id"], rows[-1]["reaction_id"]]
                cumulative = 0.0
                for row in rows:
                    hazard = float(row["rate_1au_s-1"])
                    uniforms.append((cumulative + hazard / 2) / total)
                    identities.append(row["reaction_id"])
                    cumulative += hazard
                selected = CONTEXT.batch("photo", ["select %s %d %.17g" % (key, solar, u)
                                                    for u in uniforms])
                for selection, identity in zip(selected, identities):
                    self.assertTrue(selection["ok"])
                    self.assertEqual(selection["id"], identity)
                # Deterministic quadrature of the piecewise-constant inverse
                # CDF: count error <=1 is the exact grid-rounding identity,
                # not a stochastic confidence or observational threshold.
                n = 100000
                counts = CONTEXT.batch("photo", ["grid %s %d %d" % (key, solar, n)])[0]["counts"]
                self.assertEqual(sum(counts.values()), n)
                for identity, probability in branches.items():
                    self.assertLessEqual(abs(counts.get(identity, 0) - n * probability), 1.000000001)

    def test_invalid_domain_total_absorption_and_unselected_experimental_oh_are_refused(self):
        lines = ["rate OH 0 1 1", "rate unknown 0 1 0", "rate Na 2 1 1",
                 "rate Na 0 0 1", "rate Na 0 -1 1", "rate Na 0 nan 1",
                 "rate Na 0 1 -0.1", "rate Na 0 1 1.01", "rate Na 0 1 inf",
                 "rate Na 0 1e-300 1", "channel-rate PHIDRATES:H2O:00:total 0"]
        for row in CONTEXT.batch("photo", lines):
            self.refused(row)
        # A successful query after a refusal clears the prior diagnostic;
        # callers must not mistake an old error string for the new outcome.
        recovery = CONTEXT.batch("photo", ["rate Na 9 1 1", "rate Na 0 1 1"])
        self.refused(recovery[0])
        self.assertTrue(recovery[1]["ok"])
        self.assertEqual(recovery[1]["error"], "")
        for row in CONTEXT.batch("photo", ["select H2O 0 -0.1", "select H2O 0 1",
                                           "select H2O 0 nan", "select H2O 9 0.5"]):
            self.assertFalse(row["ok"])
            self.assertEqual(row["id"], "unchanged")

    def test_valid_looking_source_changes_mixed_files_and_self_approved_hash_are_refused_transactionally(self):
        raw = (CONTEXT.photo / "photo_rates_1au.csv").read_bytes()
        changed = raw.replace(b"5.64000000000000E-8", b"5.65000000000000E-8", 1)
        self.assertNotEqual(raw, changed)
        self.assertEqual(len(raw), len(changed))
        path = CONTEXT.mutation("valid changed rates.csv", changed)
        deleted = CONTEXT.mutation("missing-rate-row.csv", b"\n".join(raw.split(b"\n")[:-2]) + b"\n")
        lines = ["reload rates %s - -" % json.dumps(str(path)),
                 "reload rates %s %s -" % (json.dumps(str(path)), sha(path)),
                 "reload rates %s - -" % json.dumps(str(deleted)),
                 "reload channels %s - -" % json.dumps(str(CONTEXT.photo / "photo_rates_1au.csv")),
                 "reload rates %s - srcMoon-D04-phidrates-v2" % json.dumps(str(CONTEXT.photo / "photo_rates_1au.csv"))]
        for row in CONTEXT.batch("photo", lines):
            self.assertFalse(row["ok"])
            self.assertEqual(row["value"], 7.26e-6)  # Retained previously loaded source.
            self.assertTrue(row["error"])


class TotalCapture(NativeCase):
    @classmethod
    def setUpClass(cls):
        missing = [str(CONTEXT.capture / r["filename"]) for r in reference.PUBLICATIONS.values()
                   if not (CONTEXT.capture / r["filename"]).is_file()]
        if missing:
            raise unittest.SkipTest("original CEX publications unavailable: " + repr(missing))
        reference.write_reference(CONTEXT.out, CONTEXT.capture)

    def test_original_he_fit_matches_independent_decimal_polynomials(self):
        energies = [float(e) for e, unused in reference.HE_PRINTED_TABLE_CM2]
        energies += [99 * (11000000 / 99) ** (i / 60) for i in range(1, 60)]
        responses = CONTEXT.batch("cex", ["query 0 %.17g" % e for e in energies])
        for energy, row in zip(energies, responses):
            self.assertTrue(row["ok"], row)
            truth = reference.helium_fit_m2(format(energy, ".17g"))
            self.close(row["value"], truth)
            CONTEXT.comparisons.append({"channel": "He_original_fit", "energy_eV_u": energy,
                                        "production_m2": row["value"], "independent_m2": str(truth),
                                        "relative_error": abs(row["value"] / float(truth) - 1)})
        # Source fit-vs-table residuals are retained as DIAGNOSTICS. Its printed
        # 10%RMS/18.7%maximum deviation is not replaced with a fictitious1%
        # experimental-agreement criterion or called observational validation.
        with (CONTEXT.out / "helium_original_fit_vs_printed_table.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(("energy_eV_u", "printed_table_cm2", "source_fit_m2", "relative_fit_table_difference"))
            for energy, printed in reference.HE_PRINTED_TABLE_CM2:
                fitted = reference.helium_fit_m2(energy)
                writer.writerow((energy, printed, str(fitted), str(fitted / (Decimal(printed) * Decimal("1e-4")) - 1)))

    def test_sodium_exact_all_state_knots_and_named_interpolation(self):
        native = CONTEXT.batch("cex", ["knots"])[0]["knots"]
        expected = [(int(e) * 1000, float(Decimal(s) * Decimal("1e-20")))
                    for e, s in reference.NA_PRINTED_TABLE_KEV_U]
        self.assertEqual(native, [list(x) for x in expected])
        energies = [e for e, unused in expected]
        energies += [math.sqrt(a[0] * b[0]) for a, b in zip(expected, expected[1:])]
        rows = CONTEXT.batch("cex", ["query 1 %.17g" % e for e in energies])
        for energy, row in zip(energies, rows):
            self.assertTrue(row["ok"])
            truth = reference.sodium_loglog_m2(format(energy, ".17g"))
            self.close(row["value"], truth)
            CONTEXT.comparisons.append({"channel": "Na_named_loglog", "energy_eV_u": energy,
                                        "production_m2": row["value"], "independent_m2": str(truth),
                                        "relative_error": abs(row["value"] / float(truth) - 1)})

    def test_metadata_domain_and_cold_beam_density_units(self):
        metadata = CONTEXT.batch("cex", ["metadata 0", "metadata 1"])
        for row, expected in zip(metadata, ([99, 11000000], [1000, 200000])):
            self.assertTrue(row["ok"])
            self.assertEqual(row["domain"], expected)
            self.assertFalse(row["daughter_kinetics"])
            self.assertFalse(row["particle_execution"])
            self.assertEqual(row["energy_convention"], "target-rest kinetic energy per projectile u")
        # 1000eV/u,4.39e7cm/s are an explicit rounded E/g pair in ORNL A-32;
        # the density here defines a mathematical fixture, never measured wind.
        rows = CONTEXT.batch("cex", ["beam 0 1000 439000 1000000",
                                     "beam 0 1000 439000 2000000", "beam 0 1000 439000 0"])
        truth = reference.helium_fit_m2(1000) * Decimal(439000) * 1000000
        self.close(rows[0]["value"], truth)
        self.close(rows[1]["value"], 2 * truth)
        self.assertEqual(rows[2]["value"], 0)

    def test_no_extrapolation_unsupported_channels_or_invalid_baths(self):
        lines = ["query 0 98.999", "query 0 11000001", "query 0 nan",
                 "query 1 999.999", "query 1 200001", "query 1 829",
                 "query 1 inf", "query 2 1000", "query -1 1000",
                 "beam 1 1000 439000 -1", "beam 1 1000 0 1",
                 "beam 1 1000 inf 1", "beam 1 999 439000 0"]
        for row in CONTEXT.batch("cex", lines):
            self.refused(row)

    def test_changed_publication_or_self_approved_metadata_cannot_replace_loaded_law(self):
        path = CONTEXT.capture / reference.PUBLICATIONS["na"]["filename"]
        altered = CONTEXT.mutation("changed original.pdf", path.read_bytes() + b"\n")
        lines = ["reload na %s - -" % json.dumps(str(altered)),
                 "reload na %s %s -" % (json.dumps(str(altered)), sha(altered)),
                 "reload he %s - -" % json.dumps(str(path)),
                 "reload na %s - unreviewed-v2" % json.dumps(str(path))]
        for row in CONTEXT.batch("cex", lines):
            self.assertFalse(row["ok"])
            self.assertEqual(row["value"], 43e-20)
            self.assertTrue(row["error"])


def main():
    global CONTEXT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--data-root", type=Path, default=Path("/data/vtenishe/moon_validation_data"))
    parser.add_argument("--cex-reference-dir", type=Path,
                        default=REPO / "test_output/srcMoon/m3-charge-exchange-research-20261117")
    parser.add_argument("--cxx", default="/usr/bin/g++")
    args = parser.parse_args()
    CONTEXT = Context(args)
    started = time.monotonic()
    try:
        CONTEXT.compile()
        suite = unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(c)
                                   for c in (HashIdentity, PhotonRegistry, TotalCapture))
        with (CONTEXT.out / "checks.log").open("w") as log:
            result = unittest.TextTestRunner(stream=log, verbosity=2, resultclass=RecordedResult).run(suite)
        counts = Counter(r["status"] for r in result.records)
        failed, errors, skipped = counts["FAIL"], counts["ERROR"], counts["SKIPPED"]
        status = "ERROR" if errors else ("FAIL" if failed else ("SKIPPED" if skipped else "PASS"))
        report = {"status": status, "scope": "supporting PIC-independent production rate/registry verification only",
                  "executed": counts["PASS"] + failed + errors, "passed": counts["PASS"],
                  "failed": failed, "errors": errors, "skipped": skipped,
                  "unit_test_methods_started": result.testsRun, "checks": result.records,
                  "formal_full_product_status": {"U17": "SKIPPED", "I22": "SKIPPED", "I24": "SKIPPED"},
                  "reason": "Daughter kinetics/production activation/formal predecessors are separate gates",
                  "checks_log": str(CONTEXT.out / "checks.log"),
                  "native_call_count": len(CONTEXT.calls), "elapsed_s": time.monotonic() - started}
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        report = {"status": "ERROR", "scope": "supporting rate/registry test infrastructure",
                  "executed": 0, "passed": 0, "failed": 0, "skipped": 0, "errors": 1,
                  "reason": str(error), "elapsed_s": time.monotonic() - started}
    (CONTEXT.out / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    (CONTEXT.out / "independent_comparisons.json").write_text(json.dumps(CONTEXT.comparisons, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return {"PASS": 0, "FAIL": 1, "ERROR": 2, "SKIPPED": 77}[report["status"]]


if __name__ == "__main__":
    raise SystemExit(main())
