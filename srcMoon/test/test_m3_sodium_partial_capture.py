#!/usr/bin/env python3
"""Check production P07 partial-capture rates/products against independent SI math.

Real native cross sections and independently acquired isotope/constants data
are evidence inputs. The forward-angle closure is explicitly experimental;
these checks do not turn that assumed direction into a measurement or a
complete U17/I24 qualification. No application configuration is changed.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from chemical_database_build import shared_sources, shared_dependencies

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REFERENCE = ROOT / "test_output/srcMoon/m3-atomic-online-databases-20261118T083400Z"
DEFAULT_INVENTORY = ROOT / "test_output/srcMoon/m3-approved-products-20261118T163142Z/source-acquisition-v2/reference_mass_inventory.json"
SOURCE_CSV_SHA256 = "2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93"
DEPENDENCIES = ("srcMoon/species/SodiumPartialCapture.h", "srcMoon/species/SodiumPartialCapture.cpp",
    "srcMoon/species/ReducedKinematics.h", "srcMoon/species/ReducedKinematics.cpp",
    "srcMoon/species/Photochemistry.h", "srcMoon/species/Photochemistry.cpp",
    "srcMoon/test/sodium_partial_capture/probe.cpp", "srcMoon/test/sodium_partial_capture/README.md",
    "srcMoon/test/sodium_partial_capture/prepare_source.py", "srcMoon/test/physical_state_sources/prepare.py",
    "srcMoon/test/test_m3_sodium_partial_capture.py")
DEPENDENCIES += tuple(str(p.relative_to(ROOT)) for p in shared_dependencies(ROOT,'photochemistry'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_inputs(reference, inventory_path):
    """Verify source-bound research records before any numerical comparison."""
    reference, inventory_path = Path(reference), Path(inventory_path)
    native = reference / "collisiondb_119990_embedded_arrays.csv"
    if not native.is_file() or not inventory_path.is_file():
        raise FileNotFoundError("explicit D119990 native CSV or NIST reference-mass inventory unavailable")
    if sha(native) != SOURCE_CSV_SHA256:
        raise ValueError("D119990 native partial-channel bytes changed")
    records = json.loads((reference / "native_inventory.json").read_text())
    source = [r for r in records if r.get("dataset_id") == 119990 and r.get("output") == native.name]
    if len(source) != 1 or source[0]["output_sha256"] != SOURCE_CSV_SHA256:
        raise ValueError("native IAEA source provenance is ambiguous or changed")
    source_html = reference / source[0]["native_file"]
    if sha(source_html) != source[0]["native_sha256"]:
        raise ValueError("retained native IAEA response does not match its manifest")
    extractor = ROOT / "srcMoon/test/sodium_partial_capture/prepare_source.py"
    spec = importlib.util.spec_from_file_location("moon_partial_source_extractor", extractor)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    native_pairs = module.extract_native_rows(source_html)
    parsed_pairs = [(float(r["energy_eV_per_u"]), float(r["cross_section_cm2"])) for r in csv.DictReader(native.open())]
    if parsed_pairs != native_pairs:
        raise ValueError("native partial-channel CSV differs from maintained source extraction")
    inventory = json.loads(inventory_path.read_text())
    packet = inventory_path.parent
    manifest = json.loads((packet / "run_manifest.json").read_text())
    extractor = ROOT / "srcMoon/test/physical_state_sources/prepare.py"
    if manifest.get("script_sha256") != sha(extractor):
        raise ValueError("mass inventory extraction script differs from its source receipt; reproduce offline with the maintained script")
    spec = importlib.util.spec_from_file_location("moon_reference_source_extractor", extractor)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    identities = []
    for row in manifest["sources"]:
        file = Path(row["path"])
        if sha(file) != row["sha256"]:
            raise ValueError("NIST raw source bytes changed: " + str(file))
        identities.append({"path": str(file), "sha256": sha(file)})
    # Bind DERIVED inventory values as well as source bytes. Otherwise editing
    # one derived mass/constant would fool both native and Decimal consumers.
    if module.extract(packet / "raw") != inventory:
        raise ValueError("derived NIST inventory disagrees with reproducible native-source extraction")
    source_result = json.loads((packet / "result.json").read_text())
    if source_result.get("status") != "PASS" or inventory.get("schema") != "srcMoon-isotope-reference-mass-evidence-v1":
        raise ValueError("reference source packet has no successful extraction/schema receipt")
    isotopes = {r["identity"]: r for r in inventory["isotopes"]}
    constants = inventory["constants"]
    ground = {r["element"]: r for r in inventory["atomic_ground_thresholds"]}
    # Validate dimensional/source roles instead of interpreting printed digits
    # as an uncertainty. Missing isotope-shift qualification stays explicit.
    assert constants["electron mass"]["unit"] == "kg"
    assert constants["atomic mass constant"]["unit"] == "kg"
    assert constants["elementary charge"]["unit"] == "C"
    return native, inventory, isotopes, constants, ground, identities


def run(output, reference=DEFAULT_REFERENCE, mass_inventory=DEFAULT_INVENTORY, compiler="g++"):
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("output must be new or empty")
    output.mkdir(parents=True, exist_ok=True)
    result = {"schema": "srcMoon-P07-partial-capture-controls-v1", "status": "ERROR",
        "scope": "source-bounded experimental partial-CEX component, not full particle chemistry",
        "formal_M3_status": "SKIPPED / NOT QUALIFIED", "results": [], "counts": {}, "commands": []}
    manifest = {"schema": "srcMoon-P07-partial-capture-manifest-v1", "source_revision":
        subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_hashes": {}, "MPI": "none; production pure-kernel checks", "random_seed": None,
        "physical_approximation": "forward COM direction, source-native partial channel, explicit target-rest beam",
        "limitations": ["not an outgoing angular measurement", "not an all-state total rate",
            "isotope shifts/threshold uncertainty not qualified", "no PIC particle/mover validation"]}
    try:
        resolved = shutil.which(compiler)
        if not resolved:
            raise FileNotFoundError("C++17 compiler unavailable: " + compiler)
        native, inventory, isotopes, constants, ground, identities = source_inputs(reference, mass_inventory)
        manifest["source_hashes"] = {name: sha(ROOT / name) for name in DEPENDENCIES}
        manifest["external_data"] = identities + [{"path": str(native), "sha256": sha(native)},
            {"path": str(mass_inventory), "sha256": sha(mass_inventory)}]
        # An edited DERIVED number must be refused even if all original raw
        # source hashes remain intact. This is an explicit corruption control,
        # never a replacement scientific product or observation.
        corrupt = output / "derived-inventory-corruption-control"
        corrupt.mkdir()
        packet = Path(mass_inventory).resolve().parent
        for filename in ("run_manifest.json", "result.json"):
            shutil.copy2(packet / filename, corrupt / filename)
        (corrupt / "raw").symlink_to(packet / "raw", target_is_directory=True)
        changed_inventory = json.loads(Path(mass_inventory).read_text())
        changed_inventory["isotopes"][0]["reference_neutral_mass_kg"] = "9.999e-27"
        (corrupt / "reference_mass_inventory.json").write_text(json.dumps(changed_inventory))
        refused = False
        try:
            source_inputs(reference, corrupt / "reference_mass_inventory.json")
        except ValueError as error:
            refused = "derived NIST inventory disagrees" in str(error)
        result["results"].append({"name": "derived_inventory_mutation_rejected", "status": "PASS" if refused else "FAIL"})
        executable = output / "partial-capture-probe"
        command = [resolved, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror", "-pedantic", "-I", str(ROOT),
            str(ROOT / "srcMoon/species/SodiumPartialCapture.cpp"), str(ROOT / "srcMoon/species/ReducedKinematics.cpp"),
            str(ROOT / "srcMoon/species/Photochemistry.cpp"), str(ROOT / "srcMoon/test/sodium_partial_capture/probe.cpp"),
            *(str(p) for p in shared_sources(ROOT,'photochemistry')), "-o", str(executable)]
        with (output / "compile.log").open("w") as log:
            built = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, timeout=120)
        result["commands"].append({"argv": command, "returncode": built.returncode, "log": str(output / "compile.log")})
        if built.returncode:
            raise RuntimeError("production partial-capture compilation failed; see compile.log")
        values = [isotopes["1H"]["reference_neutral_mass_kg"], isotopes["23Na"]["reference_neutral_mass_kg"],
            constants["electron mass"]["value"], constants["elementary charge"]["value"], constants["atomic mass constant"]["value"],
            ground["H"]["ionization_energy_eV"], ground["Na"]["ionization_energy_eV"]]
        command = [str(executable), str(native), *values, str(output / "comparison.csv")]
        completed = subprocess.run(command, cwd=output, capture_output=True, text=True, timeout=60)
        (output / "stdout.log").write_text(completed.stdout)
        (output / "stderr.log").write_text(completed.stderr)
        result["commands"].append({"argv": command, "returncode": completed.returncode})
        rows = re.findall(r"^(PASS|FAIL) (.+)$", completed.stdout, re.M)
        summary = re.findall(r"^SUMMARY checks=(\d+) failed=(\d+)$", completed.stdout, re.M)
        if len(summary) != 1 or not rows:
            raise RuntimeError("missing/malformed native component receipt")
        total, failed = map(int, summary[0])
        if total != len(rows) or failed != sum(s == "FAIL" for s, _ in rows) or completed.returncode != (1 if failed else 0):
            raise RuntimeError("native component rows/summary/exit contradict each other")
        required = Counter({name: 1 for name in (
            "load_exact_published_partial_source", "all_twelve_native_nodes_retained",
            "initial_target_state_remains_unqualified",
            "nonnode_gap_refused_unchanged", "low_energy_extrapolation_refused",
            "high_energy_extrapolation_refused", "negative_density_refused", "zero_density_exact_zero_rate",
            "missing_constant_source_refused", "failed_reload_preserves_authenticated_registry",
            "untransformed_frame_label_refused", "changed_constant_value_refused", "changed_constant_identity_refused",
            "wrong_outgoing_state_refused", "wrong_projectile_isotope_refused",
            "undefined_forward_direction_refused", "null_output_refused")})
        required.update({name: 12 for name in ("native_node_rate", "source_bound_forward_products",
            "partial_model_label_preserved", "forged_rate_source_refused_transactionally")})
        if Counter(name for _, name in rows) != required:
            raise RuntimeError("native component omitted/repeated/changed a frozen required control")
        names = Counter()
        for status, name in rows:
            names[name] += 1
            result["results"].append({"name": name + ":" + str(names[name]), "status": status})
        # Decimal reference arithmetic is independent of the compiled solver.
        # It verifies source conversions and full product COM/energy accounting,
        # while retaining the reduced model's experimental physics label.
        comparisons = list(csv.DictReader((output / "comparison.csv").open()))
        native_rows = list(csv.DictReader(native.open()))
        if len(comparisons) != 12 or len(native_rows) != 12:
            raise RuntimeError("incomplete source/node product output")
        with localcontext() as context:
            context.prec = 60
            ev = Decimal(constants["elementary charge"]["value"])
            mass_unit = Decimal(constants["atomic mass constant"]["value"])
            q = (Decimal(ground["H"]["ionization_energy_eV"]) - Decimal(ground["Na"]["ionization_energy_eV"])) * ev
            for row, source in zip(comparisons, native_rows):
                energy = Decimal(source["energy_eV_per_u"])
                sigma = Decimal(source["cross_section_cm2"]) * Decimal("1e-4")
                speed = (2 * energy * ev / mass_unit).sqrt()
                expected = {"energy_eV_per_u": energy, "cross_section_m2": sigma,
                    "relative_speed_m_s": speed, "rate_s-1": Decimal("2.5e6") * sigma * speed, "Q_J": q}
                for field, oracle in expected.items():
                    observed = Decimal(row[field])
                    ok = abs(observed - oracle) <= abs(oracle) * Decimal("1e-12")
                    result["results"].append({"name": str(energy) + ":" + field, "status": "PASS" if ok else "FAIL"})
                for field, oracle in (("px_out", Decimal(row["px_in"])), ("kinetic_out_J", Decimal(row["kinetic_in_J"]) + q)):
                    ok = abs(Decimal(row[field]) - oracle) <= abs(oracle) * Decimal("1e-12")
                    result["results"].append({"name": str(energy) + ":" + field, "status": "PASS" if ok else "FAIL"})
                # Check the small chemical gain on its OWN scale; a large
                # incident beam energy must not hide a lost Q in subtraction.
                gain_residual = (Decimal(row["kinetic_out_J"]) - Decimal(row["kinetic_in_J"]) - q) / abs(q)
                row["relative_Q_closure_residual"] = str(gain_residual)
                result["results"].append({"name": str(energy) + ":chemical_Q_gain", "status": "PASS" if abs(gain_residual) <= Decimal("1e-12") else "FAIL"})
        result["counts"] = {s: sum(r["status"] == s for r in result["results"]) for s in ("PASS", "FAIL", "ERROR", "SKIPPED")}
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        figure, axes = plt.subplots(1, 2, figsize=(10, 4))
        energy = [float(r["energy_eV_per_u"]) / 1000 for r in comparisons]
        axes[0].plot(energy, [float(r["cross_section_cm2"]) * 1e-4 for r in native_rows], "o", label="Native measured partial channel")
        axes[0].plot(energy, [float(r["cross_section_m2"]) for r in comparisons], "+", label="Production parser")
        axes[0].set(xlabel="Projectile energy (keV/u)", ylabel="Partial cross section (m²)")
        axes[1].plot(energy, [float(r["relative_Q_closure_residual"]) for r in comparisons], "o", label="Production gain minus source Q")
        axes[1].axhline(0, color="black", label="Exact closure")
        axes[1].set(xlabel="Projectile energy (keV/u)", ylabel="Relative chemical-energy closure residual")
        for axis in axes:
            axis.legend(fontsize=8, framealpha=1)
        figure.suptitle("P07 partial channel — experimental forward closure")
        figure.tight_layout()
        for suffix in ("png", "eps"):
            figure.savefig(output / ("comparison." + suffix), dpi=150)
        plt.close(figure)
        if {name: sha(ROOT / name) for name in DEPENDENCIES} != manifest["source_hashes"]:
            raise RuntimeError("component sources changed during verification")
        manifest["executable"] = {"path": str(executable), "sha256": sha(executable)}
        result["counts"] = {s: sum(r["status"] == s for r in result["results"]) for s in ("PASS", "FAIL", "ERROR", "SKIPPED")}
        result["status"] = "FAIL" if result["counts"]["FAIL"] else "PASS"
    except (FileNotFoundError, ModuleNotFoundError) as error:
        # Missing plotting/source dependencies must not hide an already
        # evaluated mismatch in available native/independent controls.
        result.update(status="FAIL" if any(r["status"] == "FAIL" for r in result["results"]) else "SKIPPED", reason=str(error))
    except (OSError, ValueError, RuntimeError, ArithmeticError, AssertionError, KeyError, TypeError, subprocess.SubprocessError) as error:
        result.update(status="ERROR", reason=str(error))
    finally:
        result["counts"] = {s: sum(r["status"] == s for r in result["results"]) for s in ("PASS", "FAIL", "ERROR", "SKIPPED")}
        (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "test_output/srcMoon" / ("m3-na-partial-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")))
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--mass-inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--cxx", default="g++")
    args = parser.parse_args()
    try:
        result = run(args.output_dir, args.reference_dir, args.mass_inventory, args.cxx)
    except (OSError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2))
    return {"PASS": 0, "FAIL": 1, "ERROR": 2, "SKIPPED": 77}[result["status"]]


if __name__ == "__main__":
    raise SystemExit(main())
