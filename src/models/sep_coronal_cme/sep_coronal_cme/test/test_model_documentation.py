#!/usr/bin/env python3
"""DOCSCCM01 tests for the generated shared-model specification.

The positive test proves that generation is byte deterministic.  Each
negative test copies the small documentation package to a temporary directory,
changes one independently meaningful contract, and executes the copied
generator.  Testing the command in an isolated copy is important: it exercises
the same path discovery and committed-file comparison used in a clean checkout
without ever modifying the developer's source tree.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


MODEL_ROOT = Path(__file__).resolve().parents[1]
GENERATOR = MODEL_ROOT / "tools" / "generate_model.py"


def write_text_lf(path: Path, content: str) -> None:
    """Write deterministic UTF-8/LF text on every supported Python version.

    ``Path.write_text`` did not gain its ``newline`` keyword until after
    Python 3.8.  NASA HEC systems still provide Python 3.8, while
    ``Path.open`` has accepted the corresponding argument for much longer.
    Keeping this compatibility shim in the test prevents a negative fixture
    from failing before the documentation validator receives the mutation.
    """

    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(content)


class DocumentationGenerationTests(unittest.TestCase):
    """Exercise the complete model-neutral ``DOCSCCM01`` release gate."""

    maxDiff = 4000

    def run_generator(
            self, root: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
        """Run a package-local generator and capture one diagnostic stream."""

        return subprocess.run(
            [sys.executable, str(root / "tools" / "generate_model.py"), *arguments],
            cwd=root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )

    def isolated_copy(self, temporary_root: Path) -> Path:
        """Copy only this shared model package for a destructive negative case."""

        destination = temporary_root / "sep_coronal_cme"
        shutil.copytree(
            MODEL_ROOT,
            destination,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
        )
        return destination

    def assert_rejected(
            self, mutate, expected_fragment: str) -> None:  # type: ignore[no-untyped-def]
        """Apply one mutation and require an informative nonzero check result."""

        with tempfile.TemporaryDirectory(prefix="docsccm01-negative-") as directory:
            root = self.isolated_copy(Path(directory))
            mutate(root)
            completed = self.run_generator(root, "--check")
            self.assertNotEqual(completed.returncode, 0, completed.stdout)
            self.assertIn(expected_fragment, completed.stdout)

    def test_clean_check_and_output_are_byte_identical(self) -> None:
        """A clean registry/modules snapshot must reproduce committed bytes."""

        checked = self.run_generator(MODEL_ROOT, "--check")
        self.assertEqual(checked.returncode, 0, checked.stdout)
        self.assertIn("DOCSCCM01 PASS", checked.stdout)
        with tempfile.TemporaryDirectory(prefix="docsccm01-output-") as directory:
            output = Path(directory) / "model.generated.md"
            generated = self.run_generator(MODEL_ROOT, "--output", str(output))
            self.assertEqual(generated.returncode, 0, generated.stdout)
            self.assertEqual(output.read_bytes(), (MODEL_ROOT / "model.md").read_bytes())

    def test_dirty_canonical_document_is_rejected(self) -> None:
        """Direct edits to generated ``model.md`` must never become authority."""

        def mutate(root: Path) -> None:
            path = root / "model.md"
            write_text_lf(path, path.read_text(encoding="utf-8") + "direct edit\n")

        self.assert_rejected(mutate, "not generated-clean")

    def test_duplicate_numbered_section_is_rejected(self) -> None:
        """Moving or copying normative sections across owners must be detected."""

        def mutate(root: Path) -> None:
            path = root / "model" / "architecture_exchange.md"
            write_text_lf(path, 
                path.read_text(encoding="utf-8")
                + "\n## 3. Illegally duplicated ownership\n\nnegative fixture\n")

        self.assert_rejected(mutate, "duplicate numbered sections")

    def test_duplicate_requirement_id_is_rejected(self) -> None:
        """The structured trace registry cannot assign one ID twice."""

        def mutate(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["requirements"][1]["id"] = payload["requirements"][0]["id"]
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(mutate, "duplicate requirement ID")

    def test_requirement_status_is_a_closed_schema_value(self) -> None:
        """A status typo must not pass merely because it is nonempty text."""

        def mutate(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["requirements"][0]["status"] = "specified-ish"
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(mutate, "is not a schema-1 status")

    def test_review_sources_are_typed_primary_and_reciprocal(self) -> None:
        """Review provenance must be canonical and agree with disposition links."""

        def malformed(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["requirements"][0]["review_sources"] = ["review-two:R1"]
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(malformed, "review_sources has invalid values")

        def missing_primary(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["requirements"][0]["review_sources"] = ["review-3:N2"]
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(missing_primary, "must include 'review-2:R1'")

        def stale_backlink(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["requirements"][0]["review_sources"].append("review-3:N2")
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(stale_backlink, "does not link back to the requirement")

    def test_unknown_requirement_and_review_members_are_rejected(self) -> None:
        """Schema-version 1 records cannot silently absorb misspelled fields."""

        def requirement_member(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["requirements"][0]["statuz"] = "specified"
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(requirement_member, "unknown members ['statuz']")

        def review_member(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["review_findings"][0]["requirement_id"] = []
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(review_member, "unknown members ['requirement_id']")

    def test_duplicate_requirement_and_review_links_are_rejected(self) -> None:
        """Duplicate array entries cannot overstate trace coverage."""

        def requirement_link(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            links = payload["requirements"][0]["config_keys"]
            links.append(links[0])
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(requirement_link, "contains duplicate values")

        def review_link(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            finding = next(
                item for item in payload["review_findings"] if item["id"] == "R1")
            finding["requirement_ids"].append(finding["requirement_ids"][0])
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(review_link, "contains duplicate values")

    def test_requirement_id_series_must_match_review_item(self) -> None:
        """An R1 record cannot carry a syntactically valid R2-prefixed ID."""

        def mutate(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["requirements"][0]["id"] = "SCCM-R2-SCS-RADIALIZATION"
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(mutate, "ID series does not match review_item")

    def test_latest_review_series_cannot_be_omitted(self) -> None:
        """Every accepted P-series review item must remain in the registry."""

        def mutate(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["requirements"] = [
                item for item in payload["requirements"]
                if item["review_item"] != "P8"
            ]
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(
            mutate,
            "requirements must contain SCCM R1--R10, N1--N7, and P1--P8",
        )

    def test_unresolved_test_link_is_rejected(self) -> None:
        """A requirement may reference only a canonical Section-17 definition."""

        def mutate(root: Path) -> None:
            path = root / "model" / "requirements.yaml"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["requirements"][0]["test_ids"] = ["NO-SUCH-TEST99"]
            write_text_lf(path, json.dumps(payload, indent=2) + "\n")

        self.assert_rejected(mutate, "references unknown canonical test")

    def test_test_definition_outside_section_17_2_is_rejected(self) -> None:
        """A definition-shaped bullet elsewhere in the module is not canonical."""

        def mutate(root: Path) -> None:
            path = root / "model" / "testing_validation.md"
            text = path.read_text(encoding="utf-8")
            text = text.replace("- `SCS3D09`:", "- **SCS3D09**:", 1)
            text = text.replace(
                "## 17. Testing and validation campaign",
                "- `SCS3D09`: misplaced definition outside Section 17.2.\n\n"
                "## 17. Testing and validation campaign",
                1,
            )
            write_text_lf(path, text)

        self.assert_rejected(mutate, "roadmap references unknown canonical tests")

    def test_roadmap_range_must_expand_to_canonical_tests(self) -> None:
        """Every member of a compact roadmap range must be registered."""

        def mutate(root: Path) -> None:
            path = root / "model" / "testing_validation.md"
            text = path.read_text(encoding="utf-8").replace(
                "| 1 | `PFSS3D01--09` |",
                "| 1 | `PFSS3D01--99` |",
                1,
            )
            write_text_lf(path, text)

        self.assert_rejected(mutate, "roadmap references unknown canonical tests")

    def test_roadmap_requires_every_canonical_test(self) -> None:
        """Removing a test from its stage row must leave a detectable gap."""

        def mutate(root: Path) -> None:
            path = root / "model" / "testing_validation.md"
            text = path.read_text(encoding="utf-8").replace(
                ", `OFX3D01`, `LOS3D01`",
                ", `LOS3D01`",
                1,
            )
            write_text_lf(path, text)

        self.assert_rejected(mutate, "canonical tests have no roadmap stage")

    def test_stage_13_retains_the_prior_stage_release_dependency(self) -> None:
        """The aggregate release stage cannot drop its Stage-0--12 dependency."""

        def mutate(root: Path) -> None:
            path = root / "model" / "testing_validation.md"
            text = path.read_text(encoding="utf-8").replace(
                " plus every mandatory record from Stages 0--12",
                "",
                1,
            )
            write_text_lf(path, text)

        self.assert_rejected(
            mutate,
            "Stage 13 must include every mandatory Stage 0--12 record",
        )

    def test_cross_stage_test_reuse_must_be_explicit(self) -> None:
        """The one deliberate Stage-14 reuse cannot become an implicit duplicate."""

        def mutate(root: Path) -> None:
            path = root / "model" / "testing_validation.md"
            text = path.read_text(encoding="utf-8").replace(
                "Stage 14A reuses `SRC3D19`",
                "Stage 14A includes `SRC3D19`",
                1,
            )
            write_text_lf(path, text)

        self.assert_rejected(mutate, "must be explicit")

    def test_unknown_generated_placeholder_is_rejected(self) -> None:
        """Misspelled or private generated regions cannot vanish silently."""

        def mutate(root: Path) -> None:
            path = root / "model" / "physics.md"
            write_text_lf(path, path.read_text(encoding="utf-8")
                            + "{{GENERATED:UNREGISTERED_TABLE}}\n")

        self.assert_rejected(mutate, "unknown generated placeholders")

    def test_noncanonical_newlines_and_missing_final_newline_are_rejected(self) -> None:
        """Platform newline conversion cannot create a second canonical byte form."""

        def crlf(root: Path) -> None:
            path = root / "model" / "configuration_validation.md"
            path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))

        self.assert_rejected(crlf, "is not LF-only")

        def no_final_lf(root: Path) -> None:
            path = root / "model" / "testing_validation.md"
            path.write_bytes(path.read_bytes().rstrip(b"\n"))

        self.assert_rejected(no_final_lf, "must end with one LF newline")


if __name__ == "__main__":
    unittest.main(verbosity=2)
