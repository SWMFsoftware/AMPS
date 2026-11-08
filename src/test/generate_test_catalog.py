#!/usr/bin/env python3
"""Generate the human-readable catalog of every registered AMPS core case.

The manifest and out-of-line case registries remain execution authorities.
This program turns that exact expanded registry into documentation and enriches
case rows with descriptions already maintained beside C++/Python fixtures,
model-local individual READMEs, JSON campaigns, and Markdown requirement
tables. When no explicit prose exists, it produces a conservative label from
the stable case ID rather than inventing a physical claim.

The output is deterministic: it contains no generation timestamp or host path.
``--check`` is therefore suitable for harness/CI drift detection. New tests
must be registered first, then the catalog must be regenerated in the same
change.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


TEST_ROOT = Path(__file__).resolve().parent
REPO_ROOT = TEST_ROOT.parents[1]
CATALOG_PATH = TEST_ROOT / "TEST_CATALOG.md"
sys.path.insert(0, str(TEST_ROOT))

from lib.manifest import ManifestError, load_manifest  # noqa: E402


CASE_ID = re.compile(r"^[A-Z][A-Z0-9-]*$")
SOURCE_SUFFIXES = {".cc", ".cpp", ".cxx", ".h", ".hpp", ".inc", ".json", ".md", ".py"}
SKIP_PARTS = {
    ".git", "__pycache__", "build", "output", "test_output", "stage1",
    # Large empirical/model data are not catalog prose authorities. Their
    # integrity/provenance remain owned by the individual test fixtures.
    "assets", "figures",
}


@dataclass(frozen=True)
class Description:
    """One prose candidate with deterministic precedence and provenance."""

    text: str
    source: str
    priority: int
    order: int


def normalize_id(value: str) -> str:
    """Normalize model-local punctuation exactly as global adapters do."""

    return re.sub(r"-+", "-", re.sub(r"[^A-Za-z0-9]+", "-", value)).strip("-").upper()


def one_line(value: str, limit: int = 240) -> str:
    """Return compact table prose without silently producing giant rows."""

    text = " ".join(str(value).replace("|", " ").split())
    if len(text) <= limit:
        return text
    return text[:limit - 3].rstrip() + "..."


def markdown_cell(value: str) -> str:
    return one_line(value).replace("|", "\\|")


def repository_link(path: str, label: str | None = None) -> str:
    """Make a link relative to src/test/ so the catalog remains relocatable."""

    target = (REPO_ROOT / path).resolve()
    relative = os.path.relpath(str(target), str(TEST_ROOT)).replace(os.sep, "/")
    return "[{}]({})".format(label or "`{}`".format(path), relative)


def first_readme_paragraph(text: str) -> str:
    """Extract the first explanatory paragraph after a Markdown heading."""

    paragraphs: List[str] = []
    current: List[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            if current:
                paragraphs.append(" ".join(current))
                current = []
            continue
        if line.startswith("#") or line.startswith("```"):
            continue
        current.append(line)
    if current:
        paragraphs.append(" ".join(current))
    return one_line(paragraphs[0]) if paragraphs else ""


def iter_catalog_sources(repo_root: Path, test_root: Path) -> Iterable[Path]:
    """Yield maintained prose/test sources without scanning generated data."""

    roots = [
        test_root / "tests",
        test_root / "native",
        repo_root / "test" / "srcCoreUnitTests",
        repo_root / "src" / "models" / "sep_common",
        repo_root / "src" / "models" / "swcme" / "test",
        repo_root / "src" / "models" / "swcme" / "docs",
        repo_root / "src" / "models" / "sep_coronal_cme" / "test",
        repo_root / "src" / "models" / "sep_coronal_cme" / "docs",
        repo_root / "src" / "models" / "sep_corona_swcme" / "test",
        repo_root / "src" / "models" / "sep_corona_swcme" / "docs",
        repo_root / "srcSEP3D" / "test",
        # Application suites keep their independent scientific descriptions
        # beside the applications.  Scanning those maintained test trees lets
        # the global catalog describe imported Moon/Enceladus identities
        # without copying or inventing application-level acceptance prose.
        repo_root / "srcMoon" / "test",
        repo_root / "srcEnceladus" / "test",
    ]
    seen = set()
    for root in roots:
        if not root.exists():
            continue
        paths = [root] if root.is_file() else root.rglob("*")
        for path in paths:
            if (not path.is_file() or path.suffix.lower() not in SOURCE_SUFFIXES or
                    any(part in SKIP_PARTS for part in path.parts)):
                continue
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            # Catalog extraction should remain cheap and deterministic. Very
            # large data/documents fall back to the stable-ID description.
            try:
                if path.stat().st_size > 4 * 1024 * 1024:
                    continue
            except OSError:
                continue
            yield path


def add_description(
    target: MutableMapping[str, List[Description]],
    identifier: str,
    text: str,
    source: Path,
    priority: int,
    order: int,
    relevant_keys: set[str],
) -> None:
    normalized = normalize_id(identifier)
    cleaned = one_line(text)
    if (not CASE_ID.match(normalized) or normalized not in relevant_keys or
            len(cleaned) < 4):
        return
    relative = source.resolve().relative_to(REPO_ROOT).as_posix()
    candidate = Description(cleaned, relative, priority, order)
    rows = target.setdefault(normalized, [])
    if not any(row.text == candidate.text and row.source == candidate.source for row in rows):
        rows.append(candidate)


def json_descriptions(
    value: object,
    source: Path,
    target: MutableMapping[str, List[Description]],
    order: List[int],
    relevant_keys: set[str],
) -> None:
    """Collect explicit id/title or id/description records recursively."""

    if isinstance(value, dict):
        identifier = value.get("id")
        description = value.get("title", value.get("description"))
        if isinstance(identifier, str) and isinstance(description, str):
            add_description(
                target, identifier, description, source, 90, order[0],
                relevant_keys)
            order[0] += 1
        for child in value.values():
            json_descriptions(child, source, target, order, relevant_keys)
    elif isinstance(value, list):
        for child in value:
            json_descriptions(child, source, target, order, relevant_keys)


def discover_descriptions(
    repo_root: Path,
    test_root: Path,
    case_ids: Sequence[str],
) -> Dict[str, List[Description]]:
    """Discover maintained descriptions; never infer scientific outcomes."""

    found: Dict[str, List[Description]] = {}
    order = [0]
    # Source documents use model-local IDs while the global registry adds
    # scope prefixes. Precompute every contiguous token subsequence that can
    # actually describe a registered case. This both rejects unrelated prose
    # IDs and keeps lookup linear-sized even for the large CME documents.
    relevant_keys: set[str] = set()
    for case_id in case_ids:
        tokens = normalize_id(case_id).split("-")
        for begin in range(len(tokens)):
            for end in range(begin + 1, len(tokens) + 1):
                candidate = "-".join(tokens[begin:end])
                if len(candidate) >= 4:
                    relevant_keys.add(candidate)
    macro = re.compile(
        r'AMPS_TEST_CASE\s*\(\s*"([A-Z][A-Z0-9-]+)"\s*,\s*"([^"\n]+)"')
    python_tuple = re.compile(
        r'\(\s*"([A-Z][A-Z0-9-]+)"\s*,\s*"([^"\n]+)"\s*,')
    # Reduced-front checks use one local RSH identity several times. Preserving
    # source order lets global occurrence suffixes (-01, -02, ...) select the
    # matching assertion prose without executing production code.
    check_call = re.compile(
        r',\s*"([A-Z][A-Z0-9-]+)"\s*,\s*"([^"\n]+)"\s*\)')
    markdown_row = re.compile(
        r'^\|\s*`?([A-Z][A-Z0-9-]+)`?\s*\|\s*([^|\n]+?)\s*\|', re.MULTILINE)

    for path in iter_catalog_sources(repo_root, test_root):
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue

        # Model-local individual README files are the strongest descriptive
        # authority because they document one exact case and its limitations.
        if (path.name == "README.md" and path.parent.parent.name == "individual" and
                CASE_ID.match(path.parent.name)):
            paragraph = first_readme_paragraph(text)
            if paragraph:
                add_description(
                    found, path.parent.name, paragraph, path, 100, order[0],
                    relevant_keys)
                order[0] += 1

        for pattern, priority in (
                (macro, 95), (check_call, 94), (python_tuple, 85),
                (markdown_row, 70)):
            for identifier, description in pattern.findall(text):
                if set(description.strip()) <= {"-", ":"}:
                    continue
                add_description(
                    found, identifier, description, path, priority, order[0],
                    relevant_keys)
                order[0] += 1

        if path.suffix.lower() == ".json":
            try:
                json_descriptions(
                    json.loads(text), path, found, order, relevant_keys)
            except json.JSONDecodeError:
                # JSON validity belongs to the owning test. Catalog generation
                # ignores a malformed optional prose source rather than
                # turning it into a second data validator.
                pass

    for rows in found.values():
        rows.sort(key=lambda row: (-row.priority, row.order, row.source, row.text))
    return found


ABBREVIATIONS = {
    "AMR", "API", "ASAN", "ABI", "BC", "CFL", "CME", "CO2", "CUDA",
    "ECSIM", "GNU", "GPU", "H2O", "IDF", "IMF", "INMS", "IO", "MPI",
    "NVCC", "PIC", "PFSS", "RNG", "SEP", "SI", "SPICE", "SWCME", "TSAN",
    "UBSAN", "UVIS", "VHS", "VSS",
}


def humanize_case(identifier: str, common_prefix: Sequence[str]) -> str:
    """Create conservative prose from an ID when no maintained text exists."""

    tokens = identifier.split("-")
    prefix = list(common_prefix)
    if prefix and tokens[:len(prefix)] == prefix and len(tokens) > len(prefix):
        tokens = tokens[len(prefix):]
    words: List[str] = []
    for token in tokens:
        if token in ABBREVIATIONS:
            words.append(token)
        elif token.isdigit():
            words.append("case {}".format(token))
        else:
            match = re.fullmatch(r"([A-Z]+)([0-9]+)D([0-9]+)", token)
            if match:
                words.append("{} {}-D {}".format(
                    match.group(1), match.group(2), match.group(3)))
            else:
                words.append(token.replace("_", " ").lower())
    text = " ".join(words).strip()
    return (text[:1].upper() + text[1:]) if text else identifier


def common_case_prefix(cases: Sequence[str]) -> List[str]:
    split = [identifier.split("-") for identifier in cases]
    if not split:
        return []
    prefix: List[str] = []
    for values in zip(*split):
        if len(set(values)) != 1:
            break
        prefix.append(values[0])
    # Never strip the complete identity or leave only one numeric token.
    if any(len(parts) - len(prefix) < 2 for parts in split):
        return []
    return prefix


def resolve_description(
    case_id: str,
    discovered: Mapping[str, List[Description]],
    common_prefix: Sequence[str],
) -> Tuple[str, str | None]:
    normalized = normalize_id(case_id)
    exact = discovered.get(normalized)
    if exact:
        return exact[0].text, exact[0].source

    matches: List[Tuple[int, str, List[Description]]] = []
    # Look up only contiguous token subsequences of this case. Iterating every
    # discovered prose ID for every one of 1,028 cases made the documentation
    # drift check unnecessarily expensive on large model trees.
    tokens = normalized.split("-")
    for begin in range(len(tokens)):
        for end in range(begin + 1, len(tokens) + 1):
            local_id = "-".join(tokens[begin:end])
            rows = discovered.get(local_id)
            if rows is not None and len(local_id) >= 4:
                matches.append((len(local_id), local_id, rows))
    if matches:
        _, local_id, rows = max(matches, key=lambda item: (item[0], item[2][0].priority))
        best_priority = rows[0].priority
        best_rows = [row for row in rows if row.priority == best_priority]
        occurrence = re.search(
            r"(?:^|-){}-(\d{{2}})$".format(re.escape(local_id)), normalized)
        if occurrence:
            index = int(occurrence.group(1)) - 1
            if 0 <= index < len(best_rows):
                return best_rows[index].text, best_rows[index].source
        return best_rows[0].text, best_rows[0].source

    return humanize_case(case_id, common_prefix), None


def suite_how(suite: Mapping[str, object]) -> str:
    parts = [str(suite["level"])]
    if bool(suite.get("exclusive", False)):
        parts.append("serialized configuration owner")
    elif suite.get("mpi_ranks"):
        parts.append("{}-rank adapter".format(suite["mpi_ranks"]))
    else:
        parts.append("isolated suite adapter")
    tools = list(suite.get("requires", {}).get("tools", []))  # type: ignore[union-attr]
    if tools:
        parts.append("tools {}".format(", ".join(str(value) for value in tools)))
    return "; ".join(parts)


def registry_fingerprint(tests: Sequence[Mapping[str, object]]) -> str:
    payload = [
        {"id": suite["id"], "cases": list(suite["expected_cases"])}
        for suite in tests
    ]
    encoded = json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def render_catalog(
    repo_root: Path,
    test_root: Path,
    tests: Sequence[Mapping[str, object]],
    inventory: Mapping[str, object],
) -> str:
    """Return the complete deterministic Markdown catalog."""

    all_case_ids = [
        str(case_id)
        for suite in tests for case_id in suite["expected_cases"]
    ]
    discovered = discover_descriptions(repo_root, test_root, all_case_ids)
    inventory_by_component = {
        str(row.get("component")): row
        for row in inventory.get("components", [])  # type: ignore[union-attr]
        if isinstance(row, dict)
    }
    total_cases = sum(len(suite["expected_cases"]) for suite in tests)
    resolved_details: Dict[str, Tuple[str, str | None]] = {}
    maintained_descriptions = 0
    for suite in tests:
        cases = [str(value) for value in suite["expected_cases"]]
        prefix = common_case_prefix(cases)
        for case_id in cases:
            detail = resolve_description(case_id, discovered, prefix)
            resolved_details[case_id] = detail
            maintained_descriptions += detail[1] is not None
    lines = [
        "# AMPS core-test catalog",
        "",
        "<!-- Generated by src/test/generate_test_catalog.py; do not edit by hand. -->",
        "",
        "This catalog describes every suite and individual case in the expanded",
        "global registry. `manifest.json` and `case_registries.json` remain the",
        "execution authorities; an actual run's `amps-core-tests.json` remains the",
        "authority for executed status, messages, logs, and evidence.",
        "",
        "Regenerate after any suite/case/ownership change:",
        "",
        "```sh",
        "python3 src/test/generate_test_catalog.py",
        "python3 src/test/generate_test_catalog.py --check",
        "```",
        "",
        "Registry summary: **{} suites / {} cases**. Registry SHA-256: `{}`.".format(
            len(tests), total_cases, registry_fingerprint(tests)),
        "Description provenance: **{} maintained source descriptions / {} "
        "stable-ID fallbacks**.".format(
            maintained_descriptions, total_cases - maintained_descriptions),
        "",
        "## Suite index",
        "",
        "| Suite | Component | Level | Campaigns | Cases | Isolation | What | Where |",
        "| --- | --- | --- | --- | ---: | --- | --- | --- |",
    ]

    for suite in tests:
        runner = str(suite["runner"])
        isolation = "exclusive" if suite.get("exclusive", False) else "nonexclusive"
        lines.append("| `{}` | `{}` | {} | {} | {} | {} | {} | {} |".format(
            suite["id"], suite["component"], suite["level"],
            ", ".join(str(value) for value in suite["tiers"]),
            len(suite["expected_cases"]), isolation,
            markdown_cell(str(suite["description"])),
            repository_link("src/test/" + runner, "runner")))

    lines.extend(["", "## Individual cases", ""])
    for suite_index, suite in enumerate(tests, 1):
        suite_id = str(suite["id"])
        component = str(suite["component"])
        runner = str(suite["runner"])
        cases = list(suite["expected_cases"])
        inventory_row = inventory_by_component.get(component, {})
        coverage = str(inventory_row.get("status", "unclassified"))
        how = suite_how(suite)
        why = "{} coverage evidence for `{}`; preserves this contract as an exact global result".format(
            coverage, component)
        component_target = component
        where_base = "{}; {}".format(
            repository_link(component_target, component),
            repository_link("src/test/" + runner, "adapter"))

        lines.extend([
            "### {}. `{}`".format(suite_index, suite_id),
            "",
            "{} **{} cases.** {}".format(
                str(suite["description"]).rstrip("."), len(cases),
                "Runs behind an exclusive checkout barrier." if suite.get("exclusive", False)
                else "Eligible for suite-level concurrency."),
            "",
            "| Case | What | How | Why | Where |",
            "| --- | --- | --- | --- | --- |",
        ])
        for case_id in cases:
            description, source = resolved_details[str(case_id)]
            where = where_base
            if source and source not in {component, "src/test/" + runner}:
                where += "; " + repository_link(source, "case source")
            lines.append("| `{}` | {} | {} | {} | {} |".format(
                case_id, markdown_cell(description), markdown_cell(how),
                markdown_cell(why), where))
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def load_inputs() -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    tests = load_manifest(TEST_ROOT / "manifest.json", REPO_ROOT, TEST_ROOT)
    inventory = json.loads(
        (TEST_ROOT / "component_inventory.json").read_text(encoding="utf-8"))
    if not isinstance(inventory, dict) or inventory.get("schema_version") != 1:
        raise ValueError("component_inventory.json must be a schema_version 1 object")
    return tests, inventory


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true",
        help="fail if TEST_CATALOG.md differs from the generated registry")
    parser.add_argument(
        "--output", type=Path, default=CATALOG_PATH,
        help="catalog output path (default: src/test/TEST_CATALOG.md)")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        tests, inventory = load_inputs()
        rendered = render_catalog(REPO_ROOT, TEST_ROOT, tests, inventory)
    except (OSError, ValueError, ManifestError) as error:
        print("catalog generation ERROR: {}".format(error), file=sys.stderr)
        return 2

    output = args.output.expanduser().resolve()
    if args.check:
        try:
            current = output.read_text(encoding="utf-8")
        except OSError as error:
            print("catalog check ERROR: {}".format(error), file=sys.stderr)
            return 1
        if current != rendered:
            print(
                "catalog check FAIL: regenerate with "
                "python3 src/test/generate_test_catalog.py", file=sys.stderr)
            return 1
        print("PASS test catalog: {} suites / {} cases".format(
            len(tests), sum(len(suite["expected_cases"]) for suite in tests)))
        return 0

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(rendered)
    os.replace(str(temporary), str(output))
    print("wrote {} ({} suites / {} cases)".format(
        output, len(tests), sum(len(suite["expected_cases"]) for suite in tests)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
