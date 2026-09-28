#!/usr/bin/env python3
"""Build and validate the canonical ``sep_coronal_cme/model.md`` document.

The maintained prose is intentionally split by semantic ownership.  This
program performs the mechanical operation that humans should not perform by
hand: it validates the structured trace registry, expands the two generated
tables, orders numbered sections, and writes the public review artifact.

The registry file is JSON syntax stored in ``requirements.yaml``.  JSON is a
strict subset of YAML 1.2, so the file remains valid YAML while this generator
uses only Python's standard library.  Keeping the documentation gate free of a
PyYAML dependency matters on AMPS build hosts with deliberately small Python
installations.
"""

from __future__ import annotations

import argparse
import difflib
import json
import os
from pathlib import Path
import re
import sys
import tempfile
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


MODEL_ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = MODEL_ROOT / "model" / "requirements.yaml"
CANONICAL_PATH = MODEL_ROOT / "model.md"

# These ownership boundaries are part of the documentation architecture, not
# configurable presentation preferences.  Mirroring them here lets the gate
# detect a registry mutation that moves normative prose to the wrong owner.
EXPECTED_OWNERSHIP: Mapping[str, Tuple[int, ...]] = {
    "physics": (1, 2, *range(4, 14), 18, 19),
    "architecture_exchange": (3, 15),
    "configuration_validation": (14,),
    "testing_validation": (16, 17, 20),
}
EXPECTED_PATHS: Mapping[str, str] = {
    "physics": "model/physics.md",
    "architecture_exchange": "model/architecture_exchange.md",
    "configuration_validation": "model/configuration_validation.md",
    "testing_validation": "model/testing_validation.md",
}
EXPECTED_PLACEHOLDERS: Mapping[str, str] = {
    "REVIEW_DISPOSITION": "physics",
    "REQUIREMENT_REGISTRY": "architecture_exchange",
}

# ``schema_version`` is meaningful only if misspelled/new members cannot be
# silently ignored.  Keep these sets next to the other version-1 constants so
# adding a registry field is an explicit schema change rather than an
# accidental no-op in the generator.
EXPECTED_REGISTRY_KEYS = frozenset({
    "schema_version", "document", "requirements", "review_findings",
})
EXPECTED_DOCUMENT_KEYS = frozenset({
    "canonical_output", "section_order", "modules", "generated_placeholders",
})
EXPECTED_MODULE_KEYS = frozenset({"id", "path", "owns_preamble", "sections"})
EXPECTED_REQUIREMENT_KEYS = frozenset({
    "id", "review_item", "title", "subject", "status", "review_sources",
    "sections", "config_keys", "api_records", "test_ids",
})
EXPECTED_REVIEW_FINDING_KEYS = frozenset({
    "id", "label", "status", "requirement_ids", "resolution",
})

# Requirement status is release metadata, not unrestricted explanatory prose.
# These are the complete schema-1 states currently used by the stable R/N/P
# records.  A new lifecycle state therefore requires a deliberate schema/gate
# update and cannot enter through a typo that is invisible in generated prose.
EXPECTED_REQUIREMENT_STATUSES = frozenset({
    "specified",
    "implemented",
    "specified-diagnostic; future-options-not-implemented",
    "documented-limitation; future-capability-not-implemented",
    "specified-validation-protocol; evidence-not-yet-produced",
    "specified-offline-benchmark; evidence-not-yet-produced",
    "documented-optional-extension; not-implemented",
    "current-gates-specified; continuous-envelope-not-implemented",
})
EXPECTED_REVIEW_STATUSES = frozenset({"resolved", "superseded", "documented"})

# The three review series are deliberately closed sets.  Treating this as a
# generator contract prevents a newly accepted review finding from appearing
# only in narrative prose without a stable requirement, links, and tests.
EXPECTED_REVIEW_ITEMS: Tuple[str, ...] = tuple(
    [f"R{number}" for number in range(1, 11)]
    + [f"N{number}" for number in range(1, 8)]
    + [f"P{number}" for number in range(1, 9)])
REQUIREMENT_ID = re.compile(
    r"^SCCM-(?P<review_item>R(?:[1-9]|10)|N[1-7]|P[1-8])-"
    r"[A-Z0-9]+(?:-[A-Z0-9]+)*$")
REVIEW_SOURCE = re.compile(
    r"^(?:review-2:R(?:[1-9]|10)|review-3:N[1-7]|review-4:P[1-8])$")
PRIMARY_REVIEW_SOURCE = {"R": "review-2", "N": "review-3", "P": "review-4"}
TEST_ID = re.compile(r"^[A-Z][A-Z0-9-]*[0-9](?:[A-Z0-9-]*)$")
# Match the complete placeholder namespace, including malformed/lower-case
# names, so a typo cannot pass through as ordinary prose.
PLACEHOLDER = re.compile(r"\{\{GENERATED:([^{}\r\n]+)\}\}")
# Capture every definition-shaped bullet first, then validate its ID.  Using a
# permissive capture prevents a malformed/lower-case ID from disappearing from
# the canonical registry merely because it did not match the ID grammar.
TEST_DEFINITION = re.compile(r"^- `([^`\r\n]+)`:", re.MULTILINE)
ROADMAP_TEST_TOKEN = re.compile(
    r"^(?P<prefix>[A-Z][A-Z0-9-]*?)(?P<first>[0-9]+)"
    r"(?:--(?P<last>[0-9]+))?$")
EXPECTED_ROADMAP_STAGES: Tuple[str, ...] = tuple(
    [str(number) for number in range(0, 11)]
    + ["11A", "11B", "12", "13", "14 (future/campaign; non-release)"])

# Stage 13's prose dependency on all prior release stages is not duplicate
# ownership.  The sole literal cross-row reuse is SRC3D19: Stage 9 owns its
# baseline gate and Stage 14A explicitly reuses it.  Any later reuse must be
# named here and annotated in the table, making the exception reviewable.
EXPECTED_ROADMAP_REUSE: Mapping[str, Tuple[str, ...]] = {
    "SRC3D19": ("9", "14 (future/campaign; non-release)"),
}


class ValidationError(RuntimeError):
    """A deterministic source, registry, or generated-output contract failed."""


def read_text_strict(path: Path) -> str:
    """Read canonical UTF-8/LF text without silently normalizing bad bytes."""

    try:
        data = path.read_bytes()
    except OSError as exc:
        raise ValidationError(f"cannot read {path}: {exc}") from exc
    if data.startswith(b"\xef\xbb\xbf"):
        raise ValidationError(f"{path} has a forbidden UTF-8 byte-order mark")
    if b"\r" in data:
        raise ValidationError(f"{path} is not LF-only")
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValidationError(f"{path} is not valid UTF-8: {exc}") from exc
    if not text.endswith("\n"):
        raise ValidationError(f"{path} must end with one LF newline")
    return text


def _fence_token(line: str) -> Tuple[str, int, str] | None:
    """Return a CommonMark fence's character, length, and trailing text."""

    match = re.match(r"^ {0,3}(`{3,}|~{3,})(.*)$", line.rstrip("\n"))
    if not match:
        return None
    token = match.group(1)
    return token[0], len(token), match.group(2)


def split_numbered_sections(
        text: str, source_name: str) -> Tuple[str, Dict[int, str]]:
    """Split real level-2 numeric sections while ignoring fenced examples.

    Byte offsets are accumulated from the original strings.  No newline,
    whitespace, or Unicode normalization occurs, which is what makes the
    generated canonical file a byte-reproducible review artifact.
    """

    starts: List[Tuple[int, int]] = []
    offset = 0
    open_character: str | None = None
    open_length = 0
    for line in text.splitlines(keepends=True):
        fence = _fence_token(line)
        if fence is not None:
            character, length, trailing = fence
            if open_character is None:
                open_character = character
                open_length = length
            elif (character == open_character and length >= open_length
                  and not trailing.strip()):
                open_character = None
                open_length = 0
        elif open_character is None:
            heading = re.match(r"^## ([1-9][0-9]*)\. ", line)
            if heading:
                starts.append((int(heading.group(1)), offset))
        offset += len(line)
    if open_character is not None:
        raise ValidationError(f"{source_name} has an unclosed Markdown fence")

    numbers = [number for number, _ in starts]
    duplicates = sorted({number for number in numbers if numbers.count(number) > 1})
    if duplicates:
        raise ValidationError(
            f"{source_name} has duplicate numbered sections: {duplicates}")

    sections: Dict[int, str] = {}
    for index, (number, start) in enumerate(starts):
        stop = starts[index + 1][1] if index + 1 < len(starts) else len(text)
        sections[number] = text[start:stop]
    preamble = text[:starts[0][1]] if starts else text
    return preamble, sections


def _require_mapping(value: object, context: str) -> Mapping[str, object]:
    if not isinstance(value, dict):
        raise ValidationError(f"{context} must be an object")
    return value


def _require_list(value: object, context: str) -> List[object]:
    if not isinstance(value, list):
        raise ValidationError(f"{context} must be an array")
    return value


def _require_exact_keys(
        value: Mapping[str, object], expected: frozenset[str], context: str) -> None:
    """Reject missing and additional schema members with stable diagnostics."""

    observed = set(value)
    missing = sorted(expected - observed)
    unknown = sorted(observed - expected)
    if missing or unknown:
        details: List[str] = []
        if missing:
            details.append(f"missing members {missing}")
        if unknown:
            details.append(f"unknown members {unknown}")
        raise ValidationError(f"{context} has " + "; ".join(details))


def _require_string_list(
        value: object, context: str, *, nonempty: bool = False) -> List[str]:
    """Return a duplicate-free string array under the registry schema."""

    raw = _require_list(value, context)
    if nonempty and not raw:
        raise ValidationError(f"{context} must be nonempty")
    if not all(isinstance(item, str) for item in raw):
        raise ValidationError(f"{context} must contain only strings")
    strings = [item for item in raw if isinstance(item, str)]
    duplicate_values = sorted(
        {item for item in strings if strings.count(item) > 1})
    if duplicate_values:
        raise ValidationError(
            f"{context} contains duplicate values: {duplicate_values}")
    return strings


def load_registry(path: Path = REGISTRY_PATH) -> Mapping[str, object]:
    """Load the JSON-compatible YAML registry with a useful parse error."""

    text = read_text_strict(path)
    try:
        registry = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValidationError(
            f"{path} is not valid JSON-compatible YAML: line {exc.lineno}, "
            f"column {exc.colno}: {exc.msg}") from exc
    root = _require_mapping(registry, "registry root")
    _require_exact_keys(root, EXPECTED_REGISTRY_KEYS, "registry root")
    return root


def render_review_disposition(registry: Mapping[str, object]) -> str:
    """Render the review table from structured rows in registry order."""

    lines = [
        "| Review item | Resolution in this revision |",
        "|---|---|",
    ]
    for index, raw in enumerate(_require_list(
            registry.get("review_findings"), "review_findings")):
        row = _require_mapping(raw, f"review_findings[{index}]")
        label = row.get("label")
        resolution = row.get("resolution")
        if not isinstance(label, str) or not isinstance(resolution, str):
            raise ValidationError(
                f"review_findings[{index}] label/resolution must be strings")
        if "\n" in label or "\n" in resolution or "|" in label or "|" in resolution:
            raise ValidationError(
                f"review finding {row.get('id')!r} cannot be rendered as one table row")
        lines.append(f"| {label} | {resolution} |")
    return "\n".join(lines)


def render_requirement_registry(registry: Mapping[str, object]) -> str:
    """Render the stable requirement table from its sole structured authority."""

    lines = [
        "| Requirement ID | Review item | Normative subject |",
        "|---|---:|---|",
    ]
    for index, raw in enumerate(_require_list(
            registry.get("requirements"), "requirements")):
        requirement = _require_mapping(raw, f"requirements[{index}]")
        requirement_id = requirement.get("id")
        review_item = requirement.get("review_item")
        subject = requirement.get("subject")
        if not all(isinstance(item, str)
                   for item in (requirement_id, review_item, subject)):
            raise ValidationError(
                f"requirements[{index}] id/review_item/subject must be strings")
        if any("\n" in item or "|" in item
               for item in (requirement_id, review_item, subject)):
            raise ValidationError(
                f"requirement {requirement_id!r} cannot be rendered as one table row")
        lines.append(f"| `{requirement_id}` | {review_item} | {subject} |")
    return "\n".join(lines)


def _section_reference_exists(reference: str, module_text: str) -> bool:
    if re.fullmatch(r"[1-9][0-9]*", reference):
        pattern = rf"^## {re.escape(reference)}\. "
    elif re.fullmatch(r"[1-9][0-9]*\.[1-9][0-9]*", reference):
        pattern = rf"^### {re.escape(reference)}(?:\s|\.)"
    else:
        return False
    return re.search(pattern, module_text, re.MULTILINE) is not None


def _configuration_key_exists(reference: str, configuration_text: str) -> bool:
    """Resolve a ``section.key`` trace against the Section-14 input grammar."""

    if reference.count(".") != 1:
        return False
    section, key = reference.split(".")
    header = re.search(rf"^\[{re.escape(section)}\]\s*$", configuration_text,
                       re.MULTILINE)
    if header is None:
        return False
    next_header = re.search(r"^\[[^]]+\]\s*$", configuration_text[header.end():],
                            re.MULTILINE)
    stop = header.end() + next_header.start() if next_header else len(configuration_text)
    block = configuration_text[header.end():stop]
    return re.search(rf"^{re.escape(key)}\s*=", block, re.MULTILINE) is not None


def _api_record_exists(name: str, architecture_text: str) -> bool:
    return re.search(rf"\b(?:enum class|class|struct)\s+{re.escape(name)}\b",
                     architecture_text) is not None


def _extract_level_three_subsection(
        text: str, reference: str, source_name: str) -> str:
    """Return one real level-three subsection, ignoring fenced examples."""

    start: int | None = None
    offset = 0
    open_character: str | None = None
    open_length = 0
    for line in text.splitlines(keepends=True):
        fence = _fence_token(line)
        if fence is not None:
            character, length, trailing = fence
            if open_character is None:
                open_character = character
                open_length = length
            elif (character == open_character and length >= open_length
                  and not trailing.strip()):
                open_character = None
                open_length = 0
        elif open_character is None:
            heading = re.match(
                r"^### ([1-9][0-9]*\.[1-9][0-9]*)(?:\s|\.)", line)
            if heading:
                if start is not None:
                    return text[start:offset]
                if heading.group(1) == reference:
                    start = offset
        offset += len(line)
    if start is None:
        raise ValidationError(
            f"{source_name} does not contain subsection {reference}")
    return text[start:]


def _expand_roadmap_test_token(token: str, context: str) -> List[str]:
    """Expand one exact ID or same-family inclusive ``NN--MM`` range."""

    match = ROADMAP_TEST_TOKEN.fullmatch(token)
    if match is None:
        raise ValidationError(f"{context} has malformed test token {token!r}")
    prefix = match.group("prefix")
    first_text = match.group("first")
    last_text = match.group("last")
    if last_text is None:
        if not TEST_ID.fullmatch(token):
            raise ValidationError(f"{context} has invalid test ID {token!r}")
        return [token]
    if len(first_text) != len(last_text):
        raise ValidationError(
            f"{context} range {token!r} must use equal-width endpoints")
    first = int(first_text)
    last = int(last_text)
    if last < first:
        raise ValidationError(f"{context} range {token!r} is descending")
    return [
        f"{prefix}{number:0{len(first_text)}d}"
        for number in range(first, last + 1)
    ]


def _validate_roadmap_test_trace(
        section_17_2: str, canonical_tests: set[str]) -> None:
    """Validate the Stage-to-test table, including deliberate test reuse.

    The prose says range completeness is a release obligation.  Parse the one
    canonical table here so a range typo, an unstaged definition, or an
    unannounced cross-stage reuse cannot survive by remaining ordinary prose.
    """

    lines = section_17_2.splitlines()
    header = "| Stage | Canonical required test records |"
    header_positions = [index for index, line in enumerate(lines) if line == header]
    if len(header_positions) != 1:
        raise ValidationError(
            "Section 17.2 must contain exactly one roadmap-to-test table")
    header_index = header_positions[0]
    if (header_index + 1 >= len(lines)
            or lines[header_index + 1] != "|---:|---|"):
        raise ValidationError("roadmap-to-test table has a malformed separator")

    rows: List[Tuple[str, str]] = []
    row_index = header_index + 2
    while row_index < len(lines) and lines[row_index].startswith("|"):
        match = re.fullmatch(r"\| ([^|]+?) \| (.*) \|", lines[row_index])
        if match is None:
            raise ValidationError(
                f"malformed roadmap-to-test row {lines[row_index]!r}")
        rows.append((match.group(1), match.group(2)))
        row_index += 1

    observed_stages = tuple(stage for stage, _ in rows)
    if observed_stages != EXPECTED_ROADMAP_STAGES:
        raise ValidationError(
            f"roadmap stage order is {observed_stages}, "
            f"expected {EXPECTED_ROADMAP_STAGES}")

    assignments: Dict[str, List[str]] = {}
    row_bodies: Dict[str, str] = {}
    for stage, body in rows:
        row_bodies[stage] = body
        tokens = re.findall(r"`([^`\r\n]+)`", body)
        if not tokens:
            raise ValidationError(f"roadmap stage {stage} has no test records")
        row_tests: set[str] = set()
        for token in tokens:
            for test_id in _expand_roadmap_test_token(
                    token, f"roadmap stage {stage}"):
                if test_id in row_tests:
                    raise ValidationError(
                        f"roadmap stage {stage} repeats test {test_id}")
                row_tests.add(test_id)
                assignments.setdefault(test_id, []).append(stage)

    # Stage 13 is the aggregate release gate.  Its dependency is prose rather
    # than 200+ duplicate literal IDs, so preserve that edge explicitly while
    # keeping those earlier tests owned by their original stages.
    if "plus every mandatory record from Stages 0--12" not in row_bodies["13"]:
        raise ValidationError(
            "roadmap Stage 13 must include every mandatory Stage 0--12 record")

    unknown = sorted(set(assignments) - canonical_tests)
    if unknown:
        raise ValidationError(
            f"roadmap references unknown canonical tests: {unknown}")
    unstaged = sorted(canonical_tests - set(assignments))
    if unstaged:
        raise ValidationError(
            f"canonical tests have no roadmap stage: {unstaged}")

    repeated = {
        test_id: tuple(stages)
        for test_id, stages in assignments.items() if len(stages) > 1
    }
    if set(repeated) != set(EXPECTED_ROADMAP_REUSE):
        unexpected = sorted(set(repeated) - set(EXPECTED_ROADMAP_REUSE))
        missing = sorted(set(EXPECTED_ROADMAP_REUSE) - set(repeated))
        raise ValidationError(
            "roadmap test reuse differs from the explicit contract: "
            f"unexpected={unexpected}, missing={missing}")
    for test_id, expected_stages in EXPECTED_ROADMAP_REUSE.items():
        observed = repeated[test_id]
        if observed != expected_stages:
            raise ValidationError(
                f"roadmap reuse for {test_id} is {observed}, "
                f"expected {expected_stages}")
        reuse_stage = expected_stages[-1]
        if f"reuses `{test_id}`" not in row_bodies[reuse_stage]:
            raise ValidationError(
                f"roadmap reuse of {test_id} in {reuse_stage} must be explicit")


def _validate_document_contract(registry: Mapping[str, object]) -> None:
    if registry.get("schema_version") != 1:
        raise ValidationError("requirements.yaml schema_version must be 1")
    document = _require_mapping(registry.get("document"), "document")
    _require_exact_keys(document, EXPECTED_DOCUMENT_KEYS, "document")
    if document.get("canonical_output") != "model.md":
        raise ValidationError("document.canonical_output must be 'model.md'")
    if document.get("section_order") != list(range(1, 21)):
        raise ValidationError("document.section_order must be the integers 1 through 20")

    raw_modules = _require_list(document.get("modules"), "document.modules")
    if len(raw_modules) != 4:
        raise ValidationError("document.modules must contain exactly four owners")
    observed: Dict[str, Tuple[str, Tuple[int, ...], bool]] = {}
    for index, raw in enumerate(raw_modules):
        module = _require_mapping(raw, f"document.modules[{index}]")
        _require_exact_keys(
            module, EXPECTED_MODULE_KEYS, f"document.modules[{index}]")
        module_id = module.get("id")
        path = module.get("path")
        sections = module.get("sections")
        preamble = module.get("owns_preamble")
        if not isinstance(module_id, str) or not isinstance(path, str):
            raise ValidationError(f"document.modules[{index}] id/path must be strings")
        if not isinstance(sections, list) or not all(isinstance(n, int) for n in sections):
            raise ValidationError(f"document.modules[{index}].sections must be integers")
        if not isinstance(preamble, bool):
            raise ValidationError(
                f"document.modules[{index}].owns_preamble must be Boolean")
        if module_id in observed:
            raise ValidationError(f"duplicate module ID {module_id!r}")
        observed[module_id] = (path, tuple(sections), preamble)

    if set(observed) != set(EXPECTED_OWNERSHIP):
        raise ValidationError("document.modules does not name the four canonical owners")
    for module_id, expected_sections in EXPECTED_OWNERSHIP.items():
        path, sections, owns_preamble = observed[module_id]
        if path != EXPECTED_PATHS[module_id]:
            raise ValidationError(f"module {module_id} must use {EXPECTED_PATHS[module_id]}")
        if sections != expected_sections:
            raise ValidationError(
                f"module {module_id} owns {sections}, expected {expected_sections}")
        if owns_preamble != (module_id == "physics"):
            raise ValidationError("only the physics module may own the preamble")

    placeholders = _require_mapping(
        document.get("generated_placeholders"),
        "document.generated_placeholders")
    if placeholders != EXPECTED_PLACEHOLDERS:
        raise ValidationError(
            "document.generated_placeholders does not match the canonical mapping")


def _load_and_validate_modules(
        registry: Mapping[str, object]) -> Tuple[Dict[str, str], str, Dict[int, str]]:
    """Load all owners, enforce placement, and return canonical source slices."""

    _validate_document_contract(registry)
    module_texts: Dict[str, str] = {}
    preamble = ""
    all_sections: Dict[int, str] = {}
    placeholder_locations: Dict[str, List[str]] = {}

    for module_id, expected_sections in EXPECTED_OWNERSHIP.items():
        path = MODEL_ROOT / EXPECTED_PATHS[module_id]
        text = read_text_strict(path)
        module_texts[module_id] = text
        module_preamble, sections = split_numbered_sections(text, str(path))
        if tuple(sections) != expected_sections:
            raise ValidationError(
                f"{path} contains sections {tuple(sections)}, expected {expected_sections}")
        if module_id == "physics":
            if not module_preamble.startswith("# "):
                raise ValidationError("physics module must own the level-1 title/preamble")
            preamble = module_preamble
        elif module_preamble:
            raise ValidationError(f"{path} has text before its first owned section")
        for number, section in sections.items():
            if number in all_sections:
                raise ValidationError(f"numbered section {number} has multiple owners")
            all_sections[number] = section
        for name in PLACEHOLDER.findall(text):
            placeholder_locations.setdefault(name, []).append(module_id)

    unknown = sorted(set(placeholder_locations) - set(EXPECTED_PLACEHOLDERS))
    if unknown:
        raise ValidationError(f"unknown generated placeholders: {unknown}")
    for name, owner in EXPECTED_PLACEHOLDERS.items():
        locations = placeholder_locations.get(name, [])
        if locations != [owner]:
            raise ValidationError(
                f"placeholder {name} must occur exactly once in {owner}; got {locations}")
    if set(all_sections) != set(range(1, 21)):
        raise ValidationError(
            f"numbered section ownership is incomplete: {sorted(all_sections)}")
    return module_texts, preamble, all_sections


def _validate_requirements_and_reviews(
        registry: Mapping[str, object], module_texts: Mapping[str, str]) -> None:
    requirements = _require_list(registry.get("requirements"), "requirements")
    expected_items = list(EXPECTED_REVIEW_ITEMS)
    if len(requirements) != len(expected_items):
        raise ValidationError(
            "requirements must contain SCCM R1--R10, N1--N7, and P1--P8")

    requirement_by_id: Dict[str, Mapping[str, object]] = {}
    observed_items: List[str] = []
    for index, raw in enumerate(requirements):
        requirement = _require_mapping(raw, f"requirements[{index}]")
        _require_exact_keys(
            requirement, EXPECTED_REQUIREMENT_KEYS, f"requirements[{index}]")
        requirement_id = requirement.get("id")
        review_item = requirement.get("review_item")
        requirement_match = (
            REQUIREMENT_ID.fullmatch(requirement_id)
            if isinstance(requirement_id, str) else None)
        if requirement_match is None:
            raise ValidationError(f"invalid requirement ID {requirement_id!r}")
        if requirement_id in requirement_by_id:
            raise ValidationError(f"duplicate requirement ID {requirement_id}")
        if not isinstance(review_item, str):
            raise ValidationError(f"{requirement_id}.review_item must be a string")
        if requirement_match.group("review_item") != review_item:
            raise ValidationError(
                f"{requirement_id} ID series does not match review_item "
                f"{review_item!r}")
        requirement_by_id[requirement_id] = requirement
        observed_items.append(review_item)

        for scalar_name in ("title", "subject"):
            if not isinstance(requirement.get(scalar_name), str) or not requirement[scalar_name]:
                raise ValidationError(f"{requirement_id}.{scalar_name} must be nonempty")
        status = requirement.get("status")
        if status not in EXPECTED_REQUIREMENT_STATUSES:
            raise ValidationError(
                f"{requirement_id}.status {status!r} is not a schema-1 status")
        review_sources = _require_string_list(
            requirement.get("review_sources"),
            f"{requirement_id}.review_sources", nonempty=True)
        invalid_sources = sorted(
            source for source in review_sources
            if REVIEW_SOURCE.fullmatch(source) is None)
        if invalid_sources:
            raise ValidationError(
                f"{requirement_id}.review_sources has invalid values: "
                f"{invalid_sources}")
        primary_source = f"{PRIMARY_REVIEW_SOURCE[review_item[0]]}:{review_item}"
        if primary_source not in review_sources:
            raise ValidationError(
                f"{requirement_id}.review_sources must include "
                f"{primary_source!r}")

        section_refs = _require_string_list(
            requirement.get("sections"),
            f"{requirement_id}.sections", nonempty=True)
        for reference in section_refs:
            top = int(reference.split(".", 1)[0]) if reference.split(".", 1)[0].isdigit() else -1
            owner = next((module_id for module_id, numbers in EXPECTED_OWNERSHIP.items()
                          if top in numbers), None)
            if owner is None or not _section_reference_exists(reference, module_texts[owner]):
                raise ValidationError(
                    f"{requirement_id} has unresolved section link {reference!r}")

        config_keys = _require_string_list(
            requirement.get("config_keys"), f"{requirement_id}.config_keys")
        for reference in config_keys:
            if not _configuration_key_exists(
                    reference, module_texts["configuration_validation"]):
                raise ValidationError(
                    f"{requirement_id} has unresolved configuration key {reference!r}")

        api_records = _require_string_list(
            requirement.get("api_records"), f"{requirement_id}.api_records")
        for record in api_records:
            if not _api_record_exists(record, module_texts["architecture_exchange"]):
                raise ValidationError(
                    f"{requirement_id} has unresolved API record {record!r}")

    if observed_items != expected_items:
        raise ValidationError(
            f"requirement review_item order is {observed_items}, expected {expected_items}")

    _, testing_sections = split_numbered_sections(
        module_texts["testing_validation"], EXPECTED_PATHS["testing_validation"])
    section_17_2 = _extract_level_three_subsection(
        testing_sections[17], "17.2", EXPECTED_PATHS["testing_validation"])
    test_definitions = TEST_DEFINITION.findall(section_17_2)
    invalid_test_definitions = sorted(
        test_id for test_id in test_definitions
        if TEST_ID.fullmatch(test_id) is None)
    if invalid_test_definitions:
        raise ValidationError(
            f"invalid canonical test definitions: {invalid_test_definitions}")
    duplicate_tests = sorted(
        {test_id for test_id in test_definitions if test_definitions.count(test_id) > 1})
    if duplicate_tests:
        raise ValidationError(f"duplicate canonical test definitions: {duplicate_tests}")
    canonical_tests = set(test_definitions)
    if "DOCSCCM01" not in canonical_tests:
        raise ValidationError("canonical test registry does not define DOCSCCM01")
    _validate_roadmap_test_trace(section_17_2, canonical_tests)

    for requirement_id, requirement in requirement_by_id.items():
        tests = _require_string_list(
            requirement.get("test_ids"),
            f"{requirement_id}.test_ids", nonempty=True)
        for test_id in tests:
            if not TEST_ID.fullmatch(test_id):
                raise ValidationError(f"{requirement_id} has invalid test ID {test_id!r}")
            if test_id not in canonical_tests:
                raise ValidationError(
                    f"{requirement_id} references unknown canonical test {test_id}")

    reviews = _require_list(registry.get("review_findings"), "review_findings")
    seen_review_ids: set[str] = set()
    referenced_requirements: set[str] = set()
    review_by_id: Dict[str, Mapping[str, object]] = {}
    for index, raw in enumerate(reviews):
        finding = _require_mapping(raw, f"review_findings[{index}]")
        _require_exact_keys(
            finding, EXPECTED_REVIEW_FINDING_KEYS,
            f"review_findings[{index}]")
        finding_id = finding.get("id")
        label = finding.get("label")
        if not isinstance(finding_id, str) or not finding_id:
            raise ValidationError(f"review_findings[{index}].id must be nonempty")
        if finding_id in seen_review_ids:
            raise ValidationError(f"duplicate review finding ID {finding_id}")
        seen_review_ids.add(finding_id)
        review_by_id[finding_id] = finding
        if not isinstance(label, str) or not label.startswith(finding_id + ":"):
            raise ValidationError(f"review finding {finding_id} has a mismatched label")
        if finding.get("status") not in EXPECTED_REVIEW_STATUSES:
            raise ValidationError(f"review finding {finding_id} has invalid status")
        linked = _require_string_list(
            finding.get("requirement_ids"),
            f"review finding {finding_id}.requirement_ids")
        for requirement_id in linked:
            if requirement_id not in requirement_by_id:
                raise ValidationError(
                    f"review finding {finding_id} links unknown requirement "
                    f"{requirement_id}")
            referenced_requirements.add(requirement_id)

    missing_review_rows = sorted(set(EXPECTED_REVIEW_ITEMS) - seen_review_ids)
    if missing_review_rows:
        raise ValidationError(f"missing review disposition rows: {missing_review_rows}")

    # A source tag is stronger than a general related-requirement link: it says
    # that the named review finding created or revised this exact requirement.
    # Require that relationship to be reciprocal so a well-formed but stale
    # source string cannot silently survive a review-table edit.
    for requirement_id, requirement in requirement_by_id.items():
        for source in requirement["review_sources"]:  # type: ignore[index]
            finding_id = str(source).split(":", 1)[1]
            finding = review_by_id.get(finding_id)
            if finding is None:
                raise ValidationError(
                    f"{requirement_id} references missing review source "
                    f"{source!r}")
            if requirement_id not in finding["requirement_ids"]:  # type: ignore[operator]
                raise ValidationError(
                    f"{requirement_id} review source {source!r} does not link "
                    "back to the requirement")
    orphan_requirements = sorted(set(requirement_by_id) - referenced_requirements)
    if orphan_requirements:
        raise ValidationError(
            f"requirements have no review-disposition link: {orphan_requirements}")


def _replace_placeholders(
        text: str, replacements: Mapping[str, str], context: str) -> str:
    for name, replacement in replacements.items():
        text = text.replace(f"{{{{GENERATED:{name}}}}}", replacement)
    remaining = PLACEHOLDER.findall(text)
    if remaining:
        raise ValidationError(f"{context} retains generated placeholders: {remaining}")
    return text


def _validate_bibliography(section_19: str) -> None:
    numbers = [int(value) for value in re.findall(r"^([1-9][0-9]*)\. ",
                                                  section_19, re.MULTILINE)]
    if not numbers:
        raise ValidationError("Section 19 contains no numbered references")
    expected = list(range(1, numbers[-1] + 1))
    if numbers != expected:
        raise ValidationError(
            f"Section 19 bibliography sequence is {numbers}, expected {expected}")


def generate(registry_path: Path = REGISTRY_PATH) -> str:
    """Validate every source authority and return deterministic canonical text."""

    registry = load_registry(registry_path)
    module_texts, preamble, sections = _load_and_validate_modules(registry)
    _validate_requirements_and_reviews(registry, module_texts)
    replacements = {
        "REVIEW_DISPOSITION": render_review_disposition(registry),
        "REQUIREMENT_REGISTRY": render_requirement_registry(registry),
    }
    rendered_preamble = _replace_placeholders(preamble, replacements, "preamble")
    rendered_sections = {
        number: _replace_placeholders(section, replacements, f"Section {number}")
        for number, section in sections.items()
    }
    _validate_bibliography(rendered_sections[19])
    output = rendered_preamble + "".join(
        rendered_sections[number] for number in range(1, 21))
    if not output.endswith("\n") or "\r" in output:
        raise ValidationError("generated output violated the LF/final-newline contract")
    if PLACEHOLDER.search(output):
        raise ValidationError("generated output contains an unresolved placeholder")
    return output


def atomic_write(path: Path, text: str) -> None:
    """Publish generated bytes atomically without temporary-file leakage."""

    path.parent.mkdir(parents=True, exist_ok=True)
    mode = path.stat().st_mode & 0o777 if path.exists() else 0o644
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", newline="\n",
                dir=path.parent, prefix=f".{path.name}.", suffix=".tmp",
                delete=False) as temporary:
            temporary.write(text)
            temporary.flush()
            os.fsync(temporary.fileno())
            temporary_name = temporary.name
        os.chmod(temporary_name, mode)
        os.replace(temporary_name, path)
        temporary_name = None
    finally:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass


def _diff(expected: str, actual: str, path: Path) -> str:
    lines = list(difflib.unified_diff(
        actual.splitlines(), expected.splitlines(),
        fromfile=str(path), tofile="generated:model.md", lineterm=""))
    limit = 200
    body = "\n".join(lines[:limit])
    if len(lines) > limit:
        body += f"\n... diff truncated ({len(lines) - limit} additional lines)"
    return body


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate or verify the shared sep_coronal_cme model specification")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--check", action="store_true",
        help="fail if committed model.md differs from deterministic generation")
    mode.add_argument(
        "--output", type=Path,
        help="write generated text to this path instead of canonical model.md")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        generated = generate()
        if args.check:
            current = read_text_strict(CANONICAL_PATH)
            if current != generated:
                print("DOCSCCM01 FAIL: committed model.md is not generated-clean",
                      file=sys.stderr)
                print(_diff(generated, current, CANONICAL_PATH), file=sys.stderr)
                return 1
            print("DOCSCCM01 PASS: registry, modules, links, and model.md are clean")
            return 0
        destination = args.output.resolve() if args.output else CANONICAL_PATH
        atomic_write(destination, generated)
        print(f"generated {destination}")
        return 0
    except (ValidationError, OSError) as exc:
        print(f"DOCSCCM01 ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
