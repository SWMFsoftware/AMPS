"""Shared, fail-closed support for importing CME model-local registries.

The global harness owns scheduling and evidence publication, but it must not
become a second physics oracle.  These helpers therefore copy the exact local
source trees into the suite work directory, execute the maintained local
runner, and translate only identity/status/message fields.  A missing,
duplicate, reordered, or foreign local row is an infrastructure ERROR rather
than a smaller green suite.

The isolated copy is important for portability: the model Makefiles place
objects beside their sources.  Running them in the checkout would make a
nominally read-only nightly depend on stale objects and would mutate a source
tree that may be mounted read-only on another machine.
"""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


STATUSES = {"PASS", "FAIL", "SKIP", "ERROR"}
SEVERITY = {"PASS": 0, "SKIP": 1, "FAIL": 2, "ERROR": 3}


def copy_tree(source: Path, target: Path) -> None:
    """Copy sources while excluding generated objects and cached evidence."""

    shutil.copytree(
        str(source), str(target),
        ignore=shutil.ignore_patterns(
            # ``output`` is deliberately *not* ignored: srcSEP3D/output is a
            # production source directory containing point-sampling headers.
            # Only the explicitly generated ``test_output`` tree is excluded.
            "build", "test_output", "__pycache__", "*.o", "*.d", "*.a",
            "stage1", "*.pyc"),
    )


def prepare_model_repo(repo: Path, work: Path, include_sep3d: bool = False) -> Path:
    """Create the smallest repository-shaped tree accepted by local Makefiles.

    Corona--SWCME includes its coronal and SEP-common siblings by relative
    paths, and both use the generic runtime-input header.  The SCCM build also
    compiles its two thin field-line adapters through repository-relative
    paths.  Copying those *real* adapter sources (rather than substituting a
    mock interface) preserves the production compile boundary without copying
    either complete application.  The reduced adapter additionally resolves
    all shared model paths from ``srcSEP3D/..`` and therefore requests the
    complete application tree.  Together these choices preserve actual
    include/link ownership; no altered compiler flags can make a local gate
    pass merely because an application dependency was hidden.
    """

    root = work / "isolated-repo"
    (root / "src" / "models").mkdir(parents=True)
    for name in ("sep_common", "sep_coronal_cme", "sep_corona_swcme", "swcme"):
        copy_tree(repo / "src" / "models" / name,
                  root / "src" / "models" / name)
    copy_tree(repo / "src" / "general" / "input",
              root / "src" / "general" / "input")
    copy_tree(repo / "srcSEP" / "adapters", root / "srcSEP" / "adapters")
    if include_sep3d:
        # The maintained reduced-front C++ harness identifies the repository
        # by walking upward to this normative task marker.  Copying the real
        # marker gives the disposable tree the same declared identity;
        # inventing an empty sentinel would let a malformed copy pass its
        # path-resolution checks.  It is inert root-discovery input, not a
        # physics oracle.  Non-reduced adapters deliberately do not acquire
        # this unrelated campaign dependency.
        shutil.copy2(str(repo / "CODEX_REDUCED_SHOCK_TASK.txt"), str(root))
        copy_tree(repo / "srcSEP3D", root / "srcSEP3D")
    else:
        copy_tree(repo / "srcSEP3D" / "adapters",
                  root / "srcSEP3D" / "adapters")
    return root


def load_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise RuntimeError("cannot read authoritative JSON {}: {}".format(path, error))
    if not isinstance(value, dict):
        raise RuntimeError("authoritative JSON {} is not an object".format(path))
    return value


def exact_rows(rows: Iterable[Mapping[str, object]], expected: Sequence[str],
               id_key: str = "id", status_key: str = "status"
               ) -> List[Tuple[str, str, str]]:
    """Validate a complete ordered local registry and return normalized rows."""

    normalized: List[Tuple[str, str, str]] = []
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise RuntimeError("local result row {} is not an object".format(index))
        local_id = row.get(id_key)
        status = str(row.get(status_key, "")).upper()
        if not isinstance(local_id, str) or status not in STATUSES:
            raise RuntimeError("malformed local result row {}".format(index))
        normalized.append((local_id, status, str(row.get("message", ""))))
    observed = [row[0] for row in normalized]
    if observed != list(expected):
        raise RuntimeError(
            "local registry drift: expected {} rows {}, observed {} rows {}".format(
                len(expected), list(expected), len(observed), observed))
    return normalized


def global_id(prefix: str, local: str) -> str:
    """Map punctuation-bearing local labels into manifest-safe stable IDs."""

    token = re.sub(r"[^A-Z0-9]+", "-", local.upper()).strip("-")
    if not token:
        raise RuntimeError("empty normalized local identity {!r}".format(local))
    return prefix + "-" + token


def selected_rows(context: Mapping[str, object], rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    """Require the complete central registry, then apply an exact selection.

    Local validation alone cannot prevent a central manifest edit from
    deleting a case: the adapter would still know about the local row but the
    global protocol would never request or report it.  Compare the complete
    ordered globalization against ``expected_cases`` *before* honoring
    ``selected_cases``.  This makes additions, deletions and reorderings a
    deliberate synchronized registry update while retaining efficient
    one-case selection after that full contract has been established.
    """

    configured = list(context["test"]["expected_cases"])  # type: ignore[index]
    observed = [row["id"] for row in rows]
    if configured != observed:
        raise RuntimeError(
            "central/local complete registry drift: configured {} rows {}, "
            "observed {} rows {}".format(
                len(configured), configured, len(observed), observed))
    wanted = list(context["test"].get("selected_cases", configured))  # type: ignore[index]
    by_id = {row["id"]: row for row in rows}
    if len(by_id) != len(rows):
        raise RuntimeError("globalized local registry contains duplicate IDs")
    missing = [case_id for case_id in wanted if case_id not in by_id]
    if missing:
        raise RuntimeError("requested cases absent from imported registry: {}".format(missing))
    return [by_id[case_id] for case_id in wanted]


def suite_status(rows: Sequence[Mapping[str, str]]) -> str:
    if not rows:
        return "ERROR"
    return max((row["status"] for row in rows), key=lambda value: SEVERITY[value])

