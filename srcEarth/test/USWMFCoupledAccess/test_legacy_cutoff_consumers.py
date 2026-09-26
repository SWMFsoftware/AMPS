#!/usr/bin/env python3
"""Regression tests for Step-10 AUXDATA consumed by the C9/C10 validators.

The production writer and live/replay comparator were tested when Step 10 was
introduced, but the older observational consumers were not.  This test loads
the real C9 and C10 parser functions and supplies the exact ordering emitted by
the production writer: TITLE, one or more AUXDATA records, VARIABLES, ZONE, and
numeric rows.  It also verifies that accepting AUXDATA does not turn the parser
into a permissive text filter: malformed metadata and inconsistent access flags
must still fail closed.
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path
from types import ModuleType
from typing import Callable


HERE = Path(__file__).resolve().parent
EARTH = HERE.parents[1]


def load_runner(name: str, path: Path) -> ModuleType:
    """Load a runner under a unique name without changing its source directory."""
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise AssertionError("cannot load validation runner: %s" % path)
    module = importlib.util.module_from_spec(specification)
    # Dataclass construction consults sys.modules while the file is executing.
    # Register the unique module name first, exactly as Python's import machinery does.
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


def require_value_error(action: Callable[[], object], label: str) -> None:
    """Require a fail-closed parser outcome for deliberately damaged input."""
    try:
        action()
    except ValueError:
        return
    raise AssertionError("parser unexpectedly accepted %s" % label)


ACCESS_PRODUCT = '''TITLE="manufactured direct-access product"
AUXDATA SNAPSHOT_ID="field-v1-reference"
AUXDATA SNAPSHOT_EPOCH_UTC="2006-12-14T09:49:00.000000000Z"
AUXDATA SNAPSHOT_MESH_REVISION="mesh-v1-reference"
AUXDATA SNAPSHOT_CONTENT_FINGERPRINT="state-v1-reference"
AUXDATA OUTER_BOUNDARY_POLICY="BOX"
VARIABLES="lon_deg" "lat_deg" "rigidity_gv" "access_state" "allowed" "unresolved"
ZONE T="shell" I=2 F=POINT
0 50 0.2 0 0 0
0 70 0.2 1 1 0
'''


# C9 needs the full penumbra diagnostic set, whereas C10 consumes the five
# cutoff coordinates.  One superset fixture exercises both real parser schemas.
PENUMBRA_PRODUCT = '''TITLE="manufactured penumbra product"
AUXDATA SNAPSHOT_ID="field-v1-reference"
AUXDATA SNAPSHOT_EPOCH_UTC="2006-12-14T09:49:00.000000000Z"
AUXDATA OUTER_BOUNDARY_POLICY="BOX"
VARIABLES="lon_deg" "lat_deg" "Rc_lower_GV" "Rc_effective_GV" "Rc_upper_GV" "n_allowed_intervals" "n_transitions" "n_unresolved" "lower_bracket_unresolved" "upper_bracket_unresolved" "lower_below_range" "lower_above_range" "upper_below_range" "upper_above_range"
ZONE T="shell" I=1 F=POINT
0 60 0.4 0.5 0.6 1 2 0 0 0 0 0 0 0
'''


def main() -> int:
    c9 = load_runner("step10_regression_c9", EARTH / "test" / "C9" / "run_C9.py")
    c10 = load_runner("step10_regression_c10", EARTH / "test" / "C10" / "run_C10.py")

    with tempfile.TemporaryDirectory(prefix="step10_legacy_consumers_") as directory:
        root = Path(directory)
        access = root / "cutoff_3d_shells_access.dat"
        penumbra = root / "cutoff_3d_shells_penumbra.dat"
        access.write_text(ACCESS_PRODUCT, encoding="utf-8")
        penumbra.write_text(PENUMBRA_PRODUCT, encoding="utf-8")

        for label, module in (("C9", c9), ("C10", c10)):
            access_rows = module.parse_tecplot_shell_access(access)
            shell_rows = module.parse_tecplot_shell_penumbra(penumbra)
            if [row.access_state for row in access_rows] != [0, 1]:
                raise AssertionError("%s changed the manufactured access states" % label)
            if len(shell_rows) != 1 or shell_rows[0].rc_effective_gv != 0.5:
                raise AssertionError("%s changed the manufactured cutoff row" % label)

        print("PASS S10-U09 C9/C10 consume validated Step-10 AUXDATA")

        malformed = root / "malformed_auxdata.dat"
        malformed.write_text(
            ACCESS_PRODUCT.replace(
                'AUXDATA SNAPSHOT_ID="field-v1-reference"',
                'AUXDATA SNAPSHOT_ID=field-v1-reference',
            ),
            encoding="utf-8",
        )
        inconsistent = root / "inconsistent_access_state.dat"
        inconsistent.write_text(
            ACCESS_PRODUCT.replace("0 70 0.2 1 1 0", "0 70 0.2 1 0 0"),
            encoding="utf-8",
        )
        for label, module in (("C9", c9), ("C10", c10)):
            require_value_error(
                lambda module=module: module.parse_tecplot_shell_access(malformed),
                "%s malformed AUXDATA" % label,
            )
            require_value_error(
                lambda module=module: module.parse_tecplot_shell_access(inconsistent),
                "%s inconsistent state flags" % label,
            )

        print("PASS S10-U10 malformed metadata/state still fail closed")

    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
