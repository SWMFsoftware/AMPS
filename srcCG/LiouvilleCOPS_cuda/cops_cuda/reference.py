"""Explicit binding to the CPU LiouvilleCOPS reference application.

The CUDA application reuses the reference's readers, selection profiles and
physical definitions instead of copying them, so there is exactly one
definition to validate against.  The reference location is chosen
explicitly (argument, then ``LIOUVILLE_COPS_REFERENCE``, then the sibling
``../LiouvilleCOPS`` directory) and is never inferred from the current working
directory.  After import, the module files are checked to lie inside that
tree, which rejects a same-named module that happens to shadow it on
``sys.path``.  Callers record ``tree_identity`` of the bound tree in every
product so a later edit of the reference cannot silently change a comparison.
"""

import importlib
import os
from pathlib import Path
import sys


ENVIRONMENT_VARIABLE = "LIOUVILLE_COPS_REFERENCE"
DEFAULT_ROOT = Path(__file__).resolve().parents[2] / "LiouvilleCOPS"
REQUIRED_FILES = ("tools/fit_cops.py", "tools/cops_io.py",
                  "tools/cops_selection.py", "tools/gravity_cache.py",
                  "tools/build_physical_activity_templates.py",
                  "tools/activity_templates.py", "tools/activity_basis.py")
# Subdirectories of the reference whose content defines the reference
# behaviour.  Generated build products and Python caches are excluded by
# provenance.tree_identity; executables are identified separately.
SOURCE_EXCLUDES = ("build", "examples")

_bound_root = None


def reference_root(path=None):
    """Resolve and verify the reference tree without importing it."""
    if path is None:
        path = os.environ.get(ENVIRONMENT_VARIABLE) or DEFAULT_ROOT
    root = Path(path).expanduser().resolve()
    missing = [name for name in REQUIRED_FILES if not (root / name).is_file()]
    if missing:
        raise FileNotFoundError("not a LiouvilleCOPS reference tree (%s missing): %s" %
                                (", ".join(missing), root))
    return root


def bind(path=None):
    """Put ``<root>/tools`` first on ``sys.path`` exactly once and return root.

    Binding a second, different root in one interpreter is refused: Python's
    module cache would otherwise mix two reference versions.
    """
    global _bound_root
    root = reference_root(path)
    if _bound_root is not None:
        if _bound_root != root:
            raise RuntimeError("reference already bound to %s; cannot rebind to %s" %
                               (_bound_root, root))
        return root
    tools = str(root / "tools")
    if tools in sys.path:
        sys.path.remove(tools)
    sys.path.insert(0, tools)
    _bound_root = root
    return root


def load(name, path=None):
    """Import one reference tool module and prove where it came from."""
    root = bind(path)
    module = importlib.import_module(name)
    source = Path(module.__file__).resolve()
    if root / "tools" not in source.parents:
        raise ImportError("module %s resolved to %s, outside the reference tree %s" %
                          (name, source, root))
    return module
