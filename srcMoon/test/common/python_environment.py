"""Python interpreter requirements and provenance for srcMoon test runners.

The maintained interpreter on this machine is /usr/bin/python3, selected by
test/qualification-machine/passion.json and documented in srcMoon/RUN.md
("Python environment for this machine"). A shell may put another python3
first on PATH (for example an unrelated virtual environment).

A missing third-party package is a missing dependency, so the status is
SKIPPED (AGENTS.md status semantics), never a late ERROR. Runners call
``require`` before any compilation, MPI launch or other native work; the
first failure names the interpreter and the maintained command. ``describe``
returns the interpreter identity recorded in run manifests.

Standard library only; importing this module never changes global state.
"""
import importlib
import platform
import sys

MAINTAINED_INTERPRETER = "/usr/bin/python3"
RUN_GUIDE = "srcMoon/RUN.md, section 'Python environment for this machine'"
# Packages named in RUN.md for the numerical references and PNG/EPS figures.
NUMERICAL_PACKAGES = ("numpy", "scipy", "matplotlib", "mpmath")


def absent(modules):
    """Return the requested top-level modules that cannot be imported."""
    missing = []
    for name in modules:
        try:
            importlib.import_module(name)
        except ImportError:
            missing.append(name)
    return missing


def message(missing):
    return ("Python package(s) %s unavailable to %s (Python %s); run with %s (see %s)"
            % (", ".join(missing), sys.executable, platform.python_version(), MAINTAINED_INTERPRETER, RUN_GUIDE))


def require(modules, error_type):
    """Raise ``error_type`` (the caller's SKIPPED exception) if any module is missing."""
    missing = absent(modules)
    if missing:
        raise error_type(message(missing))


def describe(modules=NUMERICAL_PACKAGES):
    """Interpreter path, version and package versions (None when unavailable)."""
    packages = {}
    for name in modules:
        try:
            packages[name] = str(getattr(importlib.import_module(name), "__version__", "unknown"))
        except ImportError:
            packages[name] = None
    return {"executable": sys.executable, "version": platform.python_version(),
            "implementation": platform.python_implementation(), "packages": packages}
