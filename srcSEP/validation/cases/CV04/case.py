"""CV04 entrypoint: linked adiabatic cooling versus exact characteristics."""

from pathlib import Path
import sys

# Keep the case compatible with validation/run_case.py's path-based loader
# without requiring users to install srcSEP as a Python package.
CASES = Path(__file__).resolve().parents[1]
if str(CASES) not in sys.path:
    sys.path.insert(0, str(CASES))
from controlled_case_runner import run_controlled_case


def run_case(**kwargs):
    """Execute CV04 through the supplied linked application and score it."""
    return run_controlled_case("CV04", **kwargs)
