"""CV02 entrypoint: linked Parker diffusion versus the exact Green function."""

from pathlib import Path
import sys

# Dynamic case loading does not install validation/cases as a package.  Add the
# reviewed shared-helper directory explicitly; no production module is loaded.
CASES = Path(__file__).resolve().parents[1]
if str(CASES) not in sys.path:
    sys.path.insert(0, str(CASES))
from controlled_case_runner import run_controlled_case


def run_case(**kwargs):
    """Execute CV02 through the supplied linked application and score it."""
    return run_controlled_case("CV02", **kwargs)
