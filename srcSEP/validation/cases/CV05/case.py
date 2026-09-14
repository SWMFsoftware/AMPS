"""CV05 entrypoint: linked magnetic focusing versus exact characteristics."""

from pathlib import Path
import sys

# The shared layer standardizes execution/evidence only. CV05's expected
# physics stays in its case-local independent reference_solution.py process.
CASES = Path(__file__).resolve().parents[1]
if str(CASES) not in sys.path:
    sys.path.insert(0, str(CASES))
from controlled_case_runner import run_controlled_case


def run_case(**kwargs):
    """Execute CV05 through the supplied linked application and score it."""
    return run_controlled_case("CV05", **kwargs)
