"""CV03 entrypoint: nonuniform Parker diffusion versus finite volume."""

from pathlib import Path
import sys

# Only campaign plumbing is imported; the numerical reference remains a
# separate process and never calls srcSEP's production coefficient provider.
CASES = Path(__file__).resolve().parents[1]
if str(CASES) not in sys.path:
    sys.path.insert(0, str(CASES))
from controlled_case_runner import run_controlled_case


def run_case(**kwargs):
    """Execute CV03 through the supplied linked application and score it."""
    return run_controlled_case("CV03", **kwargs)
