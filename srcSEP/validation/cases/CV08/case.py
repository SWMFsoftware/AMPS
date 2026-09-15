"""CV08 entrypoint: linked first passage versus inverse-Gaussian theory."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from advanced_case_runner import run_advanced_case
def run_case(**kwargs):
    """Execute the production-linked model, independent reference, and gates."""
    return run_advanced_case("CV08", **kwargs)
