"""XM01 linked cross-solver focused-transport entrypoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cross_model_case_runner import run_cross_model_case

def run_case(**kwargs):
    """Run the linked stochastic model and independent finite-volume solver."""
    return run_cross_model_case("XM01", **kwargs)
