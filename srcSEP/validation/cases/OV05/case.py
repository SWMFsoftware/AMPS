"""OV05 linked September-2017 compound-event diagnostic entrypoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from observational_case_runner import run_observational_case

def run_case(**kwargs):
    return run_observational_case("OV05", **kwargs)
