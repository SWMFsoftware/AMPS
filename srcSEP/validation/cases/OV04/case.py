"""OV04 linked 2014-01-06 connectivity-sensitivity entrypoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from observational_case_runner import run_observational_case

def run_case(**kwargs):
    return run_observational_case("OV04", **kwargs)
