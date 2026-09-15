"""XM02 linked M-FLAMPA Parker-spiral comparison entrypoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cross_model_case_runner import run_cross_model_case

def run_case(**kwargs):
    return run_cross_model_case("XM02", **kwargs)
