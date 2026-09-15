"""IV01 linked Parker-spiral focusing entrypoint."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from integrated_case_runner import run_integrated_case
def run_case(**kwargs): return run_integrated_case("IV01",**kwargs)
