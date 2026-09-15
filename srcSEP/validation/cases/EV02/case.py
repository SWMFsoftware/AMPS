"""EV02 sealed held-out validation campaign entrypoint."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from campaign_evidence_runner import run_campaign_case
def run_case(**kwargs): return run_campaign_case("EV02",**kwargs)
