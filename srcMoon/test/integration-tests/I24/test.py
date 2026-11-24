#!/usr/bin/env python3
"""Run the maintained I24 formal predecessor/capability owner.

Supporting scalar/storage checks do not replace the physical particle gate.
Missing kinetic policies remain SKIPPED; no reference or event is fabricated.
"""
import sys
from pathlib import Path
# The aggregate invokes this entry point as a separate program. Its own imports
# and both child CLIs keep bytecode out of the versioned test directories.
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from integration_m3 import main
if __name__=='__main__':raise SystemExit(main())
