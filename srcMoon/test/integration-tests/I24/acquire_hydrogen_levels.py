#!/usr/bin/env python3
"""Retrieve the NIST Atomic Spectra Database H I energy levels for I24.

The approved total-capture case needs the excitation energy of each H(n)
shell. This script stores the exact NIST ASD response, in eV with J and
statistical weights, together with its retrieval receipt. It never edits
the response. `prepare_total_capture.py` parses it later. Each acquisition
goes to a new directory; the raw response is immutable once written.
"""
import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
import urllib.request

# NIST ASD levels query: spectrum H I, energies in eV (units=1), CSV output
# (format=2), all configurations, with term, J, g, uncertainty and references.
URL = ("https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=H+I&submit=Retrieve+Data"
       "&units=1&format=2&output=0&page_size=15&multiplet_ordered=0&conf_out=on&term_out=on"
       "&level_out=on&unc_out=1&j_out=on&g_out=on&lande_out=on&perc_out=on&biblio=on&temp=")
RAW = "nist-asd-h-i-levels.csv"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path, help="new directory for the raw response")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    request = urllib.request.Request(URL, headers={"User-Agent": "AMPS-srcMoon-source-acquisition"})
    retrieved = datetime.now(timezone.utc).isoformat()
    with urllib.request.urlopen(request, timeout=60) as response:
        body = response.read()
        status, content_type = response.status, response.headers.get("Content-Type", "")
    if status != 200 or not body.startswith(b"Configuration,Term,J,g,"):
        raise SystemExit("unexpected NIST ASD response; nothing was retained as source")
    (output / RAW).write_bytes(body)
    receipt = dict(url=URL, retrieved_utc=retrieved, http_status=status, content_type=content_type,
                   file=RAW, bytes=len(body), sha256=sha256(body).hexdigest(),
                   source="NIST Atomic Spectra Database, H I energy levels (Kramida et al.)",
                   note="raw response; immutable after retrieval")
    (output / "retrieval.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    sys.exit(main())
