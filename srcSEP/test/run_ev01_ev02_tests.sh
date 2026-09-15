#!/bin/sh
# Dependency-light EV01/EV02 registration, observational-provenance, and syntax gate.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
out=${SEP_TEST_OUTPUT_DIR:-/tmp/srcsep-ev01-ev02-test}
mkdir -p "$out"
export PYTHONPYCACHEPREFIX="$out/pycache"
python3 -m py_compile "$root/validation/cases/campaign_evidence_runner.py" \
  "$root/validation/cases/EV01/case.py" "$root/validation/cases/EV02/case.py"
python3 -m json.tool "$root/validation/case_registry.json" >/dev/null
for case_id in EV01 EV02; do
  python3 -m json.tool "$root/validation/cases/$case_id/input.json" >/dev/null
  python3 -m json.tool "$root/validation/cases/$case_id/reference/provenance.json" >/dev/null
  grep -q '"synthetic_reference"[[:space:]]*:[[:space:]]*false' "$root/validation/cases/$case_id/input.json"
done
python3 "$root/validation/run_case.py" --list | grep -q 'EV01 | campaign-evidence'
python3 "$root/validation/run_case.py" --list | grep -q 'EV02 | campaign-evidence'

# Regression guard for the native family dispatcher.  EV01/EV02 share the
# observational-profile implementation in RunCrossModelValidationModel().  A
# missing EV branch here makes the registry callbacks fall through to the
# CV06+ advanced-validation dispatcher and the linked executable exits 2.
python3 - <<'PY_DISPATCH' "$root"
from pathlib import Path
import sys
text=(Path(sys.argv[1])/'component_tests.cpp').read_text(encoding='utf-8')
needle='identifier.substr(0, 2) == "EV"'
assert needle in text, 'EV native dispatcher branch is missing from component_tests.cpp'
assert 'RunCrossModelValidationModel' in text, 'cross-model/observational dispatcher is missing'
print('PASS EV01-EV02-DISPATCH: EV cases route to RunCrossModelValidationModel')
PY_DISPATCH
python3 - <<'PY' "$root"
import csv,sys
from pathlib import Path
root=Path(sys.argv[1]); rows=list(csv.DictReader(open(root/'validation/cases/EV01/reference/ccmc_2021_goes_observations.csv')))
assert len({r['event_id'] for r in rows})==9
assert any(r['threshold_mev']=='100' and not r['crossing_utc'] for r in rows)
assert all(r['source']=='CCMC-SEPVAL-2021' for r in rows)
print('PASS EV01-EV02-OBS: nine real CCMC events; observed non-crossings retained')
PY

# 2012-07-12 is a real CCMC edge case: flare peak and flare end are both
# 16:49 UTC.  The old writer emitted duplicate source times, which the native
# parser correctly rejected with "source times must be strictly increasing".
# Generate every event source and prove the serialized times are monotone.
python3 - <<'PY_SOURCE' "$root" "$out"
import ast,csv,sys,tempfile
from pathlib import Path
root=Path(sys.argv[1]); out=Path(sys.argv[2])
sys.path.insert(0,str(root/'validation/cases'))
import campaign_evidence_runner as ev
# Python 3.8 and older do not provide math.nextafter.  The EV source writer
# must therefore use the local IEEE-754 compatibility helper rather than the
# version-dependent math API.
source_text=(root/'validation/cases/campaign_evidence_runner.py').read_text(encoding='utf-8')
tree=ast.parse(source_text)
assert not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id=='math' and n.func.attr=='nextafter' for n in ast.walk(tree))
assert ev._next_float_toward_positive_infinity(1.2) > 1.2
assert 0.0 < ev._next_float_toward_positive_infinity(0.0) < 1.0e-300
rows=list(csv.DictReader(open(root/'validation/cases/EV01/reference/ccmc_2021_goes_observations.csv')))
grouped={}
for row in rows: grouped.setdefault(row['event_id'],[]).append(row)
for event_id,event_rows in grouped.items():
    path=out/(event_id+'_source.csv')
    ev._write_source(event_rows,path)
    generated=list(csv.DictReader(open(path)))
    times=[float(r['elapsed_hours']) for r in generated]
    assert all(b>a for a,b in zip(times,times[1:])), (event_id,times)
# Specifically retain the observed equality while serializing a right-hand
# shutdown at the next representable float rather than inventing a finite tail.
r2012=grouped['20120712'][0]
assert r2012['flare_peak_utc']==r2012['flare_end_utc']
generated=list(csv.DictReader(open(out/'20120712_source.csv')))
times=[float(r['elapsed_hours']) for r in generated]
assert times[2] > times[1] and (times[2]-times[1]) < 1.0e-12
print('PASS EV01-EV02-SOURCE: portable source timestamps strictly increase for all events, including 2012-07-12 peak=end')
PY_SOURCE
echo "PASS EV01-EV02-UNIT: registry, fixed observation provenance, source serialization, and campaign runner syntax"
