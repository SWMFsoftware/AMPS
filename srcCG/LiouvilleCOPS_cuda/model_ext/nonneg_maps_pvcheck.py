# ==== version-controlled copy: header added 2026-11-12, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_maps_pvcheck.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_maps_pvcheck.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 3c328fc63eddc51b4c649dac074ac6398ede6535fda206179fa3f049910d06bf
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool. Read-back check of the written maps with ParaView's own readers (pvpython): every .vtu and the .pvd open and give
#   the expected arrays and values (the cluster's ParaView 5.11.1 cannot read appended VTK XML data; the maps use
#   inline base64).
# USED BY: nonneg_maps.sbatch (verify step).
# ==== end of header ====
# nonneg_maps_pvcheck.py -- read a nonneg_maps directory with ParaView's own readers (run with /usr/bin/pvpython):
# the static .vtu, every per-time .vtu through the .pvd collection, and compare cell / point counts and the
# per-array min / max / sum with maps_manifest.json. Prints one PASS/FAIL line per check; exit 1 on a failure.
# usage: pvpython nonneg_maps_pvcheck.py <maps dir>
import json
import sys
from pathlib import Path

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader

d = Path(sys.argv[1])
man = json.load(open(d / "maps_manifest.json"))
model = man["model"]
fails = 0


def rd(path):
    r = vtkXMLUnstructuredGridReader()
    r.SetFileName(str(path))
    r.Update()
    return r.GetOutput()


def cmp(tag, arr, ref):
    global fails
    a = vtk_to_numpy(arr).astype(np.float64)
    got = (float(a.min()), float(a.max()), float(a.sum()))
    ok = got[0] == ref["min"] and got[1] == ref["max"] and abs(got[2] - ref["sum"]) <= 1e-12 * max(abs(ref["sum"]), 1e-300)
    fails += 0 if ok else 1
    print("%s %s min %.17g max %.17g sum %.17g (manifest %.17g %.17g %.17g)" % ("PASS" if ok else "FAIL", tag, got[0],
                                                                               got[1], got[2], ref["min"], ref["max"],
                                                                               ref["sum"]))


g = rd(d / ("%s_static.vtu" % model))
n = man["static_stats"]["A_H2O"]["n"]
ok = g.GetNumberOfCells() == n and g.GetCellType(0) == 5
fails += 0 if ok else 1
print("%s static: %d points, %d cells" % ("PASS" if ok else "FAIL", g.GetNumberOfPoints(), g.GetNumberOfCells()))
for name, key in (("A_H2O", "A_H2O"), ("CO2_activity_kappa_ACO2", "CO2_activity_kappa_ACO2"), ("seen", "seen"),
                  ("slope_bin", "slope_bin"), ("area_m2", "area_m2")):
    cmp("static " + name, g.GetCellData().GetArray(name), man["static_stats"][key])
for t in man["times"]:
    g = rd(d / t["vtu"])
    tv = g.GetFieldData().GetArray("TimeValue").GetValue(0)
    ok = g.GetNumberOfCells() == n and tv == t["TimeValue_s"]
    fails += 0 if ok else 1
    print("%s %s %s: %d cells, TimeValue %g" % ("PASS" if ok else "FAIL", t["label"], t["utc"], g.GetNumberOfCells(), tv))
    for name, key in (("cos_SZA", "cos_SZA"), ("F_H2O_per_m2_s", "F_H2O"), ("F_CO2_per_m2_s", "F_CO2")):
        cmp("%s %s" % (t["label"], name), g.GetCellData().GetArray(name), t["stats"][key])
from paraview import simple  # noqa: E402
s = simple.PVDReader(FileName=str(d / ("%s_flux.pvd" % model)))
steps = list(s.TimestepValues)
ok = steps == [t["TimeValue_s"] for t in man["times"]]
fails += 0 if ok else 1
print("%s pvd: %d time steps %s" % ("PASS" if ok else "FAIL", len(steps), steps))
print("nonneg_maps_pvcheck: %d failures" % fails)
sys.exit(1 if fails else 0)
