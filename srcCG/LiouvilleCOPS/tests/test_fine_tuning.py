"""First-visible interval groups, reduced solver, and actual CPU CLI validation.

Synthetic fixtures deliberately isolate support history from gauge dominance.
They exercise the production kernel reader, physical pressure evaluator and
C++ facet-field reader, without mission data or importing the CUDA branch.
"""
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

import numpy as np

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP/"tools"))
import fit_cops
import fit_restart as rs
import fine_tuning as ft
import fine_tune_activity as cli
from cops_io import MODEL_COLUMNS, file_identity, validate_fit_base


def response(faces, ng, rg, structural=None):
    faces = np.array(faces, dtype=int)
    return (faces, np.ones(len(faces), bool) if structural is None else np.array(structural),
            np.array(ng, float), np.array(rg, float), np.zeros(2))


class GroupTests(unittest.TestCase):
    def test_first_seen_not_just_previous_interval_and_gauge_categories(self):
        rows = [response([0], [1], [1]),
                response([1, 2, 3], [10, .001, 1], [.001, 10, 1]),
                response([0, 1, 4], [1, 1, 20], [1, 5, .001])]
        first, atom, ng, rg = ft.prepare_support(6, [-1, 0, 1], [0, 60, 60], rows, 1e-3)
        np.testing.assert_array_equal(first, [-1, 0, 0, 0, 1, -2])
        np.testing.assert_array_equal(atom, [-1, 0, 1, 2, 3, -1])
        # A facet classified NG-dominated in its first bin can have RG
        # response later. Its membership and its later response both persist.
        baseline, matrices = ft.pressure_matrices(atom, 6, np.ones(6), 1, rows)
        self.assertEqual(matrices[2, 1, 0], 5)
        self.assertFalse(np.any(matrices[0]))
        self.assertFalse(np.any(matrices[:2, :, 3:]))

    def test_tiny_early_response_excludes_facet(self):
        rows = [response([0], [1e-100], [0]), response([0, 1], [10, 1], [0, .001])]
        _, atom, _, _ = ft.prepare_support(2, [-1, 0], [0, 1], rows, 1e-3)
        np.testing.assert_array_equal(atom, [-1, 0])
        # epsilon changes the type, never the strict early-zero condition.
        _, changed, _, _ = ft.prepare_support(2, [-1, 0], [0, 1], rows, 1e-4)
        np.testing.assert_array_equal(changed, [-1, 2])

    def test_boundary_weights_and_empty_intervals(self):
        utc = ["2016-09-30T00:%02d:00" % j for j in range(7)]
        bins, weights, edges = ft.intervals(utc, 3, 3)
        np.testing.assert_array_equal(bins, [-1, -1, -1, 0, 1, 2, 2])
        self.assertEqual(weights.sum(), 180)
        for j in range(3):
            self.assertEqual(weights[bins == j].sum(), 60)
        with self.assertRaisesRegex(ValueError, "no modeled"):
            ft.intervals(utc, 3, 20)
        with self.assertRaisesRegex(ValueError, "earlier"):
            ft.intervals(utc, 8, 2)

    def test_all_overlap_redundancy_and_empty_groups(self):
        b, names, omitted, sizes, area = ft.group_basis(np.array([0, 1, 3]), np.ones(3), 2)
        self.assertEqual(names, ["late_000_ng", "late_000_rg", "late_001_ng"])
        self.assertEqual(len(omitted), 3)
        b, names, _, _, _ = ft.group_basis(np.array([0, 1, 2]), np.ones(3), 1)
        np.testing.assert_array_equal(b, [[1, 0, 1], [0, 1, 1], [0, 0, 1]])
        mixed, _, _, _, _ = ft.group_basis(np.array([0, 1, 2]), np.ones(3), 1, "mixed")
        np.testing.assert_array_equal(mixed, np.eye(3))

    def test_recovery_early_invariance_penalty_and_zero_response(self):
        rng = np.random.default_rng(3026)
        design = rng.uniform(.1, 2, (20, 2, 3)); design[:4] = 0
        base = np.ones((20, 2))
        truth = np.array([.2, .5, .3])
        refs = base+design@truth
        late = np.arange(20) >= 4
        x, pred, status = ft.fit_corrections(base, design, refs, np.ones_like(refs, bool), late, np.eye(3))
        np.testing.assert_allclose(x, truth, atol=1e-6)
        np.testing.assert_array_equal(pred[:4], base[:4])
        self.assertEqual(status["data_rank"], 3)
        regularized, _, _ = ft.fit_corrections(base, design, refs, np.ones_like(refs, bool), late, np.eye(3), 100.)
        self.assertLess(np.linalg.norm(regularized), np.linalg.norm(x))
        design[:, :, 2] = 0
        x, _, status = ft.fit_corrections(base, design, refs, np.ones_like(refs, bool), late, np.eye(3))
        self.assertEqual(x[2], 0)
        self.assertEqual(status["zero_response_columns"], [2])
        empty = np.empty((20, 2, 0))
        x, pred, status = ft.fit_corrections(base, empty, refs, np.ones_like(refs, bool), late, np.empty((0, 0)))
        self.assertEqual(len(x), 0)
        np.testing.assert_array_equal(pred, base)


def binary_fixture(root):
    """A small manifested base and frozen C12 state with newly visible faces."""
    run = root/"run"; (run/"kernels").mkdir(parents=True)
    nf = 7
    coef = np.zeros((3, 25)); coef[0, 0] = 2; coef[1, 0] = .2
    coefficients = root/"coefficients.txt"; np.savetxt(coefficients, coef, fmt="%.17e")
    table = np.zeros((nf, 34)); table[:, 0] = np.arange(nf)+1
    table[:, 3] = 1; table[:, 6:8] = 1; table[:, 32] = 2; table[:, 33] = .2
    (run/"elements.bin").write_bytes(struct.pack("<iii", 1, nf, 25)+table.astype("<f8").tobytes())
    header = ["synthetic fine tuning fixture",
              "FlightTimeRotation=none TemperatureModel=analytical SZAMode=continuous",
              "ActivityModel=file ActivityFile=%s HeliocentricModel=shap5" % coefficients,
              "AnalyticalTmax=200 AnalyticalT75=170 AnalyticalTmin=100",
              "TemperatureNG=293 RGPressureModel=equilibrium_sphere SpeciesSensitivity=on (H2O 0.893, CO2 0.704)"]
    # Facet zero contributes throughout; 1,2,3 first occur in bin zero;
    # 4,5,6 in bin one. Early zeros make strict protection testable.
    rows=[]; names=list(MODEL_COLUMNS); cols={n:i for i,n in enumerate(names)}
    for j in range(7):
        face = [0] if j < 3 else ([0, 1, 2, 3] if j == 3 else list(range(nf)))
        records = np.zeros(len(face), dtype=fit_cops.REC); records["face"] = face
        records["cos"] = .9
        for k, f in enumerate(face):
            records["Kn"][k] = 0 if f in (2, 5) else 1e-6*(1+.2*j*f)
            records["KPhi"][k] = 0 if f in (1, 4) else 1e-6*(1+.1*j*(7-f))
        geometry = [float(j), 0., 100., 0., 0., 1., 0., 0., 3*fit_cops.AU, 50.]
        (run/"kernels"/("point=%04d.bin" % j)).write_bytes(
            struct.pack("<ii10di", 1, j, *geometry, len(face))+records.tobytes())
        values = {n:"0" for n in names}
        values.update(iPoint=str(j), UTC="2016-09-30T00:%02d:00" % j,
                      et=str(j), x="100", altitude="50", rh_AU="3", nVisible=str(len(face)),
                      T_RG="293", useNG="1", useRG="1", pNG="1", pRG="1", refNG="1", refRG="1")
        rows.append([values[n] for n in names])
    snapshot = run/"input.snapshot.in"
    snapshot.write_text("ActivityModel = file\nActivityFile = %s\n" % coefficients)
    el = fit_cops.read_elements(run/"elements.bin")
    model = fit_cops.Model(str(run), header, cols, rows, el, coef)
    # Freeze a fit different from the base prior, including T and kappa.
    a=coef[0].copy(); a[0]=3
    temperature=np.array([210., 180., 110.]); kappa=1.3
    ng, rg = model.evaluate(a, temperature, kappa)
    original_ng, original_rg = model.evaluate(coef[0], model.T0, 1.)
    for j, row in enumerate(rows):
        for g, pred, orig in (("NG", ng[j], original_ng[j]), ("RG", rg[j], original_rg[j])):
            row[cols["p"+g]]="%.17e" % orig
            row[cols["ref"+g]]="%.17e" % (pred*(1.1 if j>=3 else 1))
    (run/"COPS_model.dat").write_text("".join("# "+h+"\n" for h in header)+
        "#"+" ".join(names)+"\n"+"".join(" ".join(row)+"\n" for row in rows))
    audit=validate_fit_base(run)
    manifest=dict(schema="liouville_cops.run_manifest", schema_version=1,
                  states={"process":"complete"}, effective_settings={"BodyFrame":["67P/C-G_CK"]},
                  outputs={k:audit[k] for k in ("model_table", "elements", "kernels")})
    (run/"run_manifest.json").write_text(json.dumps(manifest))
    arrays=dict(a=a, coef_prior=coef, T=temperature, kappa=np.array([kappa]),
                harmonic_values=np.empty(0), template_values=np.empty(0), template_columns=np.empty((nf, 0)),
                activity_linear=el["Y"]@a, centroid=el["centroid"], area=el["area"],
                point_ids=np.arange(7), utc=np.array([r[cols["UTC"]] for r in rows]),
                altitude_m=np.full(7,50.), ref_ng=np.array([float(r[cols["refNG"]]) for r in rows]),
                ref_rg=np.array([float(r[cols["refRG"]]) for r in rows]), use_ng=np.ones(7,int), use_rg=np.ones(7,int),
                saturated=np.zeros(7,bool), pred_ng=ng, pred_rg=rg,
                x=np.empty(0), x0=np.empty(0), lo=np.empty(0), hi=np.empty(0))
    meta=dict(mesh={"identity":__import__("activity_templates").mesh_identity(el["centroid"],el["area"])},
              enhancement_scale=3., templates={"names":[],"sources":[]}, harmonics={"modes":[]},
              base={k:file_identity(path) for k,path in (("base_elements",run/"elements.bin"),
                    ("base_model_table",run/"COPS_model.dat"),("base_input",snapshot))})
    state=rs.Restart(meta, arrays); restart=root/"fit_restart.npz"; rs.write(restart,state)
    return run, restart, state, el, model


class PipelineTests(unittest.TestCase):
    def test_cli_reuse_identity_replay_and_cpp_export(self):
        with tempfile.TemporaryDirectory(prefix="cops_fine_tuning_") as tmp:
            root=Path(tmp)
            run, restart, state, el, model=binary_fixture(root)
            args=[str(run),"--restart-from",str(restart),"--minutes","3","--intervals","3"]
            first=root/"first"
            cli.main(args+["--out",str(first)])
            report=json.loads((first/"fine_tuning.json").read_text())
            self.assertEqual(report["state"],"complete")
            self.assertEqual(report["earlier_pressure_max_absolute_change"],0)
            self.assertLess(report["optimizer"]["cost"],report["optimizer"]["initial_cost"])
            with np.load(first/"fine_tuning_state.npz") as z:
                final=np.array(z["activity_final"]); predicted=np.array(z["prediction"])
            # Independent full production Model replay with the same frozen
            # physical parameters, including non-prior temperature and kappa.
            delta=final-el["Y"]@state.a
            ng,rg=model.evaluate(state.a,state.T,float(state.kappa[0]),b=np.ones(1),
                                 enhancement_basis=delta[:,None],enhancement_scale=1.)
            np.testing.assert_allclose(predicted,np.column_stack([ng,rg]),rtol=5e-14)
            self.assertEqual(final[0],3.)
            again=root/"again"
            cli.main(args+["--out",str(again),"--prepared",str(first/"fine_tuning_prepared.npz")])
            with np.load(again/"fine_tuning_state.npz") as z:
                np.testing.assert_array_equal(final,z["activity_final"])
            with self.assertRaisesRegex(ValueError,"identity/settings"):
                cli.main(args+["--out",str(root/"bad"),"--prepared",str(first/"fine_tuning_prepared.npz"),
                               "--dominance-ratio",".01"])
            self.assertFalse((root/"bad").exists())
            prepared=root/"prepared"
            cli.main(args+["--out",str(prepared),"--prepare-only"])
            self.assertFalse((prepared/"fine_tuning_state.npz").exists())
            # Check the existing independent C++ strict facet reader with one
            # exact facet from the generated full export, keeping its identity.
            compiler=shutil.which("g++")
            if compiler:
                probe=root/"probe"
                subprocess.run([compiler,"-std=c++17","-O2","-I",str(APP/"src"),
                                str(APP/"tests/activity_facet_enhancement_probe.cpp"),
                                str(APP/"src/ActivityEnhancement.cpp"),str(APP/"src/ActivityEnhancementBasis.cpp"),
                                "-o",str(probe)],check=True,capture_output=True)
                lines=(first/"ActivityEnhancementFit.txt").read_text().splitlines()
                # Exercise enhanced facets too, not only the unchanged early
                # control. The probe has one local face, so remap its ordinal
                # while retaining the exported centroid, area and residual.
                for f in (0, 1, 3, 6):
                    one=root/("one_%d.txt" % f)
                    selected=[]
                    for line in lines:
                        if line.startswith("FacetCount"):
                            selected.append("FacetCount = 1")
                        elif line.startswith("Facet ="):
                            fields=line.split()
                            if int(fields[2]) == f:
                                fields[2]="0"; selected.append(" ".join(fields))
                        else:
                            selected.append(line)
                    one.write_text("\n".join(selected)+"\n")
                    command=[str(probe),str(one),"67P/C-G_CK","sha256:"+file_identity(first/"ActivityFit.txt")["sha256"]]
                    command += [repr(float(v)) for v in el["centroid"][f]]+["1.0",repr(float(el["Y"][f]@state.a))]
                    proc=subprocess.run(command,check=True,capture_output=True,text=True)
                    self.assertAlmostEqual(float(proc.stdout.split()[-1]),final[f])

    def test_pressure_weights_match_model_for_inactive_faces(self):
        with tempfile.TemporaryDirectory() as tmp:
            run, restart, state, el, model=binary_fixture(Path(tmp))
            header,names,rows=fit_cops.read_model_file(run/"COPS_model.dat")
            cols={n:i for i,n in enumerate(names)}
            responses=list(cli.response_iterator(run,header,cols,rows,el,state))
            activity=np.zeros(len(el["area"]))  # zero baseline must not zero sensitivities
            _,matrix=ft.pressure_matrices(np.arange(len(activity)),len(activity),activity,0,responses)
            self.assertTrue(np.any(matrix[3,:,1]))
            with self.assertRaisesRegex(ValueError,"corrupt"):
                path=Path(tmp)/"bad.npz"
                cli.save_arrays(path,{"x":np.ones(1)},{"definition":{}})
                with np.load(path) as z: data={k:np.array(z[k]) for k in z.files}
                data["x"][0]=2; np.savez_compressed(path,**data)
                cli.read_prepared(path,{})


if __name__ == "__main__":
    unittest.main()
