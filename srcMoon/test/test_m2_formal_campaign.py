#!/usr/bin/env python3
"""I19 dispatch/provenance refusal tests, NOT native or observational evidence.

All artificial JSON/CSV below exists only in TemporaryDirectory and is labelled
as a control-plane fixture. No fixture supplies a physical force or observed
value; no test asserts formal scientific PASS for those files. In particular,
the false-green tests stop before a native execution or independent science
score could be substituted by a cached summary.
"""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import m2_formal_campaign as formal
from m2_protocol import cases


class FormalCampaignContract(unittest.TestCase):
    def fixture(self, root, action="import"):
        """Build self-consistent metadata solely for individual parser checks."""
        exe, application = root / "not-an-executable", root / "input.in"
        exe.write_text("PROTOCOL FIXTURE - NEVER EXECUTE\n")
        contract = root / "contract.json"
        application.write_text("m2_formal_contract_file = " + str(contract) + "\n")
        record = lambda p: {"path": str(p), "sha256": formal._sha(p)}
        value = {"schema": formal.SCHEMA, "status": "APPROVED", "action": action,
            "review": {"approved": True, "reviewer": "PROTOCOL FIXTURE",
                       "approved_utc": "2026-01-01T00:00:00Z", "frozen_before_execution": True},
            "protocol_sha256": formal._sha(formal.PROTOCOL),
            "qualification_root": str(root), "input": record(application), "executable": record(exe),
            "production_source_signature_sha256": hashlib.sha256(b"{}").hexdigest(),
            "execution": {"mpi_ranks": 8, "workers": 1, "allocated_cores": 8, "core_hours": 1,
                          "memory_limit_bytes": 1000, "artifact_limit_bytes": 1000,
                          "host_free_floor_bytes": 0, "disk_free_floor_bytes": 0}}
        for key in ("sites_file", "forecast_file", "normalizations_file", "statistical_support_file"):
            p = root / (key + ".json")
            p.write_text('{"fixture":"not scientific evidence"}\n')
            value[key] = record(p)
        contract.write_text(json.dumps(value))
        return exe, application, contract, value

    def test_missing_contract_skips_before_preflight_or_model_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe, application = root / "not-an-executable", root / "input.in"
            exe.write_text("CONTROL FIXTURE\n")
            application.write_text("m2_environment = enabled\n")
            with patch("run_m2_cluster.preflight", side_effect=AssertionError("must not inspect later evidence")):
                result = formal.run_formal_i19(exe, application, root)
            self.assertEqual(result["status"], "SKIPPED")
            self.assertFalse((root / "I19/comparison.csv").exists())
            self.assertFalse((root / "I19/receipt.txt").exists())
            manifest = json.loads((root / "I19/run_manifest.json").read_text())
            self.assertEqual(manifest["execution"], "not-run")

    def test_duplicate_input_contract_setting_is_error_not_missing_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe, application = root / "not-an-executable", root / "input.in"
            exe.write_text("CONTROL FIXTURE\n")
            application.write_text("m2_formal_contract_file = first.json\nm2_formal_contract_file = second.json\n")
            result = formal.run_formal_i19(exe, application, root)
            self.assertEqual(result["status"], "ERROR")
            self.assertIn("duplicate", result["reason"])

    def test_unreviewed_contract_stays_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe, application, path, value = self.fixture(root)
            value["review"]["approved"] = False
            path.write_text(json.dumps(value))
            with self.assertRaises(formal.Unavailable):
                formal.load_contract(path, exe, application, root, 8)

    def test_changed_input_bytes_are_error_before_predecessor_access(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe, application, path, value = self.fixture(root)
            application.write_text(application.read_text() + "! changed control fixture\n")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                formal.load_contract(path, exe, application, root, 8)

    def test_duplicate_and_nonfinite_json_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "contradictory.json"
            for content in ('{"status":"APPROVED","status":"DRAFT"}', '{"limit":NaN}'):
                path.write_text(content)
                with self.assertRaises(ValueError):
                    formal._json(path)

    def test_same_campaign_identity_is_mandatory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe, application, path, value = self.fixture(root)
            with self.assertRaisesRegex(ValueError, "campaign identity"):
                formal.load_contract(path, exe, application, root / "another", 8)

    def test_layout_and_resource_controls_are_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe, application, path, value = self.fixture(root)
            with self.assertRaisesRegex(ValueError, "MPI layout"):
                formal.load_contract(path, exe, application, root, 4)
            value["execution"]["allocated_cores"] = 1
            path.write_text(json.dumps(value))
            with self.assertRaisesRegex(ValueError, "simultaneous MPI"):
                formal.load_contract(path, exe, application, root, 8)

    def test_native_source_signature_is_required(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe, application, path, value = self.fixture(root)
            with patch("build_m2_qualification.production_source_identity", return_value={"changed.cpp": "x"}):
                with self.assertRaisesRegex(ValueError, "source signature"):
                    formal.load_contract(path, exe, application, root, 8)

    def test_execute_action_does_not_implicitly_authorize_long_job(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe, application, path, value = self.fixture(root, "execute")
            gate = {"status": "PASS", "gaps": [], "errors": [], "evidence": {}, "documents": {}}
            with patch("build_m2_qualification.production_source_identity", return_value={}), \
                 patch("run_m2_cluster.preflight", return_value=gate), \
                 patch("campaign_budget.CampaignBudget.run", side_effect=AssertionError("long execution forbidden")):
                result = formal.run_formal_i19(exe, application, root)
            self.assertEqual(result["status"], "SKIPPED")
            self.assertIn("explicit CLI", result["reason"])
            self.assertFalse((root / "I19/comparison.csv").exists())

    def test_real_predecessor_gap_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe, application, path, value = self.fixture(root)
            gate = {"status": "SKIPPED", "gaps": ["I13 unavailable"], "errors": [], "evidence": {}, "documents": {}}
            with patch("build_m2_qualification.production_source_identity", return_value={}), \
                 patch("run_m2_cluster.preflight", return_value=gate):
                result = formal.run_formal_i19(exe, application, root)
            self.assertEqual(result["status"], "SKIPPED")
            self.assertIn("I13 unavailable", result["reason"])

    def test_cached_pass_without_native_cases_is_not_qualification(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "cluster_summary.json").write_text('{"execution":"executed","status":"PASS"}')
            (root / "job_plan.json").write_text(json.dumps({"cases": cases("cluster"),
                                                        "protocol_sha256": formal._sha(formal.PROTOCOL)}))
            (root / "preflight.json").write_text('{"status":"PASS","evidence":{}}')
            with self.assertRaisesRegex(formal.Unavailable, "selected long I19 cases"):
                formal.evaluate_completed(root, {}, {}, {"evidence": {}}, root / "no-native", 8)

    def test_import_rejects_changed_predecessor_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "cluster_summary.json").write_text('{"execution":"executed"}')
            (root / "job_plan.json").write_text(json.dumps({"cases": cases("cluster"),
                                                        "protocol_sha256": formal._sha(formal.PROTOCOL)}))
            (root / "preflight.json").write_text('{"status":"PASS","evidence":{"I13":{"result_sha256":"old"}}}')
            with self.assertRaisesRegex(ValueError, "evidence differs: I13"):
                formal.evaluate_completed(root, {}, {}, {"evidence": {"I13": {"result_sha256": "new"}}}, root / "no-native", 8)

    def test_missing_case_cannot_reduce_the_predeclared_family(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "cluster_summary.json").write_text('{"execution":"executed"}')
            plan = cases("cluster")[:-1]
            (root / "job_plan.json").write_text(json.dumps({"cases": plan,
                                                        "protocol_sha256": formal._sha(formal.PROTOCOL)}))
            (root / "preflight.json").write_text('{"status":"PASS","evidence":{}}')
            with self.assertRaisesRegex(ValueError, "288-case"):
                formal.evaluate_completed(root, {}, {}, {"evidence": {}}, root / "no-native", 8)

    def test_existing_formal_result_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe, application, path, value = self.fixture(root)
            old = root / "I19/result.json"
            old.parent.mkdir()
            old.write_text('{"preserve":"original diagnostic"}')
            digest = formal._sha(old)
            with self.assertRaisesRegex(ValueError, "overwrite"):
                formal.run_formal_i19(exe, application, root)
            self.assertEqual(formal._sha(old), digest)

    def test_native_seed_cannot_be_relabelled_by_a_cached_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe = root / "not-an-executable"
            exe.write_text("PROTOCOL FIXTURE - NEVER EXECUTE")
            case = cases("cluster")[0]
            for name in ("result.json", "run_manifest.json", "receipt.txt", "application.in",
                         "stdout.log", "stderr.log", "m2_budget.csv", "window_budget.csv"):
                (root / name).write_text("PROTOCOL FIXTURE - NOT NATIVE SCIENCE\n")
            (root / "result.json").write_text('{"status":"PASS","fixture":"control plane only"}')
            (root / "run_manifest.json").write_text(json.dumps({
                "execution_scope": "cluster-numerical-design", "case": case,
                "test_id": "I19", "mpi_ranks": 8, "openmp_threads": 1,
                "random_seed": case["seed"], "executable": {"sha256": formal._sha(exe)}}))
            (root / "receipt.txt").write_text("test_id=I19\nmpi_ranks=8\nseed=1\n")
            with self.assertRaisesRegex(ValueError, "stochastic seed"):
                formal._case_bytes(root, case, exe, 8, set())

    def test_comparison_writer_preserves_arithmetic_rows_without_plotting(self):
        # Parser/serialization fixture only. These values are mathematical
        # interval endpoints, not observations or native AMPS output.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            formal.write_comparisons({"comparisons": [{"mean": .002, "half_width": .001}]}, root, render=False)
            import csv
            with (root / "comparison.csv").open() as stream:
                row = next(csv.DictReader(stream))
            self.assertEqual(float(row["reference_value"]), 0)
            self.assertEqual(float(row["model_value"]), .002)
            self.assertEqual(float(row["uncertainty_value"]), .001)
            self.assertFalse((root / "comparison.png").exists())


class AppliedInitialState(unittest.TestCase):
    """Tiny labelled state-grid/metadata fixtures, NEVER thermal validation."""

    def fixture(self, root, mode="declared-uniform"):
        import m2_initial_state_evidence as state
        from m2_initial_state import physics_identity_sha256
        root=Path(root);pilot=root/"PROTOCOL-PILOT-NOT-A-NATIVE-RUN";pilot.mkdir()
        from integration_common import SPICE_FILES
        spice=root/'SPICE-METADATA-FIXTURE';external={}
        for name in SPICE_FILES:
            path=spice/'Kernels'/name;path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text('CONTROL METADATA ONLY: NOT A VALID SPICE KERNEL\n');external[name]=state.sha256(path)
        for name in ('NOT-A-DEM.IMG','NOT-A-LABEL.LBL'):
            path=root/name;path.write_text('CONTROL METADATA ONLY: NOT LUNAR SCIENCE\n');external[str(path)]=state.sha256(path)
        application=pilot/"application.in"
        application.write_text('#begin application: moon\n'
            'spice_path = '+str(spice)+'\nsurface_geometry = lola\n'
            'lola_image_file = '+str(root/'NOT-A-DEM.IMG')+'\nlola_label_file = '+str(root/'NOT-A-LABEL.LBL')+'\n'
            'surface_temperature_mode = thermal-model\n'
            'thermal_initialization_mode = '+mode+'\nthermal_initial_temperature_k = 40\n'
            'simulation_start_utc = 2000-01-01T12:00:00Z\n'
            'thermal_spinup_temperature_tolerance_k = 0.01\nthermal_spinup_energy_tolerance_j_m2 = 10000\n'
            'm2_verification_record_initial_thermal_state = true\n'
            'm2_verification_cohort_packets = 32\nm2_verification_source_integral_tolerance = 1e-8\n'
            'm2_verification_surface_elements = 0,1\n#end application\n')
        executable=root/"PROTOCOL-NOT-AN-EXECUTABLE";executable.write_text('CONTROL FIXTURE ONLY\n')
        grid=pilot/("receipt.txt"+state.STATE_SUFFIX)
        grid.write_text(','.join(state.STATE_HEADER)+'\n'+''.join(
            f'{element},{layer},{depth},40\n' for element in range(2)
            for layer,depth in ((-1,0),(0,.1),(1,.3))))
        receipt=pilot/"receipt.txt"
        receipt.write_text('test_id=I19\nmpi_ranks=8\ninitial_epoch_tdb_s=123\n'
            'thermal_initialization_mode='+mode+'\nthermal_initialization_et_s=123\n'
            'initial_thermal_state_csv='+str(grid)+'\ninitial_thermal_state_rows=6\n'
            'initial_thermal_state_element_count=2\ninitial_thermal_state_layer_count=2\n'
            'initial_na_cohort_atoms_per_site=64\ninitial_ar_cohort_atoms_per_site=96\n'
            'thermal_initialization_spinup_completed_cycles=2\n'
            'thermal_initialization_spinup_maximum_node_temperature_change_k=0.001\n'
            'thermal_initialization_spinup_maximum_column_enthalpy_change_j_m2=1\n'
            'thermal_initialization_spinup_cycle_start_et_s=122\n'
            'thermal_initialization_spinup_cycle_end_et_s=123\n')
        with receipt.open('a') as stream:
            stream.write('thermal_initialization_spice_kernel_count='+str(len(SPICE_FILES))+'\n'+''.join(
                'thermal_initialization_spice_kernel_'+str(index)+'_file='+str(spice/'Kernels'/name)+'\n'
                for index,name in enumerate(SPICE_FILES)))
        manifest=pilot/"run_manifest.json"
        manifest.write_text(json.dumps({'test_id':'I19','mpi_ranks':8,'openmp_threads':1,
            'started_utc':'2026-10-31T00:00:00Z','executable':{'sha256':state.sha256(executable)},
            'inputs':{str(application):state.sha256(application)},
            'external_data':external,
            'physics_identity_sha256':physics_identity_sha256(application)}))
        evidence=root/"METHOD-SOURCE-CONTROL-NOT-SCIENCE.txt"
        evidence.write_text('PROTOCOL FIXTURE ONLY: no native evolution, observations or qualification\n')
        descriptor=lambda path:{'path':str(path),'sha256':state.sha256(path)}
        q={'schema':state.SCHEMA,'status':'READY',
           'review':{'approved':True,'reviewer':'PROTOCOL FIXTURE ONLY',
                     'approved_utc':'2026-11-01T00:00:00Z','frozen_before_numerical_campaign':True},
           'initialization_mode':mode,'simulation_start_utc':'2000-01-01T12:00:00Z',
           'ephemeris_time_tdb_s':123,'physics_identity_sha256':physics_identity_sha256(application),
           'executable_sha256':state.sha256(executable),'production_source_signature_sha256':'b'*64,
           'pilot_input':descriptor(application),'pilot_manifest':descriptor(manifest),
           'pilot_receipt':descriptor(receipt),'native_state':descriptor(grid),
           'state_shape':{'logical_surface_element_count':2,'layer_count':2},
           'method_qualification':{'status':'READY','criteria':{'fixture_only':'not scientific evidence'},
                                   'evidence_files':[descriptor(evidence)]},
           'source_evidence_files':[descriptor(evidence)]}
        return {'actual_epoch_thermal_qualification':q},application,executable,pilot

    def qualify(self, root, fixture):
        import m2_initial_state_evidence as state
        support,application,executable,pilot=fixture
        with patch.object(state,'_source_signature',return_value='b'*64):
            return state.bind_frozen_sites(state.load_qualification(support,application,executable,base=root),
                                          [{'logical_surface_element':0},{'logical_surface_element':1}])

    def weights(self, application, receipt):
        import m2_initial_state_evidence as state
        return {'method':'BOOKKEEPING PROTOCOL FIXTURE ONLY','cohort_packets':32,
                'initial_atoms_per_site':{'Na':64,'Ar':96},'reference_weight_atoms':{'Na':2,'Ar':3},
                'actual_input_sha256':state.sha256(application),'native_receipt_sha256':state.sha256(receipt)}

    def case(self, root, fixture):
        import shutil
        import m2_initial_state_evidence as state
        from m2_initial_state import physics_identity_sha256
        support,application,executable,pilot=fixture
        case=root/'PROTOCOL-CASE-NOT-A-NATIVE-RUN';case.mkdir()
        for name in ('application.in','receipt.txt','run_manifest.json','receipt.txt'+state.STATE_SUFFIX):
            shutil.copy2(pilot/name,case/name)
        receipt=case/'receipt.txt'
        receipt.write_text(receipt.read_text().replace(str(pilot/('receipt.txt'+state.STATE_SUFFIX)),str(case/('receipt.txt'+state.STATE_SUFFIX))))
        manifest=json.loads((case/'run_manifest.json').read_text())
        manifest.update(started_utc='2026-11-02T00:00:00Z',inputs={str(case/'application.in'):state.sha256(case/'application.in')},
                        physics_identity_sha256=physics_identity_sha256(case/'application.in'))
        (case/'run_manifest.json').write_text(json.dumps(manifest))
        return case

    def test_old_generic_ready_proof_is_unavailable(self):
        import m2_initial_state_evidence as state
        with self.assertRaises(state.Unavailable):
            state.load_qualification({'actual_epoch_thermal_qualification':{'status':'READY','evidence_files':['generic']}},
                                     Path('not-used'),Path('not-used'),base=Path('.'))

    def test_explicit_uniform_and_periodic_states_bind_actual_grid(self):
        import m2_initial_state_evidence as state
        for mode in ('declared-uniform','periodic-spinup'):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory(prefix='initial-state-protocol-') as tmp:
                root=Path(tmp);fixture=self.fixture(root,mode);qualification=self.qualify(root,fixture)
                actual=state.verify_case_initial_state(self.case(root,fixture),qualification)
                self.assertEqual(actual['initialization_mode'],mode)
                self.assertEqual(actual['shape']['rows'],6)
                self.assertEqual(actual['sha256'],qualification['state_sha256'])

    def test_complete_grid_not_selected_site_average(self):
        import m2_initial_state_evidence as state
        with tempfile.TemporaryDirectory(prefix='initial-grid-protocol-') as tmp:
            root=Path(tmp);fixture=self.fixture(root);grid=fixture[3]/('receipt.txt'+state.STATE_SUFFIX)
            original=grid.read_text()
            for content in ('\n'.join(original.splitlines()[:-1])+'\n',
                            original.replace('1,1,0.3,40','0,1,0.3,40'),
                            original.replace('1,1,0.3,40','1,1,0.4,40'),
                            original.replace('0,0,0.1,40','0,0,0.1,nan'),
                            original.replace('0,-1,0,40','0,-1,0.1,40')):
                grid.write_text(content)
                with self.subTest(content=content),self.assertRaises(ValueError):state.read_state(grid,2,2)

    def test_unapplied_periodic_claim_does_not_qualify_uniform_receipt(self):
        import m2_initial_state_evidence as state
        with tempfile.TemporaryDirectory(prefix='initial-method-protocol-') as tmp:
            root=Path(tmp);fixture=self.fixture(root,'periodic-spinup')
            receipt=fixture[3]/'receipt.txt'
            receipt.write_text(receipt.read_text().replace('thermal_initialization_mode=periodic-spinup',
                                                         'thermal_initialization_mode=declared-uniform'))
            fixture[0]['actual_epoch_thermal_qualification']['pilot_receipt']['sha256']=state.sha256(receipt)
            with self.assertRaisesRegex(ValueError,'mode/epoch'):self.qualify(root,fixture)

    def test_actual_native_state_change_and_physics_change_are_errors(self):
        import m2_initial_state_evidence as state
        for changed in ('state','physics','epoch','review'):
            with self.subTest(changed=changed),tempfile.TemporaryDirectory(prefix='initial-case-protocol-') as tmp:
                root=Path(tmp);fixture=self.fixture(root);qualification=self.qualify(root,fixture);case=self.case(root,fixture)
                if changed=='state':
                    path=case/('receipt.txt'+state.STATE_SUFFIX);path.write_text(path.read_text().replace(',40',',41',1))
                elif changed=='physics':
                    path=case/'application.in';path.write_text(path.read_text().replace('thermal_initial_temperature_k = 40','thermal_initial_temperature_k = 41'))
                elif changed=='epoch':
                    path=case/'receipt.txt';path.write_text(path.read_text().replace('thermal_initialization_et_s=123','thermal_initialization_et_s=124'))
                else:
                    path=case/'run_manifest.json';obj=json.loads(path.read_text());obj['started_utc']='2026-10-01T00:00:00Z';path.write_text(json.dumps(obj))
                with self.assertRaises(ValueError):state.verify_case_initial_state(case,qualification)

    def test_missing_native_export_capability_is_unavailable(self):
        import m2_initial_state_evidence as state
        with tempfile.TemporaryDirectory(prefix='initial-capability-protocol-') as tmp:
            root=Path(tmp);fixture=self.fixture(root);qualification=self.qualify(root,fixture);case=self.case(root,fixture)
            (case/'receipt.txt').write_text('test_id=I19\n')
            with self.assertRaises(state.Unavailable):state.verify_case_initial_state(case,qualification)

    def test_identical_uniform_state_does_not_hide_retargeted_dem_or_kernel(self):
        import m2_initial_state_evidence as state
        for relative in ('NOT-A-DEM.IMG','SPICE-METADATA-FIXTURE/Kernels/NAIF/de421.bsp'):
            with self.subTest(file=relative),tempfile.TemporaryDirectory(prefix='physical-data-protocol-') as tmp:
                root=Path(tmp);fixture=self.fixture(root);q=self.qualify(root,fixture);case=self.case(root,fixture)
                before=state.sha256(case/('receipt.txt'+state.STATE_SUFFIX))
                physical=root/relative;physical.write_text('CHANGED PHYSICAL METADATA FIXTURE; NOT LUNAR DATA\n')
                manifest=case/'run_manifest.json';record=json.loads(manifest.read_text())
                name=str(physical) if relative=='NOT-A-DEM.IMG' else 'NAIF/de421.bsp'
                record['external_data'][name]=state.sha256(physical);manifest.write_text(json.dumps(record))
                self.assertEqual(before,state.sha256(case/('receipt.txt'+state.STATE_SUFFIX)))
                with self.assertRaisesRegex(ValueError,'identities differ'):state.verify_case_initial_state(case,q)

    def test_unrelated_i20_observation_entries_do_not_change_initial_physics(self):
        import m2_initial_state_evidence as state
        with tempfile.TemporaryDirectory(prefix='physical-extra-protocol-') as tmp:
            root=Path(tmp);fixture=self.fixture(root);q=self.qualify(root,fixture);case=self.case(root,fixture)
            unrelated=root/'OBSERVATION-TRANSPORT-FIXTURE-NOT-DATA.txt';unrelated.write_text('CONTROL FILE ONLY\n')
            manifest=case/'run_manifest.json';record=json.loads(manifest.read_text())
            record['external_data'][str(unrelated)]=state.sha256(unrelated);manifest.write_text(json.dumps(record))
            self.assertEqual(state.verify_case_initial_state(case,q)['sha256'],q['state_sha256'])

    def test_actual_kernel_pool_order_and_old_missing_provenance_are_distinct(self):
        import m2_initial_state_evidence as state
        for missing in (False,True):
            with self.subTest(missing=missing),tempfile.TemporaryDirectory(prefix='actual-kernel-protocol-') as tmp:
                root=Path(tmp);fixture=self.fixture(root);q=self.qualify(root,fixture);case=self.case(root,fixture)
                receipt=case/'receipt.txt';text=receipt.read_text()
                if missing:
                    text='\n'.join(line for line in text.splitlines() if not line.startswith('thermal_initialization_spice_kernel_'))+'\n'
                else:
                    first=next(line for line in text.splitlines() if line.startswith('thermal_initialization_spice_kernel_0_file='))
                    second=next(line for line in text.splitlines() if line.startswith('thermal_initialization_spice_kernel_1_file='))
                    text=text.replace(first,'ORDER-SWAP-FIXTURE').replace(second,second.split('=',1)[0]+'='+first.split('=',1)[1]).replace(
                        'ORDER-SWAP-FIXTURE',first.split('=',1)[0]+'='+second.split('=',1)[1])
                receipt.write_text(text)
                with self.assertRaises(state.Unavailable if missing else ValueError):state.verify_case_initial_state(case,q)

    def test_actual_kernel_filename_controls_even_outside_predicted_root(self):
        import shutil
        import m2_initial_state_evidence as state
        with tempfile.TemporaryDirectory(prefix='actual-root-protocol-') as tmp:
            root=Path(tmp);fixture=self.fixture(root);support,application,executable,pilot=fixture
            old=root/'SPICE-METADATA-FIXTURE/Kernels/NAIF/de421.bsp'
            actual=root/'ACTUAL-NATIVE-POOL-FIXTURE/de421.bsp';actual.parent.mkdir();shutil.copy2(old,actual)
            receipt=pilot/'receipt.txt';receipt.write_text(receipt.read_text().replace(str(old),str(actual)))
            support['actual_epoch_thermal_qualification']['pilot_receipt']['sha256']=state.sha256(receipt)
            q=self.qualify(root,fixture);case=self.case(root,fixture)
            self.assertIn(str(actual),[row['path'] for row in q['physical_external_data']['ordered_actual_spice_kernels']])
            self.assertEqual(state.verify_case_initial_state(case,q)['sha256'],q['state_sha256'])

    def test_numerical_recipe_and_frozen_site_mapping_cannot_change(self):
        import m2_initial_state_evidence as state
        for before,after in (('m2_verification_cohort_packets = 32','m2_verification_cohort_packets = 64'),
                             ('m2_verification_source_integral_tolerance = 1e-8','m2_verification_source_integral_tolerance = 1e-7'),
                             ('m2_verification_surface_elements = 0,1','m2_verification_surface_elements = 1,0')):
            with self.subTest(changed=before),tempfile.TemporaryDirectory(prefix='initial-recipe-protocol-') as tmp:
                root=Path(tmp);fixture=self.fixture(root);q=self.qualify(root,fixture);case=self.case(root,fixture)
                application=case/'application.in';application.write_text(application.read_text().replace(before,after))
                with self.assertRaisesRegex(ValueError,'recipe'):state.verify_case_initial_state(case,q)

    def test_normalization_requires_actual_preserved_physics_native_pilot(self):
        import m2_initial_state_evidence as state
        from m2_initial_state import physics_identity_sha256
        with tempfile.TemporaryDirectory(prefix='normalization-physics-protocol-') as tmp:
            root=Path(tmp);fixture=self.fixture(root);support,application,executable,pilot=fixture
            with self.assertRaises(state.Unavailable):state.verify_normalization_physics({},application,executable,base=root,mpi_ranks=8)
            descriptor=lambda path:{'path':str(path),'sha256':state.sha256(path)}
            norms={'production_input_tier':'cluster','physics_identity_sha256':physics_identity_sha256(application),
                'template_physics_identity_sha256':physics_identity_sha256(application),
                'actual_input':descriptor(application),'native_manifest':descriptor(pilot/'run_manifest.json'),
                'native_receipt':descriptor(pilot/'receipt.txt'),'native_weight_provenance':self.weights(application,pilot/'receipt.txt')}
            self.assertEqual(state.verify_normalization_physics(norms,application,executable,base=root,mpi_ranks=8)['physics_identity_sha256'],norms['physics_identity_sha256'])
            norms['physics_identity_sha256']='c'*64
            with self.assertRaises(ValueError):state.verify_normalization_physics(norms,application,executable,base=root,mpi_ranks=8)

    def test_native_normalization_changed_bytes_are_errors(self):
        import m2_initial_state_evidence as state
        from m2_initial_state import physics_identity_sha256
        with tempfile.TemporaryDirectory(prefix='normalization-bytes-protocol-') as tmp:
            root=Path(tmp);fixture=self.fixture(root);support,application,executable,pilot=fixture
            descriptor=lambda path:{'path':str(path),'sha256':state.sha256(path)}
            norms={'production_input_tier':'cluster','physics_identity_sha256':physics_identity_sha256(application),
                'template_physics_identity_sha256':physics_identity_sha256(application),
                'actual_input':descriptor(application),'native_manifest':descriptor(pilot/'run_manifest.json'),
                'native_receipt':descriptor(pilot/'receipt.txt'),'native_weight_provenance':self.weights(application,pilot/'receipt.txt')}
            (pilot/'receipt.txt').write_text('CHANGED CONTROL FIXTURE ONLY')
            with self.assertRaisesRegex(ValueError,'SHA-256 mismatch'):state.verify_normalization_physics(norms,application,executable,base=root,mpi_ranks=8)

    def test_i20_own_phase_is_not_implicitly_the_cohort_initial_epoch(self):
        import m2_initial_state_evidence as state
        from integration_common import parse_receipt
        with tempfile.TemporaryDirectory(prefix='i20-phase-protocol-') as tmp:
            root=Path(tmp);fixture=self.fixture(root,'periodic-spinup');cohort=self.qualify(root,fixture)
            support,application,executable,pilot=fixture
            receipt=parse_receipt(pilot/'receipt.txt');manifest=json.loads((pilot/'run_manifest.json').read_text())
            receipt.update(initial_thermal_state_file=str(pilot/('receipt.txt'+state.STATE_SUFFIX)),
                initial_thermal_surface_elements=2,initial_thermal_layer_count=2,
                initial_thermal_state_tdb_s=456,thermal_initialization_et_s=456,
                thermal_initialization_spinup_cycle_start_et_s=455,
                thermal_initialization_spinup_cycle_end_et_s=456)
            own={**cohort,'ephemeris_time_tdb_s':456,'qualification_sha256':'c'*64}
            result=state.verify_surface_initial_state(receipt,manifest,application,
                pilot/('receipt.txt'+state.STATE_SUFFIX),own)
            self.assertEqual(result['qualification_sha256'],'c'*64)
            with self.assertRaises(ValueError):state.verify_surface_initial_state(receipt,manifest,application,
                pilot/('receipt.txt'+state.STATE_SUFFIX),cohort)

    def test_i20_real_pilot_schema_requires_explicit_surface_permission(self):
        import m2_initial_state_evidence as state
        with tempfile.TemporaryDirectory(prefix='i20-pilot-protocol-') as tmp:
            root=Path(tmp);fixture=self.fixture(root,'periodic-spinup')
            support,application,executable,pilot=fixture;q=support['actual_epoch_thermal_qualification']
            grid=pilot/'receipt.txt.I20_initial_thermal_state.csv'
            (pilot/('receipt.txt'+state.STATE_SUFFIX)).rename(grid)
            receipt=pilot/'receipt.txt';receipt.write_text(receipt.read_text().replace('test_id=I19','test_id=I20')+
                'initial_thermal_state_file='+str(grid)+'\ninitial_thermal_surface_elements=2\n'
                'initial_thermal_layer_count=2\ninitial_thermal_state_tdb_s=123\n')
            manifest=pilot/'run_manifest.json';value=json.loads(manifest.read_text());value['test_id']='I20';manifest.write_text(json.dumps(value))
            q['native_state']={'path':str(grid),'sha256':state.sha256(grid)}
            q['pilot_receipt']['sha256']=state.sha256(receipt);q['pilot_manifest']['sha256']=state.sha256(manifest)
            with patch.object(state,'_source_signature',return_value='b'*64):
                with self.assertRaises(ValueError):state.load_qualification(support,application,executable,base=root)
                qualified=state.load_qualification(support,application,executable,base=root,allow_surface_pilot=True)
            self.assertEqual(qualified['pilot_test_id'],'I20')

    def test_all_288_states_gate_before_any_numerical_difference(self):
        import m2_initial_state_evidence as state
        from m2_equivalence import score_campaign
        selected=cases('cluster')
        def state_check(directory,q):
            if directory.name==selected[-1]['key']:raise ValueError('last applied state differs')
            return {'fixture':'MODEL-IDENTITY-CONTROL-ONLY'}
        support={'status':'READY','marginal_intervals_justified':True,'rng_stream_policy_reviewed':True}
        with patch.object(state,'load_qualification',return_value={'shape':{'logical_surface_element_count':1}}), \
             patch.object(state,'verify_case_initial_state',side_effect=state_check) as checked, \
             patch('m2_equivalence.load_primary',side_effect=AssertionError('numerical values must remain unopened')):
            result=score_campaign(Path('PROTOCOL-NOT-NATIVE'),{},support,application_input=Path('unused'),
                                  executable=Path('unused'),support_base=Path('.'),sites=[{'logical_surface_element':0}])
        self.assertEqual(result['status'],'ERROR')
        self.assertEqual(checked.call_count,288)
        self.assertNotIn('comparisons',result)

    def test_initialization_cost_is_not_a_free_cached_state(self):
        import m2_initial_state_evidence as state
        with tempfile.TemporaryDirectory(prefix='initial-cost-protocol-') as tmp:
            root=Path(tmp);fixture=self.fixture(root);q=self.qualify(root,fixture)
            with self.assertRaises(state.Unavailable):state.verify_initialization_forecast({},q,mpi_ranks=8,case_count=288)
            evidence=q['qualification']['source_evidence_files'][0]
            forecast={'allocated_core_hours':10,'initialization_cost':{
                'included_in_allocated_core_hours':True,'reexecuted_for_every_case':True,
                'qualification_sha256':q['qualification_sha256'],'native_state_sha256':q['state_sha256'],
                'initialization_mode':q['initialization_mode'],'case_count':288,'measurement_mpi_ranks':8,
                'measured_per_case_wall_seconds':1,'allocated_initialization_core_hours':1,
                'evidence_files':[evidence]}}
            self.assertEqual(state.verify_initialization_forecast(forecast,q,mpi_ranks=8,case_count=288)['case_count'],288)
            forecast['initialization_cost']['allocated_initialization_core_hours']=0
            with self.assertRaises(ValueError):state.verify_initialization_forecast(forecast,q,mpi_ranks=8,case_count=288)


if __name__ == "__main__":
    unittest.main(verbosity=2)
