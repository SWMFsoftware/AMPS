#!/usr/bin/env python3
"""Harness-only failures: fake subprocesses never constitute native evidence."""
from pathlib import Path
import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest

RUNNER = Path(__file__).with_name("run_native_host_step_tests.py")
FIXTURE = Path(__file__).with_name("fixtures") / "host_step_smoke.json"


class HostRunner(unittest.TestCase):
    def setUp(self):
        self.scope = tempfile.TemporaryDirectory(prefix="srcsep-host-runner ")
        self.root = Path(self.scope.name)
        self.launcher = self.root / "launcher"
        self.launcher.write_text("#!/usr/bin/env python3\nimport os,sys\nos.execv(sys.argv[3],sys.argv[3:])\n")
        self.launcher.chmod(0o755)
        self.native = self.root / "native"

    def tearDown(self):
        self.scope.cleanup()

    def run_fixture(self, defect=""):
        # The only defect selector lives in this disposable test executable;
        # it is not an environment hook or failing production registry case.
        program = textwrap.dedent('''\
            #!/usr/bin/env python3
            from pathlib import Path
            import json,sys
            DEFECT=DEFECT_LITERAL
            args=sys.argv[1:]
            value=lambda flag:args[args.index(flag)+1]
            if '--test-input' in args and '--test-output-dir' not in args:
                print('--test-input and --test-output-dir must be provided together')
                raise SystemExit(1)
            if args.count('--test')>1:
                print('The mesh is generated' if DEFECT=='guard' else
                      'end-to-end test input requires exactly one explicit --test ID')
                raise SystemExit(0 if DEFECT=='guard' else 1)
            if '--test-repeat' in args:
                print('The mesh is generated' if DEFECT=='guard' else 'requires one case in a fresh process')
                raise SystemExit(0 if DEFECT=='guard' else 2)
            case=value('--test')
            target=Path(value('--test-json'));target.parent.mkdir(parents=True,exist_ok=True)
            artifact=Path(value('--test-output-dir'))
            mover={'HOST01':'parker','HOST02':'fte-dmumu','HOST03':'fte-mfp'}[case]
            for rank in range(2):
                if DEFECT=='missing-trace' and case=='HOST02' and rank==1:continue
                path=artifact if rank==0 else artifact/('rank-'+str(rank));path.mkdir(parents=True,exist_ok=True)
                events=[]
                for step in range(2):
                    event={'step':step,'particle':100+rank+1000*step,'phase':'source','clock_s':0.5*step,'parent_start_s':0.5*step,
                           'parent_end_s':0.5*(step+1),'first_duration_s':0.5,'requested_s':0.0,'source_radius_m':1.4e10,
                           'pending_first':True,'cartesian_matches':not(DEFECT in ('bad-x','native-fail') and case=='HOST02'),
                           'canonical_amr_block_local':True}
                    events.extend([event,dict(event,phase='mover',requested_s=0.5)])
                if DEFECT=='duplicate-first' and case=='HOST02':events.append(events[-1].copy())
                for sequence in range(2):
                    if DEFECT=='incomplete' and case=='HOST02' and sequence==1:continue
                    document={'schema':'srcsep-native-host-step-trace-v2','case_id':case,'mover':mover,'rank':rank,
                              'snapshot_sequence':sequence,'complete':sequence==1,'requested_host_steps':2,
                              'completed_host_steps':2 if sequence else 0,'counter_at_end':True,
                              'events':events if sequence else events[:1]}
                    if DEFECT=='rewrite-prefix' and case=='HOST02' and sequence==1:
                        document=json.loads(json.dumps(document));document['events'][0]['source_radius_m']+=100
                    (path/('host-step-trace-'+format(sequence,'06d')+'.json')).write_text(json.dumps(document))
            metrics=[{'name':'assertion_failures','value':1 if DEFECT in ('native-fail','fail-exit') and case=='HOST02' else 0},{'name':'completed_host_steps','value':2},
                     {'name':'global_first_birth_moves','value':2}]
            status='FAIL' if DEFECT in ('native-fail','fail-exit') and case=='HOST02' else ('SKIP' if DEFECT=='skip' and case=='HOST02' else 'PASS')
            if DEFECT!='missing-report' or case!='HOST02':
                target.write_text(json.dumps({'results':[{'id':case,'status':status,'metrics':metrics}]}))
            (target.parent/'actual-cwd.txt').write_text(str(Path.cwd()))
            if DEFECT=='binary' and case=='HOST01':
                executable=Path(sys.argv[0]);executable.write_text(executable.read_text()+'\\n# mutation\\n')
            if DEFECT=='protocol' and case=='HOST01':
                protocol=Path(value('--test-input'));protocol.write_text(protocol.read_text()+'source_mu 0.25\\n')
            raise SystemExit(7 if DEFECT=='exit' and case=='HOST02' else (1 if DEFECT=='native-fail' and case=='HOST02' else 0))
            ''').replace("DEFECT_LITERAL", repr(defect))
        self.native.write_text(program);self.native.chmod(0o755)
        output = self.root / ("output-" + (defect or "pass"))
        process = subprocess.run([sys.executable, "-B", str(RUNNER), "--amps", str(self.native),
                                  "--fixture", str(FIXTURE), "--mpi-exec", str(self.launcher),
                                  "--output", str(output)], text=True, capture_output=True,
                                 env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        return process, output, json.loads((output / "summary.json").read_text())

    def test_success_preserves_separate_native_jobs_and_all_child_evidence(self):
        process, output, report = self.run_fixture()
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        self.assertEqual([row["status"] for row in report["results"]], ["PASS"] * 3)
        self.assertEqual([row["status"] for row in report["preflight_checks"]], ["PASS"] * 2)
        self.assertEqual([row["exit_code"] for row in report["preflight_checks"]], [1, 2])
        self.assertTrue(all('--test-output-dir' in row['command'] for row in report['preflight_checks']))
        working = [(output / case / "actual-cwd.txt").read_text() for case in ("HOST01", "HOST02", "HOST03")]
        self.assertEqual(len(set(working)), 3)
        self.assertTrue(all(str(output) in path for path in working))
        self.assertEqual(len(report["executable_sha256"]), 64)

    def test_native_pass_cannot_hide_bad_exit_skip_missing_evidence_or_phase(self):
        for defect in ("exit", "skip", "missing-report", "missing-trace", "bad-x", "duplicate-first", "incomplete", "rewrite-prefix", "fail-exit"):
            with self.subTest(defect=defect):
                process, _, report = self.run_fixture(defect)
                self.assertEqual(process.returncode, 1, process.stdout + process.stderr)
                rows = {row["id"]: row for row in report["results"]}
                self.assertEqual(rows["HOST02"]["status"], "ERROR")
                self.assertTrue(rows["HOST02"]["message"])
                self.assertEqual(rows["HOST03"]["status"], "PASS")
                self.assertTrue(Path(rows["HOST02"]["log"]).is_file())

    def test_real_native_fail_is_retained_with_complete_immutable_trace(self):
        process, _, report = self.run_fixture("native-fail")
        self.assertEqual(process.returncode, 1)
        self.assertEqual([row["status"] for row in report["results"]], ["PASS", "FAIL", "PASS"])
        self.assertEqual(report["results"][1]["exit_code"], 1)

    def test_binary_change_or_missing_exclusivity_blocks_green_results(self):
        for defect in ("binary", "protocol", "guard"):
            with self.subTest(defect=defect):
                process, _, report = self.run_fixture(defect)
                self.assertEqual(process.returncode, 1)
                self.assertTrue(all(row["status"] == "ERROR" for row in report["results"]))
                if defect == "guard":
                    self.assertTrue(all(row["exit_code"] is None for row in report["results"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
