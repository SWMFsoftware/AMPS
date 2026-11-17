#!/usr/bin/env python3
"""Harness-only smoke gates; fake roles never count as native physics evidence."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile
import textwrap
import unittest

RUNNER = Path(__file__).with_name("run_fresh_process_checkpoint_tests.py")


class PairedRunner(unittest.TestCase):
    def setUp(self) -> None:
        self.scope = tempfile.TemporaryDirectory(prefix="srcsep paired harness ")
        self.root = Path(self.scope.name)
        self.launcher = self.root / "launcher"
        self.launcher.write_text("#!/usr/bin/env python3\nimport os,sys\nos.execv(sys.argv[3],sys.argv[3:])\n")
        self.launcher.chmod(0o755)
        self.native = self.root / "native"
        self.native.write_text(textwrap.dedent('''\
            #!/usr/bin/env python3
            from pathlib import Path
            import json,shlex,sys
            args=sys.argv[1:]
            get=lambda key:args[args.index(key)+1]
            case=get('--test')
            tokens=shlex.split(Path(get('--test-input')).read_text())
            role,prefix=tokens[3],Path(tokens[5])
            Path(str(prefix)+'.work-'+role).write_text(str(Path.cwd()))
            if role=='write':
                Path(str(prefix)+'.checkpoint').write_bytes(b'state-one')
                Path(str(prefix)+'.expected').write_bytes(b'state-three')
            else:
                assert Path(str(prefix)+'.checkpoint').read_bytes()==b'state-one'
                Path(str(prefix)+'.resumed').write_bytes(b'wrong' if case=='CP05F03' else b'state-three')
            status='SKIP' if case=='CP05F02' and role=='resume' else 'PASS'
            Path(get('--test-json')).write_text(json.dumps({'results':[{'id':case,'status':status}]}))
            '''))
        self.native.chmod(0o755)

    def tearDown(self) -> None:
        self.scope.cleanup()

    def run_case(self, case: str):
        output = self.root / case
        result = subprocess.run([sys.executable, "-B", str(RUNNER), "--amps", str(self.native),
                                 "--mpi-exec", str(self.launcher), "--cases", case,
                                 "--output", str(output)], text=True, capture_output=True)
        return result, output, json.loads((output / "summary.json").read_text())

    def test_separate_process_scratch_and_complete_artifact_match(self) -> None:
        result, output, summary = self.run_case("CP05F01")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(summary["results"][0]["status"], "PASS")
        prefix = output / "CP05F01" / "physics"
        self.assertNotEqual(Path(str(prefix)+".work-write").read_text(),
                            Path(str(prefix)+".work-resume").read_text())
        self.assertEqual(len(summary["executable_sha256"]), 64)

    def test_zero_exit_skip_cannot_become_native_pass(self) -> None:
        result, _, summary = self.run_case("CP05F02")
        self.assertEqual(result.returncode, 1)
        self.assertEqual(summary["results"][0]["status"], "FAIL")

    def test_two_native_pass_reports_do_not_hide_artifact_mismatch(self) -> None:
        result, _, summary = self.run_case("CP05F03")
        self.assertEqual(result.returncode, 1)
        self.assertEqual(summary["results"][0]["status"], "FAIL")
        self.assertIn("differ", summary["results"][0]["message"])


if __name__ == "__main__":
    unittest.main()
