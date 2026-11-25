"""Administrative scope controls; these fixtures contain no scientific data."""
import unittest
import run_all_tests as runner


class RunnerSelection(unittest.TestCase):
    def test_protocol_scope_excludes_declared_external_sources_before_run(self):
        class Administrative(unittest.TestCase):
            def test_control(self):
                pass

        class ActualSource(Administrative):
            requires_authoritative_source = True

        loaded = unittest.TestSuite([unittest.TestSuite([Administrative('test_control')]),
                                    ActualSource('test_control')])
        excluded = []
        selected = runner.select_additional_tests(loaded, True, excluded)
        self.assertEqual(selected.countTestCases(), 1)
        self.assertEqual(excluded, [ActualSource('test_control').id()])

    def test_full_scope_keeps_external_source_case_for_honest_skip_or_run(self):
        class ActualSource(unittest.TestCase):
            requires_authoritative_source = True

            def test_control(self):
                pass

        excluded = []
        selected = runner.select_additional_tests(
            unittest.TestSuite([ActualSource('test_control')]), False, excluded)
        self.assertEqual(selected.countTestCases(), 1)
        self.assertEqual(excluded, [])

    def test_protocol_scope_records_unselected_production_owner(self):
        class Native(unittest.TestCase):
            def test_native(self):
                pass

        class PythonApi(unittest.TestCase):
            def test_api(self):
                pass

        class Production(unittest.TestCase):
            def test_release_a(self):
                pass

            def test_release_b(self):
                pass

        owners = (Native, PythonApi, Production)
        excluded = []
        selected = runner.select_owner_suites(owners, True, excluded)
        self.assertEqual(selected.countTestCases(), 2)
        self.assertEqual(excluded, [Production('test_release_a').id(), Production('test_release_b').id()])
        excluded = []
        selected = runner.select_owner_suites(owners, False, excluded)
        self.assertEqual(selected.countTestCases(), 4)
        self.assertEqual(excluded, [])

    def test_declared_source_classes_are_marked_consistently(self):
        from test_native_cea import ActualSource
        from test_native_exomol import ExoMolRealSource
        from test_native_surface_sources import NativeSurfaceSource
        for owner in (ActualSource, ExoMolRealSource, NativeSurfaceSource):
            self.assertIs(owner.requires_authoritative_source, True)
