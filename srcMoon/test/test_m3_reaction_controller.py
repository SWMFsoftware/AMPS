#!/usr/bin/env python3
"""Protocol-only resource-refusal controls; these fixtures are not observations.

The actual native runner uses the same predicate. Merely retaining expected
fatal text cannot turn a signal/launcher/parser failure into a passing control.
"""
import unittest
from run_m3_reaction_controller import expected_capacity_refusal


class CapacityRefusalProtocol(unittest.TestCase):
    text=('complete reaction packet cannot fit available PIC particle slots; '
          'Increase particle-buffer capacity; without discarding products')

    def test_only_documented_native_exit_is_accepted(self):
        self.assertTrue(expected_capacity_refusal(1,self.text))

    def test_success_signal_and_launcher_failures_remain_errors(self):
        for code in (0,-11,139,-9,137,2,127,True,False,None,'1'):
            with self.subTest(code=code):
                self.assertFalse(expected_capacity_refusal(code,self.text))

    def test_missing_actionable_diagnostic_is_not_the_named_refusal(self):
        for required in ('complete reaction packet cannot fit available PIC particle slots',
                         'Increase particle-buffer capacity','without discarding products'):
            with self.subTest(required=required):
                self.assertFalse(expected_capacity_refusal(1,self.text.replace(required,'')))


if __name__=='__main__':
    unittest.main()
