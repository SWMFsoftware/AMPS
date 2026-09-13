import unittest
from reference_solution import intersection_time_h
class ReferenceTests(unittest.TestCase):
    def test_apex(self): self.assertEqual(intersection_time_h(0,45,42),42)
    def test_outside_cap(self): self.assertIsNone(intersection_time_h(46,45,42))
if __name__ == "__main__": unittest.main(verbosity=2)

