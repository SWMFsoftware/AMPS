import unittest
from reference_solution import gaussian_spread_ratio
class ReferenceTests(unittest.TestCase):
    def test_peak(self): self.assertEqual(gaussian_spread_ratio(0,0.08),1)
    def test_symmetry(self): self.assertEqual(gaussian_spread_ratio(-30,.08),gaussian_spread_ratio(30,.08))
if __name__ == "__main__": unittest.main(verbosity=2)

