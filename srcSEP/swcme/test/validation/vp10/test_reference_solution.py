import unittest
from reference_solution import parker_metrics,SOURCE,AU_M
class ReferenceTests(unittest.TestCase):
    def test_radial_limit(self):
        p,f=parker_metrics(0,1); self.assertAlmostEqual(p,1-SOURCE/AU_M); self.assertAlmostEqual(f,.5)
if __name__ == "__main__": unittest.main(verbosity=2)

