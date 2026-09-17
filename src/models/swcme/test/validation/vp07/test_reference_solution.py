import unittest
from reference_solution import theta_bn_deg,dsa_q
class ReferenceTests(unittest.TestCase):
    def test_geometry_is_acute(self): self.assertTrue(0.0 <= theta_bn_deg(20.0) <= 90.0)
    def test_strong_shock_slope(self): self.assertAlmostEqual(dsa_q(4.0),4.0)
if __name__ == "__main__": unittest.main(verbosity=2)

