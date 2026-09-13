import unittest
from reference_solution import intensity_ratio
class ReferenceTests(unittest.TestCase):
    def test_normalization(self): self.assertEqual(intensity_ratio(20,20,4.2),1.0)
    def test_falling_spectrum(self): self.assertGreater(intensity_ratio(10,20,4.2),1.0)
if __name__ == "__main__": unittest.main(verbosity=2)

