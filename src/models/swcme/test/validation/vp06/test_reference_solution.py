import unittest
from reference_solution import expected_residual
class ReferenceTests(unittest.TestCase):
    def test_exact_target_is_zero(self): self.assertEqual(expected_residual(), 0.0)
if __name__ == "__main__": unittest.main(verbosity=2)

