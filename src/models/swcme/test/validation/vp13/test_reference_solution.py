import unittest
from reference_solution import relative_difference
class ReferenceTests(unittest.TestCase):
    def test_symmetric(self): self.assertEqual(relative_difference(2,3),relative_difference(3,2))
    def test_identity(self): self.assertEqual(relative_difference(7,7),0)
if __name__ == "__main__": unittest.main(verbosity=2)

