import unittest
from reference_solution import skill_score,median
class ReferenceTests(unittest.TestCase):
    def test_improvement_is_positive(self): self.assertGreater(skill_score(4,8),0)
    def test_median(self): self.assertEqual(median([3,1,2]),2)
if __name__ == "__main__": unittest.main(verbosity=2)

