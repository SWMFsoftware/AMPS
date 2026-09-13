import unittest
from reference_solution import profile,profile_log_form
class ReferenceTests(unittest.TestCase):
    def test_paths_agree(self): self.assertAlmostEqual(profile(3,20,4.2),profile_log_form(3,20,4.2),places=14)
    def test_positive(self): self.assertGreater(profile(1,100,4.2),0)
if __name__ == "__main__": unittest.main(verbosity=2)

