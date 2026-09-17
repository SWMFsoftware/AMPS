import unittest
from reference_solution import expected_history
class ReferenceTests(unittest.TestCase):
    def test_sta_persistent(self): self.assertTrue(all(expected_history("STEREO-A")))
    def test_stb_loses(self):
        h=expected_history("STEREO-B"); self.assertTrue(h[0] and not h[-1])
if __name__ == "__main__": unittest.main(verbosity=2)

