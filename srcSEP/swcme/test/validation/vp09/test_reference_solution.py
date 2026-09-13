import unittest
from reference_solution import classify
class ReferenceTests(unittest.TestCase):
    def test_boundaries(self): self.assertEqual([classify(x) for x in (-1,0,.999,1,2)],[0,1,1,2,0])
if __name__ == "__main__": unittest.main(verbosity=2)

