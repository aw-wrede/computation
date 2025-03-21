import unittest
import numpy as np

import matrix as mtrx

class MatrixTest(unittest.TestCase):

    a = np.random.random((3, 5))
    b = np.random.random((5, 3))
    c = np.random.random((4, 4))

    def test_dimension(self):
        self.assertTrue(np.allclose(mtrx.dot(self.a, self.b), np.dot(self.a, self.b)))

        with self.assertRaises(ValueError):
            mtrx.dot(self.a, self.c)

    def test_dot(self):
        self.assertTrue(np.allclose(mtrx.dot(self.c, self.c), np.dot(self.c, self.c)))

        big = np.random.random((100, 100))
        self.assertTrue(np.allclose(mtrx.dot(big, big), np.dot(big, big)))

    def test_close(self):
        a1 = np.random.random((3, 5))
        a2 = np.random.random((3, 5))
        self.assertTrue(np.allclose(mtrx.close(a1, a2, 1e-05, 1e-08), np.isclose(a1, a2, 1e-05, 1e-08)))
        with self.assertRaises(ValueError):
            mtrx.close(self.a, self.b)

if __name__ == '__main__':
    unittest.main()
