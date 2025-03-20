import unittest
import numpy as np

import linalg


class MyTestCase(unittest.TestCase):

    a = np.array([[2, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 1, 0, 0]], dtype=np.float64)
    b = np.array([1, 2, 3, 4], dtype=np.float64)

    def test_gaussian_elimination(self):
        reference_a = np.array([[2, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -0.5]], dtype=np.float64)
        reference_b = np.array([1, 3, 2, 0.5], dtype=np.float64)

        g_a = self.a.copy()
        g_b = self.b.copy()
        linalg.gaussian_elimination(g_a, g_b)

        self.assertTrue(np.allclose(reference_a, g_a))
        self.assertTrue(np.allclose(reference_b, g_b))
        self.assertTrue(np.allclose(np.triu(g_a), g_a))


    def test_back_substitution(self):
        reference_b = np.array([1, 3, 2, -1], dtype=np.float64)

        g_a = self.a.copy()
        g_b = self.b.copy()

        linalg.gaussian_elimination(g_a, g_b)

        x = linalg.back_substitution(g_a, g_b)

        self.assertTrue(np.allclose(x, reference_b))

    def test_solve(self):
        a = np.random.randn(5, 5)
        reference_x = np.random.rand(5)
        b = np.dot(a, reference_x)

        x = linalg.solve(a, b)

        self.assertTrue(np.allclose(reference_x, x))


if __name__ == '__main__':
    unittest.main()
