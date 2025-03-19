import unittest
import numpy as np

import linalg


class MyTestCase(unittest.TestCase):

    a = np.array([[2, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 1, 0, 0]], dtype=np.float64)
    b = np.array([1, 2, 3, 4], dtype=np.float64)

    def test_gaussian_elimination(self):
        reference_a = np.array([[2, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -0.5]], dtype=np.float64)
        reference_b = np.array([1, 3, 2, 0.5], dtype=np.float64)

        g_a, g_b = linalg.gaussian_elimination(self.a, self.b)

        self.assertTrue(np.allclose(reference_a, g_a) and np.allclose(reference_b, g_b))


if __name__ == '__main__':
    unittest.main()
