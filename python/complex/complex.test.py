import unittest

import numpy as np

import complex as c


class MyTestCase(unittest.TestCase):
    f_4 = c.dft_matrix(4)
    f_5 = c.dft_matrix(5)
    f_10 = c.dft_matrix(10)

    eye_4 = np.eye(4)

    rnd_matrix = np.random.randn(4, 3)

    def test_dft(self):
        # TODO add test cases
        self.assertEqual(True, True)  # add assertion here

    def test_unitary(self):
        self.assertTrue(c.is_unitary(self.f_4))
        self.assertTrue(c.is_unitary(self.f_5))
        self.assertTrue(c.is_unitary(self.f_10))

        self.assertTrue(c.is_unitary(self.eye_4))
        self.assertFalse(c.is_unitary(self.rnd_matrix))

if __name__ == '__main__':
    unittest.main()
